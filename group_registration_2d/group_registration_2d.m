function [new_im,res] = group_registration_2d(im,main,optim)
% GROUP_REGISTRATION_2D: Main function for the deformable 2d groupwise
% registration method for the MOLLI image sequences based on dictionary
% matching
%
%   Input
%   ------------------
%   im      - MOLLI image sequence with size Ny x Nx x Nt
%   main    - Structure of main options:
%             .D (no default value): Dictionary for matching
%             .rho (default 0.7): Coefficient to balance the dictionary- 
%              matching (DM) loss and low-rankness (LR) loss, which ranges
%              from 0 to 1. 0: only the DM loss; 1: only the LR loss 
%             .okno (default 6): Space between adjacent control points
%             .subdivide (default 1): Number of resolution levels
%             .regularizer(default 'second_order'): Regularization type.
%              Choose from 'first-order' (diffusion) and 'second-order' (
%              bending energy of a thin-plate metal)
%             .lambda1 (default 1e-4): Spatial regularization weight
%             .lambda2 (default 0): Temporal regularization weight
%             .disp_flag (default 0): Flag for the display of control meshes.
%              Set disp_flag to anything nonzero to display the deformation
%              of control meshes during registration
%   optim   - Structure of optimization options:
%             .max_alt_iter (default [16,8,4]): Maximum number of alternation
%              iterations at each resolution level
%             .max_reg_iter (default 81): Maximum number of registration
%              iterations at each alternation
%             .alt_tol (default 1e-6): Tolerance for alternation termination
%             .reg_tol (default 1e-6): Tolerance for registration termination
%             .gamma (default 10): Initial step size for registration
%             .anneal (default 0.9): Decay rate for step size for registration
%
%   Output
%   ------------------
%   new_im   - Registered MOLLI image sequence with size Ny x Nx x Nt
%   res      - Structure of resulting parameters:
%             .mesh_x/y: X/Y-coordinates of final control points with size
%              Gy x Gx x Nt
%             .okno: Space between two adjacent control points, which is
%              equal to that value in main options
%
%
% Copyright: Haiyang Chen, Chenxi Hu, SJTU 2023

% Check the options and set the defaults
if nargin < 1
    error('Error: The images must be provided!')
elseif nargin < 2
    error('Error: The main options must be provided!')
end
if ~exist('optim','var')
    optim = struct();
end
[main,optim] = check_options(main,optim,nargin);

% tic
fprintf('Non-rigid group registration with following options:\n');
main, optim

% Image size and frame number
[Ny,Nx,Nt] = size(im);

% Small image size at the lowest resolution level
Mx = ceil(Nx/2^(main.subdivide-1));
My = ceil(Ny/2^(main.subdivide-1));

% Initialize the B-splines meshes of control points equally spaced with
% mian.okno (control meshes) at the lowest resolution level
x = (1-main.okno):main.okno:(Mx+2*main.okno);
y = (1-main.okno):main.okno:(My+2*main.okno);
[X,Y] = meshgrid(x,y);
main.mesh_x = repmat(X,[1,1,Nt]);
main.mesh_y = repmat(Y,[1,1,Nt]);
% Also save the initial mesh (used for regularization)
main.init_mesh_x = X;
main.init_mesh_y = Y;

% Size of the control meshes at the lowest resolution level
[main.Gy,main.Gx] = size(X);

% The new small image size at the lowest resolution level is equal or bigger
% than the previous one
Mx = (main.Gx-3)*main.okno;
My = (main.Gy-3)*main.okno;

% Generate the B-splines coefficient matrix F
main.F = generate_F(main.okno); % size: main.okno^2 x 4^2

% At the original large image size, calculate the size described by control
% points, which will be equal to or larger than the original image size.
% Create a zero matrix of that size, and patch it with a original image. Do
% this for each orginal image, and we will get the new original images with
% border of zeros.
tmp = zeros([[My,Mx]*2^(main.subdivide-1),Nt]);
tmp(1:Ny,1:Nx,:) = im; im = tmp;
clear tmp;

% Display the deformation of control meshes
if main.disp_flag
    disp_num = min(Nt,5);  % Only show the first five control meshes
    fig = figure('Position',[150,300,250*disp_num,550]);
    main.fig = fig;
    main.ha = tight_subplot(2,disp_num,[0.01,0.01],[0.01,0.05],[0.01,0.01]);
end

% Go across all the resolution levels
D = main.D;  % Dictionary
for level = 1:main.subdivide

    main.level = level;

    % Update the size of control meshes to twice bigger for the 2nd or higher
    % resolution levels
    if level > 1
        main.Gx = 2*main.Gx-3;
        main.Gy = 2*main.Gy-3;
    end

    % Current image size
    Mx = (main.Gx-3)*main.okno;
    My = (main.Gy-3)*main.okno;
    main.Mx = Mx;
    main.My = My;

    % Resize the moving images
    im_small = imresize(im,[My,Mx],'cubic');
    % if level < main.subdivide  % Gaussian filtering
    %     im_small = imgaussfilt(im_small,0.75);
    % end
    % im_small = imgaussfilt(im_small,0.75);
    main.im_small = im_small;

    % Start the alternating optimization
    % Calculate the deformed image
    [field_x,field_y] = mesh2field(main.mesh_x,main.mesh_y,main.F,main.okno);
    def_im_small = zeros(My,Mx,Nt);
    for t = 1:Nt
        def_im_small(:,:,t) = bilinear_interp(im_small(:,:,t),field_x(:,:,t),field_y(:,:,t));
    end

    % Do while the relative function difference is above the given tolerance
    % and the meximum number of iterations has not been reached
    iter = 0;
    delta_f = inf;
    while  (abs(delta_f) > optim.alt_tol) && (iter < optim.max_alt_iter(level))

        fprintf('\nlevel=%d, alt_iter=%d:\n',level,iter+1)

        % Dictionary matching
        X = reshape(def_im_small,[My*Mx,Nt]);  % Casorati matrix

        % Nonparallel
        % [coeff,ind] = max(X*D,[],2);

        % Parallel
        coeff = zeros(My*Mx,1);
        ind = zeros(My*Mx,1);
        parfor i = 1:My*Mx
            [coeff(i),ind(i)] = max(X(i,:)*D);
        end

        ref_im_small = (D(:,ind).').*coeff;
        ref_im_small = reshape(ref_im_small,[My,Mx,Nt]);
        % ref_im_small = imgaussfilt(ref_im_small,0.75);
        main.ref_im_small = ref_im_small;

        % Perform a single level 2d non-rigid group registration
        [main.mesh_x,main.mesh_y,def_im_small,new_f] = sub_group_registration_2d(main.mesh_x,main.mesh_y,main,optim);

        iter = iter+1;

        % Calculate the relative change of cost function
        if iter == 1
            f = new_f;
        else
            delta_f = (new_f-f)/abs(f);
            f = new_f;
        end

    end

    % If the resolution level is not the last one, subdivide the control
    % meshes for the next level
    if level < main.subdivide
        [main.mesh_x,main.mesh_y] = subdivide_control_mesh(main.mesh_x,main.mesh_y);
        [main.init_mesh_x,main.init_mesh_y] = subdivide_control_mesh(main.init_mesh_x,main.init_mesh_y);
    end

end

% Prepare the outputs
res.mesh_x = main.mesh_x;
res.mesh_y = main.mesh_y;
res.okno = main.okno;

% Since we have expanded the original images with zero borders, we now need
% to remove those borders to get the new images with original image size
new_im = def_im_small(1:Ny,1:Nx,:);

fprintf('Non-rigid group registration has been completed succesfully.\n')
% toc

end

