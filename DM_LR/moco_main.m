%%
addpath(genpath(pwd))

%% Load the multimapping data
clc
clear
close all

sub_ind = 1;
slice_ind = 2;

data_file_path = fullfile('./data',['subject',num2str(sub_ind)],['slice',num2str(slice_ind)],'im.mat');
load(data_file_path);

im = abs(im);
im = im./max(im(:));

% Resample to the specified in-plane resolution 1.80*1.80
[Ny,Nx,Nt] = size(im);
for t = 1:Nt
    im(:,:,t) = interp2(vox_size(1)*linspace(-Nx,Nx,Nx),vox_size(2)*linspace(-Ny,Ny,Ny).',...
        im(:,:,t),1.8*linspace(-Nx,Nx,Nx),1.8*linspace(-Ny,Ny,Ny).','cubic',0);
end

% Flip and rotate
im = flip(im,2);
im = imrotate(im,90);

tgt_size = [144,144];  % Target image size

% Crop the image sequence to the target size
fprintf('Specify the center of the ROI.\n')

f = figure('Position',[100,400,500,500]);
ha = tight_subplot(1,1,[0,0],[0.01,0.01],[0.01,0.01]);
axes(ha(1))
imagesc(im(:,:,1),prctile(im(:),[5,95]))
colormap('gray'); axis equal off
[cen_x,cen_y] = getpts(f);
cen_x = round(cen_x(1));
cen_y = round(cen_y(1));
close(f)

[Ny,Nx] = size(im,[1,2]);
min_x = cen_x-floor(tgt_size(1)/2);
max_x = min_x+tgt_size(1)-1;
min_x = max(min_x,1);
max_x = min(max_x,Nx);
min_y = cen_y-floor(tgt_size(2)/2);
max_y = min_y+tgt_size(2)-1;
min_y = max(min_y,1);
max_y = min(max_y,Ny);
im = im(min_y:max_y,min_x:max_x,:);

climits = prctile(im(:),[5,95]);  % For display

%% Show the multimapping image sequence
close all

Nt = size(im,3);

f = figure('Position',[100,400,500,500]);
set(f,'color','white')
ha = tight_subplot(1,1,[0,0.01],[0.01,0.01],[0.01,0.01]);
rep_num = 1;
for rep_ind = 1:rep_num
    for t = 1:Nt

        axes(ha(1))
        imagesc(im(:,:,t),climits)
        colormap('gray'); axis equal off

        pause(0.5)

    end
end

%% Generate the dictionary
T1_record = 200:10:2200;
T2_record = 20:5:170;
B1_record = 0.5:0.1:1;
[D,new_T1_record,new_T2_record,new_B1_record] = generate_multimapping_dictionary(...
    T1_record,T2_record,B1_record,xDoc,userOutput);

%% Perform MoCo
% Main options
main.D = D;                         % Signal dictionary
main.rho = 0.7;                     % Coefficient to balance the DM loss and LR loss
main.okno = 12;                     % Space between adjacent control points
main.subdivide = 3;                 % Number of resolution levels
main.regularizer = 'second_order';  % Regularization type
main.lambda1 = 1e-4;                % Spatial regularization weight
main.lambda2 = 0;                   % Temporal regularization weight
main.disp_flag = 0;                 % Flag for the display of control meshes

% Optimization options
optim.max_alt_iter =...
    [10,5,5];                       % Maximum numbers of alternation iterations
optim.max_reg_iter = 51;            % Maximum number of registration iterations
optim.alt_tol = 1e-5;               % Tolerance for alternation termination
optim.reg_tol = 1e-5;               % Tolerance for registration termination
optim.gamma = 1e-2;                 % Initial step size for registration
optim.anneal = 0.9;                 % Decay rate for step size for registration

% Start registration
[~,res] = group_registration_2d(im,main,optim);
new_im = group_transform(im,res);

%% Show the registration results
% close all

[Ny,Nx,Nt] = size(im);

file_path = './results/registration.gif';

f = figure('Position',[100,400,1100,400]);
f.Name = 'Registration';
set(f,'color','white')
ha = tight_subplot(1,3,[0,0.01],[0.01,0.08],[0.01,0.01]);
rep_num = 1;
for rep_ind = 1:rep_num
    for t = 1:Nt

        axes(ha(1))
        imagesc(im(:,:,t),climits)
        colormap('gray'); axis equal off
        xlim([1,Nx]); ylim([1,Ny])
        title('Original','FontSize',18)

        axes(ha(2))
        imagesc(new_im(:,:,t),climits)
        colormap('gray'); axis equal off
        xlim([1,Nx]); ylim([1,Ny])
        title('Registered','FontSize',18)

        axes(ha(3))
        imagesc(im(:,:,t),climits); hold on
        mesh_plot(res.mesh_x(:,:,t),res.mesh_y(:,:,t),'b.-'); hold off
        colormap('gray'); axis equal off
        xlim([1,Nx]); ylim([1,Ny])
        title('Control Mesh','FontSize',18)

        f.Name = ['Frame ',num2str(t)];

        pause(0.5)

        frame = getframe(f);
        frame = frame2im(frame);
        [frame,map] = rgb2ind(frame,256);
        if t == 1
            imwrite(frame,map,file_path,'gif','LoopCount',Inf,'DelayTime',0.75);
        else
            imwrite(frame,map,file_path,'gif','WriteMode','append','DelayTime',0.75);
        end

    end
end

%% Compute the parametric maps by dictionary matching
% Uncorrected maps
[~,T1_map,T2_map,B1_map] = multimapping_dictionary_matching(im,D,new_T1_record,...
    new_T2_record,new_B1_record);

% Corrected maps
[~,new_T1_map,new_T2_map,new_B1_map] = multimapping_dictionary_matching(new_im,D,new_T1_record,...
    new_T2_record,new_B1_record);

%% Show the uncorrected T1 map, T2 map and B1 map
% close all

[Ny,Nx] = size(im,[1,2]);

file_path = './results/uncorrected_maps.png';

f = figure('Position',[100,400,1300,400]);
f.Name = 'Uncorrected Maps';
set(f,'color','white')
ha = tight_subplot(1,3,[0,0.01],[0.01,0.08],[0.01,0.01]);

axes(ha(1))
imagesc(T1_map,[0,2200])
colormap(gca,'fire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('T1','FontSize',18)

axes(ha(2))
imagesc(T2_map,[0,170])
colormap(gca,'greenfire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('T2','FontSize',18)

axes(ha(3))
imagesc(B1_map,[0,1])
colormap(gca,'parula'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('B1','FontSize',18)

frame = getframe(f);
frame = frame2im(frame);
[frame,map] = rgb2ind(frame,256);
imwrite(frame,map,file_path,'png');

%% Show the corrected T1 map, T2 map and B1 map
% close all

[Ny,Nx] = size(im,[1,2]);

file_path = './results/corrected_maps.png';

f = figure('Position',[100,400,1300,400]);
f.Name = 'Corrected Maps';
set(f,'color','white')
ha = tight_subplot(1,3,[0,0.01],[0.01,0.08],[0.01,0.01]);

axes(ha(1))
imagesc(new_T1_map,[0,2200])
colormap(gca,'fire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('T1','FontSize',18)

axes(ha(2))
imagesc(new_T2_map,[0,170])
colormap(gca,'greenfire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('T2','FontSize',18)

axes(ha(3))
imagesc(new_B1_map,[0,1])
colormap(gca,'parula'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('B1','FontSize',18)

frame = getframe(f);
frame = frame2im(frame);
[frame,map] = rgb2ind(frame,256);
imwrite(frame,map,file_path,'png');

