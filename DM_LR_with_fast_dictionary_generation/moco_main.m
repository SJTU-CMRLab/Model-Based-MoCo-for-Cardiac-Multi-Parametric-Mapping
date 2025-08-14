%%
addpath(genpath(pwd))

%% Load the multimapping data
clc
clear
close all

sub_ind = 1;
slice_ind = 2;

dir_path = sprintf('./data/in_vivo/subject%d/slice%d',sub_ind,slice_ind);

data_path = fullfile(dir_path,'data.mat');
data = load(data_path);
im = data.im;

im = abs(im);
im = im./max(im(:));

tgt_size = [144,144];  % Target image size

climits = prctile(im(:),[1,99]);  % For display

% Crop the image sequence to the target size
fprintf('Specify the center of the ROI.\n')

f = figure('Position',[400,400,500,500]);
ha = tight_subplot(1,1,[0,0],[0.01,0.01],[0.01,0.01]);
axes(ha(1))
imagesc(im(:,:,1),climits)
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

% Show the multimapping image sequence
close all

Nt = size(im,3);

f = figure('Position',[400,400,500,500]);
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

%% Perform MoCo
% Load the dictionary for MoCo
dict_path = fullfile(dir_path,'dict.mat');
dict = load(dict_path);
D = double(dict.D);

% Main options
main.D = D;                         % Signal dictionary
main.rho = 0.7;                     % Coefficient to balance the DM loss and LR loss
main.okno = 10;                     % Space between adjacent control points
main.subdivide = 3;                 % Number of resolution levels
main.regularizer = 'second_order';  % Regularization type
main.lambda1 = 1e-4;                % Spatial regularization weight
main.lambda2 = 0;                   % Temporal regularization weight
main.disp_flag = 0;                 % Flag for the display of control meshes

% Optimization options
optim.max_alt_iter =...
    [12,6,6];                       % Maximum numbers of alternation iterations
optim.max_reg_iter = 51;            % Maximum number of registration iterations
optim.alt_tol = 1e-6;               % Tolerance for alternation termination
optim.reg_tol = 1e-6;               % Tolerance for registration termination
optim.gamma = 1e-2;                 % Initial step size for registration
optim.anneal = 0.9;                 % Decay rate for step size for registration

% Start registration
[~,res] = group_registration_2d(im,main,optim);
[new_im,trans_field_x,trans_field_y] = group_transform(im,res);

%% Show the registration results
close all

[Ny,Nx,Nt] = size(im);

file_path = './results/registration.gif';

f = figure('Position',[100,400,1100,400]);
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

%% Compute the parametric maps
new_T1_record = double(dict.new_T1_record);
new_T2_record = double(dict.new_T2_record);

% Uncorrected maps
[~,T1_map,T2_map] = multimapping_dictionary_matching(im,D,new_T1_record,new_T2_record);

% Corrected maps
[~,new_T1_map,new_T2_map] = multimapping_dictionary_matching(new_im,D,new_T1_record,new_T2_record);

%% Show the T1 and T2 maps
close all

[Ny,Nx] = size(im,[1,2]);

file_path = './results/maps.png';

f = figure('Position',[100,400,1300,300]);
set(f,'color','white')
ha = tight_subplot(1,4,[0,0.01],[0.01,0.08],[0.01,0.01]);

axes(ha(1))
imagesc(T1_map,[0,2200])
colormap(gca,'fire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('Uncorrected T1','FontSize',18)

axes(ha(2))
imagesc(new_T1_map,[0,2200])
colormap(gca,'fire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('Corrected T1','FontSize',18)

axes(ha(3))
imagesc(T2_map,[0,170])
colormap(gca,'greenfire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('Uncorrected T2','FontSize',18)

axes(ha(4))
imagesc(new_T2_map,[0,170])
colormap(gca,'greenfire'); colorbar
axis equal off
xlim([1,Nx]); ylim([1,Ny])
title('Corrected T2','FontSize',18)

frame = getframe(f);
frame = frame2im(frame);
[frame,map] = rgb2ind(frame,256);
imwrite(frame,map,file_path,'png');

