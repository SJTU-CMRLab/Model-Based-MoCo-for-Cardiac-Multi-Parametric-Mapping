function [mesh_x,mesh_y,def_im_small,f] = sub_group_registration_2d(mesh_x,mesh_y,main,optim)
% SUB_GROUP_REGISTRATION_2D: Subfunction for the 2d group registration method
% based on free-form deformations model, low-rank similarity, diffusion or
% bending energy regularization and gradient descent optimization at a given
% resolution level

% Normalize the initial step size by dividing the standard deviation of the
% gradient of similarity about the derformation fields
[field_x,field_y] = mesh2field(mesh_x,mesh_y,main.F,main.okno);
[~,~,~,field_grad_x,field_grad_y] = similarity(field_x,field_y,main);
optim.gamma = optim.gamma/std([field_grad_x(:); field_grad_y(:)]);
clear field_x field_y field_grad_x field_grad_y

% Start the registration
[f,~,~,~,~,~,mesh_grad_x,mesh_grad_y] = calculate_gradient(mesh_x,mesh_y,main);

% Do while the relative function difference is above the given tolerance
% and the meximum number of iterations has not been reached
iter = 0;
delta_f = inf;
while (abs(delta_f) > optim.reg_tol) && (iter < optim.max_reg_iter)

    % Update the control meshes. Constrain the average deformation to be
    % identity transform by subtracting the mean from each update vector
    new_mesh_x = mesh_x-optim.gamma*(mesh_grad_x-mean(mesh_grad_x,3));
    new_mesh_y = mesh_y-optim.gamma*(mesh_grad_y-mean(mesh_grad_y,3));

    % Calculate the cost function and its gradient about the control meshes
    [new_f,fs,fs1,fs2,fr1,fr2,new_mesh_grad_x,new_mesh_grad_y,def_im_small] = ...
        calculate_gradient(new_mesh_x,new_mesh_y,main);

    % Calculate the relative change of cost function
    delta_f = (new_f-f)/abs(f);

    % Check if the step size is appropriate
    if delta_f > 0
        % If the new cost function value does not decrease, reduce the step
        % size and slightly increase the value of the cost function, which
        % is an optimization heuristic to avoid some local minima
        optim.gamma = optim.anneal*optim.gamma;
        f = f+0.01*abs(f);
    else
        iter = iter+1;

        % If the new cost function value does decrease, accept all of the
        % new results as true
        mesh_x = new_mesh_x;
        mesh_y = new_mesh_y;
        mesh_grad_x = new_mesh_grad_x;
        mesh_grad_y = new_mesh_grad_y;
        f = new_f;

        % Show the progress
        % Deforming images & control meshes
        climits = prctile(def_im_small(:),[10,90]);
        if main.disp_flag && (mod(iter,10) == 1)
            disp_num = numel(main.ha)/2;
            for i = 1:disp_num
                axes(main.ha(i)); imagesc(def_im_small(:,:,i),climits); colormap('gray'); axis equal off; title(['Deforming ',num2str(i)],'FontSize',16); drawnow
                axes(main.ha(i)); imagesc(main.ref_im_small(:,:,i),climits); colormap('gray'); axis equal off; title(['Deforming ',num2str(i)],'FontSize',16); drawnow
                axes(main.ha(i+disp_num)); mesh_plot(mesh_x(:,:,i),mesh_y(:,:,i)); drawnow
            end
        end

        if mod(iter,10) == 1
            fprintf('reg_iter=%d, f=%.2f, fs=%.2f, fs1=%.3f, fs2=%.2f, fr1=%.2f, fr2=%.2f, gamma=%.3f\n',iter,f,fs,fs1,fs2,fr1,fr2,optim.gamma)
        end

    end

end

end
