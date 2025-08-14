function [fr1,fr2,mesh_grad_x1,mesh_grad_y1,mesh_grad_x2,mesh_grad_y2] = regularizer(mesh_x,mesh_y,main)
% REGULARIZER: Calculate the specified regularizer and its gradient about
% the given control meshes. The regularization is applied to both spatial
% and temporal dimensions

% Displacements of control meshes
mesh_x = mesh_x-main.init_mesh_x;
mesh_y = mesh_y-main.init_mesh_y;

N = numel(mesh_x)/1000;

switch main.regularizer
    
    case 'first_order'  % Diffusion
        % Spatial regularization
        term1_x = diff1(mesh_x,2);
        term1_y = diff1(mesh_y,2);
        
        term2_x = diff1(mesh_x,1);
        term2_y = diff1(mesh_y,1);
        
        fr1 = (sum(term1_x(:).^2)+sum(term1_y(:).^2)+sum(term2_x(:).^2)+...
            sum(term2_y(:).^2))/2/N;
        
        mesh_grad_x1 = (trans_diff1(term1_x,2)+trans_diff1(term2_x,1))/N;
        mesh_grad_y1 = (trans_diff1(term1_y,2)+trans_diff1(term2_y,1))/N;
        
        % Temporal regularization
        term_x = circular_diff1(mesh_x,3);
        term_y = circular_diff1(mesh_y,3);
%         term_x = diff1(mesh_x,3);
%         term_y = diff1(mesh_y,3);
        
        fr2 = (sum(term_x(:).^2)+sum(term_y(:).^2))/2/N;
        
        mesh_grad_x2 = trans_circular_diff1(term_x,3)/N;
        mesh_grad_y2 = trans_circular_diff1(term_y,3)/N;
%         mesh_grad_x2 = trans_diff1(term_x,3)/N;
%         mesh_grad_y2 = trans_diff1(term_y,3)/N;
        
    case 'second_order'  % Bending energy of a thin-plate of metal
        % Spatial regularization
        term1_x = diff2(mesh_x,2);
        term1_y = diff2(mesh_y,2);
        
        term2_x = diff1(diff1(mesh_x,2),1);
        term2_y = diff1(diff1(mesh_y,2),1);
        
        term3_x = diff2(mesh_x,1);
        term3_y = diff2(mesh_y,1);
        
        fr1 = (sum(term1_x(:).^2)+sum(term1_y(:).^2)+2*sum(term2_x(:).^2)+...
            2*sum(term2_y(:).^2)+sum(term3_x(:).^2)+sum(term3_y(:).^2))/2/N;
        
        mesh_grad_x1 = (trans_diff2(term1_x,2)+2*trans_diff1(trans_diff1(term2_x,2),1)+...
            trans_diff2(term3_x,1))/N;
        mesh_grad_y1 = (trans_diff2(term1_y,2)+2*trans_diff1(trans_diff1(term2_y,2),1)+...
            trans_diff2(term3_y,1))/N;
        
        % Temporal regularization
        term_x = circular_diff2(mesh_x,3);
        term_y = circular_diff2(mesh_y,3);
%         term_x = diff2(mesh_x,3);
%         term_y = diff2(mesh_y,3);
        
        fr2 = (sum(term_x(:).^2)+sum(term_y(:).^2))/2/N;
        
        mesh_grad_x2 = trans_circular_diff2(term_x,3)/N;
        mesh_grad_y2 = trans_circular_diff2(term_y,3)/N;
%         mesh_grad_x2 = trans_diff2(term_x,3)/N;
%         mesh_grad_y2 = trans_diff2(term_y,3)/N;
        
    otherwise
        error('Error: The regularizer is unsupported! Choose from ''first_order'' and ''second_order''.')
        
end

end

%% Helper functions
function Y = diff1(X,dim)
% DIFF1: Calculate the first order difference of input 3d tensor along the
% specified dimension
% X: Ny x Nx x Nt --> Y: (Ny-1) x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x (Nx-1) x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x (Nt-1), if dim = 3

if dim == 1
    Y = X(2:end,:,:)-X(1:end-1,:,:);
elseif dim == 2
    Y = X(:,2:end,:)-X(:,1:end-1,:);
elseif dim == 3
    Y = X(:,:,2:end)-X(:,:,1:end-1);
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = diff2(X,dim)
% DIFF2: Calculate the second order difference of input 3d tensor along the
% specified dimension
% X: Ny x Nx x Nt --> Y: (Ny-2) x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x (Nx-2) x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x (Nt-2), if dim = 3

if dim == 1
    Y = diff1(diff1(X,1),1);
elseif dim == 2
    Y = diff1(diff1(X,2),2);
elseif dim == 3
    Y = diff1(diff1(X,3),3);
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = trans_diff1(X,dim)
% TRANS_DIFF1: Calculate the first order transposed difference of input 3d
% tensor along the specified dimension
% X: Ny x Nx x Nt --> Y: (Ny+1) x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x (Nx+1) x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x (Nt+1), if dim = 3

if dim == 1
    Y = padarray(X,[1,0,0],'pre')-padarray(X,[1,0,0],'post');
elseif dim == 2
    Y = padarray(X,[0,1,0],'pre')-padarray(X,[0,1,0],'post');
elseif dim == 3
    Y = padarray(X,[0,0,1],'pre')-padarray(X,[0,0,1],'post');
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = trans_diff2(X,dim)
% TRANS_DIFF2: Calculate the second order transposed difference of input 3d
% tensor along the specified dimension
% X: Ny x Nx x Nt --> Y: (Ny+2) x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x (Nx+2) x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x (Nt+2), if dim = 3

if dim == 1
    Y = trans_diff1(trans_diff1(X,1),1);
elseif dim == 2
    Y = trans_diff1(trans_diff1(X,2),2);
elseif dim == 3
    Y = trans_diff1(trans_diff1(X,3),3);
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = circular_diff1(X,dim)
% CIRCULAR_DIFF1: Calculate the first order circular difference of input 3d
% tensor along the specified dimension
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 3

if dim == 1
    Y = X([2:end,1],:,:)-X;
elseif dim == 2
    Y = X(:,[2:end,1],:)-X;
elseif dim == 3
    Y = X(:,:,[2:end,1])-X;
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = circular_diff2(X,dim)
% CIRCULAR_DIFF2: Calculate the second order circular difference of input
% 3d tensor along the specified dimension
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 3

if dim == 1
    Y = circular_diff1(circular_diff1(X,1),1);
elseif dim == 2
    Y = circular_diff1(circular_diff1(X,2),2);
elseif dim == 3
    Y = circular_diff1(circular_diff1(X,3),3);
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = trans_circular_diff1(X,dim)
% TRANS_CIRCULAR_DIFF1: Calculate the first order transposed circular
% difference of input 3d tensor along the specified dimension
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 3

if dim == 1
    Y = X([end,1:end-1],:,:)-X;
elseif dim == 2
    Y = X(:,[end,1:end-1],:)-X;
elseif dim == 3
    Y = X(:,:,[end,1:end-1])-X;
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end

function Y = trans_circular_diff2(X,dim)
% TRANS_CIRCULAR_DIFF2: Calculate the second order transposed circular
% difference of input 3d tensor along the specified dimension
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 1
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 2
% X: Ny x Nx x Nt --> Y: Ny x Nx x Nt, if dim = 3

if dim == 1
    Y = trans_circular_diff1(trans_circular_diff1(X,1),1);
elseif dim == 2
    Y = trans_circular_diff1(trans_circular_diff1(X,2),2);
elseif dim == 3
    Y = trans_circular_diff1(trans_circular_diff1(X,3),3);
else
    error('Error: The difference of 3d tensor should be performed along the first, second or third dimension!')
end

end



