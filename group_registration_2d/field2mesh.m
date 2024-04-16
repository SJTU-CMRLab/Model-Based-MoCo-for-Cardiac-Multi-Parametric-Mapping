function [mesh_grad_x,mesh_grad_y] = field2mesh(field_grad_x,field_grad_y,F,okno,Gx,Gy)
% FIELD2MESH: Transform the dense gradient about deformation fields to the
% gradient about control meshes

Nt = size(field_grad_x,3);
mesh_grad_x = zeros(Gy,Gx,Nt);
mesh_grad_y = mesh_grad_x;

% Go across all the 4 x 4 local control meshes
for i = 1:Gy-3
    for j = 1:Gx-3
        
        % Indices of the pixels corresponding to the given local meshes
        ind_x = (j-1)*okno+1:j*okno;
        ind_y = (i-1)*okno+1:i*okno;
        
        % Extract the gradient about the deformation fields in x-direction
        % that correspond to the given local meshes
        tmp = field_grad_x(ind_y,ind_x,:);
        
        % Calculate the gradient about control meshes in x-direction. 
        % Accumulation is because most B-splines control points are shared 
        % by multiple image patches
        tmp = reshape(tmp,[okno^2,Nt]);
        tmp = reshape(F'*tmp,[4,4,Nt]);
        mesh_grad_x(i:i+3,j:j+3,:) = mesh_grad_x(i:i+3,j:j+3,:)+tmp;
        
        % Do the same thing in y-direction
        tmp = field_grad_y(ind_y,ind_x,:);
        tmp = reshape(tmp,[okno^2,Nt]);
        tmp = reshape(F'*tmp,[4,4,Nt]);
        mesh_grad_y(i:i+3,j:j+3,:) = mesh_grad_y(i:i+3,j:j+3,:)+tmp;
        
    end
end

end


