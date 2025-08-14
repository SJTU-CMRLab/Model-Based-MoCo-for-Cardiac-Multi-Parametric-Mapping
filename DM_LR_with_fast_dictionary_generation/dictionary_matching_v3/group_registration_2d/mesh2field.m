function [field_x,field_y] = mesh2field(mesh_x,mesh_y,F,okno)
% MESH2FIELD: Calculate the deformation fields (new positions of all image 
% pixels) from the given control meshes

[Gy,Gx,Nt] = size(mesh_x);

% Go across all the 4 x 4 local control meshes
field_x = zeros((Gy-3)*okno,(Gx-3)*okno,Nt);
field_y = field_x;
for i = 1:Gy-3
    for j = 1:Gx-3
        
        % Indices of the pixels corresponding to the given local meshes
        ind_x = (j-1)*okno+1:j*okno;
        ind_y = (i-1)*okno+1:i*okno;
        
        % Take the x-coordinates of the control points in each 4 x 4 local
        % control mesh. Rearrange them into vectors and multiply by the
        % matrix of B-splines coefficients to get the new x-coordinates
        % of pixels within the corresponding okno x okno image patch
        tmp = mesh_x(i:i+3,j:j+3,:);
        tmp = reshape(tmp,[16,Nt]);
        field_x(ind_y,ind_x,:) = reshape(F*tmp,[okno,okno,Nt]);
        
        % Repeat for the y-coordinates
        tmp = mesh_y(i:i+3,j:j+3,:);
        tmp = reshape(tmp,[16,Nt]);
        field_y(ind_y,ind_x,:) = reshape(F*tmp,[okno,okno,Nt]);
        
    end
end

end

