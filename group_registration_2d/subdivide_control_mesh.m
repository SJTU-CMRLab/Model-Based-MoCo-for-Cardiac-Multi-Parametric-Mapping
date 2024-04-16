function [new_mesh_x,new_mesh_y] = subdivide_control_mesh(mesh_x,mesh_y)
% SUBDIVIDE_CONTROL_MESH: Subdivide the current control meshes to one level
% up (twise denser). The new mesh size will be (2*Gy-3) x (2*Gx-3) x Nt.
% The main idea is to take each 2 x 2 local control mesh and upsample it to
% form a new 3 x 3 local control mesh. After upsampling, we need to remove
% the control points around the border, because they are useless at higher
% resolution level, and will only increase the computational time

[Gy,Gx,Nt] = size(mesh_x);

new_mesh_x = zeros(Gy,2*Gx-2,Nt);
new_mesh_y = zeros(Gy,2*Gx-2,Nt);

mesh_x_fill = (mesh_x(:,1:end-1,:)+mesh_x(:,2:end,:))/2;
mesh_y_fill = (mesh_y(:,1:end-1,:)+mesh_y(:,2:end,:))/2;
for i = 1:Gx-1
    new_mesh_x(:,2*i-1:2*i,:) = cat(2,mesh_x(:,i,:),mesh_x_fill(:,i,:));
    new_mesh_y(:,2*i-1:2*i,:) = cat(2,mesh_y(:,i,:),mesh_y_fill(:,i,:));
end
new_mesh_x = new_mesh_x(:,2:end,:);
new_mesh_y = new_mesh_y(:,2:end,:);

mesh_x = new_mesh_x;
mesh_y = new_mesh_y;

new_mesh_x = zeros(2*Gy-2,2*Gx-3,Nt);
new_mesh_y = zeros(2*Gy-2,2*Gx-3,Nt);

mesh_x_fill = (mesh_x(1:end-1,:,:)+mesh_x(2:end,:,:))/2;
mesh_y_fill = (mesh_y(1:end-1,:,:)+mesh_y(2:end,:,:))/2;
for i = 1:Gy-1
    new_mesh_y(2*i-1:2*i,:,:) = cat(1,mesh_y(i,:,:),mesh_y_fill(i,:,:));
    new_mesh_x(2*i-1:2*i,:,:) = cat(1,mesh_x(i,:,:),mesh_x_fill(i,:,:));
end
new_mesh_x = new_mesh_x(2:end,:,:);
new_mesh_y = new_mesh_y(2:end,:,:);

new_mesh_x = 2*new_mesh_x-1;
new_mesh_y = 2*new_mesh_y-1;

end



