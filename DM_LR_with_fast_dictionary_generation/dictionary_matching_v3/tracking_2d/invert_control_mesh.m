function [mesh_x2,mesh_y2,okno2] = invert_control_mesh(mesh_x1,mesh_y1,okno1)
% INVERT_CONTROL_MESH: Invert the given control mesh. A smaller okno can
% be used for the inverse mesh to yield a more accurate approximation, but
% the increasing computational cost should be considered

[Gy1,Gx1,Nt] = size(mesh_x1);

Nx = (Gx1-3)*okno1;
Ny = (Gy1-3)*okno1;

[X,Y] = meshgrid(1:2:Nx,1:2:Ny);
points1 = [X(:),Y(:)];

okno2 = ceil(okno1);  % Okno of the inverse mesh
Gx2 = ceil(Nx/okno2)+3;  % Size of the inverse mesh
Gy2 = ceil(Ny/okno2)+3;
mesh_x2 = zeros(Gy2,Gx2,Nt);
mesh_y2 = zeros(Gy2,Gx2,Nt);
for i = 1:Nt
    
    points2 = move_points(points1,mesh_x1(:,:,i),mesh_y1(:,:,i),okno1);
    [T,ind] = mesh2points_matrix(points2,Gx2,Gy2,okno2);
    tmp = points1;
    tmp(ind,:) = [];
    
    pinv_T = pinv(T);
    mesh_x2(:,:,i) = reshape(pinv_T*tmp(:,1),[Gy2,Gx2]);
    mesh_y2(:,:,i) = reshape(pinv_T*tmp(:,2),[Gy2,Gx2]);
    
%     disp(mean(abs(T*pinv_T*tmp(:,1)-tmp(:,1))))
%     disp(mean(abs(T*pinv_T*tmp(:,2)-tmp(:,2))))
    
end

% Replace the unstable control points at borders by removing and padding
x0 = (1:Gx2)*okno2-okno2;
y0 = (1:Gy2)*okno2-okno2;
[mesh_x0,mesh_y0] = meshgrid(x0,y0);
disp_x2 = mesh_x2-mesh_x0;
disp_y2 = mesh_y2-mesh_y0;
disp_x2 = disp_x2(3:end-2,3:end-2,:);
disp_y2 = disp_y2(3:end-2,3:end-2,:);
disp_x2 = padarray(disp_x2,[2,2,0],'both','symmetric');
disp_y2 = padarray(disp_y2,[2,2,0],'both','symmetric');
mesh_x2 = mesh_x0+disp_x2;
mesh_y2 = mesh_y0+disp_y2;

end

