function [points,ind] = move_points(points,mesh_x,mesh_y,okno)
% MOVE_POINTS: Move the given points according to the control mesh. The
% mesh has no time dimension here. The points out of the mesh range will
% not be moved, and their indices will be returned

[Gy,Gx] = size(mesh_x);

Nx = (Gx-3)*okno;
Ny = (Gy-3)*okno;

points_x = points(:,1);
points_y = points(:,2);

ind  = find(points_x < 1 | points_x >= Nx+1 | points_y < 1 | points_y >= Ny+1);

for i = 1:numel(points_x)
    
    if ~ismember(i,ind)
        px = points_x(i);
        py = points_y(i);
        
        ind_x = floor((px-1)/okno)+1;
        ind_y = floor((py-1)/okno)+1;
        
        Cx = mesh_x(ind_y:ind_y+3,ind_x:ind_x+3);
        Cy = mesh_y(ind_y:ind_y+3,ind_x:ind_x+3);
        
        B = bases((py-(ind_y-1)*okno-1)/okno).'*bases((px-(ind_x-1)*okno-1)/okno);
        
        points_x(i) = sum(Cx(:).*B(:));
        points_y(i) = sum(Cy(:).*B(:));
    end
    
end
points = [points_x,points_y];

end

%% Helper function
function b = bases(u)
% BASES: Evaluate the cubic B-splines basis functions at the given points.
% The points should be between 0 and 1 here

u = u(:);

b = zeros(numel(u),4);
b(:,1) = (1-u).^3/6;
b(:,2) = (3*u.^3-6*u.^2+4)/6;
b(:,3) = (-3*u.^3+3*u.^2+3*u+1)/6;
b(:,4) = u.^3/6;

end


