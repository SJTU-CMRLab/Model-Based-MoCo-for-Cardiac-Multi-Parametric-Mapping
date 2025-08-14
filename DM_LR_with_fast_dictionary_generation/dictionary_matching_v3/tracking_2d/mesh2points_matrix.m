function [T,ind] = mesh2points_matrix(points,Gx,Gy,okno)
% MASH2POINTS_MATRIX: Generate the transformation matrix from the control
% mesh to the moved points based on the given points. The given points out
% of the mesh range will not be considered when calculating the matrix, and
% their indices will be returned

Nx = (Gx-3)*okno;
Ny = (Gy-3)*okno;

points_x = points(:,1);
points_y = points(:,2);

ind  = find(points_x < 1 | points_x >= Nx+1 | points_y < 1 | points_y >= Ny+1);

points_x(ind) = [];
points_y(ind) = [];

Np = numel(points_x);
T = zeros(Np,Gx*Gy);
tmp1 = zeros(Gy,Gx);
for i = 1:Np
    
    px = points_x(i);
    py = points_y(i);
    
    ind_x = floor((px-1)/okno)+1;
    ind_y = floor((py-1)/okno)+1;
    
    B = bases((py-(ind_y-1)*okno-1)/okno).'*bases((px-(ind_x-1)*okno-1)/okno);
    
    tmp2 = tmp1;
    tmp2(ind_y:ind_y+3,ind_x:ind_x+3) = B;
    T(i,:) = tmp2(:);
    
end

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


