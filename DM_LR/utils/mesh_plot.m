function mesh_plot(mesh_x,mesh_y,line_spec)
% MESH_PLOT: Plot a 2d control mesh. Use line_spec to specify the color,
% line type and marker symbol to use in the plot (e.g. line_spec = 'r-*')

[Gy,Gx] = size(mesh_x);

XY = [mesh_x(:),mesh_y(:)];
A = spdiags(ones(Gy*Gx,2),[1,Gy],Gy*Gx,Gy*Gx);
a = Gy:Gy:Gy*Gx-1;
A(sub2ind([Gy*Gx,Gy*Gx],a,a+1)) = 0;
if ~exist('line_spec','var') || isempty(line_spec)
    gplot(A,XY)
else
    gplot(A,XY,line_spec)
end
axis ij equal off;

end

