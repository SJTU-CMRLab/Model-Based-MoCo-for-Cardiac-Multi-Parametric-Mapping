function field_plot(field_x,field_y,s,r,line_spec)
% FIELD_PLOT: Plot a 2d displacement field. Use line_spec to specify the
% color, line type and marker symbol to use in the plot (e.g. line_spec =
% 'r-*')
% s: scaling factor, r: downsampling rate

[My,Mx] = size(field_x);

[X,Y] = meshgrid(1:Mx,1:My);
X = X(1:r:end,1:r:end);
Y = Y(1:r:end,1:r:end);
disp_x = s*(field_x(1:r:end,1:r:end)-X);
disp_y = s*(field_y(1:r:end,1:r:end)-Y);
if ~exist('line_spec','var') || isempty(line_spec)
    quiver(X,Y,disp_x,disp_y,'off')
else
    quiver(X,Y,disp_x,disp_y,'off',line_spec)
end
axis ij equal off;

end



