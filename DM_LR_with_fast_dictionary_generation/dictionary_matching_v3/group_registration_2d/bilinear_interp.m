function [im_interp,deriv_x,deriv_y] = bilinear_interp(im,field_x,field_y)
% BILINEAR_INTERP: Interpolate the input 2d image linearly at every query
% points in given deformation field (field_x,field_y), and calculate the
% x,y-derivatives of continuous interpolated image at every query points

[Ny,Nx] = size(im);
N = Nx*Ny;

Ax1 = floor(field_x); Ax2 = Ax1+1;
Ay1 = floor(field_y); Ay2 = Ay1+1;

Bx1 = Ax2-field_x; Bx2 = field_x-Ax1;
By1 = Ay2-field_y; By2 = field_y-Ay1;

ind = (Ax1-1)*Ny+Ay1;
ind((ind < 1) | (ind > N)) = 1;
V11 = im(ind);

ind = (Ax1-1)*Ny+Ay2;
ind((ind < 1) | (ind > N)) = 1;
V12 = im(ind);

ind = (Ax2-1)*Ny+Ay1;
ind((ind < 1) | (ind > N)) = 1;
V21 = im(ind);

ind = (Ax2-1)*Ny+Ay2;
ind((ind < 1) | (ind > N)) = 1;
V22 = im(ind);

im_interp = V11.*Bx1.*By1+V12.*Bx1.*By2+V21.*Bx2.*By1+V22.*Bx2.*By2;
deriv_x = -V11.*By1-V12.*By2+V21.*By1+V22.*By2;
deriv_y = -V11.*Bx1+V12.*Bx1-V21.*Bx2+V22.*Bx2;

ind = (Ax1 < 1) | (Ax2 > Nx) | (Ay1 < 1) | (Ay2 > Ny);
im_interp(ind) = 0;
deriv_x(ind) = 0;
deriv_y(ind) = 0;

end

