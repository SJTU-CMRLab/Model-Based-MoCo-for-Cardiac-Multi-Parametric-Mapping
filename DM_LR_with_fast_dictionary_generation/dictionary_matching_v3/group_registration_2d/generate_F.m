function F = generate_F(okno)
% GENERATE_F: Generate a matrix of cubic B-splines coefficients for a 2d
% patch. F times the coordinate vectors of given control points will give
% the new spatial locations of image pixels within corresponding patch
%
% 2D patch: okno x okno, Local control mesh: 4 x 4, Matrix F: okno^2 x 4^2

B = [-1  3 -3  1;
      3 -6  3  0;
     -3  0  3  0;
      1  4  1  0]/6;

u = linspace(0,1,okno+1).';
u = u(1:end-1);
T = [u.^3,u.^2,u,ones(okno,1)];

B = T*B;
F = kron(B,B);

end