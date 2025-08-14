function [fs,fs1,fs2,field_grad_x,field_grad_y,def_im_small] = similarity(field_x,field_y,main)
% SIMILARITY: Calculate the specified similarity and its dense gradient
% about the given deformation fields

[My,Mx,Nt] = size(main.im_small);

% Interpolate each 2d image and calculate the x,y-derivatives of continuous
% interpolated image at each query point in deformation field
def_im_small = zeros(My,Mx,Nt);
deriv_x = def_im_small;
deriv_y = def_im_small;
for i = 1:Nt
    [def_im_small(:,:,i),deriv_x(:,:,i),deriv_y(:,:,i)] = ...
        bilinear_interp(main.im_small(:,:,i),field_x(:,:,i),field_y(:,:,i));
end

% Calculate the similarity and its dense gradient
% N = My*Mx*Nt/1000;
% field_grad = def_im_small-main.ref_im_small;
% fs = sum(field_grad(:).^2)/2/N;
% field_grad = field_grad/N;
% field_grad_x = field_grad.*deriv_x;
% field_grad_y = field_grad.*deriv_y;

% Dictionary-matching loss
N1 = My*Mx*Nt/1000;
field_grad1 = def_im_small-main.ref_im_small;
fs1 = sum(field_grad1(:).^2)/2/N1;
field_grad1 = field_grad1/N1;

% Low-rankness loss
N2 = Nt*2^(main.level-main.subdivide);
L = 2;  % Expected rank
X = reshape(def_im_small,[My*Mx,Nt]);  % Casorati matrix
[U,S,V] = svd(X,'econ');
s = diag(S);
fs2 = sum(s(L+1:end))/N2;
U(:,1:L) = 0;
field_grad2 = U*V.';
field_grad2 = reshape(field_grad2,[My,Mx,Nt])/N2;

mu = 20;
rho = main.rho;
% rho = 0.7;
fs = (1-rho)*mu*fs1+rho*fs2;
field_grad = (1-rho)*mu*field_grad1+rho*field_grad2;

field_grad_x = field_grad.*deriv_x;
field_grad_y = field_grad.*deriv_y;

end

