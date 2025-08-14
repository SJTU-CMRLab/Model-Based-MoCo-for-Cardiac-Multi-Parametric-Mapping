function [im,field_x,field_y] = group_transform(im,res)
% GROUP_TRANSFORM: Apply the transformations stored in the resulting parameter
% structure (output of group_registration_2d) to the given image sequence

[Ny,Nx,Nt] = size(im);

% Generate the B-splines coefficient matrix F
F = generate_F(res.okno); % size: res.okno^2 x 4^2

% Calculate the transformation fields from the given control meshes
[field_x,field_y] = mesh2field(res.mesh_x,res.mesh_y,F,res.okno);

% Remove the borders so that the deformed images will be of the same size
% as the input ones
field_x = field_x(1:Ny,1:Nx,:);
field_y = field_y(1:Ny,1:Nx,:);

% Go across all the 2d images
for i = 1:Nt
    % Interpolate to get the deformed image
    im(:,:,i) = interp2(im(:,:,i),field_x(:,:,i),field_y(:,:,i),'cubic',0);
end

end