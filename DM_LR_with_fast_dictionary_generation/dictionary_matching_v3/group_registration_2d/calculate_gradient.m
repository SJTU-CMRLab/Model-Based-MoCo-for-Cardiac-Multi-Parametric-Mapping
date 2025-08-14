function  [f,fs,fs1,fs2,fr1,fr2,mesh_grad_x,mesh_grad_y,def_im_small] = calculate_gradient(mesh_x,mesh_y,main)
% CALCULATE_GRADIENT: Calculate the cost function (similarity + regularizer)
% and its gradient about the given control meshes

% Calculate the deformation fields from the given control meshes
[field_x,field_y] = mesh2field(mesh_x,mesh_y,main.F,main.okno);

% Calculate the specified similarity and its dense gradient about the given
% deformation fields
[fs,fs1,fs2,field_grad_x,field_grad_y,def_im_small] = similarity(field_x,field_y,main);

% Calculate the gradient of specified similarity about the given control
% meshes
[mesh_grad_x,mesh_grad_y] = field2mesh(field_grad_x,field_grad_y,main.F,main.okno,main.Gx,main.Gy);

% Calculate the specified regularizer and its gradient about the given
% control meshes
[fr1,fr2,mesh_grad_x1,mesh_grad_y1,mesh_grad_x2,mesh_grad_y2] = regularizer(mesh_x,mesh_y,main);

% Calculate the cost function and its gradient about the given control
% meahes
f = fs+main.lambda1*fr1+main.lambda2*fr2;
mesh_grad_x = mesh_grad_x+main.lambda1*mesh_grad_x1+main.lambda2*mesh_grad_x2;
mesh_grad_y = mesh_grad_y+main.lambda1*mesh_grad_y1+main.lambda2*mesh_grad_y2;

end



