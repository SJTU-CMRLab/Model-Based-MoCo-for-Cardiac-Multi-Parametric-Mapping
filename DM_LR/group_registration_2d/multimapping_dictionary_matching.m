function [signal_map,T1_map,T2_map,B1_map] = multimapping_dictionary_matching(im,...
    D,T1_record,T2_record,B1_record)

[Ny,Nx,Nt] = size(im);

X = reshape(im,[Ny*Nx,Nt]);  % Casorati matrix

% Nonparallel
% [coeff,ind] = max(X*D,[],2);

% Parallel
coeff = zeros(Ny*Nx,1);
ind = zeros(Ny*Nx,1);
parfor i = 1:Ny*Nx
    [coeff(i),ind(i)] = max(X(i,:)*D);
end

signal_map = (D(:,ind).').*coeff;
signal_map = reshape(signal_map,[Ny,Nx,Nt]);

T1_map = T1_record(ind);
T1_map = reshape(T1_map,[Ny,Nx]);

T2_map = T2_record(ind);
T2_map = reshape(T2_map,[Ny,Nx]);

B1_map = B1_record(ind);
B1_map = reshape(B1_map,[Ny,Nx]);

end

