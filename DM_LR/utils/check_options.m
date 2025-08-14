function [main,optim] = check_options(main,optim,n)
% CHECK_OPTIONS: Check the main options and optimization options, and set
% the defaults

% Default main options
def_main.rho = 0.7;                     % Coefficient to balance the DM loss and LR loss 
def_main.okno = 6;                      % Space between adjacent control points
def_main.subdivide = 1;                 % Number of resolution levels
def_main.regularizer = 'second_order';  % Regularization type
def_main.lambda1 = 1e-4;                % Spatial regularization weight
def_main.lambda2 = 0;                   % Temporal regularization weight
def_main.disp_flag = 0;                 % Flag for the display of control meshes

% Default optimization options
def_optim.max_alt_iter =...
    [16,8,4];                           % Maximum number of alternation iterations
def_optim.max_reg_iter = 81;            % Maximum number of registration iterations
def_optim.alt_tol = 1e-6;               % Tolerance for alternation termination
def_optim.reg_tol = 1e-6;               % Tolerance for registration termination
def_optim.gamma = 10;                   % Initial step size for registration
def_optim.anneal = 0.9;                 % Decay rate for step size for registration

% Check the options and set the defaults
if ~isfield(main,'D') || isempty(main.D)
    error('Error: The dictionary must be provided!')
end

if n < 3
    optim.max_alt_iter = def_optim.max_alt_iter;
    optim.max_reg_iter = def_optim.max_reg_iter;
    optim.alt_tol = def_optim.alt_tol;
    optim.reg_tol = def_optim.reg_tol;
    optim.gamma = def_optim.gamma;
    optim.anneal = def_optim.anneal;
end

if ~isfield(main,'rho') || isempty(main.rho); main.rho = def_main.rho; end
if ~isfield(main,'okno') || isempty(main.okno); main.okno = def_main.okno; end
if ~isfield(main,'subdivide') || isempty(main.subdivide); main.subdivide = def_main.subdivide; end
if ~isfield(main,'regularizer') || isempty(main.regularizer); main.regularizer = def_main.regularizer; end
if ~isfield(main,'lambda1') || isempty(main.lambda1); main.lambda1 = def_main.lambda1; end
if ~isfield(main,'lambda2') || isempty(main.lambda2); main.lambda2 = def_main.lambda2; end
if ~isfield(main,'disp_flag') || isempty(main.disp_flag); main.disp_flag = def_main.disp_flag; end

if ~isfield(optim,'max_alt_iter') || isempty(optim.max_alt_iter); optim.max_alt_iter = def_optim.max_alt_iter; end
if ~isfield(optim,'max_reg_iter') || isempty(optim.max_reg_iter); optim.max_reg_iter = def_optim.max_reg_iter; end
if ~isfield(optim,'alt_tol') || isempty(optim.alt_tol); optim.alt_tol = def_optim.alt_tol; end
if ~isfield(optim,'reg_tol') || isempty(optim.reg_tol); optim.reg_tol = def_optim.reg_tol; end
if ~isfield(optim,'gamma') || isempty(optim.gamma); optim.gamma = def_optim.gamma; end
if ~isfield(optim,'anneal') || isempty(optim.anneal); optim.anneal = def_optim.anneal; end

main.regularizer = lower(main.regularizer);

end

