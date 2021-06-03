function [nlin, A, B] = cvode_sens24(x, u, t0, t1)

% Call the cvode integrator and calculate states, x, and
% sensitivities, sx, at the requested times.

global model

diff = cputime;
% output size
nx = size(x, 1);
nu = size(u, 1);
nt = 2;
% Number of parameters
np = nx + nu;
sx0 = eye(nx, np);

param = [x; u];
dataCV.theta = param + 1e-10;

estflag = 1:np;

% ---------------------
% CVODES initialization
% ---------------------

% Set options for integrator.

options = CVodeSetOptions('UserData', dataCV, ...
    'RelTol', model.rtol, ...
    'AbsTol', model.atol, ...
    'MaxNumSteps', model.odesteps);


% Allocate storage and set initial conditions for integrator.

%  if (~isfield(model,'cvodesfun'))
%    disp('a')
CVodeInit(@oderhs1, 'BDF', 'Newton', t0, x, options);
%CVodeMalloc (@oderhs1, t0, x, options, dataCV);
%  else
%    CVodeMalloc (model.cvodesfun, t0, x, options, dataCV);
%  end

% Set options for forward sensitivity problem.

fsa_options = CVodeSensSetOptions('method', 'Simultaneous', ...
    'ErrControl', true, ...
    'ParamField', 'theta', ...
    'ParamList', estflag);

if (~isfield(model, 'dodedu') || ~isfield(model, 'dodedx'))
    CVodeSensInit(np, [], sx0, fsa_options);
else
    CVodeSensInit(np, @sensrhs, sx0, fsa_options);
end

% Allocate storage for forward sensitivity problem.
%  CVodeSensMalloc (np, 'Simultaneous', sx0, fsa_options);

[status, t, x_step, sx_step] = CVode(t1, 'Normal');

if (status == 0)
    x(:, 2) = x_step;
else (status < 0)
    error('CVode failed with status = %d', status);
end

nlin = x(:, 2);
Sk = [];
Sk = sx_step;
A = Sk(:, 1:nx);
B = Sk(:, nx+1:nx+nu);
CVodeFree();

end %cvode_sens

function [xsd, flag, new_data] = sensrhs(t, x, xdot, xs, dataCV)
% Right hand side of forward sensitivity equation
% sx = dx/dp
% dsx / dt = df/dx * sx + df/dp
global model
nx = length(x);
p = dataCV.theta;
u = p(nx+1:end);

dfdp = [model.dodedx(x, t, u), model.dodedu(x, t, u)];;
xsd = model.dodedx(x, t, u) * xs + dfdp;
flag = 0;
new_data = [];
end

function [xdot, flag, new_data] = oderhs1(t, x, dataCV)
global model;
u = dataCV.theta(length(x)+1:end);
xdot = model.odefun(x, t, u);
flag = 0;
new_data = [];
end
