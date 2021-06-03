% LTV-ALS case study from Lima and Rawlings (2010)
% Case 1: Noise added to reaction rate constant k(2)

clear all
close all


% Load Data and Define Parameters
load xdata_case1.mat
load ydata_case1.mat

ns = 4; % number of simulations
tsim = 250; % duration of simulation
RT = 8.21e-2 * 400;
C = [RT, RT, RT];
x0 = [0.5, 0.05, 0]';
n = length(x0);

% Parameters for simulation
delta_d = 0.25; % simulation sampling time
t_d = 0:delta_d:tsim;
nt_d = length(t_d);

p = size(C, 1);

% Define augmented model for ALS
Gbar = [0, 0, 0, 1]';
Cbar = [C, 0];
x0bar = [x0; 0]; % augmented by integrating disturbance
nbar = length(x0bar);
g = size(Gbar, 2);

% Initialize ALS estimates
Q1_calc = 0;
Q1pd_calc = 0;
R1_calc = 0;
R1pd_calc = 0;

% Necessary Functions


for s = 1:ns

    % Start ALS Estimation

    data.delt = delta_d; % data sampling time

    datalen = nt_d; % simulation length

    % Initialization necessary for sensitivity calculation
    global model
    model.odefun = @cstr_sens;
    model.atol = sqrt(eps);
    model.rtol = sqrt(eps);
    model.odesteps = 1e6;

    % Define ALS covariance guesses
    R_guess = 7.45e-2; % (1/2 variance of measurement noise)
    Q_guess = 2.30e-5; % (the other 1/2 of the variance of measurement
    % noise scaled by C*C')

    % ALS estimation parameters
    data.N = 15; % horizon length
    start = 160; % ALS starting time

    xhat(:, 1) = x0bar;
    yhat(:, 1) = Cbar * (xhat(:, 1));
    Ck{1} = Cbar;
    Ain_p = eye(nbar);

    % Generate Ak, Bk, Gk, Ck and innovations
    for i = 2:datalen

        [xhat(:, i), Ak{i}, Bk{i}] = cvode_sens24(xhat(:, i - 1), [], 0, data.delt);

        Ck{i} = Cbar;

        yhat(:, i) = Ck{i} * (xhat(:, i));

        Gk{i} = Gbar; % only 1 w with a known G

        % Solve recursive Riccati equation to update estimate error covariance and
        % calculate EKF gain
        if i == 2
            [Lr, Mr, Pr, Er] = dlqe(Ak{2}, Gbar, Cbar, Q_guess, R_guess);
            Pk{i} = Mr;
        end

        Lk{i} = Pk{i} * Ck{i}' * inv(Ck{i}*Pk{i}*Ck{i}'+R_guess);

        Pk{i+1} = Ak{i} * (Pk{i} - Pk{i} * Ck{i}' * inv(Ck{i} * Pk{i} * Ck{i}' + R_guess) * Ck{i} * Pk{i}) * Ak{i}' + Gk{i} * Q_guess * Gk{i}';

        Ain{i} = Ak{i} - Ak{i} * Lk{i} * Ck{i};

        % Check Assumption in Remark 1 of Lima & Rawlings (2010)
        Ain_p = Ain{i} * Ain_p;
        lambda(:, i-1) = abs(eig(Ain_p));

        % Calculate innovations:
        inn(:, i) = y_d{s}(:, i) - yhat(:, i);
        xhat(:, i) = xhat(:, i) + Lk{i} * inn(:, i);
    end

    % ALS (Q, R) estimation
    % Check system properties before ALS calculation
    is_stable(Ak{datalen}, sqrt(eps), 1);
    is_observable(Ak{datalen}, Cbar, sqrt(eps));
    is_detectable(Ak{datalen}, Cbar, sqrt(eps), 1);
    is_stabilizable(Ak{datalen}, Gbar, sqrt(eps), 1);

    % Compute LHS and data stacks for ALS calculation in eq. (17) of Lima & Rawlings (2010)
    LHS_stack = [];
    dat_stack = [];

    for i = 2:datalen - start - data.N + 1
        rng = i + start + data.N - 1;

        [Qdet, Rdet, LHS, Eyy] = ltv_als(inn(:, i:rng), data.N, start+1, Ain(1, i:rng), Ak(1, i:rng), Ck(1, i:rng), Gk(1, i:rng), Lk(1, i:rng));

        LHS_stack = [LHS_stack; LHS];

        dat_stack = [dat_stack; Eyy];

        cond(LHS_stack);
        var(dat_stack);

    end

    % Check condition number of LHS
    condn(s) = cond(LHS_stack);

    Eyyn(:, s) = dat_stack;

    % Initial guess for SDP
    QR0 = blkdiag(Q_guess, R_guess);

    % Call SDP to calculate semipositive definite covariance estimates
    % The calculated (Q1, R1) are constrained unique estimates using SDP

    %[QR,phi,iter] = sdp_QR(LHS_stack,dat_stack,QR0,1,10,p,1);
    [QR, phi, ~, iter] = sdp_QR_mrQ(LHS_stack, dat_stack, eye(length(dat_stack)), QR0, 0, 1, p, 1, 1);
    Q1 = QR(1:g, 1:g);
    R1 = QR(g+1:end, g+1:end);

    Q1pd = (Q1 + Q1') / 2;
    R1pd = (R1 + R1') / 2;

    Q1_calc = Q1_calc + Q1;
    Q1pd_calc = Q1pd_calc + Q1pd;

    R1_calc = R1_calc + R1;
    R1pd_calc = R1pd_calc + R1pd;
end

% Calculate average of covariance estimates for ns simulations
Q1 = Q1_calc / ns;
Q1pd = Q1pd_calc / ns;

R1 = R1_calc / ns;
R1pd = R1pd_calc / ns;

save Covariances_case1.dat Q1 R1

% Process model (CSTR example)
function xdot = cstr_model(x, ~)
k = [0.5, 0.5, 0.2, 0.01]'; % nominal value of reaction rate constants
cf = [0.5, 0.05, 0]';
Qf = 1;
Q0 = 1;
Vr = 100;
xdot = [Qf / Vr * cf(1) - Q0 / Vr * x(1) - (k(1) * x(1) - (k(2) + x(4)) * x(2) * x(3)); ...
    Qf / Vr * cf(2) - Q0 / Vr * x(2) + (k(1) * x(1) - (k(2) + x(4)) * x(2) * x(3)) - 2 * (k(3) * x(2)^2 - k(4) * x(3)); ...
    Qf / Vr * cf(3) - Q0 / Vr * x(3) + (k(1) * x(1) - (k(2) + x(4)) * x(2) * x(3)) + (k(3) * x(2)^2 - k(4) * x(3)); ...
    0];
end

% Function used for sensitivity calculations
function xdot = cstr_sens(x, t, ~)
xdot = cstr_model(x, t);
end
