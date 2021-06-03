function [Qest_cell, Rest_cell, trQ, Phi, Phi_tot, phi0, bhat_mtr, bhat_mtr_thry, cov_bound, Iter, Iter_maxed, timespent] = als_sdp_mrQ(data, N, model, estimator, varargin)

% A modified ALS-SDP function.  Uses Single Column ALS form and imposes
% semidefinite constraints.
% Options allow semidefinite constraints to be removed.

% [Qest_cell,Rest_cell,trQ,Phi,Phi_tot,phi0,bhat_mtr,bhat_mtr_thry,cov_bound, ...
% Iter,Iter_maxed,timespent] = als_sdp_mrQ(data,N,model,estimator,varargin)


% Function Inputs:

% data.yk (measurements),
% data.uk (inputs),
% data.xhatk (state estimates- optional),
% data.datapts (number of data points considered, optional),
% data.start (data to be ignored in the beginning until initial condition is
% negligible - optional, default is 100),
% N (window size),
% model.A,
% model.B (optional),
% model.C,
% model.G (optional, default is identity matrix)
% estimator.Q (Qw initial guess, optional)
% estimator.R (Rv initial guess, optional),
% estimator.L (initial estimator gain - optional)

% Either (Q,R) or L must be supplied.  If L is supplied, (Q,R) are used as
% initial guess but not in calculating the estimator gain.  If L is not supplied
% then (Q,R) are used to calculate initial estimator

% Optional inputs for vargin:
% 'rho_values', [rho1 rho2 rho3 ...] (vector of penalties on trace(Q), default is
% [0])
% 'Rform', 'sym' or 'Rform', 'diag' (imposes symmetric or diagonal constraints on
% Rv, default is diagonal)
% 'Weight', 'I', or "Weight', 'data' (imposes identity matrix weighting or
% data-based weighting on least-squares term, default is data-based)
% 'Plot', 0 (turns off all plots) or 'plot', 1 (produces plots of L-innovations
% and plots of autocovariance estimates and fits), default is to produce plots
% 'sdp', 0 (solves unconstrained problem), default is 'sdp', 1 (enforces
% semidefinite constraints)
% 'tracestates', [1 1 ... 0 0] (vector indicating which states are to be
% penalized in the trace(Q) term, default is [1 1 ... 1], to penalize all states)

%
% Function outputs :

% Qest_cell - cell containing estimated Qw for each penalty on tr(Q)
% Rest_cell - cell containing estimated Rv for each penalty on tr(Q)
% trQ - vector containing trace(Q) for each penalty on tr(Q)
% Phi - vector containing least-squares portion of objective for each
% penalty on tr(Q); scaled so that phi = 1 when rho = 0
% Phi_tot - vector containing full objective function value for each
% penalty on tr(Q)
% phi0 - vector containing objective function value at rho = 0 (prior
% to scaling)
% bhat_mtr - matrix of autocovariances estimated from data, columns are
% <y_1(k),y_1(k-j)>, ..., <y_1(k),y_p(k-j)>, ..., <y_p(k),y_1(k-j)>, ...,
% <y_p(k),y_p(k-j)>
% bhat_mtr - matrix of theoretical autocovariances calculated from model
% and ALS results (for smallest penalty on trace(Q))
% cov_bound - estimated standard deviations for each auto-cross covariances,
% [sigma_11; sigma_12; ... sigma_pp], use +/- 2*sigma as 95% confidence intervals
% Iter - vector containing number of iterations for each penalty on tr(Q)
% Iter_maxed(i) = 1 if maximum number of iterations are reached for that
% penalty on trQ
% Timespent - vector of time spent in seconds of CPU time for each
% penalty on tr(Q)

if ~isfield(data, 'datapts')
    data.datapts = size(data.yk, 2);
end
datapts = data.datapts;

Aa = model.A;
Ca = model.C;

if isfield(model, 'G')
    Ga = model.G;
else
    Ga = eye(size(Aa));
end
[n, g] = size(Ga);
p = size(Ca, 1);

if (isfield(model, 'B'))
    Ba = model.B;
    m = size(Ba, 2);
else
    Ba = zeros(n);
    m = n;
    data.uk = zeros(m, datapts);
end

if (isfield(data, 'start'))
    start = data.start;
else
    start = 100;
end


na = n;
ga = g;
pa = p;

% Deal with variable arguments

% Default values:
lam_vec = 0; % don't penalize trace
trstates = ones(ga, 1); % include all states when calculating trace
Rsym = 0; % use diagonal R
dataweight = 1; % data-based weighting
plot_flag = 1; % produce plots
sdp = 1; % enforce semidefinite constraints

% Inputs
okargs = {'rho_values', 'tracestates', 'rform', 'weight', 'plot', 'sdp'};
for j = 1:2:(nargin - 4)
    pname = varargin{j};
    pval = varargin{j+1};
    k = strcmpi(pname, okargs);
    if isempty(k)
        error('als_sdp_mrQ: Unknown parameter name: %s.', pname);
    elseif length(k) > 1
        error('als_sdp_mrQ: Ambiguous parameter name: %s.', pname);
    else
        switch (k)
            case 1 % Penalty on trQ
                lam_vec = pval(:);
            case 2 % States to be penalized
                trstates = pval(:);
            case 3 % Rv symmetric or diagonal
                if strcmpi(pval, 'sym')
                    Rsym = 1;
                elseif strcmpi(pval, 'diag')
                    Rsym = 0;
                else
                    warning('als_sdp_mrQ: Unknown structure type for R; defaulting to diagonal');
                end
            case 4 % Weighting - data-based or identity?
                if strcmpi(pval, 'I')
                    dataweight = 0;
                elseif strcmpi(pval, 'data')
                    dataweight = 1;
                else
                    warning('als_sdp_mrQ: Unknown weighting type; defaulting to data-based');
                end
            case 5 % Plot
                plot_flag = pval;
            case 6 % Semidefinite constraints
                sdp = pval;
        end
    end
end


% Estimator Simulation
y = data.yk;
u = data.uk;

% Calculate estimator gain if not supplied:
if isfield(estimator, 'L')
    L = estimator.L;
else
    L = dlqe(Aa, Ga, Ca, estimator.Q, estimator.R);
end

% Estimate states
if (~isfield (data, 'xhatk'))
    xhat = zeros(na, datapts);
    xhat_ = zeros(n, datapts);
    xhat_(1:n, 1) = model.xhat0;
    for i = 1:datapts
        xhat(:, i) = xhat_(:, i) + L * (y(:, i) - Ca * xhat_(:, i));
        xhat_(:, i+1) = Aa * xhat(:, i) + Ba * u(:, i);
    end
    xhat_ = xhat_(:, 1:end-1);
else
    xhat_ = data.xhatk;
end

% Estimate L-innovations
inntrun = y(:, start+1:end) - Ca * xhat_(:, start+1:end);

% Optional - use to look at plots to check if error looks white:
if plot_flag == 1
    figure
    for i = 1:1:p
        if p > 1
            subplot(ceil(p / 2), 2, i)
        end
        plot(inntrun(i, :), 'm', 'linewidth', 2)
        xlabel('Time')
        ylabel(sprintf('Innovation of y_%i', i))
    end
end

% Calculation of Autocorrelations for one column ALS
datatrun = datapts - start;
Eyy = [];
for i = 0:N - 1
    temp = inntrun(:, i+1:end) * inntrun(:, 1:end-i)';
    temp = temp ./ (datatrun - i);
    Eyy = [Eyy; temp];
end
Eyy = Eyy(:);

if dataweight == 1
    % Calculation of Weighting matrix for one column ALS
    Nd = size(inntrun, 2);
    nt = Nd - N + 1;
    trials = 2 * N;
    datapts2 = floor(nt/trials);
    Eyy_store2 = [];
    Yysmall = [];
    covb = zeros(pa^2*N, pa^2*N);
    for i_trial = 1:1:trials
        Yysmall = zeros(N*p, datapts2);
        for i = 1:1:datapts2
            yyst = inntrun(:, (i - 1)*trials+i_trial:(i - 1)*trials+(i_trial - 1)+N); %Matlab compatible
            Yysmall(:, i) = yyst(:);
        end
        Py = cov(Yysmall');
        Px = Py(1:p, 1:p);
        Pyx = Py(:, 1:p);
        covb = covb + 1 / datapts2 * (kron(Px, Py) + comm_mat(pa, N * pa) * kron(Pyx, Pyx'));
    end
    covb = covb / trials; %/norm(covb);
    Wm = pinv(covb);
else
    % Use identity-based weighting
    Wm = eye(length(Eyy));
end

% Building the constant matrix for the LS problem

Ain = Aa - Aa * L * Ca;

OO = [];
temp = eye(n);
for i = 1:N
    OO = [OO; Ca * temp];
    temp = temp * Ain;
end

% temporary variables
M1 = zeros(na^2, ga^2);
i = 1;
for j = 1:ga
    for k = 1:ga
        II = zeros(ga);
        II(k, j) = 1;
        t1 = dlyap(Ain, Ga*II*Ga');
        M1(:, i) = t1(:);
        i = i + 1;
    end
end
M2 = zeros(na^2, pa^2);
i = 1;
for j = 1:pa
    for k = 1:pa
        II = zeros(pa);
        II(k, j) = 1;
        t2 = dlyap(Ain, Aa*L*II*L'*Aa');
        M2(:, i) = t2(:);
        i = i + 1;
    end
end

% Single column ALS method

PSI = eye(pa);
for i = 1:N - 1
    PSI = [PSI; -Ca * Ain^(i - 1) * Aa * L];
end

OOtemp = kron(Ca, OO);
PSItemp = kron(eye(pa), PSI);

LH1 = OOtemp * M1;
LH2 = OOtemp * M2 + PSItemp;


if Rsym == 1
    LHSsingc = [LH1 * symtran(ga), LH2 * symtran(pa)]; %Adds symmetric to R_v
else
    % Adds symmetric constraint to Q_w and diagonal constraint to R_v
    LHSsingc = [LH1 * symtran(ga), LH2(:, 1:pa + 1:pa^2)];
end

M = kron(Ca, eye(n)) * M1 * symtran(ga);
Mrank = rank(M, 1e-4);
[nr, nc] = size(M);
if nc > Mrank
    fprintf('The covariance estimates are not unique!\n Use trade-off curve to find lowest rank solution\n')
end
if isfield(estimator, 'Q') && isfield(estimator, 'R')
    Q0 = celldiag({estimator.Q, estimator.R});
else
    Q0 = eye(ga+pa); % Feasible initial guess
end
numax = 1;

if sdp == 0
    % Solve least-squares problem without semi-definite constraints
    if lam_vec ~= 0
        warning('als_sdp_mrQ: Cannot solve indefinite problem with penalty on trQ');
    end
    if ~isempty(ver('matlab'))
        warning('off', 'MATLAB:nearlySingularMatrix')
    elseif ~isempty(ver('octave'))
        warning("off", "Octave:singular-matrix-div")
    else
        error('Could not determine interpreter version.');
    end
    time1 = cputime;
    X = (LHSsingc' * Wm * LHSsingc) \ LHSsingc' * Wm * Eyy;
    timespent = cputime - time1;
    Qest_cell(1) = {reshape(symtran(ga) * X(1:ga * (ga + 1) / 2), ga, ga)};
    if Rsym == 1
        Rest_cell(1) = {reshape(symtran(pa) * X(ga * (ga + 1) / 2 + 1:end), pa, pa)};
    else
        Rest_cell(1) = {diag(X(end -pa + 1:end))};
    end
    Phi = (Eyy - LHSsingc * X)' * Wm * (Eyy - LHSsingc * X);
    phi0 = Phi;
    Phi_tot = Phi;
    Iter = 0;
    Iter_maxed = 0;
    trQ = trace(Qest_cell{1}*diag(trstates));
else
    % Solve ALS with semi-definite constraints
    % Scale weighting matrix and initial objective function
    [QR0, phi0] = sdp_QR_mrQ(LHSsingc, Eyy, Wm, Q0, 0, numax, pa, Rsym, trstates);
    Wm = Wm / phi0; % Now phi = 1 for rho = 0
    % Loop through and solve for each value of lam_vec
    for i = 1:1:length(lam_vec)
        lam = lam_vec(i);
        time1 = cputime;
        [QR1, phi, phi_tot, iter, nsq, iter_maxed] = sdp_QR_mrQ(LHSsingc, Eyy, Wm, Q0, lam, numax, pa, Rsym, trstates);
        timespent(i, 1) = cputime - time1;
        Qest = QR1(1:ga, 1:ga);
        Rest = QR1(ga+1:end, ga+1:end);
        Qest_cell(i) = {(Qest + Qest') / 2};
        Rest_cell(i) = {(Rest + Rest') / 2};
        Phi(i, 1) = phi;
        Phi_tot(i, 1) = phi_tot;
        Iter(i, 1) = iter;
        Iter_maxed(i, 1) = iter_maxed;
        trQ(i, 1) = trace(Qest_cell{i}*diag(trstates));
    end
end


% Form matrices of autocovariances from data and theory
if Rsym == 1
    %Eyy_thry = LHSsingc*[symtranT(ga)*vec(Qest_cell{1}); symtranT(pa)*vec(Rest_cell{1})];
    Eyy_thry = LHSsingc * [symtranT(ga) * Qest_cell{1}(:); symtranT(pa) * vec(Rest_cell{1})]; %Matlab compatible
else
    %Eyy_thry = LHSsingc*[symtranT(ga)*vec(Qest_cell{1}); diag(Rest_cell{1})];
    Eyy_thry = LHSsingc * [symtranT(ga) * Qest_cell{1}(:); diag(Rest_cell{1})]; %Matlab compatible
end
Eyy = reshape(Eyy, N*pa, pa);
Eyy_thry = reshape(Eyy_thry, N*pa, pa);
for i = 1:1:pa
    for j = 1:1:pa
        for k = 1:1:N
            bhat_mtr(k, pa*(i - 1)+j) = Eyy(pa*(k - 1)+i, j);
            bhat_mtr_thry(k, pa*(i - 1)+j) = Eyy_thry(pa*(k - 1)+i, j);
        end
    end
end

npts = datapts - start - N;
for i = 1:1:pa; for j = 1:1:pa
        cov_bound(i, j) = sqrt(bhat_mtr(1, pa * (i - 1) + i)*bhat_mtr(1, pa * (j - 1) + j)/(npts));
    end; end
cov_bound = cov_bound(:);

if plot_flag == 1
    figure;
    for i = 1:1:pa^2
        subplot(pa, pa, i);
        if ~isempty(ver('matlab'))
            plot(0:1:N-1, bhat_mtr(:, i), '-*', 0:1:N-1, bhat_mtr_thry(:, i), '-o', 'linewidth', 2);
        elseif ~isempty(ver('octave'))
            plot(0:1:N-1, bhat_mtr(:, i), '-*', 'linewidth', 2, 0:1:N-1, bhat_mtr_thry(:, i), '-o', 'linewidth', 2); % Octave
        else
            error('Could not determine interpreter version.');
        end

        xlabel('Lag')
        hold on;
        plot([0; N], 2*kron([1, -1; 1, -1], [cov_bound(i)]), 'k');
        if i == 1
            legend('Data', 'Fit');
        end
        if i <= pa
            title('Autocorrelation');
        end
        if i > pa^2 - pa
            xlabel('Lag');
        end
        hold off;
    end
end
