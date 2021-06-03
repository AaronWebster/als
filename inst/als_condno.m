function [cond_no, svals, Mcond] = als_condno(N, model, estimator, Rform)

% Returns condition number and largest and smaller singular values for
% the ALS matrix (script A) formed from N, model, estimator.  Does not
% solve ALS problem

% [cond_no svals Mcond] = als_condno(N,model,estimator,Rform)

% Function Inputs :
% N (window size),
% model.A,
% model.C,
% model.G (optional, default is identity matrix)
% estimator.L (initial estimator gain - optional)
% or estimator.Q and estimator.R (process and measurement noise
% covariances, used to calculate estimator gain),
% Rform (optional input to specify structure of R_v.  If
% Rform = "sym", R_v is symmetric, if Rform = "diag", R_v is
% diagonal.  Default is diagonal)
%
% Function outputs :
% cond_no (condition number of script A)
% svals (singular values of script A)
% Mcond (condition number of M = kron(C,I)*(I-kron(A,A))^-1*kron(G,G)*D_g


Aa = model.A;
Ca = model.C;
if isfield(model, 'G')
    Ga = model.G;
else
    Ga = eye(size(Aa));
end
[n, g] = size(Ga);
p = size(Ca, 1);

if isfield(estimator, 'L')
    L = estimator.L;
else
    L = dlqe(Aa, Ga, Ca, estimator.Q, estimator.R);
end

na = n;
ga = g;
pa = p;

if nargin == 3
    Rsym = 0;
elseif strcmpi(Rform, 'sym')
    Rsym = 1;
elseif strcmpi(Rform, 'diag')
    Rsym = 0;
else
    warning('als_sdp_mrQ: Unknown structure type for R; defaulting to diagonal');
    Rsym = 0;
end


% Building the constant matrix for the LS problem

Ain = Aa - Aa * L * Ca;

OO = [];
temp = eye(n);
for i = 1:N
    OO = [OO; Ca * temp];
    temp = temp * Ain;
end
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
    % Adds symmetric to R_v
    LHSsingc = [LH1 * symtran(ga), LH2 * symtran(pa)];
else
    % Adds symmetric constraint to Q_w and diagonal constraint to R_v
    LHSsingc = [LH1 * symtran(ga), LH2(:, 1:pa + 1:pa^2)];
end
M = kron(Ca, eye(n)) * M1 * symtran(ga);
Mcond = cond(M);

cond_no = cond(LHSsingc);
svals = svd(LHSsingc);
