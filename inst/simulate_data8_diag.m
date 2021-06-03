%% Diagonal ALS Implementation Example

clear all

rng(100);
Aa = diag([0.1, 0.2, 0.3]);
Aa(1, 3) = 0.1;

%Ga = [1;0.2;0.3];
Ga = eye(3);

%Ca = [0.1 0.2 0];
Ca = eye(3);

Q_w = diag([0.5, 0.2, 0.1]);
%Q_w = 0.5;
%R_v = 0.1;
R_v = diag([0.5, 0.2, 0.8]);

S = dlyap(Aa, Ga*Q_w*Ga');

[vec_Qw, eig_Qw] = eig(Q_w);
[vec_Rv, eig_Rv] = eig(R_v);
mult_Qw = vec_Qw * sqrt(eig_Qw);
mult_Rv = vec_Rv * sqrt(eig_Rv);

%%% initial guesses

G_hat = eye(3);
Qw_hat = diag([1, 2, 3]);
%Qw_hat= 5*Q_w;
Rv_hat = 1e-3 * R_v;

[pa, ~] = size(Ca);
[na, ga] = size(Ga);

n = na;
p = pa;
g = ga;

datapts = 5000;

L = dlqe(Aa, G_hat, Ca, Qw_hat, Rv_hat);
%[L,P]=dlqe(Aa,Ga,Ca,Q_w,R_v);
%L = zeros(n,p);

P = dlyap((Aa-Aa*L*Ca), [Ga, -Aa * L]*[Q_w, zeros(g, p); zeros(p, g), R_v]*[Ga, -Aa * L]');

xhat = zeros(na, datapts);
xhat_ = zeros(na, datapts);

x(:, 1) = 10 * ones(na, 1); %% x0

xhat_(1:na, 1) = x(:, 1); %% assume initial state perfectly known

for i = 1:datapts
    y(:, i) = Ca * x(:, i) + mult_Rv * randn(pa, 1);
    xhat(:, i) = xhat_(:, i) + L * (y(:, i) - Ca * xhat_(:, i));
    x(:, i+1) = Aa * x(:, i) + Ga * (mult_Qw * randn(ga, 1));
    xhat_(:, i+1) = Aa * xhat(:, i);

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SETUP ALS PROBLEM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

model.A = Aa;
model.C = Ca;
model.G = G_hat;
model.xhat0 = xhat_(:, 1);

data.datapts = datapts;
data.yk = y;
data.xhatk = xhat_(:, 1:end-1);
data.start = 100;

N = 15;

estimator.L = L;
%estimator.Q = Qw_hat;
%estimator.R = Rv_hat;

[Qest, Rest, Lest, As, bhat] = als_diag(data, N, model, estimator);