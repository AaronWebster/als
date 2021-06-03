%% Example of calculation of minimum rank Q using the TRADE-OFF CURVE method

clear all

close all

rng(200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLANT SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aa = [0.53565, 0.23972; 0.75473, 0.40629];

Ca = [1, 1];

[m, n] = size(Ca);
R_v = 0.1 * eye(m);
Q_w = 0.1 * eye(n);
Ga = eye(n);

[pa, ~] = size(Ca);
[na, ga] = size(Ga);

n = na;
p = pa;
g = ga;

datapts = 5000;

[L, P] = dlqe(Aa, eye(n), Ca, Q_w, R_v);

xhat = zeros(na, datapts);
xhat_ = zeros(na, datapts);

x(:, 1) = 10 * ones(na, 1); %% x0

xhat_(1:na, 1) = x(:, 1); %% assume initial state perfectly known

for i = 1:datapts

    y(:, i) = Ca * x(:, i) + sqrt(R_v) * randn(pa, 1);
    xhat(:, i) = xhat_(:, i) + L * (y(:, i) - Ca * xhat_(:, i));
    x(:, i+1) = Aa * x(:, i) + Ga * (sqrt(Q_w) * randn(ga, 1));
    xhat_(:, i+1) = Aa * xhat(:, i);

end

model.A = Aa;
model.C = Ca;
model.G = Ga;
data.datapts = datapts;
data.yk = y;
data.xhatk = xhat_(:, 1:end-1);
data.start = 100;

N = 15;

estimator.Q = Q_w;
estimator.R = R_v;
rho = logspace(-6, 6, 25)';
[Qest, Rest, trace_Q, phi_Q] = als_sdp_mrQ(data, N, model, estimator, 'rho_values', rho);

% Build Tradeoff plots using calculated data
for i = 1:length(Qest);
    rank_Q(i) = rank(Qest{i}, 1e-4);
end
figure(1)

plot(phi_Q, trace_Q)
xlabel('\phi')
ylabel('tr(Q)')

figure(2)
subplot(3, 1, 1)
semilogx(rho, phi_Q); ylabel('\phi')
subplot(3, 1, 2);
semilogx(rho, rank_Q, 'ro'); ylabel('rank(Q)')
subplot(3, 1, 3)
semilogx(rho, trace_Q, 'k'); ylabel('tr(Q)')
