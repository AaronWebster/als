% Motivating Example (Section 2.2 of Murali Rajamani's Thesis)

function [Qest_vectr, Rest_vectr] = motivating_example

clear all

rng(100);

for k = 1:200
    Aa = diag([0.1, 0.2, 0.3]);
    Aa(1, 3) = 0.1;

    Ga = [1; 2; 3];
    Ca = [0.1, 0.2, 0];

    Q_w = 0.5;
    R_v = 0.1;

    S = dlyap(Aa, Ga*Q_w*Ga');

    [vec_Qw, eig_Qw] = eig(Q_w);
    [vec_Rv, eig_Rv] = eig(R_v);
    mult_Qw = vec_Qw * sqrt(eig_Qw);
    mult_Rv = vec_Rv * sqrt(eig_Rv);

    % initial guesses

    G_hat = [1; 2; 3];
    Qw_hat = .4 * Q_w;
    Rv_hat = 4 * R_v;

    [pa, ~] = size(Ca);
    [na, ga] = size(Ga);

    n = na;
    p = pa;
    g = ga;

    datapts = 1000;

    L = dlqe(Aa, G_hat, Ca, Qw_hat, Rv_hat);
    %[L,P]=dlqe(Aa,Ga,Ca,Q_w,R_v);
    %L = zeros(n,p);

    P = dlyap((Aa-Aa*L*Ca), [Ga, -Aa * L]*[Q_w, zeros(g, p); zeros(p, g), R_v]*[Ga, -Aa * L]');

    xhat = zeros(na, datapts);
    xhat_ = zeros(na, datapts);

    x(:, 1) = 10 * ones(na, 1); % x0

    xhat_(1:na, 1) = x(:, 1); % assume initial state perfectly known

    for i = 1:datapts

        y(:, i) = Ca * x(:, i) + mult_Rv * randn(pa, 1);
        xhat(:, i) = xhat_(:, i) + L * (y(:, i) - Ca * xhat_(:, i));
        x(:, i+1) = Aa * x(:, i) + Ga * (mult_Qw * randn(ga, 1));
        xhat_(:, i+1) = Aa * xhat(:, i);

    end

    % SETUP ALS PROBLEM

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

    [Qest, Rest] = als_sdp_mrQ(data, N, model, estimator);

    Qest_vectr(k, 1) = Qest{1};
    Rest_vectr(k, 1) = Rest{1};

    k = k + 1;

end

plot(Qest_vectr, Rest_vectr, 'o')
axis([0.1, 1, 0, 0.2])
xlabel('Qw')
ylabel('Rv')
line([0.3, 0.7], [0.1, 0.1])
line([0.5, 0.5], [0.05, 0.15])
legend('estimate', 'true value')
