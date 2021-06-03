function [b_mtr, b] = autocov_calc_thry(A, C, L, Q, R, N)

%% function [b_mtr b] = autocov_calc_thry(A,C,L,Q,R,N)

%% Calculates theoretical autocovariances given system matrices

%% and noise covariances

%%

%% Function inputs :

%% System matrices A and C

%% Estimator gain L

%% Process noise covariance Q (note G is assumed identity)

%% Measurement noise covariance R

%% Number of lags N

%%

%% Function outputs :

%% b_mtr - N x p matrix of autocovariances

%% b_mtr =

%% [ E(y_1(k)y_1(k))   E(y_1(k)y_2(k))   ...   E(y_p-1(k)y_p(k))     E(y_p(k)y_p(k))

%%	...	              ...                     ...                 ...

%% E(y_1(k)y_1(k-N+1)) E(y_1(k)y_2(k-N+1)) ... E(y_p-1(k)y_p(k-N+1)) E(y_p(k)y_p(k-N+1))];

%% cov_bound - estimated standard deviations of autocovariances

%% b - vector of autocovariances as used in ALS

Abar = A - A * L * C;
P = dlyap(Abar, Q+A*L*R*L'*A');
Eyy = C * P * C' + R;
for i = 1:1:N - 1;
    Eyy = [Eyy; C * Abar^i * P * C' - C * Abar^(i - 1) * A * L * R];
end
b = Eyy(:);

p = size(C, 1);
for i = 1:1:p
    for j = 1:1:p
        for k = 1:1:N
            b_mtr(k, p*(i - 1)+j) = Eyy(p*(k - 1)+i, j);
        end
    end
end
