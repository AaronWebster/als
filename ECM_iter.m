% Function to calculate the state error covariance matrix for the DLQE
% using iterations

function [L, P] = ECM_iter(A, G, C, Q, R)

[~, n] = size(C);
Pold = zeros(n);
P = eye(n);
tol = 1e-12;
iter = 0;

while (norm(Pold - P, 'fro') > tol)
    iter = iter + 1;
    Pold = P;
    P = A * Pold * A' + G * Q * G' - A * Pold * C' * inv(C*Pold*C'+R) * C * Pold * A';
    %    norm(Pold-P,"fro")
end

L = P * C' * inv(C*P*C'+R);
