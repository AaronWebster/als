% Open source code of linear time-varying autocovariance least-squares (LTV-ALS)
% technique to estimate covariances for nonlinear and time-varying models

% This function is used to build up the stacks for the ALS
% structure of eq. (17) of Lima & Rawlings (2010). These stacks are used to
% calculate (Q,R) using semidefinite programming (SDP)
% Please see file case1_als.m for an example of the LTV-ALS implementation for 
% a simple reactor system 

function [Qdet, Rdet, LHS, Eyy] = ltv_als(yinn, N, bgn, Ain, Ak, Ck, Gk, Lk)

  % Inputs:
    % yinn: matrix with vectors of innovations
    % N: number of lags used in the autocovariance matrix
    % bgn: ALS starting time (k) 
    % Ain = Ak-Ak*Lk*Ck
    % Ak, Ck, Gk and Lk (EKF gain): time-varying system matrices

  % Outputs:
    % Qdet, Rdet: covariance estimates
    % LHS: left hand side ALS matrix
    % Eyy: data vector

  % Load innovations data from k up to k+N-1
  yinn = yinn(:,bgn:end);

  % Define p, n and g
  [p,n] = size(Ck{1});
  g = columns(Gk{1});

  % Build data vector ([Rkhat(N)]_s):
  Eyyfl = vec(yinn)*yinn(:,1)';
  Eyy = vec(Eyyfl);

  % Calculate Gamma
  Apr{1} = eye(n);
  for j = 2:N; Apr{j} = Ain{bgn+j-2}*Apr{j-1}; endfor

  Gamma = [];
  for i = bgn:-1:2
    tempG = [];
    for j = 1:N
      tempG = [tempG; Ck{bgn+j-1}*Apr{j}];
      Apr{j} = Apr{j}*Ain{i-1};
    endfor
    Gamma = [tempG Gamma];
  endfor

  % Calculate Omega1 and Omega2
  for j = 1:bgn-1; AL{j} = -Ak{j}*Lk{j}; endfor
  
  Omega1= blkdiag(Gk{1:bgn-1});
  Omega2 = blkdiag(AL{1:bgn-1});

  % Calculate PSI
  PSI = eye(p); 
  Apr1 = eye(n);

  for i = 1:N-1;
    PSI = [PSI; -Ck{bgn+i}*Apr1*Ak{bgn}*Lk{bgn}];
    Apr1 = Ain{bgn+i}*Apr1;
  endfor

  % Gamma x Omega1
  Gam1 = Gamma*Omega1;
  % Gamma1 x Omega1
  Gam11 = Gamma(1:p,:)*Omega1;
  % Gamma x Omega2
  Gam2 = Gamma*Omega2;
  % Gamma1 x Omega2
  Gam22 = Gamma(1:p,:)*Omega2;

  LHS_Q = 0;
  LHS_R = 0;

  for i = 1:bgn-1

    ee = eye(bgn-1)(:,i);

    LHS_Q = LHS_Q + kron(Gam11*kron(ee,eye(g)), Gam1*kron(ee,eye(g)));
    LHS_R = LHS_R + kron(Gam22*kron(ee,eye(p)), Gam2*kron(ee,eye(p)));

  endfor

  LHS_R = LHS_R + kron(eye(p), PSI);

  LHS = [LHS_Q*duplication_matrix(g) LHS_R*duplication_matrix(p)];

  X = ols(Eyy, LHS);

  Qdet = reshape(duplication_matrix(g)*X(1:g*(g+1)/2,1),g,g);

  Rdet = reshape(duplication_matrix(p)*X(g*(g+1)/2+1:end,1),p,p);

endfunction
