
%% The Diagonal ALS function.  Uses Diagonal ALS form to estimate only
%% the diagonal elements of Qw and Rv.
%% Function Inputs : data.yk (measurements), data.uk (inputs),
%% data.xhatk (state estimates- optional), data.datapts (number of data
%% points considered), data.start (data to be ignored in the beginning 
%% till initial condition is negligible),
%% model.A, model.B (optional), model.C, model.G , N (window size),
%% estimator.Q (Qw initial guess), estimator.R (Rv initial guess),
%% estimator.L (initial estimator gain - optional).
%% Function outputs : estimated Qw, estimated Rv, estimated filter gain
%% L, and ALS LHS matrix As and RHS vector bhat.
%% Please see file simulate_data8_diag.m for an example of the 
%% Diagonal ALS implementation for a simple system.
%% For Matlab implementation, the call of qp must be replaced by
%% quadprog - see line 133 


function [Qest,Rest,Lest,As,bhat] = als_diag(data,N,model,estimator)

if( nargin ~= 4) error('als: invalid number of arguments');
end 
  
datapts = data.datapts;
Aa = model.A;
Ca = model.C;
Ga = model.G;
[n,g] = size(Ga);
p = size(Ca,1);
if( isfield(model, 'B')) Ba = model.B; m = columns(Ba);
else Ba = zeros(n);m = n; data.uk = zeros(m,datapts);
end
start = data.start;
  
na = n;ga = g;pa = p;
%% Estimator Simulation
y = data.yk;
u = data.uk;

if( isfield(estimator,'L')) L = estimator.L;
else L = dlqe(Aa,Ga,Ca,estimator.Q,estimator.R);
end

if(~isfield (data,'xhatk'))
  xhat=zeros(na,datapts);
  xhat_ = zeros(n,datapts);
  xhat_(1:n,1) = model.xhat0;
  for i = 1:datapts
    xhat(:,i) = xhat_(:,i) + L*(y(:,i)-Ca*xhat_(:,i));
    xhat_(:,i+1) = Aa*xhat(:,i) + Ba*u(:,i);
  end
  xhat_ = xhat_(:,1:end-1);
else
  xhat_ = data.xhatk;
end

inntrun = y(:,start+1:end)-Ca*xhat_(:,start+1:end);
datatrun = datapts-start;

%% Calculation of Autocorrelations for one column ALS
Eyy=[];
for i=0:N-1
  temp=inntrun(:,i+1:end)*inntrun(:,1:end-i)';
  temp=temp./(datatrun - i);
  Eyy = [Eyy; temp];
end
Eyy= Eyy(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Building the constant matrix for the LS problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
Ain = Aa-Aa*L*Ca;

OO = [];
temp = eye(n);
for i = 1:N
  OO = [OO ; Ca*temp];
  temp = temp*Ain;
end

%% temporary variables
i=1;
for j = 1:ga
  for k = 1:ga
    II = zeros(ga);
    II(k,j) = 1;
    t1 = dlyap(Ain,Ga*II*Ga');
    M1(:,i) = t1(:);
    i = i+1;
  end
end
i=1;
for j = 1:pa
  for k = 1:pa
    II = zeros(pa);
    II(k,j) = 1;
    t2 = dlyap(Ain,Aa*L*II*L'*Aa');
    M2(:,i) = t2(:);
    i = i+1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagonal ALS method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
PSI = eye(pa);
for  i= 1:N-1
  PSI = [PSI; -Ca*Ain^(i-1)*Aa*L];
end
  
OOtemp = kron(Ca,OO);
PSItemp = kron(eye(pa),PSI);

LH1 = OOtemp*M1;   
LH2 = OOtemp*M2 + PSItemp;

As_diag = [LH1(:,1:ga+1:ga^2) LH2(:,1:pa+1:pa^2)];
  
% Testing the uniqueness of covariance estimates  
Arank = rank(As_diag, 1e-4);

[nr, nc] = size(As_diag);

if nc > Arank
  printf('Warning: Covariance estimates are not unique!\n')
  pause
end

Xest_diag = qp(ones(ga+pa,1), As_diag'*As_diag, -As_diag'*Eyy, [], [], zeros(ga+pa,1), [], [], eye(pa+ga), []); 

% For Matlab, call line 134 instead of 131:
%Xest_diag = quadprog(As_diag'*As_diag, -As_diag'*Eyy, -eye(pa+ga), zeros(ga+pa,1));

if prod(Xest_diag)==0 
  printf('Warning: Covariance estimate(s) is (are) at constraints! You may have bad data! \n')
end

Qest=diag(Xest_diag(1:ga));
Rest=diag(Xest_diag(ga+1:end));

Lest = ECM_iter(Aa,Ga,Ca,Qest,Rest);

As = As_diag;
bhat = Eyy;