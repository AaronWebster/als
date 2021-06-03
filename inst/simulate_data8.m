% ALS Implementation Example

clear all

rng(100);
Aa = diag([0.1 0.2 0.3]); Aa(1,3)=0.1; Aa(1,2) = 0.1;


%Ga = [1;0.2;0.3];
Ga = eye(3);

%Ca = [0.1 0.2 0];
Ca = eye(3);

Q_w = diag([0.5,0.2,0.1]);
%Q_w = 0.5;
%R_v = 0.1;
R_v = diag([0.5,0.2,0.8]);

S = dlyap(Aa,Ga*Q_w*Ga');

[vec_Qw,eig_Qw] = eig(Q_w);
[vec_Rv,eig_Rv] = eig(R_v);
mult_Qw = vec_Qw*sqrt(eig_Qw);
mult_Rv = vec_Rv*sqrt(eig_Rv);

% initial guesses

G_hat = eye(3);
Qw_hat = diag([1,2,3]);
%Qw_hat= 5*Q_w;
Rv_hat= 1e-3*R_v;

[pa,na]=size(Ca);
[na,ga]=size(Ga);

n=na;
p=pa;
g=ga;

datapts = 5000;

L=dlqe(Aa,G_hat,Ca,Qw_hat,Rv_hat);
%[L,P]=dlqe(Aa,Ga,Ca,Q_w,R_v);
%L = zeros(n,p);

P = dlyap((Aa-Aa*L*Ca),[Ga -Aa*L]*[Q_w zeros(g,p);zeros(p,g) R_v]*[Ga -Aa*L]');

xhat=zeros(na,datapts);
xhat_=zeros(na,datapts);

x(:,1) = 10*ones(na,1);  % x0

xhat_(1:na,1) = x(:,1); % assume initial state perfectly known

for i = 1:datapts
  
  y(:,i) = Ca*x(:,i)+mult_Rv*randn(pa,1);
  xhat(:,i) = xhat_(:,i) + L*(y(:,i)-Ca*xhat_(:,i));
  x(:,i+1) = Aa*x(:,i) +Ga*(mult_Qw*randn(ga,1));
  xhat_(:,i+1) = Aa*xhat(:,i);

end

% SETUP ALS PROBLEM

model.A = Aa;
model.C = Ca;
model.G = G_hat;
model.xhat0 = xhat_(:,1);

data.datapts = datapts;
data.yk = y;
% data.xhatk = xhat_(:,1:end-1);
data.start = 100;

N = 15;

estimator.L = L;
%estimator.Q = Qw_hat;
%estimator.R = Rv_hat;

% Using updated codes:
% This is in general how to call the ALS method:
[Qest_cell,Rest_cell] = als_sdp_mrQ(data,N,model,estimator);
Qest1 = Qest_cell{1}
Rest1 = Rest_cell{1}
% Without semidefinite constraints
% you usually should keep the semidefinite constraints, but could remove them to see how good your model is
[Qest_celln,Rest_celln] = als_sdp_mrQ(data,N,model,estimator,'sdp',0,'plot',0);
Qest_indef = Qest_celln{1}
Rest_indef = Rest_celln{1}
% With identity weighting
% Only use identity weighting if you have a good reason; the default of data-based weighting gives lower variance estimates
[Qest_celli,Rest_celli] = als_sdp_mrQ(data,N,model,estimator,'weight','I','plot',0);
QestI = Qest_celli{1}
RestI = Rest_celli{1}
% Tradeoff curve example
% This example shows what to do if you have fewer outputs than states
% here we pretend that we only know the first row of C
% we generate a tradeoff curve and look for the elbow in the curve
data2 = data;
data2.yk = y(1,:);
model2 = model;
model2.C = model.C(1,:);
rho_vec = logspace(-6,6,25);
estimator2.L = dlqe(Aa,G_hat,Ca(1,:),Qw_hat,Rv_hat(1,1));
[Qest_cell2,Rest_cell2,trQ,Phi2] = als_sdp_mrQ(data2,N,model2,estimator2,'rho_values',rho_vec);
figure; plot(Phi2,trQ,'-o',Phi2(9),trQ(9),'r*')
Qest2 = Qest_cell2{9}
Rest2 = Rest_cell2{9}
