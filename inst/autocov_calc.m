function [bhat_mtr cov_bound bhat] = autocov_calc(Yi,N)
%% function [bhat_mtr cov_bound bhat] = autocov_calc(Yi,N)
%% Given data Yi, estimates autocovariances up to lag N-1
%%
%% Function inputs :
%% Yi - Nd x p matrix of data, where Nd is number of samples
%% and p is number of measured variables
%% N - number of lags
%%
%% Function outputs :
%% bhat_mtr - N x p matrix of autocovariances
%% bhat_mtr = 
%% [ <y_1(k)y_1(k)>   <y_1(k)y_2(k)>   ...   <y_p-1(k)y_p(k)>     <y_p(k)y_p(k)>
%%	...	              ...                     ...                 ...  
%% <y_1(k)y_1(k-N+1)> <y_1(k)y_2(k-N+1)> ... <y_p-1(k)y_p(k-N+1)> <y_p(k)y_p(k-N+1)>];
%% cov_bound - estimated standard deviations of autocovariances
%% cov_bound = [sigma_11 sigma_12 ... sigma_1p sigpa_pp]'
%% bhat - vector of autocovariances as used in ALS 
Eyy=[];
p = size(Yi,1);
for i=0:N-1
  temp=Yi(:,i+1:end)*Yi(:,1:end-i)';
  temp=temp./(length(Yi)-i);
  Eyy = [Eyy; temp];
end
bhat= Eyy(:);

for i = 1:1:p
	for j = 1:1:p
		for k = 1:1:N
			bhat_mtr(k,p*(i-1)+j) = Eyy(p*(k-1)+i,j);
		end
	end
end

npts = size(Yi,2)-N;
for i = 1:1:p; for j = 1:1:p
	cov_bound(i,j) = sqrt(bhat_mtr(1,p*(i-1)+i)*bhat_mtr(1,p*(j-1)+j)/(npts));
end; end
cov_bound = cov_bound(:);
