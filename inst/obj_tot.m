% Evaluate objection function value, including penalty on tr(Q) and log-barrier term
function y = obj_tot(Q,A,b,W,nu,vlam_mtr)
[r,p] = chol(Q);
% p = 0 if Q>0 
if p == 0
  y = (A*Q(:) - b)'*W*(A*Q(:)-b) + trace(Q*vlam_mtr) - nu*2*sum(log(diag(r)));
  % Note - using Cholesky decomposition, 2*sum(log(diag(r))) is equivalent to sum(log(eig(Q)))
else
% if Q is infeasible, objective function is infinity
  y=1e100;
end

