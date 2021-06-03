%% Evaluates only least-squares portion of objective 
%% (for use in choosing Q from tradeoff curve)
function y = obj_ls(Q,A,b,W)
	y = (A*Q(:) - b)'*W*(A*Q(:) - b);
	if(imag(y)~=0) y=1e10;end
