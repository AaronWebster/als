% SemiDefinite programming subroutine.  Uses Newton steps to satisfy
% the KKT conditions with log-barrier function

function [Q, phi2, phi, itertot, nsq, iter_maxed, phi_uw] = sdp_QR_mrQ(A, b, W, Q0, lam, nu, pa, Rsym, trstates)

if ~isempty(ver('matlab'))
    warning('off', 'MATLAB:nearlySingularMatrix')
elseif ~isempty(ver('octave'))
    warning("off", "Octave:singular-matrix-div")
else
    error('Could not determine interpreter version.');
end

ntot = size(Q0, 1);

na = ntot - pa;

ind = reshape([1:ntot^2], ntot, ntot);
ind = (ind & celldiag({rand(na), rand(pa)})) .* ind;
Az = zeros(size(A, 1), ntot^2);

% Reconstructs least squares matrix for R:
if (pa ~= 0)
    LH2 = A(:, na*(na + 1)/2+1:end);
    if Rsym == 1
        % For symmetric R matrix
        mtrx = symtranT(pa);
    else
        % For diagonal R matrix
        mtrx = zeros(pa, pa^2);
        for i = 1:1:pa
            mtrx(i, (i - 1)*(pa + 1)+1) = 1;
        end
    end
    LH2 = LH2 * mtrx;
end

% Reconstructs least-squares matrix for Q:
A = [A(:, 1:na * (na + 1) / 2) * symtranT(na), LH2];
idx = setdiff(ind(:), 0);
Az(:, idx) = A;
alphmax = 1; % maximum step size

% Initialize
tol = 1e-10; % Absolute tolerance for stopping criterion
reltol = 1e-5; % Relative tolerance for stopping criterion
itertot = 0;
maxiter = 100; % Maximum iterations for each log-barrier penalty
Q = Q0; % Initial guess
Qin = inv(Q);
con = 1;
con1 = con;
con2 = con;
tol1 = tol;
tol2 = tol;
phi = 0;
M = 2 * Az' * W * Az; % Hessian part of norm fit
dtemp = 2 * Az' * W * b; % data to be used for gradient

vlam_mtr = lam * celldiag({diag(trstates), zeros(pa)});
vlam = vlam_mtr(:);
iter = 0;

phi_old = 1e3;
nu_first = 1;
while (nu > 1e-12)
    % Decrease penalty on log-barrier term (for sdp constraints)
    if iter > 0
        if iter / maxiter < 0.75;
            nu = nu * 0.75;
        else
            nu = nu * (0.95 * iter) / maxiter;
        end
    end
    iter = 1;
    phi = obj_tot(Q, Az, b, W, nu, vlam_mtr); % Objective function value
    con1 = abs(phi-phi_old); % Stopping criterion 1
    tol1 = tol + reltol * abs(phi); % Tolerance for stopping criterion 1
    % Optimize for each penalty on log-barrier term
    while ((con1 > tol1 || con2 > tol2) && iter < maxiter)
        t1 = cputime;
        itertot = itertot + 1;
        iter = iter + 1;
        LH = M + nu * kron(Qin, Qin); % Hessian
        RH = M * Q(:) - dtemp + vlam - nu * Qin(:);
        delQ = -reshape(pinv(LH)*RH, ntot, ntot); %Newton direction
        %delQ = -reshape(LH\RH,ntot,ntot); % Can't calculate this way if Hessian is poorly conditioned
        % Optimal Step Size
        alph = golden_section_Q_mrQ(0, alphmax, Q, delQ, Az, b, W, nu, vlam_mtr);
        Qold = Q;
        Qnew = Qold + 0.99 * alph * delQ;
        if (nu_first < 1 && min(eig(Qnew)) < 0);
            % Keep solution from becoming infeasible
            Qnew = Qold;
        end
        % Update Q
        Q = Qnew;
        Q = (Q + Q') / 2;
        Qin = pinv(Q);
        % Calculate stopping criteria and tolerances:
        phi_old = phi;
        phi = obj_tot(Q, Az, b, W, nu, vlam_mtr);
        con1 = abs(phi-phi_old);
        con2 = norm(Q-Qold);
        tol1 = tol + reltol * abs(phi);
        tol2 = tol + reltol * norm(Q);
        t2 = cputime;
        if t2 - t1 > 5
            % Give warning if code is too slow
            warning('ALS iterations are really slow for some reason - you might want to check what is going on')
        end
    end
    % If solution with highest penalty on log-barrier term is infeasible, increase log-barrier penalty:
    if nu_first == 1
        [r, p] = chol(Q);
        % If solution is feasible, continue normally
        if p == 0
            nu_first = 0;
        else
            % If solution is infeasible, increase nu
            nu = nu * 100;
            iter = 0;
            Q = Q0;
            Qin = inv(Q);
            con = 1;
            con1 = con;
            con2 = con;
            tol1 = tol;
            tol2 = tol;
            phi = 0;
        end
        % Make sure we can't get stuck in this loop
        if nu > 1e10
            error('sdp_QR_mrQ: Unable to find >0 solution');
        end
    end
end

phi2 = obj_ls(Q, Az, b, W); % returns only least-squares portion of objective function
phi_uw = obj_ls(Q, Az, b, eye(size(W))); % returns least-squares portion of objective function without weighting
nsq = -nu * sum(log(eig(Q))); % log-barrier term in final objective value, for debugging purposes (should be ~0)
iter_maxed = 0;
if iter == maxiter
    iter_maxed = 1; % Let user know if maximum number of iterations have been reached
end
