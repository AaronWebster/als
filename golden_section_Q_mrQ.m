function h = golden_section_Q_mrQ(x00, x33, Q, delQ, A, b, W, nu, vlam_mtr)
% Uses golden section method to calculate optimal step size for given search direction
% Objective function: y =  norm(A*Q(:) - b)^2_w + trace(Q*vlam_mtr) - nu*log(det(Q));

tol = 1e-3;
x0 = x00;
x3 = x33;
alpha = (3 - sqrt(5)) / 2;
x1 = x0 + alpha * (x3 - x0);
x2 = x3 - alpha * (x3 - x0);

Q0 = Q + x0 * delQ;
J0 = obj_tot(Q0, A, b, W, nu, vlam_mtr);
Q1 = Q + x1 * delQ;
J1 = obj_tot(Q1, A, b, W, nu, vlam_mtr);
Q2 = Q + x2 * delQ;
J2 = obj_tot(Q2, A, b, W, nu, vlam_mtr);
Q3 = Q + x3 * delQ;
J3 = obj_tot(Q3, A, b, W, nu, vlam_mtr);
while (norm(x3 - x0)) > tol
    if (J2 < J1)
        x0 = x1;
        J0 = J1;
        x1 = x2;
        J1 = J2;
        x2 = x3 - alpha * (x3 - x0);
        Q2 = Q + x2 * delQ;
        J2 = obj_tot(Q2, A, b, W, nu, vlam_mtr);
    else
        x3 = x2;
        J3 = J2;
        x2 = x1;
        J2 = J1;
        x1 = x0 + alpha * (x3 - x0);
        Q1 = Q + x1 * delQ;
        J1 = obj_tot(Q1, A, b, W, nu, vlam_mtr);
    end
end

if (J1 < J2)
    h = x1;
else
    h = x2;
end
