%% Function to give the duplication matrix that will impose symmetry constraints on

%% the covariance matrices, (Q)_s = tran*(Q)_ss
function tran = symtran(n)

r = n * (n + 1) / 2;

tran = zeros(n^2, r);
k = 1;
for i = 1:n
    for j = i:n
        temp = zeros(n);
        temp(i, j) = 1;
        if (i == j) div = 2;
        else div = 1;
        end
        t2 = (temp + temp') / div;
        tran(:, k) = t2(:);
        k = k + 1;
    end
end
