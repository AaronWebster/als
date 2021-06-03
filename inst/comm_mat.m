% Creates mn x mn commutation matrix Kmn so that for any m x n matrix A, 
% Kmn*vec(A) = vec(A')
function Kmn = comm_mat(m,n)
Kmn = zeros(m*n,m*n);
for i = 1:1:m
  for j = 1:1:n
    Kmn((i-1)*n + j, (j-1)*m + i) = 1;
  end
end
