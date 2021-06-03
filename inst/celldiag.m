function A = celldiag(F, n)
% Function for diagonalizing matrices
% function A = celldiag(F,n)
% Function inputs:
% F - cell containing submatrices
% n - number of rows of zeros to add
% if n < 0, zeros are added to upper right
% if n > 0, zeros are added to lower left
% Function outputs:
% A - block diagonal matrix containing elements of F
% (augmented with zeros if n specified)
% A = [0 0; blockdiag(F) 0] if n < 0
% A = [0 blockdiag(F); 0 0] if n > 0
if (nargin == 1) n = 0; end
len = prod(size(F));
A = F{1};
[x, y] = size(A);
for i = 2:len
    [a, b] = size(F{i});
    A = [A, zeros(x, b); zeros(a, y), F{i}];
    x = x + a;
    y = y + b;
end
up = (n < 0);
dn = (n > 0);
n = abs(n);
A = [zeros(up * n, y + n); zeros(x, dn * n), A, zeros(x, up * n); zeros(dn * n, y + n)];
