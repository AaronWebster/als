% Function to "undo" the symmetry imposing of the duplication matrix
% tran*(x)s = (x)ss
% Equivalent to pinv(symtran(n)), but much faster
function tran = symtranT(n)
tran = symtran(n)';
for i = 1:1:size(tran, 1)
    tran(i, :) = tran(i, :) / sum(tran(i, :));
end
end
