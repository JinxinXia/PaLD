function C = get_ctrib_mat(D,beta)
% D is the distance matrix
if beta < 0
    error('beta must be positive');
end

n = size(D,1);
C = zeros(n);

ux = find(dx <= beta*D(x,y));  
uy = find(dy <= beta*D(x,y));



