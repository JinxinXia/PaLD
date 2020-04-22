function C = get_contribution_matrix(D,beta)
% D is the distance matrix
if beta < 0
    error('beta must be positive');
end


n = size(D,1);
C = zeros(n);

for x = 1:(n-1)
    for y = (x+1):n
        dx = D(x,:);
        dy = D(y,:);
        uxy = [dx( (dx <= beta*D(x,y))  dy( (dy <= beta*D(x,y))];
        wx = d
