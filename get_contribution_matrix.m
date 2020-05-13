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
        uxy = [find(dx <= beta*D(x,y))  find(dy <= beta*D(x,y))];
        uxy = unique(uxy);
        wx = sum(dx(uxy) < dy(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        wy = sum(dy(uxy) < dx(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        C(x,uxy) = C(x,uxy) + wx/(size(uxy,1)); % fix this by only adding to uxy
        C(y,uxy) = C(y,uxy) + wy/(size(uxy,1));
        
    end
end
% diag(C) = diag(C)




C = C/(n-1);
end
        
  
