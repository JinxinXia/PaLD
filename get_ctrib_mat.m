function C = get_ctrib_mat(D,beta)
% D is the distance matrix
if beta < 0
    error('beta must be positive');
end

n = size(D,1);
C = zeros(n);

% sort distance matrix by row

for i = 1:n
    D(i,:) = sort(D(i,:),'descend');
end


% use sparse accumulator to get union of ux and uy


for x = 1:(n-1)
    for y = (x+1):n
        ux = find(dx <= beta*D(x,y));  
        uy = find(dy <= beta*D(x,y));
        b = zeros(1,n);
        b(1,ux) = 1;
        b(1,uy) = 1;

        uxy = find(b ~= 0);
        
        
        
        dx = D(x,:);
        dy = D(y,:);
        wx = sum(dx(uxy) < dy(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        wy = sum(dy(uxy) < dx(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        C(x,uxy) = C(x,uxy) + wx/(size(uxy,1)); % fix this by only adding to uxy
        C(y,uxy) = C(y,uxy) + wy/(size(uxy,1));
    end
end
% diag(C) = diag(C)

C = C/(n-1);
end











