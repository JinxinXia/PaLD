function C = get_ctrib_mat(D,beta)
% D is the distance matrix
if beta < 0
    error('beta must be positive');
end

if D' ~= D 
    error('distance matrix must be symmetric');
end


n = size(D,1);
C = zeros(n);


% use sparse accumulator to get union of ux and uy
% need to sort each row of the distance matrix

for x = 1:(n-1)
    for y = (x+1):n
        dx = D(x,:);
        dy = D(y,:);
        ux = find(dx <= beta*D(x,y));  
        uy = find(dy <= beta*D(x,y));
        b = zeros(1,n);
        b(1,ux) = 1;
        b(1,uy) = 1;

        uxy = find(b ~= 0);
        
        wx = sum(dx(uxy) < dy(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        wy = sum(dy(uxy) < dx(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        u_size = size(uxy,2);
        if u_size ~= 0
            % check for zeros
            C(x,uxy) = C(x,uxy) + wx/u_size; % fix this by only adding to uxy
            C(y,uxy) = C(y,uxy) + wy/u_size;
        end
      
    end
end

C = C/(n-1);

end











