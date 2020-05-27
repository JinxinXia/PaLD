function C = get_cmat(D,beta)


% D is the distance matrix
if beta < 0
    error('beta must be positive');
end

if D' ~= D 
    error('distance matrix must be symmetric');
end

n = size(D,1);
C = zeros(n);

% sort the upper triangular of distance matrix D on each row by a ascending
% order
sort_D = sort(triu(D),2);

sort_D_indices = zeros(n);
% get the sorted indices matrix 
for i = 1:n
    for j = 1:n
        % get the original index in distance matrix
        index = find(D(i,:) == sort_D(i,j));
        index_size = size(index,2);
        
        % check if there are same distance between pairs of points
        if index_size == 1
            sort_D_indices(i,j) = index;
        else
            sort_D_indices(i,j:j+index_size) = index;
        end
        
        % skip index_size of loops
        if index_size > 0
            index_size = index_size - 1;
            continue;
        end
        
    end
end





% use sparse accumulator to get union of ux and uy
% need to sort each row of the distance matrix

for x = 1:(n-1)
    for y = (x+1):n
        
        [~,uxy_size] = find(sort_D(x,:) == D(x,y));
        D(x,y)
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
            C(x,uxy) = C(x,uxy) + wx/(size(uxy,1)); % fix this by only adding to uxy
            C(y,uxy) = C(y,uxy) + wy/(size(uxy,1));
        end
      
    end
end

C = C/(n-1);



end