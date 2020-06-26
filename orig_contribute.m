function C = orig_contribute(D,beta)
% D is the distance matrix, the value of D(x,y) is the distance between x 
% and y (d(x,y)).
% beta is the parameter to control the radius of conflict focus, the output
% matrix C is the cohension matrix. 

if beta < 0
    error('beta must be positive');
end

if D' ~= D 
    error('distance matrix must be symmetric');
end


n = size(D,1);
C = zeros(n);


for x = 1:(n-1)
    for y = (x+1):n
     
        dx = D(x,:); % get the row of distance between x and all other points 
        dy = D(y,:); % get the row of distance between y and all other points
        
        % store all the indices whose element values is smaller than d(x,y)
        % either in dx or dy
        ux = find(dx <= beta*D(x,y)); 
        uy = find(dy <= beta*D(x,y));
        
        % set a zero vector b and put one to the corresponding indices from
        % ux and uy, we can get all the unique indices from this step
        b = zeros(1,n);
        b(1,ux) = 1;
        b(1,uy) = 1;
        
        % get all the unique indices(or points) to form conflict focus
        uxy = find(b ~= 0);
       s
        % calculate local depth
        wx = sum(dx(uxy) < dy(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        wy = sum(dy(uxy) < dx(uxy)) + 0.5*sum(dy(uxy) == dx(uxy));
        u_size = size(uxy,2);
        
        % assign local depth value to the corresponding position in C
        if u_size ~= 0
            C(x,uxy) = C(x,uxy) + wx/u_size; 
            C(y,uxy) = C(y,uxy) + wy/u_size;
        end
      
    end
end

% convert summed local depth to cohesion matrix
C = C/(n-1);

end











