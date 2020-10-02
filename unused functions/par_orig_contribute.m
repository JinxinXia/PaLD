function C = par_orig_contribute(D,beta)
% D is the distance matrix, the value of D(x,y) is the distance between x 
% and y (d(x,y)).
% beta is the parameter to control the radius of conflict focus, the output
% matrix C is the cohension matrix. 

% this is the paralle version of the code, need to fix parfor

if beta < 0
    error('beta must be positive');
end

if D' ~= D 
    error('distance matrix must be symmetric');
end


n = size(D,1);
C = zeros(n);


parfor x = 1:(n-1)
    for y = (x+1):n
     
        dx = D(x,:); % get the row of distance between x and all other points 
        dy = D(y,:); % get the row of distance between y and all other points   
       
        % get all the unique indices(or points) to form conflict focus
        uxy = (dx <= beta*D(x,y) | dy <= beta*D(x,y));
       
        % calculate local depth
        zx = dx < dy; % z's closer to x
        zz = dx == dy; % z's equidistant
        zy = dx > dy; % z's closer to y
        u_size = sum(uxy);
        
        % assign local depth value to the corresponding position in C
        C(x,zx & uxy) = C(x,zx & uxy) + 1/u_size; 
        C(x,zz & uxy) = C(x,zz & uxy) + .5/u_size;
        C(y,zy & uxy) = C(y,zy & uxy) + 1/u_size;
        C(y,zz & uxy) = C(y,zz & uxy) + .5/u_size;
      
    end
end

% convert summed local depth to cohesion matrix
C = C/(n-1);

end











