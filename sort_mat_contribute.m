function [C,suxy] = sort_mat_contribute(D,beta)
% compute coherence matrix C from distance matrix D using pre-sorting

% error checking
if beta < 0
    error('beta must be positive');
end
if ~issymmetric(D) 
    error('distance matrix must be symmetric');
end

n = size(D,1);
C = zeros(n);

% sorted indices matrix 
[sorted_D,sorted_D_indices] = sort(D,2);

% initialize sparse accumulator (SPA)
suxy = zeros(n);

% loop over all pairs x and y (in order of distance from x)
for x = 1:(n-1)
    for j = 1:n  
        
        % initialize sparse accumulator (SPA)
        uxy = false(1,n);
        
        %% Determine ux and uy using sorted rows
        % let ux be set of points within D(x,y) of x
        % first j entries in sorted_D_indices(x,:) are in ux
        y = sorted_D_indices(x,j);
        % x and y are two points we consider in this iteration
        
        % perform update only for x < y case
        if y <= x
            continue;
        end
        
        % let uy be set of points within D(x,y) of y
        % find index of D(x,y) in sorted row of y        
        k = binary_search(sorted_D(y,:),D(x,y));
        % first k entries in sorted_D_indices(y,:) are in uy
        
        %% Compute size_uxy = |ux U uy| using SPA
        % set the values of SPA to 1 for entries in ux
        uxy(sorted_D_indices(x,1:j)) = true;
        uxy(sorted_D_indices(y,1:k)) = true;
        size_uxy = nnz(uxy);
        
        % calculate local depth
        zx = uxy & D(x,:) < D(y,:); % z's closer to x
        zz = uxy & D(x,:) == D(y,:); % z's equidistant
        zy = uxy & D(x,:) > D(y,:); % z's closer to y
        
        % assign local depth value to the corresponding position in C
        C(x,zx) = C(x,zx) + 1/size_uxy; 
        C(x,zz) = C(x,zz) + .5/size_uxy;
        C(y,zy) = C(y,zy) + 1/size_uxy;
        C(y,zz) = C(y,zz) + .5/size_uxy;
            
        suxy(x,y) = size_uxy;
    end
end

C = C/(n-1);



end


function ind = binary_search(sorted_array,val)
% search sorted_array for val (should be in there)

    % initialize indices
    left = 1;
    right = length(sorted_array);
    ind = floor((left+right)/2);
    
    % loop until value is found
    while sorted_array(ind) ~= val
        
        % value not in array, report error
        if left == ind
            if sorted_array(right) == val
                ind = right;
                return;
            else
                error('not found')
            end
        end
        
        % perform binary search
        if sorted_array(ind) < val
            % search right
            left = ind;
        else
            % search left
            right = ind;
        end
        ind = floor((left+right)/2);
        
    end
            
end
