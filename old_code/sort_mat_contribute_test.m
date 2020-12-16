function C = sort_mat_contribute_test(D,beta)
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
b = zeros(1,n);

% loop over all pairs x and y (in order of distance from x)
for x = 1:(n-1)
    for j = 1:n        
        
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
        b(sorted_D_indices(x,1:j)) = 1;
        size_uxy = j;
        % loop through entries in uy and count if not in ux
        for i = 1:k
            if b(sorted_D_indices(y,i)) == 0
                b(sorted_D_indices(y,i)) = 1;
                size_uxy = size_uxy+1;
            end
        end
        
        %% Accumulate coherence values for every entry in uxy using SPA
        % loop through ux first
        for i = 1:j
            z = sorted_D_indices(x,i);
            if D(x,z) < D(y,z)
                C(x,z) = C(x,z) + 1 / size_uxy;
            elseif D(x,z) > D(y,z)
                C(y,z) = C(y,z) + 1 / size_uxy;
            else
                C(x,z) = C(x,z) + 0.5 / size_uxy;
                C(y,z) = C(y,z) + 0.5 / size_uxy;
            end
            b(z) = 0;
        end
        % loop through uy next
        for i = 1:k
            z = sorted_D_indices(y,i);
            % only process if not in ux
            if b(z) == 1
                if D(x,z) < D(y,z)
                    C(x,z) = C(x,z) + 1 / size_uxy;
                elseif D(x,z) > D(y,z)
                    C(y,z) = C(y,z) + 1 / size_uxy;
                else
                    C(x,z) = C(x,z) + 0.5 / size_uxy;
                    C(y,z) = C(y,z) + 0.5 / size_uxy;
                end
                b(z) = 0;
            end 
        end

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
