function C = sort_mat_contribute_eff(D,beta)


% D is the distance matrix
if beta < 0
    error('beta must be positive');
end

if ~issymmetric(D) 
    error('distance matrix must be symmetric');
end

n = size(D,1);
C = zeros(n);

% get the sorted indices matrix 
[sorted_D,sorted_D_indices] = sort(D,2);

% initialize sparse accumulator (SPA)
b = zeros(1,n);

% use sparse accumulator to get union of ux and uy
% need to sort each row of the distance matrix

for x = 1:(n-1)
    for j = 1:n
        
        
        
        % x and y are two points that we selected
        y = sorted_D_indices(x,j);
        
        % perform update only for x < y case
        if y <= x
            continue;
        end
        
        % ux contains the points that have a smaller distance with x
        % comparing with d(x,y)
        %ux = sorted_D_indices(x,1:j);
        
        
        
        % uy contains the points that have a smaller distance with y
        % comparing with d(x,y), to 
        % get uy we need find d(x,y) in the yth
        % row, k give the index where d(x,y) lives in the sorted
        % distance indices matrix
        
        k = binary_search(sorted_D(y,:),D(x,y));
        %uy = sorted_D_indices(y,1:k);
        
        %% Compute size_uxy = |ux U uy|
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
        
        %% Accumulate values for every entry in uxy
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
% search sorted_array 

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
