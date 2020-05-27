function sort_D_indices = get_sorted_indices(D)
% this function will return the indices after sorting a square matrix 

n = size(D,1);

% sort the upper triangular of distance matrix D on each row by a ascending
% order
sort_D = sort(D,2);

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

end 