function [C,U] = optimal_matmul_block(D,beta,b)
% run optimial block algorithm on PaLD method
% b <= sqrt(M+1)-1, where M is the cache size


% error checking
if beta < 0
    error('beta must be positive');
end

if ~issymmetric(D) 
    error('distance matrix must be symmetric');
end

n = size(D,1);
C = zeros(n);
U = zeros(n);

% use optimal block (k=1) to do comparsion 
% take advantage of the vectorize operation in MATLAB
% same logic as optimal_blocked.m in matmul-test folder


% loop over blocks to update conflict focus sizes
for i = 1:b:n
    Ib = i:min(i+b-1,n);
    for j = i:b:n
        Jb = j:min(j+b-1,n);
        % update size of conflict focus for block I,J of U
        Uij = zeros(size(U(Ib,Jb))); % temp buffer
        Dij = D(Ib,Jb); % temp buffer
        for k=1:n
            % outer-product-like update
            Uij = Uij + double(min(D(Ib,k),D(k,Jb)) < Dij);
        end
        U(Ib,Jb) = Uij; % copy temp buffer back to U
    end
end
 
%U = U+U';


% loop over blocks to update cohesion matrix
for i = 1:b:n
    Ib = i:min(i+b-1,n);
    for j = i:b:n
        Jb = j:min(j+b-1,n);
        
        % compute size of conflict focus for block I,J
        Uij = zeros(min(b,n-i+1),min(b,n-j+1)); % conflict focus sizes
        Dij = D(Ib,Jb); % temp buffer
        for k=1:n
            % outer-product-like update
            Uij = Uij + double(min(D(Ib,k),D(k,Jb)) < Dij);
        end
        
        % update cohesion values according to conflict foci (I,J)
        % update this function call to another k-loop with outer product
        % like computations (updates will be to subcolumns of C corr. to
        % C(Ib,k) or C(Jb,k)
        [C(Ib,:),C(Jb,:)] = update_coh(C(Ib,:),C(Jb,:),D(Ib,Jb),D(Ib,:),D(:,Jb),Uij);
    end
end


C = C/(n-1);


end




function Uij = update_cfs(Uij,Dij,Dik,Dkj,diagbool)
% update conflict focus sizes for pairs in I,J block-pair 
% based on vertices in K block

    [~,p] = size(Dik);
    
    % make k the outer loop, then outer product like computation 
    % updates the ij block
    if diagbool
        for k=1:p
            Uij = Uij + double(triu(min(Dik(:,k),Dkj(k,:)) < Dij));
        end
    else
        for k=1:p
            Uij = Uij + double(min(Dik(:,k),Dkj(k,:)) < Dij);
        end
    end
    

end

function [Cik,Cjk] = update_coh(Cik,Cjk,Dij,Dik,Dkj,Uij)
% update cohesion values based on (I,J) conflict foci
% based on vertices in K block that are in focus

    [m,n] = size(Uij);
    [~,p] = size(Dik);
    
    for k=1:p
        for i = 1:m
            for j=1:n
                % skip case where point i is point j
                if Dij(i,j) == 0
                    continue
                end
            
                % determine if point k is in (i,j)'s focus
                if min(Dik(i,k),Dkj(k,j)) <= Dij(i,j)
                    % determine where k's contribution belongs
                    if Dik(i,k) < Dkj(k,j)
                        Cik(i,k) = Cik(i,k) + 1 / Uij(i,j);
                    elseif Dik(i,k) > Dkj(k,j)
                        Cjk(j,k) = Cjk(j,k) + 1 / Uij(i,j);
                    else
                        Cik(i,k) = Cik(i,k) + .5 / Uij(i,j);
                        Cjk(j,k) = Cjk(j,k) + .5 / Uij(i,j);
                    end
                end
            end
        end
        
    end








end