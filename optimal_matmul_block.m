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
    diagbool = 1;
    for j = i:b:n
        if j > i
            diagbool = 0;
        end
        Jb = j:min(j+b-1,n);
        % update size of conflict focus for block I,J of U
        U(Ib,Jb) = update_cfs(U(Ib,Jb),D(Ib,Jb),D(Ib,:),D(:,Jb),diagbool);
    end
end
 
U = U+U';


% loop over blocks to update cohesion matrix
for i = 1:b:n
    Ib = i:min(i+b-1,n);
    for j = i:b:n
        Jb = j:min(j+b-1,n);
        % update cohesion values according to conflict foci (I,J)
        [C(Ib,:),C(Jb,:)] = update_coh(C(Ib,:),C(Jb,:),D(Ib,Jb),D(Ib,:),D(:,Jb),U(Ib,Jb));
    end
end


C = C/(n-1);


end




function Uij = update_cfs(Uij,Dij,Dik,Dkj,diagbool)
% update conflict focus sizes for pairs in I,J block-pair 
% based on vertices in K block

    [m,n] = size(Uij);
    [~,p] = size(Dik);
    
    for i = 1:m
        % if working on diagonal block, compute only upper triangle
        if diagbool == 1
            startj = i+1;
        else
            startj = 1;
        end
        for j=startj:n
            for k=1:p
                % determine if point k is in (i,j)'s focus
                if min(Dik(i,k),Dkj(k,j)) < Dij(i,j)
                    Uij(i,j) = Uij(i,j) + 1;
                end
            end
        end
    end
    

end

function [Cik,Cjk] = update_coh(Cik,Cjk,Dij,Dik,Dkj,Uij)
% update cohesion values based on (I,J) conflict foci
% based on vertices in K block that are in focus

    [m,n] = size(Uij);
    [~,p] = size(Dik);
    
    for i = 1:m
        for j=1:n
            % skip case where point i is point j
            if Dij(i,j) == 0
                continue
            end
            for k=1:p
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