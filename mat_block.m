
function [C,U] = mat_block(D,beta,b)

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
U = zeros(n);
U2 = zeros(n);

% loop over all pairs x and y to get size of local conflict foci

% U checked
for x = 1:n
    for y = x+1:n
        U2(x,y) = sum(min(D(x,:),D(y,:)) <= D(x,y));
    end
end


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
        for k = 1:b:n
            Kb = k:min(k+b-1,n);
            U(Ib,Jb) = update_cfs(U(Ib,Jb),D(Ib,Jb),D(Ib,Kb),D(Kb,Jb),diagbool);
        end
    end
end
 
assert(norm(U - U2) == 0);
U = U+U';

% loop over all blocks to update contribution matrix C
% consider three sets of points: I, J, and K
%  - consider conflict focus on (I,J)
%  - determine who K points contribute to
%  - update blocks C(I,K) and C(J,K)
%  - need access to D(I,J), D(I,K), D(K,J) and U(I,J)

% loop over blocks to update cohesion matrix
for i = 1:b:n
    Ib = i:min(i+b-1,n);
    for j = i:b:n
        Jb = j:min(j+b-1,n);
        % update cohesion values according to conflict foci (I,J)
        for k = 1:b:n
            Kb = k:min(k+b-1,n);
            [C(Ib,Kb),C(Jb,Kb)] = update_coh(C(Ib,Kb),C(Jb,Kb),D(Ib,Jb),D(Ib,Kb),D(Kb,Jb),U(Ib,Jb));
        end
    end
end

% block things messed up
% for x = 1:b:n
%     for y = 1:b:n
%         for z = 1:b:n
%             xend = min(x+b-1,n);
%             zend = min(z+b-1,n);
%             yend = min(y+b-1,n);
%             
%             if x == y && y == z
%                 continue
%             end
%             
%             % d(x,z) < d(x,y)  && d(x,z) < d(y,z), thus z contribute to x
%             C(x:xend,z:zend) = C(x:xend,z:zend) + (D(x:xend,z:zend)...
%                 < D(x:xend,y:yend) & D(x:xend,z:zend) < D(z:zend,y:yend))./U(x:xend,y:yend);
%             
%             % d(y,z) < d(x,y)  && d(y,z) < d(x,z), thus z contribute to y
%             C(y:yend,z:zend) = C(y:yend,z:zend) + (D(y:yend,z:zend)...
%                 < D(x:xend,y:yend) & D(y:yend,z:zend) < D(z:zend,x:xend))./U(x:xend,y:yend);
%         end
%     end
% end

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
