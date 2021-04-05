function [C,U] = pald_triplet_block(D,b)


if D' ~= D 
    error('distance matrix must be symmetric');
end

n = size(D,1);
U = triu(2*ones(n),1); % init each conflict focus size to include 2 points
C = zeros(n);



% consider all unique triplet blocks to compute conflict focus sizes
for x = 1:b:n
    Xb = x:min(x+b-1,n);
    
    
    for y = x+b:b:n
        Yb = y:min(y+b-1,n);
       
        % update size of conflict focus for block I,J of U
        for z = y+b:b:n
            Zb = z:min(z+b-1,n);
              
            [U(Xb,Yb),U(Xb,Zb),U(Yb,Zb)] = ...
            update_cfs(U(Xb,Yb),U(Xb,Zb),U(Yb,Zb),D(Xb,Yb),D(Xb,Zb),D(Yb,Zb));
        end
    end
end



% fill in lower triangle of U
U = U + U';

% initialize C with diagonal entries determined by row sums of 1/U
C = diag(sum(1./(U+diag(Inf*ones(n,1))),2)); 

% consider all unique triplets to compute contributions to cohesion
for x = 1:(n-1)
    for y = (x+1):n
        for z = (y+1):n
            if D(x,y) < D(x,z) && D(x,y) < D(y,z)     % xy is closest pair
                C(x,y) = C(x,y) + 1 / U(x,z);
                C(y,x) = C(y,x) + 1 / U(y,z);
            elseif D(x,z) < D(x,y) && D(x,z) < D(y,z) % xz is closest pair
                C(x,z) = C(x,z) + 1 / U(x,y);
                C(z,x) = C(z,x) + 1 / U(y,z);
            else                                      % yz is closest pair
                C(y,z) = C(y,z) + 1 / U(x,y);
                C(z,y) = C(z,y) + 1 / U(x,z);
            end
        end
    end
end      

C = C/(n-1);
end




function [Uxy,Uxz,Uyz] = update_cfs(Uxy,Uxz,Uyz,Dxy,Dxz,Dyz)
% update conflict focus sizes for pairs in I,J block-pair 
% based on vertices in K block
    
    [m,n] = size(Dxy); % m is the size of X block, n is the size of Y block
    [~,p] = size(Dxz); % p is the size of Z block
    
    
    for x = 1:m
        % since x, y, z are different there won't be equal situation
        for y=1:n
            
            for z=1:p
                if Dxy(x,y) < Dxz(x,z) && Dxy(x,y) < Dyz(y,z)     % xy is closest pair
                    Uxz(x,z) = Uxz(x,z) + 1;
                    Uyz(y,z) = Uyz(y,z) + 1;
                elseif Dxz(x,z) < Dxy(x,y) && Dxz(x,z) < Dyz(y,z) % xz is closest pair
                    Uxy(x,y) = Uxy(x,y) + 1;
                    Uyz(y,z) = Uyz(y,z) + 1;
                else                                      % yz is closest pair
                    Uxy(x,y) = Uxy(x,y) + 1;
                    Uxz(x,z) = Uxz(x,z) + 1;
                end
            end
        end
    end
    

end
