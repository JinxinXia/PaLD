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
    for y = x:b:n
        Yb = y:min(y+b-1,n);
        for z = y:b:n
            Zb = y:min(z+b-1,n);
             
            if x == y && y == z
                index = 1;
            elseif x == y
                index = 2;
            elseif x ~= y && y ~= z
                index = 3;
            else
                continue;
            end
            
            
            [U(Xb,Yb),U(Xb,Zb),U(Yb,Zb)] = ...
            update_cfs(U(Xb,Yb),U(Xb,Zb),U(Yb,Zb),D(Xb,Yb),D(Xb,Zb),D(Yb,Zb),index);
        end
    end
end



% fill in lower triangle of U
% U = U + U';

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




function [Uxy,Uxz,Uyz] = update_cfs(Uxy,Uxz,Uyz,Dxy,Dxz,Dyz,index)
% update conflict focus sizes for triplets
% if x == y == z 
%    index = 1;
% elseif x == y
%    index = 2;
% elseif x~=y & y~=z
%    index = 3;
% end
 
    [m,n] = size(Dxy); % m is the size of X block, n is the size of Y block
    [~,p] = size(Dxz); % p is the size of Z block
    
    if index == 1
        for x = 1:m
            for y=(x+1):n
                for z=(y+1):p
                    
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
    elseif index == 2
        %{
        
        % two points from x and y, one point from z
        for x = 1:m
            for y=(x+1):n
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
        
         % two points from z and one point from x or y
        
       for i = 1:p % i is the first point select from Z
          for j = (z+1):p  % j is the second point select from Z
               for k = 1:m % k is the first point select from X(Y)
                    
                    if Dxz(k,i) < Dxy(k,j) && Dxz(k,i) < Dyz(k,j) % ki is closest
                        Uxy(k,j) = Uxy(k,j) + 1;
                        Uyz(i,j) = Uyz(j,i) + 1;
                    elseif Dxy(k,j) < Dxz(k,i) && Dxy(k,j) < Dyz(k,j) % kj is closest
                        Uxy(k,i) = Uxy(k,j) + 1;
                        Uyz(i,j) = Uyz(j,i) + 1;
                    else                                      % ij is closest 
                        Uxy(k,i) = Uxy(k,i) + 1;
                        Uyz(k,j) = Uyz(k,j) + 1;
                   
                    end
                    
                end
            end
        end
       
        %}
    
    elseif index == 3
        for x = 1:m
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
        
    
        
end
