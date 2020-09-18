
function [C,U] = mat_block(D,beta)

% compute coherence matrix C from distance matrix D using pre-sorting

% error checking
if beta < 0
    error('beta must be positive');
end
if ~issymmetric(D) 
    error('distance matrix must be symmetric');
end

n = size(D,1);
b = n/2;
C = zeros(n);
U = zeros(n);

% loop over all pairs x and y to get size of local conflict foci

% U checked
for x = 1:n
    for y = x+1:n
        U(x,y) = sum(min(D(x,:),D(y,:)) <= D(x,y));
    end
end


U = U + U' + eye(n)*10^17;


% loop over all blocks to update contribution matrix C


% block things messed up
for x = 1:b:n
    for y = 1:b:n
        for z = 1:b:n
            xend = min(x+b-1,n);
            zend = min(z+b-1,n);
            yend = min(y+b-1,n);
            
            if x == y && y == z
                continue
            end
            
            % d(x,z) < d(x,y)  && d(x,z) < d(y,z), thus z contribute to x
            C(x:xend,z:zend) = C(x:xend,z:zend) + (D(x:xend,z:zend)...
                < D(x:xend,y:yend) & D(x:xend,z:zend) < D(z:zend,y:yend))./U(x:xend,y:yend);
            
            % d(y,z) < d(x,y)  && d(y,z) < d(x,z), thus z contribute to y
            C(y:yend,z:zend) = C(y:yend,z:zend) + (D(y:yend,z:zend)...
                < D(x:xend,y:yend) & D(y:yend,z:zend) < D(z:zend,x:xend))./U(x:xend,y:yend);
        end
    end
end



C = C/(2*(n-1));


end


