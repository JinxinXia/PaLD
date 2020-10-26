function C = optimal_blocked(A,B,b)
% optimal matmul

C = zeros(size(A,1),size(B,2));
          
for i = 1:b:size(C,1)
   for j = 1:b:size(C,2)
       iend = min(i+b-1,size(C,1));
       jend = min(j+b-1,size(C,2));
       C(i:iend,j:jend) = C(i:iend,j:jend) + A(i:iend,:)*B(:,j:jend);
    end
end


end
