
clear;
rng(124);

% create a random distance matrix that is symmetric with diagonal elements
% equal to zeros
d = rand(200);

D = (d+d')/2;
D = D - diag(diag(D));

profile on


disp('their method')

tic
C1 = orig_contribute(D,1);
toc

%disp('sorted method')
%tic
%C2 = sort_mat_contribute(D,1);
%toc
% check distance between the output of two methods

disp('block method')
tic
C3 = mat_block(D,1);
toc 

profile viewer

fprintf('norm(C1-C3): %f \n',norm(C1-C3))

