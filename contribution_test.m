
clear;
rng(122);

% create a random distance matrix that is symmetric with diagonal elements
% equal to zeros
d = rand(100);

D = (d+d')/2;
D = D - diag(diag(D));


%profile on


disp('their method')

tic
C1 = orig_contribute(D,1);
toc

disp('sorted method')
tic
C2 = sort_mat_contribute(D,1);
toc
% check distance between the output of two methods

%profile viewer

fprintf('norm(C1-C2): %f \n',norm(C1-C2))
