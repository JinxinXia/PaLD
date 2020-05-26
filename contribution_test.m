
clear;
rng(12);

% create a random distance matrix that is symmetric with diagonal elements
% equal to zeros
d = rand(500,500);
D = (d+d')/2;
D = D - diag(diag(D));

disp('our method')

tic
C1 = get_ctrib_mat(D,1);
toc

disp('their method')
tic
C = get_contribution_matrix(D,1);
toc

% check distance between the output of two methods
fprintf('norm(C1-C): %f \n',norm(C-C1))
