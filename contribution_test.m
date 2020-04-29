clear;
T = readtable('distance.csv');

%distance = textscan(a, '%s %s %s %s', 'delimiter', ',', 'CollectOutput',true);
cell_dist = T{:,:};
b = cell_dist(2:end,2:end);
%distance = cell2mat(b);

class(b)
size(b)

C = get_contribution_matrix(b,1);