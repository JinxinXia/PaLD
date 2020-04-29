
T = readtable('distance.csv');
cell_dist = T{:,:};
size(cell_dist)
%distance = cell2mat(1,cell_dist{:});

%class(distance)
%size(distance)