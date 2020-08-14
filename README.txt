########
This is the updated version of method PaLD

The original algorithm has an order of O({n^3}/2) because of the symmetric property of the distance matrix. We are trying to beat it using sorted algorithm with O(n^2 + \sum{|u_xy|}), where \sum{|u_xy|} is the sum of each size of local conflict focus. 

