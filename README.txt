########
This is the updated version of algorithm PaLD from Berenhaut, K. S., Moore, K. E. and Melvin, R. L.** (2021) A social perspective on perceived distances reveals deep community structure. Proceedings of the National Academy of Sciences, in press.

The original algorithm has an order of O({n^3}/2). The new algorithm beats it with O(n^2 + \sum{|u_xy|}), where \sum{|u_xy|} is the sum of each size of local conflict focus. 
