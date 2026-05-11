*version
1.0
*title
Single hex dielectric elastomer 3D test
*dimensions
8   # number of nodes
3   # number of spatial dimensions
1   # number of element sets
#  [ID]  [nel]  [nen]
1  1  8
5   # number of node sets
#  [ID]  [nnd]
1   4
2   4
3   1
4   2
5   1
0  # number of side sets
# end dimensions
*nodesets
*set
4  # number of nodes
1  2  3  4
*set
4  # number of nodes
5  6  7  8
*set
1  # number of nodes
1
*set
2  # number of nodes
2  6
*set
1  # number of nodes
3
# end node sets
*sidesets
# end side sets
*elements
*set
1  # number of elements
8  # number of element nodes
1  1  2  4  3  5  6  8  7
# end elements
*nodes
8  # number of nodes
3   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00   0.0000000e+00
    2   1.0000000e+00   0.0000000e+00   0.0000000e+00
    3   0.0000000e+00   1.0000000e+00   0.0000000e+00
    4   1.0000000e+00   1.0000000e+00   0.0000000e+00
    5   0.0000000e+00   0.0000000e+00   1.0000000e+00
    6   1.0000000e+00   0.0000000e+00   1.0000000e+00
    7   0.0000000e+00   1.0000000e+00   1.0000000e+00
    8   1.0000000e+00   1.0000000e+00   1.0000000e+00
