*version
1.0
*title
4x4 square grid
*dimensions
25	# number of nodes
2	# number of spatial dimensions
1	# number of element sets
#  [ID]  [nel]  [nen]
1  16  4
4  # number of node sets
#  [ID]  [nnd]
1	5
2	5
3	5
4	5
2  # number of side sets
#  [ID]  [associate block ID]  [ns]
1	1	4
2	1	4
# end dimensions
*nodesets
*set
5  # number of nodes
1  6  11  16  21  
*set
5  # number of nodes
21  22  23  24  25  
*set
5  # number of nodes
5  10  15  20  25  
*set
5  # number of nodes
1  2  3  4  5  
# end node sets
*sidesets
*set
4
1	4
5	4
9	4
13	4
*set
4
4	2
8	2
12	2
16	2
*elements
*set
16  # number of elements
4   # number of element nodes
	1	1	6	7	2
	2	6	11	12	7
	3	11	16	17	12
	4	16	21	22	17
	5	2	7	8	3
	6	7	12	13	8
	7	12	17	18	13
	8	17	22	23	18
	9	3	8	9	4
	10	8	13	14	9
	11	13	18	19	14
	12	18	23	24	19
	13	4	9	10	5
	14	9	14	15	10
	15	14	19	20	15
	16	19	24	25	20
# end elements
*nodes
25  # number of nodes
2   # number of spatial dimensions
	1	-5.000000e-01	-5.000000e-01
	2	-5.000000e-01	-2.500000e-01
	3	-5.000000e-01	0
	4	-5.000000e-01	2.500000e-01
	5	-5.000000e-01	5.000000e-01
	6	-2.500000e-01	-5.000000e-01
	7	-2.500000e-01	-2.500000e-01
	8	-2.500000e-01	0
	9	-2.500000e-01	2.500000e-01
	10	-2.500000e-01	5.000000e-01
	11	0	-5.000000e-01
	12	0	-2.500000e-01
	13	0	0
	14	0	2.500000e-01
	15	0	5.000000e-01
	16	2.500000e-01	-5.000000e-01
	17	2.500000e-01	-2.500000e-01
	18	2.500000e-01	0
	19	2.500000e-01	2.500000e-01
	20	2.500000e-01	5.000000e-01
	21	5.000000e-01	-5.000000e-01
	22	5.000000e-01	-2.500000e-01
	23	5.000000e-01	0
	24	5.000000e-01	2.500000e-01
	25	5.000000e-01	5.000000e-01
