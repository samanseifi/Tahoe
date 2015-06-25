#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

const long double PI = 3.1415927;
const int NODE  = 134;
const int ELEM  = 144;
const int CONNECT[ELEM][4] = {
  {1	,	1	,	2	,	3	},
  {1	,	1	,	3	,	4	},
  {1	,	1	,	4	,	5	},
  {1	,	1	,	5	,	6	},
  {1	,	1	,	6	,	7	},
  {1	,	1	,	7	,	8	},
  {1	,	1	,	8	,	9	},
  {1	,	1	,	9	,	10	},
  {1	,	1	,	10	,	11	},
  {1	,	1	,	11	,	12	},
  {1	,	1	,	12	,	13	},
  {1	,	1	,	13	,	2	},
  {3	,	2	,	14	,	15	},
  {4	,	3	,	15	,	16	},
  {5	,	4	,	16	,	17	},
  {6	,	5	,	17	,	18	},
  {7	,	6	,	18	,	19	},
  {8	,	7	,	19	,	20	},
  {9	,	8	,	20	,	21	},
  {10	,	9	,	21	,	22	},
  {11	,	10	,	22	,	23	},
  {12	,	11	,	23	,	24	},
  {13	,	12	,	24	,	25	},
  {2	,	13	,	25	,	14	},
  {15	,	14	,	26	,	27	},
  {16	,	15	,	27	,	28	},
  {17	,	16	,	28	,	29	},
  {18	,	17	,	29	,	30	},
  {19	,	18	,	30	,	31	},
  {20	,	19	,	31	,	32	},
  {21	,	20	,	32	,	33	},
  {22	,	21	,	33	,	34	},
  {23	,	22	,	34	,	35	},
  {24	,	23	,	35	,	36	},
  {25	,	24	,	36	,	37	},
  {14	,	25	,	37	,	26	},
  {27	,	26	,	38	,	39	},
  {28	,	27	,	39	,	40	},
  {29	,	28	,	40	,	41	},
  {30	,	29	,	41	,	42	},
  {31	,	30	,	42	,	43	},
  {32	,	31	,	43	,	44	},
  {33	,	32	,	44	,	45	},
  {34	,	33	,	45	,	46	},
  {35	,	34	,	46	,	47	},
  {36	,	35	,	47	,	48	},
  {37	,	36	,	48	,	49	},
  {26	,	37	,	49	,	38	},
  {39	,	38	,	50	,	51	},
  {40	,	39	,	51	,	52	},
  {41	,	40	,	52	,	53	},
  {42	,	41	,	53	,	54	},
  {43	,	42	,	54	,	55	},
  {44	,	43	,	55	,	56	},
  {45	,	44	,	56	,	57	},
  {46	,	45	,	57	,	58	},
  {47	,	46	,	58	,	59	},
  {48	,	47	,	59	,	60	},
  {49	,	48	,	60	,	61	},
  {38	,	49	,	61	,	50	},
  {51	,	50	,	62	,	63	},
  {52	,	51	,	63	,	64	},
  {53	,	52	,	64	,	65	},
  {54	,	53	,	65	,	66	},
  {55	,	54	,	66	,	67	},
  {56	,	55	,	67	,	68	},
  {57	,	56	,	68	,	69	},
  {58	,	57	,	69	,	70	},
  {59	,	58	,	70	,	71	},
  {60	,	59	,	71	,	72	},
  {61	,	60	,	72	,	73	},
  {50	,	61	,	73	,	62	},
  {63	,	62	,	74	,	75	},
  {64	,	63	,	75	,	76	},
  {65	,	64	,	76	,	77	},
  {66	,	65	,	77	,	78	},
  {67	,	66	,	78	,	79	},
  {68	,	67	,	79	,	80	},
  {69	,	68	,	80	,	81	},
  {70	,	69	,	81	,	82	},
  {71	,	70	,	82	,	83	},
  {72	,	71	,	83	,	84	},
  {73	,	72	,	84	,	85	},
  {62	,	73	,	85	,	74	},
  {75	,	74	,	86	,	87	},
  {76	,	75	,	87	,	88	},
  {77	,	76	,	88	,	89	},
  {78	,	77	,	89	,	90	},
  {79	,	78	,	90	,	91	},
  {80	,	79	,	91	,	92	},
  {81	,	80	,	92	,	93	},
  {82	,	81	,	93	,	94	},
  {83	,	82	,	94	,	95	},
  {84	,	83	,	95	,	96	},
  {85	,	84	,	96	,	97	},
  {74	,	85	,	97	,	86	},
  {87	,	86	,	98	,	99	},
  {88	,	87	,	99	,	100	},
  {89	,	88	,	100	,	101	},
  {90	,	89	,	101	,	102	},
  {91	,	90	,	102	,	103	},
  {92	,	91	,	103	,	104	},
  {93	,	92	,	104	,	105	},
  {94	,	93	,	105	,	106	},
  {95	,	94	,	106	,	107	},
  {96	,	95	,	107	,	108	},
  {97	,	96	,	108	,	109	},
  {86	,	97	,	109	,	98	},
  {99	,	98	,	110	,	111	},
  {100	,	99	,	111	,	112	},
  {101	,	100	,	112	,	113	},
  {102	,	101	,	113	,	114	},
  {103	,	102	,	114	,	115	},
  {104	,	103	,	115	,	116	},
  {105	,	104	,	116	,	117	},
  {106	,	105	,	117	,	118	},
  {107	,	106	,	118	,	119	},
  {108	,	107	,	119	,	120	},
  {109	,	108	,	120	,	121	},
  {98	,	109	,	121	,	110     },
  {111	,	110	,	122	,	123     },
  {112	,	111	,	123	,	124     },
  {113	,	112	,	124	,	125     },
  {114	,	113	,	125	,	126     },
  {115	,	114	,	126	,	127     },
  {116	,	115	,	127	,	128     },
  {117	,	116	,	128	,	129     },
  {118	,	117	,	129	,	130     },
  {119	,	118	,	130	,	131     },
  {120	,	119	,	131	,	132     },
  {121	,	120	,	132	,	133     },
  {110	,	121	,	133	,	122     },
  {123	,	122	,	134	,	134     },
  {124	,	123	,	134	,	134     },
  {125	,	124	,	134	,	134     },
  {126	,	125	,	134	,	134     },
  {127	,	126	,	134	,	134     },
  {128	,	127	,	134	,	134     },
  {129	,	128	,	134	,	134     },
  {130	,	129	,	134	,	134     },
  {131	,	130	,	134	,	134     },
  {132	,	131	,	134	,	134     },
  {133	,	132	,	134	,	134     },
  {122	,	133	,	134	,	134     }  };

int main(int argc, char *argv[])
{
  if(argc < 7) {
    cout << endl
	 << "-- Plot a Particle and its Possibly Contacting Neighbors--" << endl
	 << "Usage: plot1pp  particle_number  contact_file_prefix  particle_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: plot1pp 215 dep_contact dep_particle 0 100 5" << endl << endl;
    return -1;
  }	

  int first, last, incre, pnumber;
  first   = atoi(argv[4]);
  last    = atoi(argv[5]);
  incre   = atoi(argv[6]);
  pnumber = atoi(argv[1]);

  ifstream ifs1,ifs2;
  ofstream ofs;
  char filein1[50], filein2[50];
  char fileout[50];
  char num[4], s[20];

  int TotalCntct, ptcl_1, ptcl_2;
  long double point1_x, point1_y, point1_z, point2_x, point2_y, point2_z;
  long double radius_1, radius_2, penetration, tang_dispmt, contact_radius, R0, E0, normal_force, shear_force;

  int id, type, TotalNum;
  long double cx, cy, cz, rd, wd, lt, ht;
  long double a, b, c, x0, y0, z0, l1, l2, l3, m1, m2, m3, n1, n2, n3, tmp;
  long double x, y, z, xp, yp, zp, t, theta; // x, y, z are local coodinates, xp, yp, zp are global.
  int n, i, j, k, ic, pn, particles;

  char pnum[6];
  sprintf(pnum, "%05d", pnumber);

  // find possible neighbors from all designated snapshots (files)
  vector<int> A;
  vector<int>::iterator it,nit;
  A.push_back(pnumber);
  for(n=first; n<=last; n+=incre) {
    sprintf(num, "%03d", n);
    strcpy(filein1, argv[2]);
    strcat(filein1, "_");
    strcat(filein1, num);

    ifs1.open(filein1);
    if(!ifs1)  { cout<<"ifs1 stream error!"<<endl; exit(-1);}
    ifs1 >> TotalCntct;
    ifs1 >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;
    for (ic=0;ic<TotalCntct;++ic) {
      ifs1>>ptcl_1>>ptcl_2>>point1_x>>point1_y>>point1_z>>point2_x>>point2_y>>point2_z
	  >>radius_1>>radius_2>>penetration>>tang_dispmt>>contact_radius>>R0>>E0>>normal_force>>shear_force;
      if (ptcl_1==pnumber) // find the particle in contact file
	A.push_back(ptcl_2);
    }
    ifs1.close();
  }

  sort(A.begin(), A.end());
  nit=unique(A.begin(), A.end());
  vector<int> possNeighbor(distance(A.begin(), nit));
  copy(A.begin(), nit, possNeighbor.begin());
  //copy(possNeighbor.begin(), possNeighbor.end(), ostream_iterator<int>(cout, "\n"));
  particles=possNeighbor.size();

  // plot all possible particles from all designated snapshots (files)
  for(n=first; n<=last; n+=incre) {
    sprintf(num, "%03d", n);
    strcpy(filein2, argv[3]);
    strcat(filein2, "_");
    strcat(filein2, num);

    strcpy(fileout, "particle_");
    strcat(fileout, pnum);
    strcat(fileout, "_");
    strcat(fileout, num);
    strcat(fileout, ".dat");
    cout << "generating file " << fileout << " ......" << endl;

    // open particle file and the output file
    ifs2.open(filein2);
    if(!ifs2)  { cout<<"ifs2 stream error!"<<endl; exit(-1);}
    ofs.open(fileout);
    if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);
    ifs2 >> TotalNum;
    ifs2>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

    ofs << "ZONE N="<<NODE*particles<<", E="<<ELEM*particles<<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
	
    pn=0;
    while(pn < TotalNum) {
      ifs2 >> id >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
	   >>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;

      it=find(possNeighbor.begin(), possNeighbor.end(), id);
      if(it != possNeighbor.end()) { // the id is found from the vector
	ofs << setfill(' ') << setw(16) << cos(l1)*(-a) + x0
	    << setfill(' ') << setw(16) << cos(m1)*(-a) + y0
	    << setfill(' ') << setw(16) << cos(n1)*(-a) + z0
	    << endl;
		
	for(i=0, x = -a+0.2*a/3; i < 11; ++i) {
	  t = sqrtl(1 - powl(x/a, 2));
	  for(j=0, theta = 0; j < 12; ++j, theta += PI/6) {
	    y = b*t*sinl(theta);
	    z = c*t*cosl(theta);
			
	    xp  = cosl(l1)*x + cosl(l2)*y + cosl(l3)*z + x0;
	    yp  = cosl(m1)*x + cosl(m2)*y + cosl(m3)*z + y0;
	    zp  = cosl(n1)*x + cosl(n2)*y + cosl(n3)*z + z0;
			
	    ofs << setfill(' ') << setw(16) << xp
		<< setfill(' ') << setw(16) << yp
		<< setfill(' ') << setw(16) << zp
		<< endl;
	  }
	  if (i==0 || i==9)
	    x += 0.2*a/3*2;
	  else
	    x += 0.2*a;
	}
		
	ofs << setfill(' ') << setw(16) << cosl(l1)*a + x0
	    << setfill(' ') << setw(16) << cosl(m1)*a + y0
	    << setfill(' ') << setw(16) << cosl(n1)*a + z0
	    << endl;
	/*		
	// plot the two points
	ofs << "ZONE I=2, DATAPACKING=POINT" << endl;
	ofs << setfill(' ') << setw(16) << point1_x
	<< setfill(' ') << setw(16) << point1_y
	<< setfill(' ') << setw(16) << point1_z << endl
	<< setfill(' ') << setw(16) << point2_x
	<< setfill(' ') << setw(16) << point2_y
	<< setfill(' ') << setw(16) << point2_z << endl;
	*/		
      }
      ++pn;
    }
	
	
    for(k = 0; k < particles; ++k) {
      for(i = 0; i < ELEM; ++i)
	ofs << setw(10) << CONNECT[i][0] + NODE*k
	    << setw(10) << CONNECT[i][1] + NODE*k
	    << setw(10) << CONNECT[i][2] + NODE*k
	    << setw(10) << CONNECT[i][3] + NODE*k << endl;
    } 
	
    ifs1.close();
    ifs2.close();
    ofs.close();  
  }	

  return 0;
}
