/*
  File: LatticeOrient.h
*/

#ifndef _LATTICE_ORIENT_H_
#define _LATTICE_ORIENT_H_

#include <fstream>
#include "ArrayT.h"
#include "dArrayT.h"
#include "Array2DT.h"
#include "ofstreamT.h"
#include "StringT.h"


namespace Tahoe {

class PolyCrystalMatT;
class ifstreamT;
class dMatrixT;

class LatticeOrient
{
 public:
  // constructor
  LatticeOrient(PolyCrystalMatT& poly);

  // destructor
  ~LatticeOrient();

  // assign Euler angles to grain/elem/IPs 
  void AssignEulerAngles(int kcode, int nelem, int nint, int ngrn,
			 ArrayT<Array2DT<dArrayT> >& euler); 
  void AssignEulerAngles_block(int nelem, int start_elem, int nint, int ngrn,
			 ArrayT<Array2DT<dArrayT> >& euler); 

  // output Euler angles 
  void WriteTexture(int group, int elem, int intpt, int ngrn, int step,
		    const ArrayT<dArrayT>& angle);

  // compute rotation matrix from Euler angles
  void AnglesToRotMatrix(const dArrayT& angle, dMatrixT& C) const;

  // compute Euler angles from rotation matrix 
  void RotMatrixToAngles(const dMatrixT& C, dArrayT& angle) const;

 private:
  // flag to distribute crystal orientations
  enum ODFDist { kODF_same_all   = 1,     // same ODF in all ELs/IPs  
 		 kODF_diff_elems = 2,     // diff ODF in all ELs; same ODF in all IPs 
		 kODF_diff_inpts = 3,     // diff ODF in all ELs/IPs
		 kODF_from_file  = 4,     // ODF read from file for all Els/IPs/GRs
		 kDOF_diff_blocks = 5 };   

  // input Euler angles
  void ReadTexture(ifstreamT& in, int numgrain);

  // assign ODF at each IP/EL as read from file
  void AssignAnglesFromFile(int nelem, int nint, int ngrn,
		        ArrayT<Array2DT<dArrayT> >& euler) const;

 private:
  // number of euler angles in input file
  int fNumAngle;

  // seed for random distribution of euler angles
  long fSeed;

  // array to hold euler angles from input file
  ArrayT<dArrayT> fAngles;

  // stream to output texture data
  ofstreamT fTextOut;
  StringT   fOutFilename;
};

} // namespace Tahoe 
#endif /* _LATTICE_ORIENT_H_ */
