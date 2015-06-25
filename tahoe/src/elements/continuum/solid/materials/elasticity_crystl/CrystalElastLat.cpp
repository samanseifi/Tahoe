/* $Id: CrystalElastLat.cpp,v 1.6 2009/05/21 22:30:27 tdnguye Exp $ */
/*
  File: CrystalElastLat.cpp
*/

#include "CrystalElastLat.h"
#include "CrystalElast.h"
#include "StringT.h"
#include "Array2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

#include "Utils.h"

// for random generation of numbers
#define RANDOM_MAX 32767

// some constants

using namespace Tahoe;

const double pi = 4.0 * atan(1.0);
const double pi180 = pi / 180.;

CrystalElastLat::CrystalElastLat(CrystalElast& poly)
{
  // number of grains
  const int numgrain = poly.NumGrain();

  // check number of grains
  if (numgrain != 1)
    throwRunTimeError("CrystalElastLat::CrystalElastLat: numgrain != 1");

  // main input file
  ifstreamT& in = poly.Input_x();

  // input file for texture data
  ifstreamT tmp;
  ifstreamT& input = OpenExternal(in, tmp, "CrystalElastLat data");

  // ouput file for texture data
  StringT outfilename;
  outfilename.Root(input.filename());
  outfilename.Append(".dat");
  fTextOut.open(outfilename);
  SetStreamPrefs(fTextOut);

  // read/write initial texture
  ReadTexture(input, numgrain);
  //WriteTexture(-1, -1, fNumAngle, 0, fAngles);
}

CrystalElastLat::~CrystalElastLat() { }

void CrystalElastLat::AssignEulerAngles(int kcode, int nelem, int nint, 
                  int ngrn, ArrayT<Array2DT<dArrayT> >& euler) const
{
  switch(kcode)
    {
    case kODF_same_all: 
      // same ODF in all ELs and IPs
	for (int ig = 0; ig < ngrn; ig++)
	{
		int random = (int) (rand()/(RAND_MAX+1.0)*fNumAngle);
      for (int ie = 0; ie < nelem; ie++)
		for (int ip = 0; ip < nint; ip++)
		{
//			euler[ie](ip,ig).SetToScaled(1.0, fAngles[ig]);
			euler[ie](ip,ig).SetToScaled(1.0, fAngles[random]);
//	cout << "\neuler: "<<euler[ie](ip,ig);
		}
	}
	break;

    case kODF_diff_elems:
      // different ODF in all ELs; same ODF in all IPs
      srand(fSeed);
      for (int ie = 0; ie < nelem; ie++)
	for (int ig = 0; ig < ngrn; ig++)
	  {
	    int random = (int) (rand()/(RANDOM_MAX+1.0)*fNumAngle);
	    for (int ip = 0; ip < nint; ip++)
	      euler[ie](ip,ig).SetToScaled(1.0, fAngles[random]);
	  }
      break;

    case kODF_diff_inpts:
      // different ODF in all ELs and IPs
      srand(fSeed);
      for (int ie = 0; ie < nelem; ie++)
	for (int ip = 0; ip < nint; ip++)
	  for (int ig = 0; ig < ngrn; ig++)
	    {
	      int random = (int) (rand()/(RANDOM_MAX+1.0)*fNumAngle);
	      euler[ie](ip,ig).SetToScaled(1.0, fAngles[random]);
	    }
      break;

    case kODF_from_file:
      // ODF read from file for all ELs, IPs, and GRs
      AssignAnglesFromFile(nelem, nint, ngrn, euler);
      break;

    default:
      throwRunTimeError("CrystalElastLat::AssignEulerAngles: Bad kcode");
    }
}

void CrystalElastLat::WriteTexture(int elem, int intpt, int ngrn, int step,
				 const ArrayT<dArrayT>& angle)
{
  // output heading
  if (elem == 0 && intpt == 0) {
     fTextOut << "\nEULER ANGLES AT STEP # " << step << endl;
  }

  // print euler angles: (Kocks, radians) -> (Kocks, degree)
  double psi, the, phi;
  for (int ig = 0; ig < ngrn; ig++)
    {
      psi = angle[ig][0] / pi180;
      the = angle[ig][1] / pi180;
      phi = angle[ig][2] / pi180;
      fTextOut << psi << "  " << the << "  " << phi << "    "
               << elem << "   " << intpt << "   " << ngrn << endl;
    }
}

void CrystalElastLat::AnglesToRotMatrix(const dArrayT& angle, dMatrixT& C) const 
{
  double sps, cps, sth, cth, sph, cph;

  sps = sin(angle[0]);
  cps = cos(angle[0]);
  sth = sin(angle[1]);
  cth = cos(angle[1]);
  sph = sin(angle[2]);
  cph = cos(angle[2]);
  
  C(0,0) = -sps * sph - cps * cph * cth;
  C(1,0) =  cps * sph - sps * cph * cth;
  C(2,0) =  cph * sth;
  C(0,1) =  cph * sps - sph * cps * cth;
  C(1,1) = -cps * cph - sps * sph * cth;
  C(2,1) =  sph * sth;
  C(0,2) =  cps * sth;
  C(1,2) =  sps * sth;
  C(2,2) =  cth;
}

void CrystalElastLat::RotMatrixToAngles(const dMatrixT& C, dArrayT& angle) const
{
  double sth;

  // compute angles in radians
  angle[1] = acos(C(2,2));  
  if (fabs(C(2,2)) != 1.0)
    {
      sth = sin(angle[1]);
      angle[0] = atan2(C(1,2)/sth, C(0, 2)/sth);
      angle[2] = atan2(C(2,1)/sth, C(2, 0)/sth);
    }
  else
    {
      angle[0] = 0.;
      angle[2] = atan2(-C(0,1), -C(0,0));
    }
}

void CrystalElastLat::ReadTexture(ifstreamT& in, int numgrain)
{
  // number of euler angles in input file
  in >> fNumAngle;
  if (fNumAngle < numgrain) 
    throwRunTimeError("CrystalElastLat: ReadTexture: NumAngle < NumGrain");

  // header to output initial texture
  fTextOut << "\nINITIAL TEXTURE,  No Euler angles = " << fNumAngle << endl;

  // allocate space to read euler angles 
  fAngles.Dimension(fNumAngle);
  for (int i = 0; i < fNumAngle; i++) fAngles[i].Dimension(3);

  // flags for angle convention
  // iikc = 0 : angles input in Kocks convention :  (psi,the,phi)
  //        1 : angles input in Canova convention : (ph,th,om)
  //            ph = 90 + phi; th = the; om = 90 - psi
  // iidr = 0 : angles input in degrees
  //          : angles input in radians
  int iikc, iidr;
  in >> iikc >> iidr;

  // some constants for angle conversion
  double piby2 = 90.;
  if (iidr == 1) piby2 = pi / 2.0;
  
  // read Euler angles & storage them in : (Cocks, radians)
  for (int i = 0; i < fNumAngle; i++)
    {
      in >> fAngles[i];
      fTextOut << fAngles[i][0] << " "
               << fAngles[i][1] << " "
               << fAngles[i][2] << endl;
      if (iikc == 1)
	{
	  double ph = fAngles[i][0];
	  fAngles[i][0] = piby2 - fAngles[i][2];
	  fAngles[i][2] = ph - piby2;
	}
      if (iidr == 0) fAngles[i] *= pi180;
    }

  // seed for the case of random distribution of euler angles
  in >> fSeed;
}

void CrystalElastLat::AssignAnglesFromFile(int nelem, int nint, int ngrn, 
			      ArrayT< Array2DT<dArrayT> >& euler) const
{
  // for specific cases  
  if (nelem == fNumAngle && ngrn == 1)    // BICRYSTAL 
    {
      for (int ie = 0; ie < nelem; ie++)
	for (int ip = 0; ip < nint; ip++)
	  for (int ig = 0; ig < ngrn; ig++)      
	    euler[ie](ip,ig).SetToScaled(1.0, fAngles[ie]);
    }
  else if (nelem*nint*ngrn == fNumAngle)  // assign to each GR/IP/EL from file
    {
      int n = 0;
      for (int ie = 0; ie < nelem; ie++)
	for (int ip = 0; ip < nint; ip++)
	  for (int ig = 0; ig < ngrn; ig++)      
	    euler[ie](ip,ig).SetToScaled(1.0, fAngles[n++]);
    }
  else
    throwRunTimeError("CrystalElastLat::AssignAnglesFromFile: check input file");
}
