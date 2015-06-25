/* $Id: EAMFCC3DSym_110surf.cpp,v 1.3 2007/06/13 01:50:32 hspark Exp $ */
/* created: paklein (12/06/1996) */
#include "EAMFCC3DSym_110surf.h"

using namespace Tahoe;

/* Bond table parameters */
/* This file assumes <100> bulk with {110} surfaces */
/* Need 6 layers for {110} instead of 4 for {100} */
const int kEAMFCC3DSurfBonds        = 222;	// updated to include surface3 and surface4 atoms (186)
const int kEAMFCC3DNumBonds			= 42;	// 54
const int kEAMFCC3DSurf1Bonds       = 25;	// 33
const int kEAMFCC3DSurf2Bonds       = 33;	// 45
const int kEAMFCC3DSurf3Bonds       = 38;
const int kEAMFCC3DNumLatticeDim 	=  3;
const int kEAMFCC3DNumAtomsPerCell	=  4;
const int kEAMFCC3DNumAtomsPerArea  =  2;
const double piby2 = 4.0 * atan(1.0) / 2.0;
const double root2by2 = sqrt(2.0)/2.0;
const double root2by4 = sqrt(2.0)/4.0;
const double threeroot2by4 = 3.0*sqrt(2.0)/4.0;

/* constructor */
EAMFCC3DSym_110surf::EAMFCC3DSym_110surf(int nshells, int normal):
	EAMFCC3D_surf(nshells, normal),
	fNormalCode(normal)
{
	SetName("FCC_EAM_Cauchy-Born");
}

/**********************************************************************
 * Protected
 **********************************************************************/
	
void EAMFCC3DSym_110surf::LoadBondTable(void)
{
	/* dimension work space */
	/* THESE ARRAYS ACTUALLY DEFINED IN BondLatticeT? */
	/* NEED TO ADD FUNCTIONALITY FOR THESE IN BondLatticeT.cpp */
	fBondCounts.Dimension(kEAMFCC3DSurfBonds);
	fBulkCounts.Dimension(kEAMFCC3DNumBonds);
	fSurf1Counts.Dimension(kEAMFCC3DSurf1Bonds);
	fSurf2Counts.Dimension(kEAMFCC3DSurf2Bonds);
	fSurf3Counts.Dimension(kEAMFCC3DSurf3Bonds);
	fDefLength.Dimension(kEAMFCC3DSurfBonds);
	fDefBulk.Dimension(kEAMFCC3DNumBonds);
	fDefSurf1.Dimension(kEAMFCC3DSurf1Bonds);
	fDefSurf2.Dimension(kEAMFCC3DSurf2Bonds);
	fDefSurf3.Dimension(kEAMFCC3DSurf3Bonds);
	fBonds.Dimension(kEAMFCC3DSurfBonds, kEAMFCC3DNumLatticeDim);
	fBulkBonds.Dimension(kEAMFCC3DNumBonds,3);
	fSurf1Bonds.Dimension(kEAMFCC3DSurf1Bonds,3);
	fSurf2Bonds.Dimension(kEAMFCC3DSurf2Bonds,3);
	fSurf3Bonds.Dimension(kEAMFCC3DSurf3Bonds,3);
	fAtomType.Dimension(kEAMFCC3DSurfBonds);

	dArray2DT temp_bonds, temp_bonds2, temp_bulk, temp_surf1, temp_surf2, tempb, temps1, temps2;
	dArray2DT temp_surf3, temps3;
	temp_bonds.Dimension(kEAMFCC3DSurfBonds, 3);	// temporary bond table before rotation
	temp_bonds2.Dimension(kEAMFCC3DSurfBonds, 3);	// Currently have # of bonds for {100} surfaces
	temp_bulk.Dimension(kEAMFCC3DNumBonds, 3);
	temp_surf1.Dimension(kEAMFCC3DSurf1Bonds, 3);
	temp_surf2.Dimension(kEAMFCC3DSurf2Bonds, 3);
	temp_surf3.Dimension(kEAMFCC3DSurf3Bonds, 3);
	tempb.Dimension(kEAMFCC3DNumBonds, 3);
	temps1.Dimension(kEAMFCC3DSurf1Bonds, 3);
	temps2.Dimension(kEAMFCC3DSurf2Bonds, 3);
	temps3.Dimension(kEAMFCC3DSurf3Bonds, 3);

	/* all bonds appear once */
	fBondCounts = 1;
	fBulkCounts = 1;
	fSurf1Counts = 1;
	fSurf2Counts = 1;
	fSurf3Counts = 1;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;
	fDefBulk = 0.0;
	fDefSurf1 = 0.0;
	fDefSurf2 = 0.0;
	fDefSurf3 = 0.0;

	/* undeformed bond data for bulk atom with 4th neighbor interactions */
	double bulkbond[kEAMFCC3DNumBonds][kEAMFCC3DNumLatticeDim] = {
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors)
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},
		{0.0, 0.0, root2by2},	// New NN for third layer atoms
		{0.0, root2by2, root2by2},	// New 2NN for third layer atoms
		{0.0, -root2by2, root2by2},
		{1.0, 0.0, root2by2},	// New 3NN for third layer atoms
		{-1.0, 0.0, root2by2},
		{0.5, root2by4, threeroot2by4},	// New 3NN to complete bulk
		{0.5, -root2by4, threeroot2by4},
		{-0.5, root2by4, threeroot2by4},
		{-0.5, -root2by4, threeroot2by4}
	};

//	int bulk4[42] = {15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,13,13,13,13,13,12,12,12,12};
//	int bulk5[42] = {15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,13,13,13,13};
//	int bulk6[42] = {15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14};
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_bulk(i,j) = bulkbond[i][j];

	/* Bond table for an atom on the surface */
	/* Note that surface normal is [001] instead of [100] for {100} */
	double surf1bond[kEAMFCC3DSurf1Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors)
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4}
	};

//	int surf1[25] = {0,0,2,1,1,1,1,0,0,2,2,0,0,0,0,2,2,3,3,3,3,1,1,1,1};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_surf1(i,j) = surf1bond[i][j];

	/* Bond table for an atom 1 layer into the bulk */
	double surf2bond[kEAMFCC3DSurf2Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors)
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4}
	};

//	int surf2[33] = {4,4,6,5,5,5,5,4,4,6,6,4,4,4,4,6,6,7,7,7,7,5,5,5,5,7,7,7,7,7,7,7,7};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_surf2(i,j) = surf2bond[i][j];

	/* Bond table for an atom 2 layers into the bulk */
	double surf3bond[kEAMFCC3DSurf3Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors)
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},
		{0.0, 0.0, root2by2},	// New NN for third layer atoms
		{0.0, root2by2, root2by2},	// New 2NN for third layer atoms
		{0.0, -root2by2, root2by2},
		{1.0, 0.0, root2by2},	// New 3NN for third layer atoms
		{-1.0, 0.0, root2by2}
	};

//	int surf3[38] = {8,8,9,9,9,9,9,8,8,9,9,8,8,8,8,9,9,11,11,11,11,9,9,9,9,11,11,11,11,11,11,11,11,10,10,10,10,10};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf3Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_surf3(i,j) = surf3bond[i][j];


	/* work space arrays for storing interaction types */
	// OTHER BOND TABLES ARE USED TO CALCULATE REPRESENTATIVE ELECTRON DENSITIES FOR VARIOUS
	// SURFACE LAYERS - CORRELATE TO THOSE USED IN EAMT.CPP
	iArrayT allbonds(222);
	/* Interaction type key */
	// 0 = surface1/surface1
	// 1 = surface1/surface2
	// 2 = surface1/surface3
	// 3 = surface1/bulk
	// 4 = surface2/surface1
	// 5 = surface2/surface2
	// 6 = surface2/surface3
	// 7 = surface2/bulk
	// 8 = surface3/surface3
	// 9 = surface3/bulk
	// 10= surface3/surface1
	// 11= surface3/surface2
	// 12= bulk/surface1
	// 13= bulk/surface2
	// 14= bulk/surface3
	// 15= bulk/bulk

	/* FILL THIS IN WITH ALL BOND TABLE INTERACTIONS */
	int surf1n[222]={0,0,2,1,1,1,1,0,0,2,2,0,0,0,0,2,2,3,3,3,3,1,1,1,1,4,4,6,5,5,5,5,4,4,6,6,4,4,4,4,6,6,7,7,7,7,5,5,5,5,7,7,7,7,7,7,7,7,8,8,9,9,9,9,9,8,8,9,9,8,8,8,8,9,9,11,11,11,11,9,9,9,9,11,11,11,11,11,11,11,11,10,10,10,10,10,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,13,13,13,13,13,12,12,12,12,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,13,13,13,13,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14};
	for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		allbonds[i] = surf1n[i];
	
	fAtomType.CopyIn(0, allbonds);
	
	/* New bond table for surface clusters - change dimensions! */
	double bonddata[kEAMFCC3DSurfBonds][kEAMFCC3DNumLatticeDim] = {
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors)	// LAYER 1 STARTS
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},	// LAYER 1 ENDS
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors)	// LAYER 2 STARTS
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},	// LAYER 2 ENDS
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors) // LAYER 3 STARTS
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},
		{0.0, 0.0, root2by2},	// New NN for third layer atoms
		{0.0, root2by2, root2by2},	// New 2NN for third layer atoms
		{0.0, -root2by2, root2by2},
		{1.0, 0.0, root2by2},	// New 3NN for third layer atoms
		{-1.0, 0.0, root2by2},	// LAYER 3 ENDS
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors) // BULK 4 STARTS
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},
		{0.0, 0.0, root2by2},	// New NN for third layer atoms
		{0.0, root2by2, root2by2},	// New 2NN for third layer atoms
		{0.0, -root2by2, root2by2},
		{1.0, 0.0, root2by2},	// New 3NN for third layer atoms
		{-1.0, 0.0, root2by2},
		{0.5, root2by4, threeroot2by4},	// New 3NN to complete bulk
		{0.5, -root2by4, threeroot2by4},
		{-0.5, root2by4, threeroot2by4},
		{-0.5, -root2by4, threeroot2by4},	// BULK 4 ENDS
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors) // BULK 5 STARTS
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},
		{0.0, 0.0, root2by2},	// New NN for third layer atoms
		{0.0, root2by2, root2by2},	// New 2NN for third layer atoms
		{0.0, -root2by2, root2by2},
		{1.0, 0.0, root2by2},	// New 3NN for third layer atoms
		{-1.0, 0.0, root2by2},
		{0.5, root2by4, threeroot2by4},	// New 3NN to complete bulk
		{0.5, -root2by4, threeroot2by4},
		{-0.5, root2by4, threeroot2by4},
		{-0.5, -root2by4, threeroot2by4},	// BULK 5 ENDS
		{0.0, root2by2, 0.0}, // Surface cluster (7 nearest neighbors) // BULK 6 STARTS
		{0.0, -root2by2, 0.0},
		{0.0, 0.0, -root2by2},
		{0.5, root2by4, -root2by4},
		{0.5, -root2by4, -root2by4},
		{-0.5, root2by4, -root2by4},
		{-0.5, -root2by4, -root2by4},
		{1.0, 0.0, 0.0},	// Surface cluster (4 2nd shell neighbors)
		{-1.0, 0.0, 0.0},
		{0.0, root2by2, -root2by2},
		{0.0, -root2by2, -root2by2},
		{1.0, root2by2, 0.0},	// Surface cluster (14 3rd shell neighbors)
		{1.0, -root2by2, 0.0},
		{-1.0, root2by2, 0.0},
		{-1.0, -root2by2, 0.0},
		{1.0, 0.0, -root2by2},
		{-1.0, 0.0, -root2by2},
		{0.5, root2by4, -threeroot2by4},
		{0.5, -root2by4, -threeroot2by4},
		{-0.5, root2by4, -threeroot2by4},
		{-0.5, -root2by4, -threeroot2by4},
		{0.5, threeroot2by4, -root2by4},
		{0.5, -threeroot2by4, -root2by4},
		{-0.5, threeroot2by4, -root2by4},
		{-0.5, -threeroot2by4, -root2by4},
		{0.5, root2by4, root2by4},	// New NN for second layer atoms
		{0.5, -root2by4, root2by4},
		{-0.5, root2by4, root2by4},
		{-0.5, -root2by4, root2by4},
		{0.5, threeroot2by4, root2by4},	// New 3NN for second layer atoms
		{0.5, -threeroot2by4, root2by4},
		{-0.5, threeroot2by4, root2by4},
		{-0.5, -threeroot2by4, root2by4},
		{0.0, 0.0, root2by2},	// New NN for third layer atoms
		{0.0, root2by2, root2by2},	// New 2NN for third layer atoms
		{0.0, -root2by2, root2by2},
		{1.0, 0.0, root2by2},	// New 3NN for third layer atoms
		{-1.0, 0.0, root2by2},
		{0.5, root2by4, threeroot2by4},	// New 3NN to complete bulk
		{0.5, -root2by4, threeroot2by4},
		{-0.5, root2by4, threeroot2by4},
		{-0.5, -root2by4, threeroot2by4}	// BULK 6 ENDS
	};
	
	/* Rotate Bond Tables based on fNormalCode and rotation matrices */
	/* Create temporary bond table temp_bonds that combines bonddata */
	for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_bonds(i,j) = bonddata[i][j];
	
	/* Now manipulate temp_bonds */
	dMatrixT blah1(3);
	dArrayT asdf(3), prod(3), asdfb(3), prodb(3), asdfs1(3), prods1(3), asdfs2(3), prods2(3);
	dArrayT asdfs3(3), prods3(3);
	if (fNormalCode == 0)	// rotate [0,0,1] to [1,0,0]
	{
		temp_bonds2 = temp_bonds;
		tempb = temp_bulk;
		temps1 = temp_surf1;
		temps2 = temp_surf2;
		temps3 = temp_surf3;
		
		blah1 = RotationMatrixA(piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		{
			tempb.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tempb.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		{
			temps1.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			temps1.SetRow(i,prods1);
		}
		for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		{
			temps2.RowCopy(i,asdfs2);
			blah1.Multx(asdfs2,prods2);
			temps2.SetRow(i,prods2);
		}
		for (int i = 0; i < kEAMFCC3DSurf3Bonds; i++)
		{
			temps3.RowCopy(i,asdfs3);
			blah1.Multx(asdfs3,prods3);
			temps3.SetRow(i,prods3);
		}		
		fBonds = temp_bonds2;
		fBulkBonds = tempb;
		fSurf1Bonds = temps1;
		fSurf2Bonds = temps2;
		fSurf3Bonds = temps3;
	}
	else if (fNormalCode == 1)	// rotate [0,0,1] to [-1,0,0]
	{
		temp_bonds2 = temp_bonds;
		tempb = temp_bulk;
		temps1 = temp_surf1;
		temps2 = temp_surf2;
		temps3 = temp_surf3;
		
		blah1 = RotationMatrixA(-piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		{
			tempb.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tempb.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		{
			temps1.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			temps1.SetRow(i,prods1);
		}
		for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		{
			temps2.RowCopy(i,asdfs2);
			blah1.Multx(asdfs2,prods2);
			temps2.SetRow(i,prods2);
		}
		for (int i = 0; i < kEAMFCC3DSurf3Bonds; i++)
		{
			temps3.RowCopy(i,asdfs3);
			blah1.Multx(asdfs3,prods3);
			temps3.SetRow(i,prods3);
		}		
		fBonds = temp_bonds2;
		fBulkBonds = tempb;
		fSurf1Bonds = temps1;
		fSurf2Bonds = temps2;
		fSurf3Bonds = temps3;
	}
	else if (fNormalCode == 2)	// rotate [0,0,1] to [0,1,0]
	{
		temp_bonds2 = temp_bonds;
		tempb = temp_bulk;
		temps1 = temp_surf1;
		temps2 = temp_surf2;
		temps3 = temp_surf3;
		
		blah1 = RotationMatrixB(piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		{
			tempb.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tempb.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		{
			temps1.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			temps1.SetRow(i,prods1);
		}
		for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		{
			temps2.RowCopy(i,asdfs2);
			blah1.Multx(asdfs2,prods2);
			temps2.SetRow(i,prods2);
		}
		for (int i = 0; i < kEAMFCC3DSurf3Bonds; i++)
		{
			temps3.RowCopy(i,asdfs3);
			blah1.Multx(asdfs3,prods3);
			temps3.SetRow(i,prods3);
		}		
		fBonds = temp_bonds2;
		fBulkBonds = tempb;
		fSurf1Bonds = temps1;
		fSurf2Bonds = temps2;
		fSurf3Bonds = temps3;
	}
	else if (fNormalCode == 3)	// rotate [0,0,1] to [0,-1,0]
	{
		temp_bonds2 = temp_bonds;
		tempb = temp_bulk;
		temps1 = temp_surf1;
		temps2 = temp_surf2;
		temps3 = temp_surf3;
		
		blah1 = RotationMatrixB(-piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		{
			tempb.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tempb.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		{
			temps1.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			temps1.SetRow(i,prods1);
		}
		for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		{
			temps2.RowCopy(i,asdfs2);
			blah1.Multx(asdfs2,prods2);
			temps2.SetRow(i,prods2);
		}
		for (int i = 0; i < kEAMFCC3DSurf3Bonds; i++)
		{
			temps3.RowCopy(i,asdfs3);
			blah1.Multx(asdfs3,prods3);
			temps3.SetRow(i,prods3);
		}		
		fBonds = temp_bonds2;
		fBulkBonds = tempb;
		fSurf1Bonds = temps1;
		fSurf2Bonds = temps2;
		fSurf3Bonds = temps3;
	}
	else if (fNormalCode == 4)	// [0,0,1]
	{
		fBonds = temp_bonds;	// this table is the default, i.e. [0,0,1]
		fBulkBonds = temp_bulk;
		fSurf1Bonds = temp_surf1;
		fSurf2Bonds = temp_surf2;
		fSurf3Bonds = temp_surf3;
	}	
	else if (fNormalCode == 5)	// rotate [0,0,1] to [0,0,-1]
	{
		temp_bonds2 = temp_bonds;
		fBonds = temp_bonds2;
		fBonds *= -1.0;
		fBulkBonds = temp_bulk;
		fBulkBonds *= -1.0;
		fSurf1Bonds = temp_surf1;
		fSurf1Bonds *= -1.0;
		fSurf2Bonds = temp_surf2;
		fSurf2Bonds *= -1.0;
		fSurf3Bonds = temp_surf3;
		fSurf3Bonds *= -1.0;
	}	

	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
	fBulkBonds *= fLatticeParameter;
	fSurf1Bonds *= fLatticeParameter;
	fSurf2Bonds *= fLatticeParameter;
	fSurf3Bonds *= fLatticeParameter;
}

/*************************************************************************
 * Private
 *************************************************************************/
 
 /* Rotate bonds with [0,0,1] normal to bonds with [1,0,0] normal using pi/2 */
dMatrixT EAMFCC3DSym_110surf::RotationMatrixA(const double angle)
 {
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = cos(angle);
	rmatrix(0,2) = sin(angle);
	rmatrix(2,0) = -sin(angle);
	rmatrix(2,2) = cos(angle);
	rmatrix(1,1) = 1.0;
	
	return rmatrix;
 }
 
/* Rotate bonds with [0,0,1] normal to bonds with [0,1,0] normal using pi/2 */
dMatrixT EAMFCC3DSym_110surf::RotationMatrixB(const double angle)
{
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(1,1) = cos(angle);
	rmatrix(1,2) = sin(angle);
	rmatrix(0,0) = 1.0;
	rmatrix(2,1) = -sin(angle);
	rmatrix(2,2) = cos(angle);
	
	return rmatrix;
}
