/* $Id: VoterChenCu.h,v 1.5 2007/07/05 00:10:54 hspark Exp $ */
/* created: paklein (10/25/1998) */
#ifndef _VOTERCHEN_CU_H_
#define _VOTERCHEN_CU_H_

/* base class */
#include "EAM.h"

namespace Tahoe {

/** Voter and Chen EAM copper  */
class VoterChenCu: public EAM
{
public:

	/* constructor */
	VoterChenCu(CBLatticeT& lattice);

	/* unstressed lattice parameter */
	 virtual double LatticeParameter(void) const;

	/* atomic mass */
	 virtual double Mass(void) const;

private:

	/* set the spline data - called by the constructor */
	virtual void SetPairPotential(void);
	virtual void SetEmbeddingEnergy(void);
	virtual void SetElectronDensity(void); 
};

} // namespace Tahoe

/* specific glue functions */
#include "C1FunctionT.h"

namespace Tahoe {

class VCPairPotentialCu: public Tahoe::C1FunctionT
{
	friend class VoterChenCu;
	
public:

	/* constructor */
	VCPairPotentialCu(void);  	

	/* I/O */
	virtual void Print(ostream& out) const;     	    	   	
	virtual void PrintName(ostream& out) const;     	    	   	

	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime.  Default behavior is just to map the virtual functions
	 * above */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	
	/* return 0th, 1st, and 2nd derivative in the respective
	 * fields of the dArrayT */  	
	virtual void SetAll(double x, dArrayT& data) const;   	

private:

	/* parameters */
	double fD;
	double falpha;
	double fR;

private:

	/* non-virtual function calls */
	double function(double x) const;
	double Dfunction(double x) const;
	double DDfunction(double x) const;	
	    	   	    	
};

class VCElectronDensityCu: public Tahoe::C1FunctionT
{
	friend class VoterChenCu;
	
public:

	/* constructor */
	VCElectronDensityCu(void);  	

	/* I/O */
	virtual void Print(ostream& out) const;     	    	   	
	virtual void PrintName(ostream& out) const;     	    	   	

	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime.  Default behavior is just to map the virtual functions
	 * above */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	
	/* return 0th, 1st, and 2nd derivative in the respective
	 * fields of the dArrayT */  	
	virtual void SetAll(double x, dArrayT& data) const;

private:

	/* parameters */
	double fbeta;
	
private:

	/* non-virtual function calls */
	double function(double x) const;
	double Dfunction(double x) const;
	double DDfunction(double x) const;	
	  		   	    	
};

} // namespace Tahoe 
#endif /* _VOTERCHEN_CU_H_ */
