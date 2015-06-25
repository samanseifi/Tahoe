/* $Id: JoinResults.h,v 1.2 2004/11/11 03:57:24 paklein Exp $ */
#ifndef _TRANSLATE_JOIN_H_
#define _TRANSLATE_JOIN_H_

/* base class */
#include "TranslateIOManager.h"

namespace Tahoe {

/** join results data from multiple results files */
class JoinResults: public TranslateIOManager
{
public:

	/** constructor */
	JoinResults(ostream& message, istream& in, bool write);

	/** run */
	virtual void Translate(const StringT& program, const StringT& version, 
		const StringT& title);

protected:

	/** set input sources */
	void SetInput(void);

private:

	/** init output sets */
	void InitOutputSets(void);

private:

	/** array of input sources */
	ArrayT<StringT> fSources;
	
	/** file types */
	ArrayT<IOBaseT::FileTypeT> fFileTypes;

	/** output ID's */
	int fOutputID;
};

} /* namespace Tahoe */

#endif /* _TRANSLATE_JOIN_H_ */
