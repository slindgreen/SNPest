#ifndef __optimizexx_h
#define __optimizexx_h

#include <fstream>

#define SETUP_C_SUBSCRIPTS // allow zero-based []-indexation of newmat vectors and matrices

#include "boost/function.hpp"
#include "newmat/newmatio.h"
#include "optpp/OptQNewton.h"
#include "optpp/Constraint.h"
#include "optpp/BoundConstraint.h"
#include "optpp/OptBCQNewton.h"
#include "optpp/OptBaQNewton.h"
#include "phy/PhyDef.h"

//#include "NLF.h"
//#include "NLP.h"
//#include "OptBCFDNewton.h"
//#include "OptFDNIPS.h"
//#include "OptFDNewton.h"
//#include "OptBCEllipsoid.h"
//#include "OptFDNewton.h"
//#include "OptCG.h"


using namespace OPTPP;

/** global optimization */
phy::vector_t const minimizexx(phy::vector_t x, boost::function< double (phy::vector_t const &) > const & objFct, string const & statusFileName = "", unsigned maxFEval = 10000, unsigned maxIter = 1000);

////////////////////////////////////////////////////////////////
// useful transformations and conversions
////////////////////////////////////////////////////////////////

/** transform maps ]l;oo[ into ]-oo;oo[ */
double transform(double x, double l);
/** deTransforms maps ]-oo;oo[ into ]l;oo[ (inverse of transform for given l) */
double deTransform(double y, double l);

/** same as transform and deTransform, but acts on whole vectors. Mapping done in place. */
void transformVec(phy::vector_t & v, double l);
void deTransformVec(phy::vector_t & v, double l);

/** force x[i] to be in [min_value; max_value]. */
phy::vector_t const forceInRange(phy::vector_t const &x, double min_value, double max_value);

/** Conversion from NEWMAT::ColumnVector to phy::vector_t and visa versa. */
NEWMAT::ColumnVector const toColumnVector(phy::vector_t const & v);
phy::vector_t const toNumber(NEWMAT::ColumnVector const & v);

 
// old stuff
//NEWMAT::ColumnVector const boundMinimizexx(NEWMAT::ColumnVector x, 
//					   NEWMAT::ColumnVector lower, 
//					   NEWMAT::ColumnVector upper, 
//					   BaseFunctionObject & f_obj, 
//					   string const & statusFileName);
//
//
//NEWMAT::Matrix const estCovarianceMatrix(NEWMAT::ColumnVector x, 
//					 BaseFunctionObject & f_obj, 
//					 ostream & str = cout);
//
//NEWMAT::Matrix const estGradient(NEWMAT::ColumnVector x, 
//				 BaseFunctionObject & f_obj, 
//				 ostream & str = cout);
//


#endif /* __optimizexx_h */
