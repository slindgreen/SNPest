#include "optimizexx.h"

//global pointers used by initWrapper and fctWrapper
boost::function< double (phy::vector_t const &) > const * glbObjFct;
phy::vector_t * glbInitParams;


void initWrapper(int ndim, NEWMAT::ColumnVector & x)
{
  assert( ndim == x.Nrows() );
  assert( ndim == static_cast<signed>( glbInitParams->size() ) );
  x = toColumnVector( *glbInitParams );
}


// wrap boost function in interface required by opt++
void fctWrapper(int ndim, const NEWMAT::ColumnVector & x, double & fx, int & result)
{
  fx = (*glbObjFct)( toNumber(x) );
  result = NLPFunction;
}


phy::vector_t const minimizexx(phy::vector_t x, boost::function< double (phy::vector_t const &) > const & objFct, string const & statusFileName, unsigned maxFEval, unsigned maxIter)
{
  glbObjFct = & objFct;
  glbInitParams = & x;

  FDNLF1 nlp(x.size(), fctWrapper, initWrapper);

  OptQNewton objfcn(&nlp);
  //OptFDNewton objfcn(&nlp);
  //OptCG objfcn(&nlp);

  objfcn.setSearchStrategy(LineSearch); //alt: TrustRegion
  // set stopping criteria: see http://csmr.ca.sandia.gov/opt%2B%2B/opt++2.4_doc/html/ControlParameters.html
  objfcn.setMaxFeval(maxFEval); 
  objfcn.setMaxIter(maxIter);

  // The "0" in the second argument says to create a new file.  A "1"
  // would signify appending to an existing file.
  if (statusFileName != "")
    if (!objfcn.setOutputFile(statusFileName.c_str(), 1) )
      cerr << "main: output file open failed" << endl;

  objfcn.optimize();

  objfcn.printStatus("Solution from quasi-newton");
  objfcn.cleanup();

  return toNumber( objfcn.getXPrev() );
}


// map ]l;oo[ into ]-oo;oo[
double transform(double x, double l)
{
  return log(10*(x-l) );
}


// map ]-oo;oo[ back into ]l;oo[  
double deTransform(double y, double l)
{
  return l + exp(y)/10;
}


void transformVec(phy::vector_t & v, double l)
{
  for (unsigned i = 0; i < v.size(); i++)
    v[i] = transform(v[i], l);
}


void deTransformVec(phy::vector_t & v, double l)
{
  for (unsigned i = 0; i < v.size(); i++)
    v[i] = deTransform(v[i], l);
}

// 
// void InitObject::operator()(phy::vector_t & x) 
// {
//   assert( x.size() == x_.size() ); 
//   x = x_;
// }
// 
// 
// double FctObject::operator()(phy::vector_t const & x)
// {
//   return objFct_(x);
// }
// 


phy::vector_t const forceInRange(phy::vector_t const &x, double min_value, double max_value)
{
  int n     = x.size();
  phy::vector_t y(x); 

  for (int i = 0; i < n; i++) {
    if (x[i] < min_value)
      y[i] = min_value;
    if (x[i] > max_value)
      y[i] = max_value;
  }
  return y;
}

NEWMAT::ColumnVector const toColumnVector(phy::vector_t const & v)
{
  unsigned n = v.size();
  NEWMAT::ColumnVector u(n);
  for (unsigned i = 0; i < n; i++)
    u[i] = v[i];
  return u;
}


phy::vector_t const toNumber(NEWMAT::ColumnVector const & v)
{
  unsigned n = v.Nrows();
  phy::vector_t u(n);
  for (unsigned i = 0; i < n; i++)
    u[i] = v[i];
  return u;
}


///////// old approach

//void InitFunctionObject::operator()(int ndim, NEWMAT::ColumnVector & x) 
//{
//  if (x.Nrows() != init_x.Nrows())
//    throw string("initFunctionObject: Size of vectors do not match");
//
//  x = init_x;
//}
//
//
//void FunctionObject::operator()(int ndim, const NEWMAT::ColumnVector & x, double & fx, int & result) 
//{
//  fx     = f_obj_(x);
//  result = NLPFunction;
//}
//
//
//void init_f_wrapper(int ndim, NEWMAT::ColumnVector & x)
//{
//  (*global_init_f_obj)(ndim, x);
//}
//
//
//void f_wrapper(int ndim, const NEWMAT::ColumnVector & x, double & fx, int & result)
//{
//  (*global_f_obj)(ndim, x, fx, result);
//}
//
//
//NEWMAT::ColumnVector const minimizexx(NEWMAT::ColumnVector x, BaseFunctionObject & f_obj, string const & statusFileName) 
//{
//  int                  ndim(x.Nrows() );
//  InitFunctionObject   init_f(x);
//  FunctionObject       f(f_obj);
//
//  global_init_f_obj = & init_f; 
//  global_f_obj      = & f;
//
//  FDNLF1 nlp(ndim, f_wrapper, init_f_wrapper);
//
//  OptQNewton objfcn(&nlp);
//  //OptFDNewton objfcn(&nlp);
//  //OptCG objfcn(&nlp);
//
//  objfcn.setSearchStrategy(LineSearch); //TrustRegion
//  objfcn.setMaxFeval(5000);
//  //objfcn.setGradTol(0.001);
//
//  // objfcn.setFcnTol(1.e-5);
//  // objfcn.setStepTol(0.001);
//
//  // The "0" in the second argument says to create a new file.  A "1"
//  // would signify appending to an existing file.
//  if (statusFileName != "")
//    if (!objfcn.setOutputFile(statusFileName.c_str(), 1) )
//      cerr << "main: output file open failed" << endl;
//
//  objfcn.optimize();
//
//  objfcn.printStatus("Solution from quasi-newton");
//  objfcn.cleanup();
//
//  return objfcn.getXPrev();
//}
//
//
//NEWMAT::ColumnVector const boundMinimizexx(NEWMAT::ColumnVector x, 
//				   NEWMAT::ColumnVector lower, 
//				   NEWMAT::ColumnVector upper, 
//				   BaseFunctionObject & f_obj, 
//				   string const & statusFileName = "")
//{
//  int                  ndim(x.Nrows() );
//  InitFunctionObject   init_f(x);
//  FunctionObject       f(f_obj);
//
//  global_init_f_obj = & init_f; 
//  global_f_obj      = & f;
//
//  // Construct the compound constraint 
//  Constraint bc            = new BoundConstraint(ndim, lower, upper);
//  CompoundConstraint * cc  = new CompoundConstraint(bc);
//
//  FDNLF1 nlp(ndim, & f_wrapper, & init_f_wrapper, cc);
//
//  //OptBCEllipsoid objfcn(&nlp);
//  //OptBCFDNewton objfcn(&nlp); // does not respect bounds...
//  OptBaQNewton objfcn(&nlp);
//
//  //  objfcn.setMethod.setMethod("OptCG");
//  objfcn.setSearchStrategy(LineSearch);
//  //objfcn.setSearchStrategy(TrustRegion);
//  objfcn.setMaxFeval(1000);
//  objfcn.setFcnTol(1.e-5);
//  objfcn.setConTol(0.00001);
//
//  // The "0" in the second argument says to create a new file.  A "1"
//  // would signify appending to an existing file.
//  if (statusFileName != "")
//    if (!objfcn.setOutputFile(statusFileName.c_str(), 0) )
//      cerr << "main: output file open failed" << endl;
//
//  objfcn.optimize();
//
//  objfcn.printStatus("Solution from quasi-newton");
//  objfcn.cleanup();
//
//
//  return objfcn.getXPrev();
//}
//
//
//NEWMAT::Matrix const estCovarianceMatrix(NEWMAT::ColumnVector x, 
//				 BaseFunctionObject & f_obj, 
//				 ostream & str)
//{
//  int                  ndim(x.Nrows() );
//  FunctionObject       f(f_obj);
//  global_f_obj     = & f;
//  FDNLF1 nlp(ndim, & f_wrapper, & init_f_wrapper);
//
//  nlp.setX(x);
//  nlp.evalF();
//  NEWMAT::SymmetricMatrix H = nlp.evalH();
//  return -H.i();
//}
//
//
//NEWMAT::Matrix const estGradient(NEWMAT::ColumnVector x, 
//				 BaseFunctionObject & f_obj, 
//				 ostream & str)
//{
//  int                  ndim(x.Nrows() );
//  FunctionObject       f(f_obj);
//  global_f_obj     = & f;
//  FDNLF1 nlp(ndim, & f_wrapper, & init_f_wrapper);
//  
//  nlp.setX(x);
//  nlp.evalF();
//  return nlp.evalG();
//}
