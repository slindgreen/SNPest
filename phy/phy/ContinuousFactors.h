#ifndef __ContinuousFactors_h
#define __ContinuousFactors_h

#include "phy/PhyDef.h"
#include "phy/utils.h"
#include "phy/utilsLinAlg.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

/** Classes for specifying and optimizing discrete factor potential
    (matrices). Potential matrices are specifies in terms of
    parameters. Complete implementation for fully parameterized,
    normalized potentials are given. All classes allow sets of
    both variable and fixed parameters to be specified. */

/** To do: In some cases it would be convenient to be able to fix a
    subset of the otherwise free parameters in the optimization. This
    could be specified by a vector of bools, which the
    optimizationImpl could query when performing the
    optimization. */

namespace phy {

  using namespace std;

  class AbstractBaseContinuousFactor; // forward declaration

  /** Pointer to base factors */
  typedef boost::shared_ptr<AbstractBaseContinuousFactor> AbsBasFacPtr_t;

  class AbstractBaseContinuousFactor {

  public:

    virtual ~AbstractBaseContinuousFactor() {};

    /** Constructor for factor matrices. Parameters are defined in derived classes.*/
    AbstractBaseContinuousFactor(string const & type, string const & name, unsigned factorMatrixSize1, unsigned factorMatrixSize2) 
      : size1_(factorMatrixSize1), 
	size2_(factorMatrixSize2), 
	hasChanged_(false),
	type_(type),
	name_(name)
    {initCounts();}
    
    /** Submit the sufficient statistics needed for the optimization of the continuous distributions. These are the observed data (x) and the expectations for the neighboring discrete random variable. */
    void submitCounts(vector<number_t> const & observations, matrix_t const & counts) {counts_ += counts; observations_.append(observations); hasChanged_ = true;}

    /** return zero 0 on failure and >0 on success. Success values may
	optionally report on other aspects of optimization. NOTE:
	derived classes should implement optimizeParametersImpl and
	NOT this function. */
    int optimizeParameters() {if (hasChanged_) {hasChanged_ = false; return optimizeParametersImpl();} else return 2;}

    /** return or set matrix defining factor */
    matrix_t mkFactor(number_t const x) const {matrix_t m(size1_, size2_); mkFactor(x, m); return m;}
    virtual void mkFactor(number_t const x, matrix_t & m) const = 0;

    /** clear all submitted expectation counts */
    void clearCounts() {observations_.clear(); reset(counts_);}

    /** returns factor type */
    string const & type() {return type_;}

    /** returns factor name */
    string const & name() {return name_;}

  protected:

    /** Optimize parameters implementation. Called by optimizeParameters and must be implemented by derived classes. */
    virtual int optimizeParametersImpl() = 0;  
   
    /** intialize the count matrix to the right size and all zero entries. */
    void initCounts() {counts_.resize(size1_, size2_); clearCounts();}

    /** Defines sizes of factor matrix */
    unsigned size1_;  // rows
    unsigned size2_;  // columns 

    /** stores the observations */
    vector<number_t> observations_;

    /** stores expectation counts */
    matrix_t counts_;

    /** counts has changed -- optimization needed. */
    bool hasChanged_;

    /** Type string */
    string const type_;

    /** Name string */
    string const name_;
  };


  class AbstractParameterizedContinuousFactor : public AbstractBaseContinuousFactor {

  public:

    virtual ~AbstractParameterizedContinuousFactor() {};

    /** Constructor for parametrized continuous distributions. Note that pseudocounts need to be replaced with some other type of (continuous) prior. */
    AbstractParameterizedContinuousFactor(string const & type, string const & name, matrix_t const & m, vector<vector<string>> const & parameterNames, vector<vector<number_t>> const & parameterValues, matrix_t const & pseudoCounts = matrix_t() );

  protected:
    vector<vector<string> > parameterNames;
    vector<vector<number_t> > parameterValues;
    matrix_t const pseudoCounts_;
    
    friend   void writeAbstractParameterizedContinuousFactor(ostream & str, AbsBasFacPtr_t const & factorPtr);
  };



  /** matrix normalized to sum to one */
  class GausianContinuousFactor : public AbstractParameterizedContinuousFactor {
  public:

    virtual ~GausianContinuousFactor() {};
    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    GlobalNormFactor(string const & name, matrix_t const & m, vector<vector<string> > const & parameterNames, vector<vector<number_t> > const & parameterValues, matrix_t const & pseudoCounts = matrix_t()) : AbstractParameterizedContinuousFactor("globNorm", name, m, parameterNames, parameterValues, pseudoCounts) {};

    void mkFactor(number_t const & x, matrix_t & m) const; // this implements the various Gaussians [ P(x|z_1) = g(x | u_1, s_1) ]
    using AbstractBaseContinuousFactor::mkFactor; // bringing other mkFactor definition into this name space (hidden otherwise)

  protected:
    int optimizeParametersImpl(); // this should implement optimization of the parameters defining the gaussian distribution(s) given the observations_ and the counts_;   
  };


  /** matrix columns normalized to sum to one */
  class ColumnNormFactor : public AbstractParameterizedContinuousFactor  {
  public: 

    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    ColumnNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractParameterizedContinuousFactor("colNorm", name, m, pseudoCounts) {};
    virtual ~ColumnNormFactor() {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  /** matrix rows normalized to sum to one */
  class RowNormFactor : public AbstractParameterizedContinuousFactor  {
  public: 

    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    RowNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractParameterizedContinuousFactor("rowNorm", name, m, pseudoCounts) {};
    virtual ~RowNormFactor() {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  class AbstractBaseContinuousFactorSet {

  public:

    AbstractBaseContinuousFactorSet(unsigned facCount) : facCount_(facCount) {};
    virtual ~AbstractBaseContinuousFactorSet() {};

    /** submit expectations for matrix entries of each factor. */
    void submitCounts(vector<matrix_t> const & countsVec) {for (unsigned i = 0; i < facCount_; i++) submitCounts(countsVec[i], i);}
    virtual void submitCounts(matrix_t const & counts, unsigned idx) = 0;  // idx is the factor index within the class

    /** return value may optionally report on success or other aspects of optimization. */
    virtual int optimizeParameters() = 0;

    /** return or set matrix defining factor idx  */
    virtual matrix_t mkFactor(unsigned idx) const = 0; 
    virtual void mkFactor(matrix_t & m, unsigned idx) const = 0;

    /** return or set vector of all factor matrices */
    virtual vector<matrix_t> mkFactorVec() const;
    virtual void mkFactorVec(vector<matrix_t> & v) const;

    /**     clear all submitted expectation counts */
    virtual void clearCounts() = 0;
    
    /** Return factor count */
    unsigned facCount() const {return facCount_;}

  protected:
    
    unsigned facCount_; // number of factors in set
  };


  class CompositeFactorSet : public AbstractBaseContinuousFactorSet {
  public:

    CompositeFactorSet(vector<AbsBasFacPtr_t> const factorPtrs) : AbstractBaseContinuousFactorSet( factorPtrs.size() ), factorPtrs_(factorPtrs) {}
    virtual ~CompositeFactorSet() {};

    /** submit expectations for matrix entries of each factor. */
    virtual void submitCounts(matrix_t const & counts, unsigned idx) {factorPtrs_[idx]->submitCounts(counts);}
    using AbstractBaseContinuousFactorSet::submitCounts; // bringing other definitions into this name space (hidden otherwise)

    /** return value may optionally report on success or other aspects of optimization. */
    virtual int optimizeParameters() {int success = 1; for (unsigned i = 0; i < facCount_; i++) success *= factorPtrs_[i]->optimizeParameters(); return success;}

    /** return or set matrix defining factor idx  */
    virtual matrix_t mkFactor(unsigned idx) const {return factorPtrs_[idx]->mkFactor();}
    virtual void mkFactor(matrix_t & m, unsigned idx) const  {return factorPtrs_[idx]->mkFactor(m);}

    using AbstractBaseContinuousFactorSet::mkFactor; // bringing other mkFactor definition into this name space (hidden otherwise)

    /** clear all submitted expectation counts */
    virtual void clearCounts() {for (unsigned i = 0; i < facCount_; i++) factorPtrs_[i]->clearCounts();}

  protected:

    vector<AbsBasFacPtr_t> const factorPtrs_; // smart pointers, defined above
  };


} // end namespace phy

#endif  // __ContinuousFactors_h
