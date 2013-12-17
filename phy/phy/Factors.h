#ifndef __Factors_h
#define __Factors_h

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

  class AbstractBaseFactor; // forward declaration

  /** Pointer to base factors */
  typedef boost::shared_ptr<AbstractBaseFactor> AbsBasFacPtr_t;

  class AbstractBaseFactor {

  public:

    virtual ~AbstractBaseFactor() {};

    /** Constructor for factor matrices. Parameters are defined in derived classes.*/
    AbstractBaseFactor(string const & type, string const & name, unsigned factorMatrixSize1, unsigned factorMatrixSize2) 
      : size1_(factorMatrixSize1), 
	size2_(factorMatrixSize2), 
	hasChanged_(false),
	type_(type),
	name_(name)
    {initCounts();}
    
    /** submit expectations for matrix entries of each factor. */
    void submitCounts(matrix_t const & counts) {counts_ += counts; hasChanged_ = true;}
										
    /** return zero 0 on failure and >0 on success. Success values may
	optionally report on other aspects of optimization. NOTE:
	derived classes should implement optimizeParametersImpl and
	NOT this function. */
    int optimizeParameters() {if (hasChanged_) {hasChanged_ = false; return optimizeParametersImpl();} else return 2;}

    /** return or set matrix defining factor */
    matrix_t mkFactor() const {matrix_t m(size1_, size2_); mkFactor(m); return m;}
    virtual void mkFactor(matrix_t & m) const = 0;

    /** clear all submitted expectation counts */
    void clearCounts() {reset(counts_);}

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

    /** stores expectation counts */
    matrix_t counts_;

    /** counts has changed -- optimization needed. */
    bool hasChanged_;

    /** Type string */
    string const type_;

    /** Name string */
    string const name_;
  };


  class AbstractFullyParameterizedFactor : public AbstractBaseFactor {

  public:

    virtual ~AbstractFullyParameterizedFactor() {};

    /** Constructor for fully parametrized factor matrices. One
	parameter value is given for each matrix entry in matrix
	m. Optionally, one pseudoCount value can be given for each
	parameter in the matrix pseudoCount. The pseudoCount are added
	to the submitted counts before optimization. */
    AbstractFullyParameterizedFactor(string const & type, string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t() );

    virtual void mkFactor(matrix_t & m) const {m = m_;}
    using AbstractBaseFactor::mkFactor; // bringing other mkFactor definition into this name space (hidden otherwise)

  protected:
    matrix_t m_;
    matrix_t const pseudoCounts_;
    friend   void writeAbstractFullyParameterizedFactor(ostream & str, AbsBasFacPtr_t const & factorPtr);
  };


  /** matrix normalized to sum to one */
  class GlobalNormFactor : public AbstractFullyParameterizedFactor {
  public:

    virtual ~GlobalNormFactor() {};
    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    GlobalNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractFullyParameterizedFactor("globNorm", name, m, pseudoCounts) {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  /** matrix columns normalized to sum to one */
  class ColumnNormFactor : public AbstractFullyParameterizedFactor  {
  public: 

    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    ColumnNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractFullyParameterizedFactor("colNorm", name, m, pseudoCounts) {};
    virtual ~ColumnNormFactor() {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  /** matrix rows normalized to sum to one */
  class RowNormFactor : public AbstractFullyParameterizedFactor  {
  public: 

    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    RowNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractFullyParameterizedFactor("rowNorm", name, m, pseudoCounts) {};
    virtual ~RowNormFactor() {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  class AbstractBaseFactorSet {

  public:

    AbstractBaseFactorSet(unsigned facCount) : facCount_(facCount) {};
    virtual ~AbstractBaseFactorSet() {};

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


  class CompositeFactorSet : public AbstractBaseFactorSet {
  public:

    CompositeFactorSet(vector<AbsBasFacPtr_t> const factorPtrs) : AbstractBaseFactorSet( factorPtrs.size() ), factorPtrs_(factorPtrs) {}
    virtual ~CompositeFactorSet() {};

    /** submit expectations for matrix entries of each factor. */
    virtual void submitCounts(matrix_t const & counts, unsigned idx) {factorPtrs_[idx]->submitCounts(counts);}
    using AbstractBaseFactorSet::submitCounts; // bringing other definitions into this name space (hidden otherwise)

    /** return value may optionally report on success or other aspects of optimization. */
    virtual int optimizeParameters() {int success = 1; for (unsigned i = 0; i < facCount_; i++) success *= factorPtrs_[i]->optimizeParameters(); return success;}

    /** return or set matrix defining factor idx  */
    virtual matrix_t mkFactor(unsigned idx) const {return factorPtrs_[idx]->mkFactor();}
    virtual void mkFactor(matrix_t & m, unsigned idx) const  {return factorPtrs_[idx]->mkFactor(m);}

    using AbstractBaseFactorSet::mkFactor; // bringing other mkFactor definition into this name space (hidden otherwise)

    /** clear all submitted expectation counts */
    virtual void clearCounts() {for (unsigned i = 0; i < facCount_; i++) factorPtrs_[i]->clearCounts();}

  protected:

    vector<AbsBasFacPtr_t> const factorPtrs_; // smart pointers, defined above
  };


} // end namespace phy

#endif  // __Factors_h