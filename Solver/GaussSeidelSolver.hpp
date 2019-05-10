#ifndef GAUSSSEIDELSOLVER_HPP
#define GAUSSSEIDELSOLVER_HPP
#include "config.h"

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/helper/map.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/DiagonalMatrix.h>
#include <sofa/core/behavior/RotationMatrix.h>
#include <math.h>

typedef sofa::defaulttype::Mat3x3 Matrix;
typedef sofa::defaulttype::Vec3f Vector;


namespace sofa{
namespace core{

/**
 * @brief The GaussSeidelSolver class
 * Probably the best GaussSeidel class you'll ever see. It does solve your problems.
 * All you need is a square matrix, a vector and boom bada boom SOLVED
 * For only 99$/minutes you can live a life without any problem
 */
class GaussSeidelSolver : public objectmodel::BaseObject
{
    GaussSeidelSolver(const GaussSeidelSolver& n) ;
    GaussSeidelSolver& operator=(const GaussSeidelSolver& n) ;

public:
    SOFA_CLASS(GaussSeidelSolver, sofa::core::objectmodel::BaseObject);
    GaussSeidelSolver();
    virtual ~GaussSeidelSolver();

    void test() {std::cout << "Gauss-Seidel" << std::endl;}
    /// Solve Mx=b usgin Gauss-Seidel algorithm
    ///
    ///  ^k+1   1          i-1       ^k+1    n        ^k
    /// x_i  = --- ( b_i - sum a_ij x_j   - sum a_ij x_j ), i /in {1,...,n}
    ///        a_ii        j=1              i+1
    void solve (Matrix& M, Vector& x, Vector& b);
public:
    Data<unsigned> m_maxIter; ///< maximum number of iterations of the Conjugate Gradient solution
    ///Ca vavite dégager je pense
    Data<SReal> m_tolerance; ///< desired precision of the Conjugate Gradient Solution (ratio of current residual norm over initial residual norm)
    Data<SReal> m_smallDenominatorThreshold; ///< minimum value of the denominator

};
} //core namespace
} //sofa namespace

#endif //GAUSSSEIDELSOLVER_HPP
