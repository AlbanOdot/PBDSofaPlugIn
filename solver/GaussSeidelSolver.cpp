#include "GaussSeidelSolver.hpp"
#include <sofa/core/ObjectFactory.h>

using namespace sofa;
using namespace core;

int GaussSeidelSolverClass = RegisterObject("A simple Gausse Seidel solver")
                             .add< GaussSeidelSolver >()
                             .addAlias("GSSolver")
                             ;

void GaussSeidelSolver::solve (Matrix& M, Vector& x, Vector& b){

    //Prepare the matrix in case some elements on the diagonal are too close to 0
    uint n = M.getNbCols ();//In case the type changes
    for( uint i = 0; i < n; ++i )
    {
        M[i][i] = std::abs( M[i][i] ) - m_smallDenominatorThreshold.getValue() < 0 ? m_smallDenominatorThreshold.getValue() : M[i][i];
    }

    for( uint nbIter = 0; nbIter < m_maxIter.getValue(); ++nbIter )
    {
        for( uint i = 0; i < n; ++i )
        {
            double newX = 0.0;
            uint j = 0;
            for(    ; j<i; ++j )
            {
                newX += M[i][j] * x[j];
            }
            for( ++j; j < n; ++j )
            {
                newX += M[i][j] * x[j];
            }
            x[i] = ( b[i] - newX ) / M[i][i];
        }

        Vector dist = b-x;
        if ( dist.norm2()< m_tolerance.getValue ())
            break;
    }
}

GaussSeidelSolver::GaussSeidelSolver(){
    m_tolerance.setValue (1e-5);
    m_maxIter.setValue(25);
    m_smallDenominatorThreshold.setValue(1e-5);
}

GaussSeidelSolver::~GaussSeidelSolver (){
}
