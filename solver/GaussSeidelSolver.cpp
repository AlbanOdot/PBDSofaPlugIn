#include "GaussSeidelSolver.hpp"

using namespace sofa::core;

void GaussSeidelSolver::init(){

    if(m_maxIter.getValue() < 0)
    {
        msg_warning() << "'iterations' must be a positive value" << msgendl
                      << "default value used: 25";
        m_maxIter.setValue(25);
    }

    if(m_smallDenominatorThreshold.getValue() < 0.0)
    {
        msg_warning() << "'threshold' must be a positive value" << msgendl
                      << "default value used: 1e-5";
        m_smallDenominatorThreshold.setValue(1e-5);
    }
}

void GaussSeidelSolver::reinit(){
    init();
}

void GaussSeidelSolver::solve (Matrix& M, Vector& x, Vector& b){

    //Prepare the matrix in case some elements on the diagonal are to close to 0
    for( uint i = 0; i < M.rowSize(); ++i )
    {
        M[i][i] = std::abs( M[i][i] ) - m_smallDenominatorThreshold.getValue() < 0 ? m_smallDenominatorThreshold.getValue() : M[i][i];
    }

    uint n = M.rowSize();

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
    }
}

bool GaussSeidelSolver::insertInNode( objectmodel::BaseNode* node )
{
    node->addLinearSolver(this);
    Inherit1::insertInNode(node);
    return true;
}

bool GaussSeidelSolver::removeInNode( objectmodel::BaseNode* node )
{
    node->removeLinearSolver(this);
    Inherit1::removeInNode(node);
    return true;
}

