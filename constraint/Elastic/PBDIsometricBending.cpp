#include "PBDIsometricBending.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDIsometricBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                               .add< PBDIsometricBending >();

void PBDIsometricBending::solve(PBDObject &object, WriteCoord &x)
{

    uint pointCount = x.ref().size();
    if(m_indices.getValue().empty())
    {
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( uint a = 0; a < pointCount; ++a)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : object.topology()[a])
                {

                    //Enforce the isometric property
                    const sofa::defaulttype::Vec3& x_ij = x[a] - x[voisin.first];
                    SReal l = x_ij.norm();
                    const auto& displacement = (0.5*(l-voisin.second)/l) * x_ij;
                    x[a]            -= displacement;
                    x[voisin.first] += displacement;

                    correction (object,a,voisin.first,x);
                }
            }
        }
    }
    else
    {
        const auto& idx = m_indices.getValue ();
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( uint a : idx)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : object.topology()[a])
                {

                    //Enforce the isometric property
                    const sofa::defaulttype::Vec3& x_ij = x[a] - x[voisin.first];
                    SReal l = x_ij.norm();
                    const auto& displacement = (0.5*(l-voisin.second)/l) * x_ij;
                    x[a]            -= displacement;
                    x[voisin.first] += displacement;

                    correction (object,a,voisin.first,x);
                }
            }
        }
    }
}

