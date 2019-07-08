#include "PBDIsometricBending.hpp"

#include <sofa/core/ObjectFactory.h>

int PBDIsometricBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                               .add< PBDIsometricBending >();

void PBDIsometricBending::solve(PBDObject<sofa::defaulttype::Vec3Types> &object, WriteCoord &x)
{

    if( !object.hasDataType(PBDObject<sofa::defaulttype::Vec3Types>::BENDING) || !object.hasDataType (PBDObject<sofa::defaulttype::Vec3Types>::STRETCH) )
    {
        if( !object.hasDataType(PBDObject<sofa::defaulttype::Vec3Types>::STRETCH) )
            object.computeStretchTopology ();
        object.computeBendingTopology ();
    }
    uint pointCount = x.ref().size();
    const auto& vel = object.object()->readVelocities ();

    if(m_indices.getValue().empty())
    {
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( uint a = 0; a < pointCount; ++a)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : object.topology().data ()[a])
                {
                    //Enforce the isometric property
                    const sofa::defaulttype::Vec3& x_ij = x[a] - x[voisin.first];
                    SReal l = x_ij.norm();
                    const auto& dx = (0.5*(l-voisin.second)/l) * x_ij;
                    x[a]            -= dx;
                    x[voisin.first] += dx;
                    correction (object,a,voisin.first,x,vel);//PBDBending
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
                for( const auto& voisin : object.topology().data ()[a])
                {

                    //Enforce the isometric property
                    const sofa::defaulttype::Vec3& x_ij = x[a] - x[voisin.first];
                    SReal l = x_ij.norm();
                    const auto& displacement = (0.5*(l-voisin.second)/l) * x_ij;
                    x[a]            -= displacement;
                    x[voisin.first] += displacement;
                    correction (object,a,voisin.first,x,vel);//PBDBending
                }
            }
        }
    }
}

