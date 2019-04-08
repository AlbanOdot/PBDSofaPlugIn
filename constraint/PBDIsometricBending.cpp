#include "PBDIsometricBending.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDIsometricBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                               .add< PBDIsometricBending >();

PBDIsometricBending::PBDIsometricBending(sofa::simulation::Node *gnode)
    : PBDBaseConstraint(true),
      m_k(initData(&m_k,(SReal)1.0,"bend","Bending factor"))
{

}

PBDIsometricBending::PBDIsometricBending(unsigned int objectSize)
    : PBDBaseConstraint(true),
      m_k(initData(&m_k,(SReal)1.0,"bend","Bending factor"))
{

}

sofa::defaulttype::BaseMatrix * PBDIsometricBending::getConstraintMatrix ()
{
    return nullptr;
}

void PBDIsometricBending::bwdInit ()
{
    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    coeff = m_k.getValue () * node->getDt () * node->getDt ();
}

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


                    //Enforce the bending property
                    uint edge_ID = object.sofaTopology ()->getEdgeIndex(a,voisin.first);
                    const auto& hessian_and_idx = object.bend_topology ()[edge_ID];
                    //Apply correction on the points
                    for(const auto& data : hessian_and_idx)
                    {
                        sofa::defaulttype::Vec3 dx(0,0,0);
                        dx  =  data.second(0,0) * x[a];
                        dx +=  data.second(0,1) * x[voisin.first];
                        dx +=  data.second(0,2) * x[data.first[0]];
                        dx +=  data.second(0,3) * x[data.first[1]];

                        sofa::defaulttype::Vec3 dx1(0,0,0);
                        dx1  =  data.second(0,1) * x[a];
                        dx1 +=  data.second(1,1) * x[voisin.first];
                        dx1 +=  data.second(1,2) * x[data.first[0]];
                        dx1 +=  data.second(1,3) * x[data.first[1]];

                        sofa::defaulttype::Vec3 dx2(0,0,0);
                        dx2  =  data.second(0,2) * x[a];
                        dx2 +=  data.second(1,2) * x[voisin.first];
                        dx2 +=  data.second(2,2) * x[data.first[0]];
                        dx2 +=  data.second(2,3) * x[data.first[1]];

                        sofa::defaulttype::Vec3 dx3(0,0,0);
                        dx3  =  data.second(0,3) * x[a];
                        dx3 +=  data.second(1,3) * x[voisin.first];
                        dx3 +=  data.second(2,3) * x[data.first[0]];
                        dx3 +=  data.second(3,3) * x[data.first[1]];

                        x[a]             += coeff * dx ;
                        x[voisin.first]  += coeff * dx1;
                        x[data.first[0]] += coeff * dx2;
                        x[data.first[1]] += coeff * dx3;
                    }
                }
            }
        }
    }
    else
    {
        const auto& idx = m_indices.getValue ();
        //#pragma omp parallel for
        for( uint i = 0; i < idx.size(); ++i)
        {
            //correction(idx[i]);
        }
    }
}

