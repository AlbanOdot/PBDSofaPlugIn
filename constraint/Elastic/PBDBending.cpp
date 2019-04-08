#include "PBDBending.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                      .add< PBDBending >();

void PBDBending::bwdInit ()
{
    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    double m_K = 1.0-std::pow(1.0-m_k.getValue (),1.0 / ((double)m_nbIter.getValue ()));
    coeff = m_K * node->getDt () * node->getDt ();
}

void PBDBending::solve(PBDObject &object, WriteCoord &x)
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
                    correction(object,a,voisin.first,x);
                }
            }
        }
    }
    else
    {
        const auto& idx = m_indices.getValue ();
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( const auto& a : idx)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : object.topology()[a])
                {
                    correction(object,a,voisin.first,x);
                }
            }
        }
    }
}


void PBDBending::correction (PBDObject &object, uint a, uint b, WriteCoord &x)
{
    uint edge_ID = object.sofaTopology ()->getEdgeIndex(a,b);
    const auto& hessian_and_idx = object.bend_topology ()[edge_ID];
    //Apply correction on the points
    for(const auto& data : hessian_and_idx)
    {
        sofa::defaulttype::Vec3 dx(0,0,0);
        dx  =  data.second(0,0) * x[a];
        dx +=  data.second(0,1) * x[b];
        dx +=  data.second(0,2) * x[data.first[0]];
        dx +=  data.second(0,3) * x[data.first[1]];

        sofa::defaulttype::Vec3 dx1(0,0,0);
        dx1  =  data.second(0,1) * x[a];
        dx1 +=  data.second(1,1) * x[b];
        dx1 +=  data.second(1,2) * x[data.first[0]];
        dx1 +=  data.second(1,3) * x[data.first[1]];

        sofa::defaulttype::Vec3 dx2(0,0,0);
        dx2  =  data.second(0,2) * x[a];
        dx2 +=  data.second(1,2) * x[b];
        dx2 +=  data.second(2,2) * x[data.first[0]];
        dx2 +=  data.second(2,3) * x[data.first[1]];

        sofa::defaulttype::Vec3 dx3(0,0,0);
        dx3  =  data.second(0,3) * x[a];
        dx3 +=  data.second(1,3) * x[b];
        dx3 +=  data.second(2,3) * x[data.first[0]];
        dx3 +=  data.second(3,3) * x[data.first[1]];

        x[a]             += coeff * dx ;
        x[b]  += coeff * dx1;
        x[data.first[0]] += coeff * dx2;
        x[data.first[1]] += coeff * dx3;
    }
}