#include "PBDIsometricBending.hpp"

#include <sofa/core/ObjectFactory.h>

int PBDIsometricBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                               .add< PBDIsometricBending >();

PBDIsometricBending::PBDIsometricBending(sofa::simulation::Node *gnode)
    : PBDBaseConstraint(true),
      m_k(initData(&m_k,(SReal)1.0,"bend","Bending factor"))
{
    if(!gnode)
    {
        dt2 = 0.0001;
        return;
    }
    dt2 = gnode->getContext ()->getDt ();
    dt2 *= dt2;
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

void PBDIsometricBending::solve(PBDObject &object, WriteCoord &x)
{

    uint pointCount = x.ref().size();
    float k = m_k.getValue ();
    Eigen::Matrix4d Qt;
    std::pair<float,float> area(0,0);
    SReal coeff = k * dt2;

    if(m_indices.getValue().empty())
    {
        for( uint a = 0; a < pointCount; ++a)
        {
            //Get the edge of the corresponding neighbors
            for( const auto& voisin : object.topology()[a])
            {
                uint edge_ID = object.sofaTopology ()->getEdgeIndex(a,voisin.first);
                const auto& hessian_and_idx = object.bend_topology ()[edge_ID];
                //Apply correction on the points
                for(const auto& data : hessian_and_idx)
                {
                    const sofa::defaulttype::Vec3 *y[4] = {&(x[a]),&(x[voisin.first]),&(x[data.first[0]]),&(x[data.first[1]])};
                    object.computeQ(y,Qt,area);
                    Qt = Qt - data.second;

                    sofa::defaulttype::Vec3 dx(0,0,0);
                    dx  =  Qt(2,0) * x[a];
                    dx +=  Qt(2,1) * x[voisin.first];
                    dx +=  Qt(2,2) * x[data.first[0]];
                    dx +=  Qt(2,3) * x[data.first[1]];
                    sofa::defaulttype::Vec3 dx1(0,0,0);
                    dx1  =  Qt(3,0) * x[a];
                    dx1 +=  Qt(3,1) * x[voisin.first];
                    dx1 +=  Qt(3,2) * x[data.first[0]];
                    dx1 +=  Qt(3,3) * x[data.first[1]];
                    x[data.first[0]] += coeff * dx;
                    x[data.first[1]] += coeff * dx1;
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

////                    sofa::defaulttype::Vec3 dx1(0,0,0);
////                    dx1 +=  Qt(0,0) * x[idx[0]];
////                    dx1 +=  Qt(0,1) * x[idx[1]];
////                    dx1 +=  Qt(0,2) * x[idx[2]];
////                    dx1 +=  Qt(0,3) * x[idx[3]];
////                    sofa::defaulttype::Vec3 dx2(0,0,0);
////                    dx2 +=  Qt(1,0) * x[idx[0]];
////                    dx2 +=  Qt(1,1) * x[idx[1]];
////                    dx2 +=  Qt(1,2) * x[idx[2]];
////                    dx2 +=  Qt(1,3) * x[idx[3]];
//                    sofa::defaulttype::Vec3 dx3(0,0,0);
//                    dx3 +=  Qt(2,0) * x[idx[0]];
//                    dx3 +=  Qt(2,1) * x[idx[1]];
//                    dx3 +=  Qt(2,2) * x[idx[2]];
//                    dx3 +=  Qt(2,3) * x[idx[3]];
//                    sofa::defaulttype::Vec3 dx4(0,0,0);
//                    dx4 +=  Qt(3,0) * x[idx[0]];
//                    dx4 +=  Qt(3,1) * x[idx[1]];
//                    dx4 +=  Qt(3,2) * x[idx[2]];
//                    dx4 +=  Qt(3,3) * x[idx[3]];
//                    std::cout << coeff * dx3 << std::endl;
////                    x[idx[0]] += coeff * dx1;
////                    x[idx[1]] += coeff * dx2;
//                    x[idx[2]] += coeff * dx3;
//                    x[idx[3]] += coeff * dx4;
