#include "PBDBending.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                      .add< PBDBending >();

void PBDBending::bwdInit ()
{
    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_k = (1.0-std::pow(1.0-m_k.getValue (),1.0 / ((double)m_nbIter.getValue ()))) * node->getDt () * node->getDt ();
}

void PBDBending::solve(PBDObject<sofa::defaulttype::Vec3Types> &object, WriteCoord &x)
{
    if( !object.hasDataType (PBDObject<sofa::defaulttype::Vec3Types>::BENDING) || !object.hasDataType (PBDObject<sofa::defaulttype::Vec3Types>::STRETCH) )
    {
        if( !object.hasDataType (PBDObject<sofa::defaulttype::Vec3Types>::STRETCH) )
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
                    correction(object,a,voisin.first,x,vel);
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
                for( const auto& voisin : object.topology().data ()[a])
                {
                    correction(object,a,voisin.first,x,vel);
                }
            }
        }
    }
}


void PBDBending::correction (PBDObject<sofa::defaulttype::Vec3Types> &object, uint a, uint b, WriteCoord &p, const ReadDeriv& vel)
{
    static SReal eps = 1e-6;
    uint edge_ID = object.sofaTopology ()->getEdgeIndex(a,b);
    const auto& hessian_and_idx = object.bendTopology ().bendingData ()[edge_ID];
    //Apply correction on the points
    for(const auto& data : hessian_and_idx)
    {
        const Vec3 *x[4] = { &p[a], &p[b], &p[data.first[0]], &p[data.first[1]] };
        Real invMass[4] = { object.invMass (a),object.invMass (b), object.invMass (data.first[0]), object.invMass (data.first[1])};

        Real energy = 0.0;
        for (unsigned char k = 0; k < 4; k++)
            for (unsigned char j = k; j < 4; j++)
                energy += data.second(j, k)*(dot(*x[k],*x[j]));

        Vec3 gradC[4];
        gradC[0]  =  data.second(0,0) * *x[0];
        gradC[0] +=  data.second(0,1) * *x[1];
        gradC[0] +=  data.second(0,2) * *x[2];
        gradC[0] +=  data.second(0,3) * *x[3];

        gradC[1]  =  data.second(0,1) * *x[0];
        gradC[1] +=  data.second(1,1) * *x[1];
        gradC[1] +=  data.second(1,2) * *x[2];
        gradC[1] +=  data.second(1,3) * *x[3];

        gradC[2]  =  data.second(0,2) * *x[0];
        gradC[2] +=  data.second(1,2) * *x[1];
        gradC[2] +=  data.second(2,2) * *x[2];
        gradC[2] +=  data.second(2,3) * *x[3];

        gradC[3]  =  data.second(0,3) * *x[0];
        gradC[3] +=  data.second(1,3) * *x[1];
        gradC[3] +=  data.second(2,3) * *x[2];
        gradC[3] +=  data.second(3,3) * *x[3];

        Real sum_normGradC =  invMass[0] * gradC[0].norm2()
                              +invMass[1] * gradC[1].norm2()
                              +invMass[2] * gradC[2].norm2()
                              +invMass[3] * gradC[3].norm2();

        // exit early if required
        if (sum_normGradC > eps)
        {
            // compute impulse-based scaling factor
            const Real s = energy / sum_normGradC;

            p[a]             += (m_k.getValue() * s * invMass[0]) * gradC[0];
            p[b]             += (m_k.getValue() * s * invMass[1]) * gradC[1];
            p[data.first[0]] += (m_k.getValue() * s * invMass[2]) * gradC[2];
            p[data.first[1]] += (m_k.getValue() * s * invMass[3]) * gradC[3];


        }
    }
}
