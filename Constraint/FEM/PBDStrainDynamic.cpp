#include "PBDStrainDynamic.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDStrainDynamicClass = sofa::core::RegisterObject("Constraint that correct deformed Dynamic.")
                            .add< PBDStrainDynamic >();



void PBDStrainDynamic::bwdInit ()
{
    SReal p = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    // THIS IS NOT THE MATRIX SHOWN IN THE PAPER (PBDStrainDynamic.hpp)
    // THIS IS SAINT VENANT KIRCHOFF MODEL

    SReal lambda = (E*p) / (( 1.0 + p ) * ( 1.0 - 2.0 * p ));
    SReal mu = E / (2.0 * ( 1.0 + p ));
    //C(i,i) defines x/y/z stretching while K[3] define the intensity of the shear motion
    m_C(0,0) = lambda + 2.0 * mu  ; m_C(0,1) = lambda         ; m_C(0,2) = lambda         ;
    m_C(1,0) = lambda           ; m_C(1,1) = lambda + 2 * mu; m_C(1,2) = lambda         ;
    m_C(2,0) = lambda           ; m_C(2,1) = lambda         ; m_C(2,2) = lambda + 2 * mu;

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
}

void PBDStrainDynamic::solve(PBDObject &object, WriteCoord &x)
{

    if(!object.hasDataType(PBDObject::TETRAHEDRON))
    {
        object.computeTetrahedraBasis ();
    }

    const auto& Dm_inv = object.tetrahedraBases ().data ();
    const uint tetCount = object.sofaTopology ()->getNbTetrahedra ();
    static Eigen::Matrix3d I; I.setIdentity ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint i = 0; i < tetCount; ++i)
        {
            const auto& t = object.sofaTopology ()->getTetra (i);
            //First compute Ds to get F
            const auto& r1 = x[t[0]] - x[t[3]];
            const auto& r2 = x[t[1]] - x[t[3]];
            const auto& r3 = x[t[2]] - x[t[3]];
            Eigen::Matrix3d Ds;
            Ds << r1[0],r1[1],r1[2],
                    r2[0],r2[1],r2[2],
                    r3[0],r3[1],r3[2];
            auto F = (Ds*Dm_inv[i].second).eval();

            //From F compute Green Strain Tensor { 0.5 (F^t * F - I) }
            auto T = (F.transpose() * F).eval();
            Eigen::Matrix3d S;
            S <<std::sqrt(T(0,0)), T(0,1), T(0,2),
                    T(1,0), std::sqrt(T(1,1)), T(1,2),
                    T(2,0), T(2,1), std::sqrt(T(2,2));
            const Eigen::Matrix3d& GST = ( S - I ).eval(); //What happend if F is orthogonal ? :-> if F is orthogonal then Ds is a pure rotation of Dm hence no strain is applied

            //compute Piola-Kirchhoff stress tensor
            const Eigen::Matrix3d& P = (F*m_C*GST).eval();

            //Then the gradient and divergence (which are actually displacement since m_C is premultiplied by dt^2)
            const Eigen::Matrix3d& gradEs = Dm_inv[i].first * (P * Dm_inv[i].second.transpose());
            const Eigen::Matrix3d& displacement = (dt2*gradEs).eval();
            //Update positions
            x[t[0]][0] += displacement(0,0); x[t[0]][1] += displacement(0,1); x[t[0]][2] += displacement(0,2);
            x[t[1]][0] += displacement(1,0); x[t[1]][1] += displacement(1,1); x[t[1]][2] += displacement(1,2);
            x[t[2]][0] += displacement(2,0); x[t[2]][1] += displacement(2,1); x[t[2]][2] += displacement(2,2);
            //The displacement is minus the divergence, hence the -=
            x[t[3]][0] -= displacement(0,0) + displacement(1,0) + displacement(2,0);
            x[t[3]][1] -= displacement(0,1) + displacement(1,1) + displacement(2,1);
            x[t[3]][2] -= displacement(0,2) + displacement(1,2) + displacement(2,2);
        }
    }
}

