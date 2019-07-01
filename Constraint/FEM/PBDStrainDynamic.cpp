#include "PBDStrainDynamic.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDStrainDynamicClass = sofa::core::RegisterObject("Constraint that correct deformed Dynamic.")
                            .add< PBDStrainDynamic >();

using namespace sofa::defaulttype;

void PBDStrainDynamic::solve(PBDObject<sofa::defaulttype::Vec3Types> &object, WriteCoord &x)
{
    static SReal eps = 1e-6;
    if(!object.hasDataType(PBDObject<sofa::defaulttype::Vec3Types>::TETRAHEDRON))
    {
        object.computeTetrahedraBasis ();
    }
    const auto& Dm_inv = object.tetrahedraBases ().data ();
    const uint tetCount = object.sofaTopology ()->getNbTetrahedra ();

    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint i = 0; i < tetCount; ++i)
        {
            const auto& t = object.sofaTopology ()->getTetra (i);
            //First compute Ds to get F
            const auto& r1 = x[t[0]] - x[t[3]];
            const auto& r2 = x[t[1]] - x[t[3]];
            const auto& r3 = x[t[2]] - x[t[3]];
            Matrix3 Ds;
            Ds.x () = r1;
            Ds.y () = r2;
            Ds.z () = r3;

            const Matrix3& F = Ds*Dm_inv[i].second;

            //From F compute Green Strain Tensor { 0.5 (F^t * F - I) }
            Matrix3 epsilon = 0.5 * ( F.transposed() * F); //What happend if F is orthogonal ? :-> if F is orthogonal then Ds is a pure rotation of Dm hence no strain is applied
            epsilon(0,0) = epsilon(0,0) - 0.5;epsilon(1,1) = epsilon(1,1) - 0.5 ;epsilon(2,2) = epsilon(2,2) - 0.5;

            //compute Piola-Kirchhoff stress tensor
            const Matrix3& P = F*m_C*epsilon;

            //Then the gradient and divergence (which are actually displacement since m_C is premultiplied by dt^2)
            const Matrix3& displacement = dt2 * Dm_inv[i].first * (P * Dm_inv[i].second.transposed());

            x[t[0]] += object.invMass(t[0]) * displacement.x ();
            x[t[1]] += object.invMass(t[1]) * displacement.y ();
            x[t[2]] += object.invMass(t[2]) * displacement.z ();
            x[t[3]] -= object.invMass(t[3]) * displacement.x () + displacement.y () + displacement.z ();
        }
    }
}

