#include "PBDFixedPoint.hpp"
#include <sofa/core/ObjectFactory.h>

//int PBDFixedPointClass = sofa::core::RegisterObject("Constraint that fixe a point.")
//                         .add< PBDFixedPoint<sofa::defaulttype::Vec3Types> >(true)
//                         .add< PBDFixedPoint<sofa::defaulttype::RigidTypes>>();

//template < class T>
//void PBDFixedPoint<T>::solve(PBDObject<T> &object, WriteCoord &p)
//{
//    for(const auto& idx : PBDBaseConstraint<T>::m_indices.getValue ())
//    {
//        p[idx] = object.rest()[idx];
//        std::cout << "Point : "<<p[idx]<<std::endl;
//        std::cout << "Rest Pos : "<<object.rest()[idx]<<std::endl;
//    }
//}


int PBDFixedPointClass = sofa::core::RegisterObject("Constraint that fixes a point")
                         .add< PBDFixedPoint >();
void PBDFixedPoint::solve(PBDObject<sofa::defaulttype::Vec3Types> &object, WriteCoord &p)
{
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = object.rest()[idx];
    }
}


int PBDFixedRigidPointClass = sofa::core::RegisterObject("Constraint that fixes a rigid point")
                         .add< PBDFixedRigidPoint >();
void PBDFixedRigidPoint::solve(PBDObject<sofa::defaulttype::RigidTypes> &object, WriteCoord &p)
{
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = object.rest()[idx];
    }
}
