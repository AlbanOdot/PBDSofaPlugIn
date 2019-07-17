#ifndef PBDOBJECT_HPP
#define PBDOBJECT_HPP

#include <SofaBaseMechanics/MechanicalObject.h>

#include "./Data/PBDVertexMass.hpp"
#include "./Data/PBDOrientation.hpp"

#include "./Data/Elastic/PBDVertexTopology.hpp"
#include "./Data/Elastic/PBDBendingTopology.hpp"

#include "./Data/FEM/PBDTetrahedronBasis.hpp"
#include "./Data/FEM/PBDStiffRodData.hpp"
#include "./Data/FEM/PBDElasticRodData.hpp"
#include "./Data/FEM/PDCosseratRodData.hpp"

template < class T >
class PBDObject
{
    typedef sofa::defaulttype::Vec3 Vec3;

    typedef typename T::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef typename T::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

public:

    PBDObject(sofa::component::container::MechanicalObject<T> * obj = nullptr) : m_obj(obj)
    {
        m_freePos.setValue (obj->x.getValue());
        m_freeVel.setValue (obj->v.getValue());
    }

    inline void resetFreePosition();
    inline void resetFreeVelocity();

    inline WriteCoord getFreePosition()    { return WriteCoord(m_freePos);}
    inline WriteDeriv getFreeVelocity()    { return WriteDeriv(m_freeVel);}

private:
    sofa::component::container::MechanicalObject<T> * m_obj;
    Coordinates m_freePos;
    Derivatives m_freeVel;
};

template class PBDObject<sofa::defaulttype::Vec3Types>;
template class PBDObject<sofa::defaulttype::RigidTypes>;
#endif
