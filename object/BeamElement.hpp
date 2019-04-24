#ifndef BEAMELEMENT_HPP
#define BEAMELEMENT_HPP

#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Quat.h>
#include <Eigen/MatrixFunctions>
#include <SofaBaseMechanics/MechanicalObject.h>

class BeamElement
{
public:
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;
    typedef sofa::defaulttype::Vec3 vec3;
    typedef Eigen::Quaternion<SReal> Quaternion;
    typedef Eigen::Matrix<double,6,6> Mat6;
    typedef Eigen::Matrix<double,6,1> Vec6;

    BeamElement(const uint edge[2], const sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj);
    BeamElement(const sofa::defaulttype::Vec3& x1,  const sofa::defaulttype::Vec3& x2);
    BeamElement(const sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj, const uint edge[3]);
    BeamElement(BeamElement& b, bool ghost);

    inline SReal length() const { return m_length;}
    inline SReal invLength() const{return m_inv_length;}
    inline uint current() const { return m_orig_idx;}
    inline uint next() const { return m_end_idx;}
    inline Quaternion& q() { return m_q;}
    inline Quaternion& restDarboux() {return m_restDarboux;}

    Quaternion m_q;//<<Quaternion describing the orientation of the segment

private:
    //Const variables
    Quaternion m_restDarboux;
    uint m_orig_idx;
    uint m_end_idx;
    SReal m_length;
    SReal m_inv_length;
};

#endif // BEAMELEMENT_HPP
