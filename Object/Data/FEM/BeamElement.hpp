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

    BeamElement(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > *mobj,
                const uint edge[2],
                const uint e,
                const SReal l);
    void initMassMatrix(SReal m, vec3 stiffness);
    void initRestDarboux(const Quaternion& q);

    inline SReal length() const { return m_averageLength;}
    inline uint extremity(uint i) const { return i == 0 ? m_segmentIndicies.first : m_segmentIndicies.second;}
    inline Quaternion& q() { return m_q;}
    inline Quaternion& restDarboux() {return m_restDarboux;}
    inline Vec6& l() {return m_lambda;}
    inline Mat6& invMass() {return m_invMass;}
    inline unsigned int edge() {return m_segmentIndex;}
    inline SReal wq() {return m_wq;}
    void setwq(SReal w) {m_wq = w;}

    Quaternion m_q;//<<Quaternion describing the orientation of the segment
    Quaternion m_restDarboux;
    std::pair<unsigned int,unsigned int> m_segmentIndicies;
    unsigned int m_segmentIndex;
    SReal m_averageLength;
    Mat6 m_invMass;
    Vec6 m_lambda;
    SReal m_wq;

private:
    //Const variables

};

#endif // BEAMELEMENT_HPP
