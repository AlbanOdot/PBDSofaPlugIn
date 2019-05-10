#include "BeamElement.hpp"

BeamElement::BeamElement(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
                         const uint edge[2],
                         const uint e,
                         const SReal l)
{
    const auto& dir = mobj->readRestPositions ()[edge[1]] - mobj->readRestPositions ()[edge[0]];
    m_q.setFromTwoVectors(Eigen::Vector3d(0,0,1),Eigen::Vector3d(dir[0],dir[1],dir[2]));
    m_q.normalize ();
    m_segmentIndicies = {edge[0],edge[1]};
    m_segmentIndex = e;
    m_lambda.setZero();
    m_averageLength = l;
    m_invMass.setZero ();
    m_wq = 1.0;
}

void BeamElement::initRestDarboux(const Quaternion& q)
{
    m_restDarboux = m_q.conjugate () * q;
}

void BeamElement::initMassMatrix(SReal m, vec3 stiffness)
{

    m_invMass << 1.0/m,  0.0,  0.0,  0.0,  0.0,  0.0,
                   0.0,1.0/m,  0.0,  0.0,  0.0,  0.0,
                   0.0,  0.0,1.0/m,  0.0,  0.0,  0.0,
                   0.0,  0.0,  0.0, 1.0/(stiffness[0] * m_averageLength),  0.0, 0.0,
                   0.0,  0.0,  0.0,  0.0,1.0/(stiffness[1] * m_averageLength),  0.0,
                   0.0,  0.0,  0.0,  0.0,  0.0,1.0/(stiffness[2] * m_averageLength);


}
