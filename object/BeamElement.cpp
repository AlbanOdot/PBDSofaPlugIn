#include "BeamElement.hpp"

BeamElement::BeamElement(const uint edge[2],
                         const sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj)
{
    //Origin, Orientation, length init
    m_orig_idx = edge[0];
    m_end_idx = edge[1];
    const auto& rest = mobj->readRestPositions ();
    const vec3& a = rest[edge[0]];
    const vec3& b = rest[edge[1]];
    const vec3& q = b-a;
    m_length = q.norm();
    m_inv_length = 1.0/m_length;
    m_q.setFromTwoVectors(Eigen::Vector3d(0,0,1),Eigen::Vector3d(q[0],q[1],q[2]));
    m_q.normalize ();
    m_restDarboux = m_q.conjugate () * m_q;
}

BeamElement::BeamElement(const sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
                         const uint edge[3])
{
    //Origin, Orientation, length init
    m_orig_idx = edge[0];
    m_end_idx = edge[1];
    const auto& rest = mobj->readRestPositions ();
    const vec3& a = rest[edge[0]];
    const vec3& b = rest[edge[1]];
    const vec3& c = rest[edge[2]];
    const vec3& d = b-a;
    const vec3& e = c-b;
    m_length = 0.5 * (d.norm() + e.norm());//<< Mean length of the segment and the enxt one
    m_inv_length = 1.0/m_length;
    m_q.setFromTwoVectors(Eigen::Vector3d(0,0,1),Eigen::Vector3d(d[0],d[1],d[2]));
    m_q.normalize ();
    m_restDarboux = m_q.conjugate () * Quaternion::FromTwoVectors(Eigen::Vector3d(0,0,1),Eigen::Vector3d(e[0],e[1],e[2]));
}

BeamElement::BeamElement(BeamElement& b,bool ghostEnd)
{
    if(ghostEnd)
    {
    //Origin, Orientation, length init
    m_orig_idx = b.next();
    m_end_idx = b.next();
    }else {
        //Origin, Orientation, length init
        m_orig_idx = b.current();
        m_end_idx = b.current();
    }
    m_length = b.length ();
    m_inv_length = b.invLength ();
    m_q = b.q ();
    m_restDarboux = b.m_restDarboux;
}

