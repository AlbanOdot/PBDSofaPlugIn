#include "./PBDOrientation.hpp"

PBDOrientation::PBDOrientation(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDOrientation::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    const auto& edges = m_sofa_topology->getEdges ();
    for(const auto& e : edges)
    {
        const auto& q0 = rest[e[0]].getOrientation ();
        const auto& q1 = rest[e[1]].getOrientation ();
        Quaternion qe; qe.slerp (q0,q1,0.5);
        Quaternionr q(qe[3],qe[0],qe[1],qe[2]);
        //Orientation
        m_orientation.emplace_back(q);
        m_freeOrientation.emplace_back(q);
        //Velocity
        m_angularSpeed.emplace_back(Eigen::Vector3d(0,0,0));
        m_torque.emplace_back(Eigen::Vector3d(0,0,0));
        m_inertia.emplace_back(Eigen::Vector3d(1,1,1));
    }
    for(uint i = 0; i < m_orientation.size(); ++i)
    {
        uint a = i == edges.size() - 1 ? i : i+1;
        Quaternionr rd = m_orientation[i].conjugate () * m_orientation[a];
        Quaternionr omega_plus, omega_minus;
        omega_plus.coeffs() = rd.coeffs() + Quaternionr(1, 0, 0, 0).coeffs();
        omega_minus.coeffs() = rd.coeffs() - Quaternionr(1, 0, 0, 0).coeffs();
        if (omega_minus.squaredNorm() > omega_plus.squaredNorm())
            rd.coeffs() *= -1.0;
        m_restDarboux.emplace_back(rd);
    }
}

void PBDOrientation::update()
{
    m_orientation.clear ();
    m_freeOrientation.clear ();
    m_restDarboux.clear ();
    m_angularSpeed.clear ();
    m_torque.clear ();
    m_inertia.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

void PBDOrientation::setAngularVelocity(const std::vector<Vector3r> &as)
{
    if( as.size () == 0)
    {
        const auto& s = Vector3r(0,0,0);
        for(auto& speed : m_angularSpeed)
        {
            speed = s;
        }
    }else if (as.size () == 1)
    {
        for(auto& speed : m_angularSpeed)
        {
            speed = as[0];
        }
    }else if( as.size () == m_angularSpeed.size ()) {
        for(uint i = 0; i < as.size(); ++i)
        {
            m_angularSpeed[i] = as[i];
        }
    }
}

void PBDOrientation::setTorque(const std::vector<Vector3r> &as)
{
    if( as.size () == 0)
    {
        const auto& s = Vector3r(0,0,0);
        for(auto& speed : m_torque)
        {
            speed = s;
        }
    }else if (as.size () == 1)
    {
        for(auto& speed : m_torque)
        {
            speed = as[0];
        }
    }else if( as.size () == m_torque.size ()) {
        for(uint i = 0; i < as.size(); ++i)
        {
            m_torque[i] = as[i];
        }
    }
}

void PBDOrientation::setInertia(const std::vector<Vector3r> & as)
{
    if( as.size () == 0)
    {
        const auto& s = Vector3r(1,1,1);
        for(auto& speed : m_inertia)
        {
            speed = s;
        }
    }else if (as.size () == 1)
    {
        for(auto& speed : m_inertia)
        {
            speed = as[0];
        }
    }else if( as.size () == m_inertia.size ()) {
        for(uint i = 0; i < as.size(); ++i)
        {
            m_inertia[i] = as[i];
        }
    }
}
