#include "./PBDOrientation.hpp"

PBDOrientation::PBDOrientation(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDOrientation::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    const Eigen::Matrix3d id;
    for(uint i = 0; i < rest.size () - 1; ++i)
    {
        const auto& dir = rest[i+1] - rest[i];
        Quaternionr q;
        q.setFromTwoVectors(Vector3r(0,0,1),Vector3r(dir[0],dir[1],dir[2]));
        q.normalize ();
        m_orientation.emplace_back(q);
        m_freeOrientation.emplace_back(q);
        m_angularSpeed.emplace_back(Eigen::Vector3d(0,0,0));
        m_torque.emplace_back(Eigen::Vector3d(0,0,0));
        m_inertia.emplace_back(id.Identity());
    }
    //Last particle is set as the previous one.
    m_orientation.emplace_back(m_orientation[m_orientation.size () - 1]);
    m_freeOrientation.emplace_back(m_freeOrientation[m_freeOrientation.size () - 1]);
    m_angularSpeed.emplace_back(Eigen::Vector3d(0,0,0));
    m_torque.emplace_back(Eigen::Vector3d(0,0,0));
    m_inertia.emplace_back(id.Identity());

    for(uint i = 0; i < rest.size() - 1; ++i)
    {
        m_restDarboux.emplace_back(m_orientation[i+1].conjugate () * m_orientation[i]);
    }
    m_restDarboux.emplace_back(m_restDarboux[m_restDarboux.size() - 1] );
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
            speed = s.asDiagonal ();
        }
    }else if (as.size () == 1)
    {
        for(auto& speed : m_inertia)
        {
            speed = as[0].asDiagonal ();
        }
    }else if( as.size () == m_inertia.size ()) {
        for(uint i = 0; i < as.size(); ++i)
        {
            m_inertia[i] = as[i].asDiagonal ();
        }
    }
}
