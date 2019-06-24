#include "./PBDOrientation.hpp"
#include <sofa/defaulttype/RigidTypes.h>

using namespace sofa::defaulttype;
PBDOrientation::PBDOrientation(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDOrientation::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    auto rigid = m_mechanicalObject->writePositions ();
    for(uint i = 0; i < rest.size () - 1; ++i)
    {
        m_freeOrientation.emplace_back(rigid[i].getOrientation ());
        m_angularSpeed.emplace_back(Vec3(0,0,0));
        m_torque.emplace_back(Vec3(0,0,0));
        m_inertia.emplace_back(Vec3(1,1,1));
    }
    //Last particle is set as the previous one.
    m_freeOrientation.emplace_back(m_freeOrientation[m_freeOrientation.size () - 1]);
    m_angularSpeed.emplace_back(Vec3(0,0,0));
    m_torque.emplace_back(Vec3(0,0,0));
    m_inertia.emplace_back(Vec3(1,1,1));

    for(uint i = 0; i < rest.size(); ++i)
    {
        uint a = i+1 == rest.size() ? i : i+1;
        m_restDarboux.emplace_back(rigid[i].getOrientation ().inverse () * rigid[a].getOrientation ());
    }
}

void PBDOrientation::update()
{
    m_freeOrientation.clear ();
    m_restDarboux.clear ();
    m_angularSpeed.clear ();
    m_torque.clear ();
    m_inertia.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

void PBDOrientation::setAngularVelocity(const std::vector<Vec3> &as)
{
    if( as.size () == 0)
    {
        const auto& s = Vec3(0,0,0);
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

void PBDOrientation::setTorque(const std::vector<Vec3> &as)
{
    if( as.size () == 0)
    {
        const auto& s = Vec3(0,0,0);
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

void PBDOrientation::setInertia(const std::vector<Vec3> & as)
{
    if( as.size () == 0)
    {
        const auto& s = Vec3(1,1,1);
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
