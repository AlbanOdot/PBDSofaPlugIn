#include "PBDVertexMass.hpp"
#include <SofaBaseMechanics/UniformMass.h>

PBDVertexMass::PBDVertexMass(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDVertexMass::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    //Compute the vertice oriented Mass
    for(uint i = 0; i < rest.size(); ++i)
    {
        SReal m  = static_cast<sofa::component::mass::UniformMass<sofa::defaulttype::Vec3Types,SReal> *>(m_mechanicalObject->getContext ()->getMass())->getVertexMass ();
        SReal w = 1.0/m;
        m_mass.emplace_back(m);
        m_weight.emplace_back(w);
    }
}

void PBDVertexMass::update()
{
    m_mass.clear();
    m_weight.clear();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}
