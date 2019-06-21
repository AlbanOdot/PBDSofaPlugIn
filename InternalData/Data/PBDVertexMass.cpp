#include "PBDVertexMass.hpp"
#include <SofaBaseMechanics/UniformMass.h>

template < class T>
PBDVertexMass<T>::PBDVertexMass(Mech * m, Topo * t) : PBDBaseConstraintData<T> (m,t)
{
    if( m && t )
        init ();
}

template < class T>
void PBDVertexMass<T>::init()
{
    const auto& rest = PBDBaseConstraintData<T>::m_mechanicalObject->readRestPositions ();
    //Compute the vertice oriented Mass
    for(uint i = 0; i < rest.size(); ++i)
    {
        SReal m  = static_cast<sofa::component::mass::UniformMass<sofa::defaulttype::Vec3Types,SReal> *>(PBDBaseConstraintData<T>::m_mechanicalObject->getContext ()->getMass())->getVertexMass ();
        SReal w = 1.0/m;
        m_mass.emplace_back(m);
        m_weight.emplace_back(w);
    }
}

template < class T>
void PBDVertexMass<T>::update()
{
    m_mass.clear();
    m_weight.clear();
    if( PBDBaseConstraintData<T>::m_mechanicalObject && PBDBaseConstraintData<T>::m_sofa_topology )
        init ();
}
