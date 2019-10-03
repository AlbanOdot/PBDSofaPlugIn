#include "PBDVertexMass.hpp"
#include <SofaBaseMechanics/UniformMass.h>
#include <SofaBaseMechanics/DiagonalMass.h>

template < class T>
PBDVertexMass<T>::PBDVertexMass(Mech * m, Topo * t) : PBDBaseConstraintData<T> (m,t)
{
    if( m )
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
        if( m <= 1e-50 )
            m = static_cast<sofa::component::mass::DiagonalMass<sofa::defaulttype::Vec3Types,SReal> *>(
                        PBDBaseConstraintData<T>::m_mechanicalObject->getContext ()->getMass())->getElementMass(i);
        m_mass.emplace_back(m);
        m_weight.emplace_back(1.0/m);
    }
}

template < class T>
void PBDVertexMass<T>::update()
{
    m_mass.clear();
    m_weight.clear();
    if( PBDBaseConstraintData<T>::m_mechanicalObject  )
        init ();
}
