#include "PBDObject.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <sofa/core/objectmodel/BaseContext.h>
#include <SofaBaseMechanics/UniformMass.h>
#include <sofa/core/ObjectFactory.h>

template < class T>
void PBDObject<T>::resetFreePosition ()
{
    const ReadCoord& positions(m_obj->read(sofa::core::VecCoordId::position ()));
    uint size = positions.ref().size();
    WriteCoord freePos(m_freePos);
    #pragma omp parallel for if(size >= 1024)
    for(uint i = 0; i < size; ++i)
    {
        freePos[i] = positions[i];
    }
}

template < class T>
void PBDObject<T>::resetFreeVelocity ()
{
    const ReadDeriv& velocity(m_obj->read(sofa::core::VecDerivId::velocity ()));
    WriteDeriv freeVel(m_freeVel);
    uint size = velocity.ref().size();
    #pragma omp parallel for if(size >= 1024)
    for(uint i = 0; i < size; ++i)
    {
        freeVel[i] = velocity[i];
    }
}
