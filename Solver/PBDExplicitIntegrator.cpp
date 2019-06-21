/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "PBDExplicitIntegrator.hpp"
#include <sofa/core/ObjectFactory.h>
//#define SOFA_NO_VMULTIOP

using namespace sofa::defaulttype;
int PBDExplicitIntegratorClass = sofa::core::RegisterObject("A simple explicit time integrator")
                                 .add< PBDExplicitIntegrator<Vec3Types> >(true)
                                 .add< PBDExplicitIntegrator<Rigid3Types>>()
                                                                          .addAlias("PBDExplicit")
                                                                          ;
template < class T >
PBDExplicitIntegrator<T>::PBDExplicitIntegrator()
{
}

template < class T >
PBDExplicitIntegrator<T>::~PBDExplicitIntegrator ()
{
}

inline sofa::defaulttype::Vec3 cwProduct(const sofa::defaulttype::Vec3& diagMat, const sofa::defaulttype::Vec3& vec)
{
    return sofa::defaulttype::Vec3(diagMat[0] * vec[0], diagMat[1] * vec[1], diagMat[2] * vec[2]);
}

inline sofa::defaulttype::Vec3 inverseDiag(sofa::defaulttype::Vec3& diagMat)
{
    return sofa::defaulttype::Vec3(1.0/diagMat[0], 1.0 / diagMat[1], 1.0/diagMat[2]);
}

template < class T >
void PBDExplicitIntegrator<T>::integrateExternalForces( const sofa::simulation::Node * gnode,
                                                        const sofa::core::MechanicalParams * mparams,
                                                        Derivatives& f,
                                                        WriteCoord& p,
                                                        const WriteCoord& x,
                                                        WriteDeriv& v ,
                                                        SReal dt)
{
    //Accumulates all external forces (Gravities,interactions etc)
    for(uint i = 0;  i < gnode->forceField.size(); ++i)
    {
        auto ff = dynamic_cast<sofa::core::behavior::ForceField< T > *>(gnode->forceField.get(i));
        ff->addForce(mparams,f,x.ref(),v.ref());
    }

    const WriteDeriv& Fext = f;
    uint pointCount = v.ref().size();
    for(uint i = 0; i < pointCount; ++i)
    {
        v[i] += Fext[i] * dt;
        p[i] = x[i] + v[i] * dt;
    }
}

template <>
void PBDExplicitIntegrator<sofa::defaulttype::RigidTypes>::updatePosAndVel (PBDObject<sofa::defaulttype::RigidTypes>& object,
                                                                            const WriteCoord &p,
                                                                            WriteCoord &x,
                                                                            WriteDeriv &v,
                                                                            const SReal &inv_dt)
{

    auto pointCount = v.ref().size();
    auto& orientation = object.orientation ();
    auto& omega = orientation.angularSpeed ();
    auto& u = orientation.freeOrientation();
    auto rigid = object.position ();
    Quaternion q;

    for(uint i = 0; i < pointCount; ++i)
    {
        v[i] = sofa::defaulttype::RigidTypes::coordDifference (p[i],x[i]) * inv_dt;//(p[i].coordDifference(x[i]) * inv_dt;
        x[i] = p[i];
        q = (rigid[i].getOrientation().inverse()*u[i])* 2.0*inv_dt;
        omega[i][0] = q[0]; omega[i][1] = q[1];omega[i][1] = q[1];
        rigid[i].getOrientation () = u[i];
    }
}

template <>
void PBDExplicitIntegrator<sofa::defaulttype::Vec3Types>::updatePosAndVel (PBDObject<sofa::defaulttype::Vec3Types>& object,
                                                                            const WriteCoord &p,
                                                                            WriteCoord &x,
                                                                            WriteDeriv &v,
                                                                            const SReal &inv_dt)
{

    auto pointCount = v.ref().size();
    for(uint i = 0; i < pointCount; ++i)
    {
        v[i] = (p[i]-x[i]) * inv_dt;
        x[i] = p[i];
    }
}

template < class T >
void PBDExplicitIntegrator<T>::setUpIntegrator(sofa::simulation::Node* node, int nbIter)
{
    if(!node)
        return;

    m_constraint = node->getContext()->getObjects<PBDBaseConstraint<T>>(sofa::core::objectmodel::BaseContext::SearchDown);

    m_nbIter = nbIter;
    if(nbIter != 1){
        for(auto& constraint : m_constraint)
        {
            constraint->setIterCount(1);
        }
    }

}

template < class T >
void PBDExplicitIntegrator<T>::solveConstraint (PBDObject<T>& object, WriteCoord& p)
{
    //From here we solve all of the constraints -> solve on p
    for(int iter = 0 ; iter < m_nbIter; ++iter)
    {
        for(auto& constraint : m_constraint)
        {
            constraint->solve(object,p);
        }
    }

}

template <>
void PBDExplicitIntegrator<sofa::defaulttype::RigidTypes>::integrateAngularVelocity(PBDObject<sofa::defaulttype::RigidTypes>& object,const SReal &dt)
{

    auto& orientation = object.orientation ();
    auto& omega = orientation.angularSpeed ();
    auto& I = orientation.inertia ();
    auto& tau = orientation.torque ();
    auto& u = orientation.freeOrientation();
    const auto& rigid = object.position ();
    for(uint j = 0; j < omega.size (); ++j)
    {
        omega[j] += dt*cwProduct( inverseDiag(I[j]), tau[j] - omega[j].cross(cwProduct(I[j],omega[j])));
        u[j] += (rigid[j].getOrientation () * Quaternion(0,omega[j][0],omega[j][1],omega[j][2])) * 0.5 * dt;
        u[j].normalize ();
    }
}

template < >
void PBDExplicitIntegrator<sofa::defaulttype::Vec3Types>::integrateAngularVelocity(PBDObject<sofa::defaulttype::Vec3Types>& object,const SReal &dt)
{
}
