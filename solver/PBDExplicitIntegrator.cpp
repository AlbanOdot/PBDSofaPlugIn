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


int PBDExplicitIntegratorClass = sofa::core::RegisterObject("A simple explicit time integrator")
                                 .add< PBDExplicitIntegrator >()
                                 .addAlias("PBDExplicit")
                                 ;

PBDExplicitIntegrator::PBDExplicitIntegrator()
{
}

PBDExplicitIntegrator::~PBDExplicitIntegrator ()
{
}

void PBDExplicitIntegrator::integrateExternalForces( const sofa::simulation::Node * gnode,
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
        auto ff = dynamic_cast<sofa::core::behavior::ForceField< sofa::defaulttype::Vec3Types > *>(gnode->forceField.get(i));
        ff->addForce(mparams,f,x.ref(),v.ref());
    }

    const WriteDeriv& Fext = f;
    uint pointCount = v.ref().size();
    for(uint i = 0; i < pointCount; ++i)
    {
        v[i] = v[i] + dt * Fext[i];
        p[i] = x[i] + dt * v[i];
    }
}

void PBDExplicitIntegrator::updatePosAndVel (PBDObject& object,
                                             const WriteCoord &p,
                                             WriteCoord &x,
                                             WriteDeriv &v,
                                             const SReal &inv_dt)
{

    auto pointCount = v.ref().size();
    for(uint i = 0; i < pointCount; ++i)
    {
        v[i] = (p[i] - x[i]) * inv_dt;
        x[i] = p[i];
    }

    auto orientationCount = object.freeOrientation ().size ();
    auto& omega = object.angularSpeed ();
    auto& u = object.freeOrientation ();
    auto& beam = object.beam ();
    Eigen::Quaterniond n;
    for(uint i = 0; i < orientationCount; ++i)
    {
        n.coeffs () = 2.0*inv_dt*(beam[i].q().conjugate()*u[i]).coeffs ();
        omega[i] = Eigen::Vector3d(n.x(),n.y (),n.z ());
        beam[i].q() = u[i];
    }
}

void PBDExplicitIntegrator::setUpIntegrator(sofa::simulation::Node* node, int nbIter)
{
    if(!node)
        return;

    m_constraint = node->getContext()->getObjects<PBDBaseConstraint>(sofa::core::objectmodel::BaseContext::SearchDown);

    m_nbIter = nbIter;
    if(nbIter != 1){
        for(auto& constraint : m_constraint)
        {
            constraint->setIterCount(1);
        }
    }

}


void PBDExplicitIntegrator::solveConstraint (PBDObject& object, WriteCoord& p)
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

void PBDExplicitIntegrator::integrateAngularVelocity(PBDObject& object,const SReal &dt)
{
    if(object.beam ().empty ())
        return;
    auto& omega = object.angularSpeed ();
    auto& I = object.inertia ();
    auto& tau = object.torque ();
    auto& u = object.freeOrientation ();
    auto& beam = object.beam ();
    for(uint j = 0; j < omega.size (); ++j)
    {
        omega[j] += dt*I[j].inverse ()*(tau[j] - omega[j].cross(I[j]*omega[j])).eval ();
        u[j].coeffs() += (0.5-1e-3)*dt*(beam[j].m_q*Eigen::Quaterniond(0,omega[j].x (),omega[j].y (),omega[j].z ())).coeffs();//beam[j].m_q*
        u[j].normalize ();
    }
}
