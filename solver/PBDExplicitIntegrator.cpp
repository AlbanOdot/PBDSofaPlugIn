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

void PBDExplicitIntegrator::updatePosAndVel (const WriteCoord &p, WriteCoord &x, WriteDeriv &v, const SReal &inv_dt)
{
    uint pointCount = v.ref().size();
    for(uint i = 0; i < pointCount; ++i)
    {
        v[i] = (p[i] - x[i]) * inv_dt;
        x[i] = p[i];
    }
}

void PBDExplicitIntegrator::setUpIntegrator(sofa::simulation::Node* node)
{
    if(!node)
        return;

    m_constraint = node->getContext()->getObjects<PBDBaseConstraint>(sofa::core::objectmodel::BaseContext::SearchDown);

}


void PBDExplicitIntegrator::solveConstraint (PBDObject& object, WriteCoord& p, uint iter)
{
    //From here we solve all of the constraints -> solve on p
    m_max_iter = iter;
    for(uint nbIter = 0; nbIter < m_max_iter; ++nbIter)
    {
        for(auto& constraint : m_constraint)
        {
            constraint->solve(object,p);
        }
    }
}
