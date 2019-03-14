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
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/helper/Quater.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/MultiVec.h>
#include <math.h>
#include <iostream>
#include <sofa/helper/AdvancedTimer.h>

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
/*
void PBDExplicitIntegrator::integrate(Coordinates& x,
                                      const Coordinates& dx,
                                      const SReal dt)
{

//#pragma omp for
    for(int currCoord = 0; currCoord < nbCoord; currCoord += 3)
    {
        x[currCoord  ] += dt * dx[currCoord  ];
        x[currCoord+1] += dt * dx[currCoord+1];
        x[currCoord+2] += dt * dx[currCoord+2];
    }
}

void PBDExplicitIntegrator::integrateTmp(Coordinates& tmpposition,
                                         const Coordinates& position,
                                         const Derivatives& velocities,
                                         SReal dt)
{

    uint nbCoord = velocities.size ();
    tmpposition.resize (nbCoord);

//#pragma omp for
    for(int currCoord = 0; currCoord < nbCoord; currCoord += 3)
    {
        tmpposition[currCoord  ] = position[currCoord  ] + dt * velocities[currCoord  ];
        tmpposition[currCoord+1] = position[currCoord+1] + dt * velocities[currCoord+1];
        tmpposition[currCoord+2] = position[currCoord+2] + dt * velocities[currCoord+2];
    }
}


void PBDExplicitIntegrator::integrateUniformExternalForces(Derivatives& velocities,
                                                           const Derivatives& extForces,
                                                           const SReal dt)
{
    uint nbCoord = velocities.size ();
    for(int currFCoord = 0; currFCoord < extForces.size (); currFCoord += 3){
//#pragma omp for
        for(int currCoord = 0; currCoord < nbCoord; currCoord += 3)
        {
            velocities[currCoord  ] += dt * extForces[currFCoord];
            velocities[currCoord+1] += dt * extForces[currFCoord+1];
            velocities[currCoord+2] += dt * extForces[currFCoord+2];
        }
    }
}


void PBDExplicitIntegrator::solveDistanceConstraint(Coordinates& position,
                                                    const Coordinates& truth,
                                                    const std::vector<uint>& pointIdx)
{
    //It's highly recommended to do not give any indicies if you want to compute for all since the algorithm is twice as fast
    if(pointIdx.empty ()){
        solveDistanceConstraintAll(position, truth, pointIdx);
        return;
    }

    uint nbCoord = position.size ();
    for(auto idx: pointIdx)
    {
        uint j = 3*idx;
        Vector p1 = Vector(position[j],position[j+1],position[j+2]);
        Vector t1 = Vector(truth[j],truth[j+1],truth[j+2]);
        uint i = 0;
        for(        ; i < j; i += 3){
            Vector p2 = Vector(position[i],position[i+1],position[i+2]);
            Vector t2 = Vector(truth[i],truth[i+1],truth[i+2]);
            Vector p2p1 = p1-p2;
            SReal dist = p2p1.norm();
            //need to change 0.5 by weight1 / (weight2 + weight1) and weight1 / (weight2 + weight1) -> We suppose equal mass
            Vector disparity = (0.5 * (1.0 - (t1-t2).norm() / dist ) ) * p2p1; //parenthesis to avoid real to vec3 multiplications.
            p1 += disparity;
            p2 -= disparity;
            position[j] = p1[0]; position[j+1] = p1[1]; position[j+2] = p1[2];
            position[j] = p2[0]; position[j+1] = p2[1]; position[j+2] = p2[2];
        }
        for(i += 3; i < nbCoord; i += 3){
            Vector p2 = Vector(position[i],position[i+1],position[i+2]);
            Vector t2 = Vector(truth[i],truth[i+1],truth[i+2]);
            Vector p2p1 = p1-p2;
            SReal dist = p2p1.norm();
            //need to change 0.5 by weight1 / (weight2 + weight1) and weight1 / (weight2 + weight1) -> We suppose equal mass
            Vector disparity = (0.5 * (1.0 - (t1-t2).norm() / dist ) ) * p2p1; //parenthesis to avoid real to vec3 multiplications.
            p1 += disparity;
            p2 -= disparity;
            position[j] = p1[0]; position[j+1] = p1[1]; position[j+2] = p1[2];
            position[j] = p2[0]; position[j+1] = p2[1]; position[j+2] = p2[2];
        }
    }
}


void PBDExplicitIntegrator::solveDistanceConstraintAll(Coordinates& position,
                                                       const Coordinates& truth,
                                                       const std::vector<uint>& pointIdx)
{

    //What's the point man ?
    if(pointIdx.empty()){
        //solveFixedConstraintAll();
        return;
    }

    //No threading since usually there is few attached points
    for(uint idx: pointIdx)
    {
        std::cout << "Fixed point is :"<<idx<<std::endl;

        uint j = 3*idx;
        SReal t = position[3*j];
        /*position[3*j  ] = truth[3*j  ];
        position[3*j+1] = truth[3*j+1];
        position[3*j+2] = truth[3*j+2];

    }
}

void PBDExplicitIntegrator::solveFixedPointConstraint(Coordinates& position,
                                                      const Coordinates& truth,
                                                      const std::vector<uint>& pointIdx)
{
    uint nbCoord = position.size ();

    //We cannot //#pragma omp for here, random order might introduce instabillity
    for(uint j = 0; j < nbCoord; j += 3 )
    {
        Vector p1 = Vector(position[j],position[j+1],position[j+2]);
        Vector t1 = Vector(truth[j],truth[j+1],truth[j+2]);
        //#pragma omp for
        for(uint i = j+3; i < nbCoord; i += 3){
            Vector p2 = Vector(position[i],position[i+1],position[i+2]);
            Vector t2 = Vector(truth[i],truth[i+1],truth[i+2]);
            Vector p2p1 = p1-p2;
            SReal dist = p2p1.norm();
            //need to change 0.5 by weight1 / (weight2 + weight1) and weight1 / (weight2 + weight1) -> We suppose equal mass
            Vector disparity = (0.5 * (1.0 - (t1-t2).norm() / dist ) ) * p2p1; //parenthesis to avoid real to vec3 multiplications.
            p2 -= disparity;
            position[i] = p2[0]; position[i+1] = p2[1]; position[i+2] = p2[2];
            //#pragma omp critical
            {
                position[j] = p1[0] + disparity[0]; position[j+1] = p1[1] + disparity[1]; position[j+2] = p1[2] + disparity[2];
            }
        }
    }
}

void PBDExplicitIntegrator::PBDUpdate(const Coordinates& newPosition,
                                      Derivatives& velocity,
                                      Coordinates& position,
                                      const SReal dt)
{
    uint nbCoord = position.size ();
    SReal one_over_dt = 1.0 / dt;
    for(int i = 0; i < nbCoord; i += 3)
    {
        velocity[i  ] = (newPosition[i  ] - position[i  ]) * one_over_dt;
        velocity[i+1] = (newPosition[i+1] - position[i+1]) * one_over_dt;
        velocity[i+2] = (newPosition[i+2] - position[i+2]) * one_over_dt;

        position[i  ] = newPosition[i  ];
        position[i+1] = newPosition[i+1];
        position[i+2] = newPosition[i+2];
    }
}
*/
