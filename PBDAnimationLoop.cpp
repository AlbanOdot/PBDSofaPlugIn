#include "PBDAnimationLoop.hpp"
#include <sofa/core/ObjectFactory.h>

//Mecha stuffs
#include <SofaBaseMechanics/MechanicalObject.h>
#include<sofa/core/behavior/MultiVec.h>

//Visitors
#include <sofa/simulation/AnimateVisitor.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/PropagateEventVisitor.h>
#include <sofa/simulation/BehaviorUpdatePositionVisitor.h>
#include <sofa/simulation/SolveVisitor.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/simulation/UpdateMappingEndEvent.h>
#include <sofa/simulation/UpdateBoundingBoxVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>

#include <omp.h>



using namespace sofa::component::container;
using namespace sofa::core::objectmodel;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;

int PBDAnimationLoopClass = sofa::core::RegisterObject("Simulation loop to use in scene without constraints nor contact.")
                            .add< PBDAnimationLoop >()
                            .addDescription(R"(
                                            This loop do the following steps:
                                            - build and solve all linear systems in the scene : collision and time integration to compute the new values of the dofs
                                            - update the context (dt++)
                                            - update the mappings
                                            - update the bounding box (volume covering all objects of the scene))");

PBDAnimationLoop::PBDAnimationLoop(sofa::simulation::Node* _gnode)
    : Inherit()
    , gnode(_gnode)
{
}

PBDAnimationLoop::~PBDAnimationLoop()
{
}


void PBDAnimationLoop::init()
{
    if (!gnode)
        gnode = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_context = gnode->getContext();
    m_mechanicalObjects = m_context->getObjects< MechanicalObject< sofa::defaulttype::Vec3Types > >(BaseContext::SearchDown);
}

void PBDAnimationLoop::bwdInit ()
{
    //Get all rest positions
    for(auto& obj : m_mechanicalObjects)
    {
        m_rest.emplace_back(obj->readRestPositions());
    }

    //On récupère les topologies
    auto topologies = m_context->getObjects<sofa::core::topology::BaseMeshTopology>(BaseContext::SearchDown);

    for(uint t = 0; t < topologies.size (); ++t)
    {
        // object[currentVertex][neighbor].first == neighbors's index
        // object[currentVertex][neighbor].second == norm(rest_pos[currentVertex],rest_pos[neighbor])
        std::vector<std::vector<std::pair<uint,SReal>>> object;
        //We assume 1 topology per object and vice versa
        for(uint i = 0; i < m_rest[t].size(); ++i)
        {
            //Get the neighbors of point I
            const auto& neighbors = topologies[t]->getVerticesAroundVertex (i);
            std::vector<std::pair<uint,SReal>> neighborhood;
            for(uint j = 0; j < neighbors.size(); ++j)
            {
                if( neighbors[j] > i )//Unidirectionnal neighborhood
                {
                    SReal d = (m_rest[t][i] - m_rest[t][neighbors[j]]).norm();
                    neighborhood.emplace_back(std::pair<uint,SReal>(neighbors[j],d));
                }
            }
            object.emplace_back(neighborhood);
        }
        m_topology.emplace_back(object);
    }
}

void PBDAnimationLoop::setNode( sofa::simulation::Node* n )
{
    gnode=n;
}

void PBDAnimationLoop::step(const sofa::core::ExecParams* params,
                            SReal dt)
{

    if (dt == 0)
    {
        dt = gnode->getDt();
    }
    sofa::core::MechanicalParams mparams(*params);
    static const Vec3 zero(0,0,0);
    static const float inv_dt = 1.0/dt;

    for(uint mObj = 0; mObj < m_mechanicalObjects.size (); ++mObj)
    {
        //Object parameters
        WriteCoord x = m_mechanicalObjects[mObj]->writePositions ();
        WriteDeriv v = m_mechanicalObjects[mObj]->writeVelocities ();
        uint pointCount = x.ref ().size();

        //External forces
        Derivatives dFext(VecDeriv((int)pointCount,zero));

        //We will compute constrainst on p
        Coordinates freeCoord(x.ref());
        WriteCoord p = freeCoord;

        //Apply external forces on p
        m_integrator.integrateExternalForces(gnode,&mparams,dFext,p,x,v,dt);

        //Solve all of the constraints
        solveConstraints(mObj,p);

        //Integrate using PBD method
        m_integrator.updatePosAndVel (p,x,v,inv_dt);

    }

    gnode->execute<UpdateMappingVisitor>(params);
    {
        UpdateMappingEndEvent ev ( dt );
        PropagateEventVisitor act ( params , &ev );
        gnode->execute ( act );
    }
}

void PBDAnimationLoop::solveConstraints(const uint mID, WriteCoord& p)
{
    //From here we solve all of the constraints -> solve on p
    uint max_iter = 20;
    for(uint nbIter = 0; nbIter < max_iter; ++nbIter)
    {
        solveStretch(mID,p);
        solveFixedPoint(mID,p);
    }
}

void PBDAnimationLoop::solveStretch(const uint mID, WriteCoord& p)
{

    //Stretching constraint
    //Solving Gauss-Seidel's style
    //Use an accumulator to put it as Jacobi
    uint pointCount =  p.ref().size();
    const auto& topology = m_topology[mID];
    static const SReal k = 0.3;//std::pow(1.f-0.999f,1/20);
    for( uint i = 0; i < pointCount; ++i)
    {
        const auto& voisins = topology[i];
        for( const auto& voisin : voisins)
        {
            const Vec3& p_ij = p[i] - p[voisin.first];
            SReal l = p_ij.norm();
            SReal d = voisin.second;
            const auto& displacement = 0.5*(l-d)/l * p_ij;
            p[i] -= k * displacement;
            p[voisin.first] += k * displacement;
        }
    }

}

void PBDAnimationLoop::solveFixedPoint(const uint mID, WriteCoord& p)
{
    //Constrainst related stuff
    std::vector<uint> fixedIdx = {3,39,64};

    //Fixed point constraint
    for(const auto& idx : fixedIdx )
    {
        p[idx] = m_rest[mID][idx];
    }
}
