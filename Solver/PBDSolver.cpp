
#include "./PBDSolver.hpp"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/ObjectFactory.h>
#include <math.h>
#include <iostream>
#include <utility>
#include <sofa/helper/system/thread/CTime.h>
#include "../Constraint/FEM/PBDElasticRod.hpp"
#include <sofa/core/collision/NarrowPhaseDetection.h>
#include <sofa/core/collision/DetectionOutput.h>
#include <SofaBaseMechanics/SubsetMapping.h>

using namespace sofa;
using sofa::core::VecId;
using namespace sofa::defaulttype;
using namespace core::behavior;
using namespace  sofa::component::container;
using namespace sofa::defaulttype;
using namespace sofa::core::objectmodel;
using namespace sofa::component::topology;
/**
 *
 *  Given current and last position $u_n$ and $u_{n-1}$, the next position $u_{n+1}$ is computed as follow:
 *
 *  v_{n+1} = damping * (u_{n} - u_{n-1}) * inv_dt
 *  u_{n+1} = u_{n}
 */

void PBDSolver::integrate(SReal dt)
{
    m_damping_times_inv_dt = m_damping / dt;
    // vel = (p-x) / dt
    // pos = p
    computeVec3Integration(m_v3m,m_damping_times_inv_dt);
    computeRigidIntegration(m_node,m_damping_times_inv_dt);

}

void PBDSolver::setupSolver(sofa::simulation::Node* node, int nbIter,SReal r)
{
    if(!node)
        return;
    m_node = node;
    m_constraints = node->getContext()->getObjects<PBDBaseConstraint>(BaseContext::SearchDown);
    m_narrowPhase = m_node->getContext()->get<core::collision::NarrowPhaseDetection>();
    std::vector<PBDConstraint<Vec3Types> *> constraints = node->getContext()->getObjects<PBDConstraint<Vec3Types>>(BaseContext::SearchDown);
    std::vector<PBDConstraint<RigidTypes> *> constraintsr = node->getContext()->getObjects<PBDConstraint<RigidTypes>>(BaseContext::SearchDown);
    m_nbIter = nbIter;
    m_damping =  fabs(r) >= 1.0 ? 0.0 : 1.0-r;
    if(nbIter != 1){
        for(auto& constraint : m_constraints)
        {
            constraint->setIterCount(1);
        }
    }


    //Linking des constraintes Vec3
    for( auto& constraint : constraints)
    {
        std::pair < MechanicalObject<Vec3Types> * , PBDObject<Vec3Types> *> pair = {constraint->mechanical (),new PBDObject<Vec3Types>(constraint->mechanical())};
        m_v3m.emplace(pair);
        constraint->linkPBDObject(m_v3m[constraint->mechanical ()]);
    }

    //Linking des contraintes Rigid
    for( auto& constraint : constraintsr)
    {
        std::pair < MechanicalObject<RigidTypes> * , PBDObject<RigidTypes> *> pair = {constraint->mechanical (),new PBDObject<RigidTypes>(constraint->mechanical())};
        m_r3m.emplace(pair);
        constraint->linkPBDObject(m_r3m[constraint->mechanical ()]);
    }

}

void PBDSolver::updateFreePositionsAndVelocities()
{
    for(auto& obj : m_v3m)
    {
        obj.second->resetFreePosition ();
        obj.second->resetFreeVelocity ();
    }
    for(auto& obj : m_r3m)
    {
        obj.second->resetFreePosition ();
        obj.second->resetFreeVelocity ();
    }
}

void PBDSolver::generateCollisions()
{
    //Collision constraints generation
    if(m_narrowPhase != nullptr)
    {
        auto detectionOutputsMap = m_narrowPhase->getDetectionOutputs();
        if ( detectionOutputsMap.size() != 0 )
        {
            for (sofa::core::collision::NarrowPhaseDetection::DetectionOutputMap::iterator it = detectionOutputsMap.begin(); it!=detectionOutputsMap.end(); ++it )
            {
                //If the PBD constraint has been generated we remove so sofa doesn't compute it a second time
                generateCollisionConstraint(it->first,it->second);
            }
        }
    }
}

void PBDSolver::solvePBDConstraints (const core::ExecParams* params)
{
    bool modification = true;
    //From here we solve all of the constraints -> solve on p
    for(int iter = 0 ; modification && iter < m_nbIter ; ++iter)
    {
        for(auto& constraint : m_constraints)
        {
            modification |= constraint->solve(m_node);
        }
        if(!modification)
            break;
        //Uncomment this when the collision will be handled
        //        for(auto& collision : m_collision_constraints)
        //        {
        //            collision->solve(m_node);
        //        }
    }


    m_collision_constraints.clear ();
}

void PBDSolver::computeVec3Integration(std::unordered_map<MechanicalObject<Vec3Types> * ,PBDObject<Vec3Types> *>& mechanicalObjects, SReal damping_times_inv_dt)
{
    //We only loop over the object which are subject to at least one constraint
    for(auto& object : mechanicalObjects)
    {
        WriteCoord x = object.first->writePositions ();
        WriteCoord p = object.second->getFreePosition ();
        WriteDeriv v = object.first->writeVelocities ();
        uint size = x.size();
        for( uint i = 0 ; i < size; ++i)
        {
            v[i] += damping_times_inv_dt * (p[i] - x[i]);
            x[i] = p[i];
        }
    }
}

void PBDSolver::computeRigidIntegration(sofa::simulation::Node* node,SReal damping_times_inv_dt)
{
    std::vector<OrientedConstraint*>  orientedConstraints = node->getContext()->getObjects<OrientedConstraint>(sofa::core::objectmodel::BaseContext::SearchDown);
    for(OrientedConstraint * c : orientedConstraints)
    {
        auto& omega = c->orientation().angularSpeed ();
        auto& u = c->orientation().freeOrientation ();
        WriteCoordR x = c->mechanical ()->writePositions ();
        WriteCoordR p = c->getPBDObject ()->getFreePosition ();
        WriteDerivR v = c->mechanical ()->writeVelocities ();
        const auto& edges = c->topology ()->getEdges ();
        uint size = edges.size();
        for(uint e = 0; e < size; ++e)
        {
            omega[e] = (2.0*damping_times_inv_dt)*(c->orientation().orientation (e).conjugate()*u[e]).vec();
            c->orientation().orientation (e) = u[e];

            //La vitesse angulaire est le slerp des edges l'entourant
            //Quaternionr ov(0,omega[e][0],omega[e][1],omega[e][2]); //ov.slerp(0.5,Quaternionr(0,omega[e+1][0],omega[e+1][1],omega[e+1][2]));
            v[edges[e][1]].getVCenter () += (p[edges[e][1]].getCenter () - x[edges[e][1]].getCenter ())* damping_times_inv_dt;
            v[edges[e][1]].getAngular().set(omega[e][0],omega[e][1],omega[e][2]);

            //Le quaternion est l'interpolation linéaire associe des quat des segments entourant le sommet
            //Quaternionr qx(u[e]); //qx.slerp(0.5,u[e+1]);
            x[edges[e][1]].getCenter () = p[edges[e][1]].getCenter ();
            x[edges[e][1]].getOrientation ().set(u[e].x (),u[e].y (),u[e].z (),u[e].w ());
        }
        //Cas particulier pour le premier sommet
        v[edges[0][0]].getVCenter () += (p[edges[0][0]].getCenter () - x[edges[0][0]].getCenter ()) * damping_times_inv_dt;
        v[edges[0][0]].getAngular().set(omega[0][0],omega[0][1],omega[0][2]);
        x[edges[0][0]].getCenter () = p[0].getCenter ();
        x[edges[0][0]].getOrientation ().set(u[0].x (),u[0].y (),u[0].z (),u[0].w ());
    }

}

void PBDSolver::setupOrientations(SReal dt)
{
    std::vector<OrientedConstraint*>  orientedConstraints = m_node->getContext()->getObjects<OrientedConstraint>(sofa::core::objectmodel::BaseContext::SearchDown);
    for(OrientedConstraint * c : orientedConstraints)
    {
        auto& omega = c->orientation().angularSpeed ();
        auto& I = c->orientation().inertia ();
        auto& tau = c->orientation().torque ();
        auto& u = c->orientation().freeOrientation ();
        uint size = omega.size ();
        for(uint j = 0; j < size; ++j)
        {
            //Line 9 of the algorithm
            omega[j] += dt*I[j].asDiagonal ().inverse ()*(tau[j] - omega[j].cross(I[j].asDiagonal ()*omega[j])).eval ();
            //Line 11 of the algorithm
            u[j].coeffs() += (0.5*dt)*(c->orientation().orientation (j)*Quaternionr(0,omega[j].x (),omega[j].y (),omega[j].z ())).coeffs();//beam[j].m_q*
            u[j].normalize ();
        }
    }
}


void PBDSolver::findObjectInMaps(PBDObject<Rigid3Types>** pbdObjr,PBDObject<Vec3Types>** pbdObjv, BaseContext * cont)
{
    MechanicalObject<Vec3Types>* objv = cont->get<MechanicalObject<Vec3Types>>();

    //Test if he is in the Vec3 map or the Rigid one
    if( objv != nullptr )
    {
        const auto& correspondance = m_v3m.find (objv);
        *pbdObjv = correspondance != m_v3m.end() ? correspondance->second : nullptr;
    }else{
        MechanicalObject<RigidTypes>* objr = cont->get<MechanicalObject<Rigid3Types>>();
        const auto& correspondance = m_r3m.find(objr);
        *pbdObjr = correspondance != m_r3m.end () ? correspondance->second : nullptr;
    }
}

bool PBDSolver::generateCollisionConstraint(const std::pair<sofa::core::CollisionModel *, sofa::core::CollisionModel *> &collModels, sofa::core::collision::DetectionOutputVector *outputs)
{
    //Get the PBDobject of the mechanical object of the first collision model
    PBDObject<Rigid3Types>* pbdObj1r = nullptr;
    PBDObject<Vec3Types>* pbdObj1v = nullptr;
    findObjectInMaps (&pbdObj1r,&pbdObj1v,collModels.first->getContext ());

    //True : Object 1 is not subject to any PBD constraint so we cannot handle it the way PBD does
    if( (pbdObj1v == nullptr)  && (pbdObj1r == nullptr))
        return false;

    //Get the PBDobject of the mechanical object of the second collision model
    PBDObject<Rigid3Types>* pbdObj2r = nullptr;
    PBDObject<Vec3Types>* pbdObj2v = nullptr;
    findObjectInMaps (&pbdObj2r,&pbdObj2v,collModels.second->getContext ());

    //True : Object 2 is not subject to any PBD constraint so we cannot handle it the way PBD does
    if( (pbdObj2v == nullptr)  && (pbdObj2r == nullptr))
        return false;

    if(pbdObj1v && pbdObj2v)
    {
        PBDCollisionConstraint<Vec3dTypes,Vec3Types>* collvv = new PBDCollisionConstraint<Vec3dTypes,Vec3Types>(pbdObj1v,pbdObj2v,collModels.first,collModels.second,outputs,m_nbIter);
        m_collision_constraints.emplace_back(collvv);
        return true;
    }
    if(pbdObj1v && pbdObj2r)
    {
        PBDCollisionConstraint<Vec3dTypes,Rigid3Types>* collvr = new PBDCollisionConstraint<Vec3dTypes,Rigid3Types>(pbdObj1v,pbdObj2r,collModels.first,collModels.second,outputs,m_nbIter);
        m_collision_constraints.emplace_back(collvr);
        return true;
    }
    if(pbdObj1r && pbdObj2v)
    {
        PBDCollisionConstraint<Rigid3Types,Vec3Types>* collrv = new PBDCollisionConstraint<Rigid3Types,Vec3Types>(pbdObj1r,pbdObj2v,collModels.first,collModels.second,outputs,m_nbIter);
        m_collision_constraints.emplace_back(collrv);
        return true;
    }
    if(pbdObj1r && pbdObj2r)
    {
        PBDCollisionConstraint<Rigid3Types,Rigid3Types>* collrr = new PBDCollisionConstraint<Rigid3Types,Rigid3Types>(pbdObj1r,pbdObj2r,collModels.first,collModels.second,outputs,m_nbIter);
        m_collision_constraints.emplace_back(collrr);
        return true;
    }
    return false;
}

