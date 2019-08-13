
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

using namespace sofa;
using sofa::core::VecId;
using namespace sofa::defaulttype;
using namespace core::behavior;
using namespace  sofa::component::container;
using namespace sofa::defaulttype;
using namespace sofa::core::objectmodel;

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
    m_constraint = node->getContext()->getObjects<PBDBaseConstraint>(BaseContext::SearchDown);
    std::vector<PBDConstraint<Vec3Types> *> constraints = node->getContext()->getObjects<PBDConstraint<Vec3Types>>(BaseContext::SearchDown);
    std::vector<PBDConstraint<RigidTypes> *> constraintsr = node->getContext()->getObjects<PBDConstraint<RigidTypes>>(BaseContext::SearchDown);
    m_nbIter = nbIter;
    m_damping =  fabs(r) >= 1.0 ? 0.0 : 1.0-r;
    if(nbIter != 1){
        for(auto& constraint : m_constraint)
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

void PBDSolver::solvePBDConstraints (const core::ExecParams* params)
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
    //From here we solve all of the constraints -> solve on p
    for(int iter = 0 ; iter < m_nbIter; ++iter)
    {
        for(auto& constraint : m_constraint)
        {
            constraint->solve(m_node);
        }
    }
}

void PBDSolver::solveSofaConstraints(const core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId freePosId, sofa::core::MultiVecDerivId freeVelId)
{
    sofa::simulation::common::VectorOperations vop( params, m_node->getContext() );
    MultiVecCoord p(&vop,freePosId /*core::VecCoordId::position()*/ );
    MultiVecDeriv vf(&vop,freeVelId /*core::VecDerivId::velocity()*/ );
    sofa::simulation::common::MechanicalOperations mop( params, m_node->getContext() );
    mop->setImplicit(false); // this solver is explicit only
    mop.solveConstraint(vf,core::ConstraintParams::VEL);
    mop.solveConstraint(p,core::ConstraintParams::POS);
}

void PBDSolver::computeVec3Integration(std::unordered_map<MechanicalObject<Vec3Types> * ,PBDObject<Vec3Types> *>& mechanicalObjects, SReal damping_times_inv_dt)
{
    //We only loop over the object which are subject to at least one constraint
    for(auto& object : mechanicalObjects)
    {
        WriteCoord x = object.first->writePositions ();
        WriteCoord p = object.second->getFreePosition ();
        WriteDeriv v = object.first->writeVelocities ();
        for( uint i = 0 ; i < x.size (); ++i)
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
        uint s1 = size - 1;
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

//        //Cas particulier pour le dernier sommet
//        v[edges[s1][1]].getVCenter () += (p[edges[s1][1]].getCenter () - x[edges[s1][1]].getCenter ()) * damping_times_inv_dt;
//        v[edges[s1][1]].getAngular() = sofa::defaulttype::Vec3(omega[s1][0],omega[s1][1],omega[s1][2]);
//        x[edges[s1][1]].getCenter () = p[edges[s1][1]].getCenter ();
//        x[edges[s1][1]].getOrientation ().set(u[s1].x (),u[s1].y (),u[s1].z (),u[s1].w ());
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
        for(uint j = 0; j < omega.size (); ++j)
        {
            //Line 9 of the algorithm
            omega[j] += dt*I[j].asDiagonal ().inverse ()*(tau[j] - omega[j].cross(I[j].asDiagonal ()*omega[j])).eval ();
            //Line 11 of the algorithm
            u[j].coeffs() += (0.5*dt)*(c->orientation().orientation (j)*Quaternionr(0,omega[j].x (),omega[j].y (),omega[j].z ())).coeffs();//beam[j].m_q*
            u[j].normalize ();
        }
    }
}
