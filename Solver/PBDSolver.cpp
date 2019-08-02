
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
#include "../Constraint/FEM/PDCosseratRod.hpp"

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
    SReal damping_times_inv_dt = m_damping / dt;
    // vel = (p-x) / dt
    // pos = p
    computeVec3Integration(m_v3m,damping_times_inv_dt);
    computeRigidIntegration(m_node,dt,damping_times_inv_dt);

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

void PBDSolver::computeRigidIntegration(sofa::simulation::Node* node, SReal dt, SReal damping_times_inv_dt)
{
    // Loop over r3m ??
    std::vector<OrientedConstraint*>  orientedConstraints = node->getContext()->getObjects<OrientedConstraint>(sofa::core::objectmodel::BaseContext::SearchDown);
    for(OrientedConstraint * c : orientedConstraints)
    {
        auto& omega = c->orientation().angularSpeed ();
        auto& I = c->orientation().inertia ();
        auto& tau = c->orientation().torque ();
        auto& u = c->orientation().freeOrientation ();
        WriteCoordR x = c->mechanical ()->writePositions ();
        WriteCoordR p = c->getPBDObject ()->getFreePosition ();
        WriteDerivR v = c->mechanical ()->writeVelocities ();
        for(uint j = 0; j < omega.size (); ++j)
        {
            omega[j] += dt*I[j].asDiagonal ().inverse ()*(tau[j] - omega[j].cross(I[j].asDiagonal ()*omega[j])).eval ();
            u[j].coeffs() += (0.5*dt)*(c->orientation().orientation (j)*Quaternionr(0,omega[j].x (),omega[j].y (),omega[j].z ())).coeffs();//beam[j].m_q*
            u[j].normalize ();
            omega[j] = (2.0*damping_times_inv_dt)*(c->orientation().orientation (j).conjugate()*u[j]).vec();
            c->orientation().orientation (j) = u[j];

            //v[j] += RigidTypes::coordDifference (p[j],x[j]) * damping_times_inv_dt;
            v[j].getLinear () += damping_times_inv_dt * (p[j].getCenter () - x[j].getCenter ());
            v[j].getAngular() = sofa::defaulttype::Vec3(omega[j][0],omega[j][1],omega[j][2]);
            x[j].getCenter () = p[j].getCenter ();
            //x[j].getOrientation ().set(u[j].x (),u[j].y (),u[j].z (),u[j].w ());
        }
    }

}
