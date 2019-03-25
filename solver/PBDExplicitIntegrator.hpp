#ifndef PBDEXPLICITEINTEGRATOR_HPP
#define PBDEXPLICITEINTEGRATOR_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/simulation/Node.h>
#include <SofaBaseMechanics/MechanicalObject.h>

class PBDExplicitIntegrator : public virtual sofa::core::objectmodel::BaseObject
{
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

public:
    SOFA_CLASS(PBDExplicitIntegrator, sofa::core::objectmodel::BaseObject);
    PBDExplicitIntegrator();

    virtual ~PBDExplicitIntegrator();

private:
    PBDExplicitIntegrator(const PBDExplicitIntegrator& n) ;
    PBDExplicitIntegrator& operator=(const PBDExplicitIntegrator& n) ;

public:


    //Compute new velocity from external forces : velocity = velocity + dt * Forces_ext
    void integrateExternalForces(const sofa::simulation::Node * gnode,
                                 const sofa::core::MechanicalParams * mparams,
                                 Derivatives& f,
                                 WriteCoord& p,
                                 const WriteCoord& x,
                                 WriteDeriv& v ,
                                 SReal dt);

    void updatePosAndVel(const WriteCoord& p,
                         WriteCoord& x,
                         WriteDeriv& v,
                         const SReal& inv_dt);


protected:
    //GaussSeidelSolver m_solver;
};
#endif //PBDEXPLICITINTEGRATOR_HPP
