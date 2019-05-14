#include "PBDStiffRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDStiffRodClass = sofa::core::RegisterObject("Constraint that correct stiff rod.")
                       .add< PBDStiffRod >();

typedef sofa::defaulttype::Vec3 V3;

void PBDStiffRod::bwdInit ()
{

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    //This constant is here to make the unit of the scene in GPa -> 193 GPA = stainless steel
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*(2.0*mu));
    m_bendingAndTwistingKs /= m_nbIter.getValue ();

}


void PBDStiffRod::solve(PBDObject &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject::STIFFROD) || !object.hasDataType(PBDObject::ORIENTED))
    {
        initObject(object);
    }
    auto& rod = object.stiffRod ();
    auto& u    = object.orientation ().freeOrientation ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0; e < 1; ++e )
        {

        }
        u[u.size ()-1] = u[u.size()-2];
    }
}


void PBDStiffRod::initObject(PBDObject& object)
{
    if(!object.hasDataType(PBDObject::ORIENTED))
        object.computeOrientation ();
    object.computeStiffRod();

    const auto& m = object.mass ();
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    std::vector<Vec6> Minv;
    std::vector<Vec6> alphainv;
    SReal stiffness = 1e-8;
    const auto& l = object.stiffRod ().length ();

    for(uint i = 0; i < m.size (); ++i)
    {
        alphainv.emplace_back(Vec6());
        alphainv[i] << stiffness,stiffness,stiffness,1.0/(E*I),1.0/(E*I),1.0/(mu*2*I);
        alphainv[i] /= dt2;
        Minv.emplace_back(Vec6());
        Minv[i] << 1.0/m[i],1.0/m[i],1.0/m[i],1.0/(l[i]*I),1.0/(l[i]*I),1.0/(l[i]*2*I);
    }
    object.stiffRod().setInvAlpha(alphainv);
    object.stiffRod().setInvMass (Minv);
    object.stiffRod().applyFixedPoint(m_indices.getValue());
}
