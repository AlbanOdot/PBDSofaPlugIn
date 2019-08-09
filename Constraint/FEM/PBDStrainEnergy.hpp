#ifndef PBDSTRAINENERGY_HPP
#define PBDSTRAINENERGY_HPP

#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>
#include "PBDFEMConstraint.hpp"
#include "../../InternalData/Data/FEM/StrainEnergyData.hpp"
/**
 * @brief The PBDStrainEnergy class
 *  This class implement http://matthias-mueller-fischer.ch/publications/strainBasedDynamics.pdf
*/

class PBDStrainEnergy : public PBDFEMConstraint<sofa::defaulttype::Vec3Types>
{

public:
    PBDStrainEnergy(sofa::simulation::Node* gnode = nullptr):PBDFEMConstraint<sofa::defaulttype::Vec3Types> (gnode),
        m_normalizeShear(initData(&m_normalizeShear,false,"normalizeShear","Normalize shear properties (make the simulation really slow)")),
        m_normalizeStretch(initData(&m_normalizeStretch,false,"normalizeStretch","Normalize stretch properties")),
        m_volumeConservation(initData(&m_volumeConservation,false,"volumeConservation","enforce the constraint to keep the same volume")),
        m_shear(initData(&m_shear,Vec3(1,1,1),"shear","shear compliance")),
        m_stretch(initData(&m_stretch,Vec3(1,1,1),"stretch","stretch compliance"))
        {}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(sofa::simulation::Node * node) override;
    virtual void draw(const sofa::core::visual::VisualParams* vparams) override;
    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    static inline void computeGreenStrainAndPiolaInversion(const Matrix3 &F,Matrix3 &S,Matrix3 &P,SReal restVolume, SReal mu, SReal lambda, SReal& energy,Vec3& stretch, Vec3& shear);

    static inline void computeGreenStrainAndPiola(const Matrix3 &F,Matrix3 &S,Matrix3 &P,SReal restVolume, SReal mu, SReal lambda, SReal& energy,Vec3& stretch, Vec3& shear);

    static inline void computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &sigma, Vec3 *J);

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

protected:
    Data<bool> m_normalizeShear;
    Data<bool> m_normalizeStretch;
    Data<bool> m_volumeConservation;
    Data<Vec3> m_shear;
    Vec3 m_sshear;
    Data<Vec3> m_stretch;
    Vec3 m_sstretch;
    StrainEnergyData m_basis;
    SReal m_lambda;
    SReal m_mu;
};


#endif // PBDStrainEnergy_HPP
