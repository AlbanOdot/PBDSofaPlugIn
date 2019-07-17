#ifndef PBDSTRAINDYNAMIC_HPP
#define PBDSTRAINDYNAMIC_HPP


#include "PBDStrainShape.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

/**
 * @brief The PBDStrainDynamic class
 *  This class implement http://matthias-mueller-fischer.ch/publications/strainBasedDynamics.pdf
*/

class PBDStrainDynamic : public PBDFEMConstraint<sofa::defaulttype::Vec3Types>
{

public:
    PBDStrainDynamic(sofa::simulation::Node* gnode = NULL):PBDFEMConstraint<sofa::defaulttype::Vec3Types> (){}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(sofa::simulation::Node * node);

    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    static void computeGreenStrainAndPiolaStressInversion(const Matrix3 &F,
            const Real restVolume,
            const Real mu, const Real lambda, Matrix3 &epsilon, Matrix3 &sigma, Real &energy);

    static void computeGreenStrainAndPiolaStress(const Matrix3 &F,
            const Real restVolume,
            const Real mu, const Real lambda, Matrix3 &epsilon, Matrix3 &sigma, Real &energy);

    static void computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &sigma, Vec3 *J);

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
    PBDTetrahedronBasis m_basis;
    SReal m_lambda;
    SReal m_mu;
};


#endif // PBDSTRAINDynamic_HPP
