#ifndef PBDSTRAINSHAPE_HPP
#define PBDSTRAINSHAPE_HPP


#include "PBDFEMConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>
#include "../../InternalData/Data/FEM/PBDTetrahedronBasis.hpp"

/**
 * @brief The PBDStrainShape class
 *  This class implement https://sci-hub.tw/10.1016/j.cag.2014.07.004
*/

class PBDStrainShape : public PBDFEMConstraint<sofa::defaulttype::Vec3Types>
{
public:
    PBDStrainShape(sofa::simulation::Node* gnode = nullptr):PBDFEMConstraint<sofa::defaulttype::Vec3Types>(),
        m_shear(initData(&m_shear,Vec3(1,1,1),"shear","shear compliance")),
        m_stretch(initData(&m_stretch,Vec3(1,1,1),"stretch","stretch compliance"))
    {}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(sofa::simulation::Node * node) override;

    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    static inline void computeGreenStrainAndPiolaStressInversion(const Matrix3 &F,
                                                                 const Real restVolume,
                                                                 const Real mu, const Real lambda, Matrix3 &epsilon, Matrix3 &sigma, Real &energy);

    static inline void computeGreenStrainAndPiolaStress(const Matrix3 &F,
                                                        const Real restVolume,
                                                        const Real mu, const Real lambda, Matrix3 &epsilon, Matrix3 &sigma, Real &energy);

    static inline void computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &sigma, Vec3 *J, Vec3 stretch, Vec3 shear);

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
    Data<Vec3> m_shear;
    Vec3 m_sshear;
    Data<Vec3> m_stretch;
    PBDTetrahedronBasis m_basis;
    SReal m_lambda;
    SReal m_mu;
};

#endif // PBDSTRAINSHAPE_HPP
