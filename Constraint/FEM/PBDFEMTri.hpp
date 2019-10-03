#ifndef PBDFEMTRI_HPP
#define PBDFEMTRI_HPP


#include "PBDFEMConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>
#include "../../InternalData/Data/FEM/PBDTriBasis.hpp"
/**
 * @brief The PBDFEMTri class
 *  This class implement https://sci-hub.tw/10.1016/j.cag.2014.07.004
*/

class PBDFEMTri : public PBDFEMConstraint<sofa::defaulttype::Vec3Types>
{
public:
    PBDFEMTri(sofa::simulation::Node* gnode = nullptr):PBDFEMConstraint<sofa::defaulttype::Vec3Types>(),
        m_shear(initData(&m_shear,Vec3(1,1,1),"shear","shear compliance")),
        m_stretch(initData(&m_stretch,Vec3(1,1,1),"stretch","stretch compliance")),
        m_volumeConservation(initData(&m_volumeConservation,false,"volumeConservation","enforce the constraint to keep the same volume"))
    {}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual bool solve(sofa::simulation::Node * node) override;

    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

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
    Data<bool> m_volumeConservation;
    PBDTriBasis m_basis;
    SReal m_lambda;
    SReal m_mu;
};

#endif // PBDFEMTetra_HPP
