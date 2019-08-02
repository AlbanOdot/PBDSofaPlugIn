#ifndef PBDCloth_HPP
#define PBDCloth_HPP


#include "../PBDFEMConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>
#include "PBDStrainShapeTri.hpp"

/**
 * @brief The PBDStrainShape class
 *  This class implement https://sci-hub.tw/10.1016/j.cag.2014.07.004
*/

class PBDCloth : public PBDFEMConstraint<sofa::defaulttype::Vec3Types>
{
    typedef sofa::defaulttype::Vec2 Vec2;
public:
    PBDCloth(sofa::simulation::Node* gnode = nullptr):PBDFEMConstraint<sofa::defaulttype::Vec3Types>(),
        m_shear(initData(&m_shear,Vec2(1,1),"shear","shear compliance")),
        m_stretch(initData(&m_stretch,Vec2(1,1),"stretch","stretch compliance"))
    {}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(sofa::simulation::Node * node) override;

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
    Data<Vec2> m_shear;
    Vec3 m_sshear;
    Data<Vec2> m_stretch;
    PBDTriBasis m_basis;
};

#endif // PBDSTRAINSHAPE_HPP
