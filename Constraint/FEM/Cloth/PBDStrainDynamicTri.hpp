#ifndef PBDSTRAINDynamicTRI_HPP
#define PBDSTRAINDynamicTRI_HPP


#include "../PBDFEMConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>
#include "../../../InternalData/Data/FEM/PBDTriDynamicBasis.hpp"

/**
 * @brief The PBDStrainDynamic class
 *  This class implement https://sci-hub.tw/10.1016/j.cag.2014.07.004
*/

class PBDStrainDynamicTri : public PBDFEMConstraint<sofa::defaulttype::Vec3Types>
{
    typedef sofa::defaulttype::Vec2 Vec2;
public:
    PBDStrainDynamicTri(sofa::simulation::Node* gnode = nullptr):PBDFEMConstraint<sofa::defaulttype::Vec3Types>(),
        m_sshear(initData(&m_sshear,Vec3(1,1,1),"stretch_and_shear"," stretch along x and y axis and shear along xy axis [stretchX,stretch,Y,shearXY]"))
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
    Data<Vec3> m_sshear;
    Vec3 m_stretch_and_shear;
    PBDTriDynamicBasis m_basis;
    SReal m_lambda;
    SReal m_mu;
};

#endif // PBDSTRAINDynamic_HPP
