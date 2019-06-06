#ifndef PBDSTRAINDYNAMIC_HPP
#define PBDSTRAINDYNAMIC_HPP


#include "PBDFEMConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

/**
 * @brief The PBDStrainDynamic class
 *  This class implement http://matthias-mueller-fischer.ch/publications/strainBasedDynamics.pdf
*/

class PBDStrainDynamic : public PBDFEMConstraint
{
public:
    PBDStrainDynamic(sofa::simulation::Node* gnode = NULL):PBDFEMConstraint(){}
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(PBDObject& object, WriteCoord& p) override;

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
};

#endif // PBDSTRAINDynamic_HPP
