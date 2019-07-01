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

class PBDStrainDynamic : public PBDStrainShape
{

public:
    PBDStrainDynamic(sofa::simulation::Node* gnode = NULL):PBDStrainShape(){}
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(PBDObject<sofa::defaulttype::Vec3Types>& object, WriteCoord& p) override;

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
