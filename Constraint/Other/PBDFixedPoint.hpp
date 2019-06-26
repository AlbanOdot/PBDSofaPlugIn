#ifndef PBDFIXEDPOINT_HPP
#define PBDFIXEDPOINT_HPP
#include "../PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

//template < class DataTypes>
//class PBDFixedPoint : public PBDBaseConstraint<DataTypes>
//{
//    typedef typename DataTypes::Coord       Coord;
//    typedef sofa::helper::vector<Coord>               VecCoord;
//    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
//    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
//    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

//    typedef typename DataTypes::Deriv       Deriv;
//    typedef sofa::helper::vector<Deriv>               VecDeriv;
//    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
//    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
//    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;
//public:
//    SOFA_CLASS(SOFA_TEMPLATE(PBDFixedPoint, DataTypes),SOFA_TEMPLATE(PBDBaseConstraint, DataTypes));
//    PBDFixedPoint(sofa::simulation::Node* gnode = NULL){}
//    PBDFixedPoint(uint objectSize);
//    /*
//     * Inputs : PBDObject   -> Object on wich we will solve the constraint
//     *          WriteCoord  -> Free positions on wich we apply the dispalcement
//     *
//     * Output : Solve the constraint adding in WriteCoord the computed displacement
//     */

//    virtual void solve( PBDObject<DataTypes>& object, WriteCoord& p);

//    /// Construction method called by ObjectFactory.
//    template<class T>
//    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
//    {
//        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
//        if (context) context->addObject(obj);
//        if (arg) obj->parse(arg);
//        return obj;
//    }
//};

//template class PBDFixedPoint<sofa::defaulttype::Vec3Types>;
//template class PBDFixedPoint<sofa::defaulttype::RigidTypes>;

class PBDFixedPoint : public PBDBaseConstraint<sofa::defaulttype::Vec3Types>
{
    typedef typename sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef typename sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;
public:
    PBDFixedPoint(sofa::simulation::Node* gnode = NULL){}
    PBDFixedPoint(uint objectSize);
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */

    virtual void solve( PBDObject<sofa::defaulttype::Vec3Types>& object, WriteCoord& p);

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


class PBDFixedRigidPoint : public PBDBaseConstraint<sofa::defaulttype::Rigid3Types>
{
    typedef typename sofa::defaulttype::RigidTypes::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef typename sofa::defaulttype::RigidTypes::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;
public:
    PBDFixedRigidPoint(sofa::simulation::Node* gnode = NULL){}
    PBDFixedRigidPoint(uint objectSize);
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */

    virtual void solve( PBDObject<sofa::defaulttype::RigidTypes>& object, WriteCoord& p);

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
#endif // PBDFIXEDPOINT_HPP
