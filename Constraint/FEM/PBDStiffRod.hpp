#ifndef PBDSTIFFROD_HPP
#define PBDSTIFFROD_HPP

//#include "./PBDFEMConstraint.hpp"

///*
// * This class implement this paper
// * https://animation.rwth-aachen.de/media/papers/2018-CGF-Rods.pdf
// */

//class PBDStiffRod : public PBDFEMConstraint
//{
//    typedef sofa::defaulttype::Vec3 vec3;
//    typedef Eigen::Matrix<SReal,6,1> Vec6;
//public:
//    PBDStiffRod() : PBDFEMConstraint(),
//        m_radius(initData(&m_radius,0.1,"radius","Diameter of the rod")),
//        m_maxError(initData(&m_maxError,1e-6,"maxError","Max error allowed before goign to the next step"))
//    {}
//    /*
//     * Inputs : PBDObject   -> Object on wich we will solve the constraint
//     *          WriteCoord  -> Free positions on wich we apply the dispalcement
//     *
//     * Output : Solve the constraint adding in WriteCoord the computed displacement
//     */
//    virtual void solve(PBDObject& object, WriteCoord& p) override;

//    /*
//     * Init function of sofa. It's called after the first init of the tree.
//     */
//    virtual void bwdInit () override;

//    /*
//     * Inputs : PBDObject   -> Object on wich we will solve the constraint
//     *
//     * Output : Init all of the datii needed by the constraint (Orientation and FixedPoint)
//     */
//            void initObject(PBDObject& object);

//    /// Construction method called by ObjectFactory.
//    template<class T>
//    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
//    {
//        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
//        if (context) context->addObject(obj);
//        if (arg) obj->parse(arg);
//        return obj;
//    }
//private:
//    /*
//     * Inputs : Quaternionr         -> Orientation of the edge
//     *          Matrix<SReal,4,3>   -> The matrix in wich it will write the data
//     *
//     * Output : The matrix G is initialized as said in the paper.
//     */
//    void computeG(const Quaternionr& q, Eigen::Matrix<SReal,4,3>& G);

//    /*
//     * Inputs : Quaternionr         -> Orientation of the edge
//     *          Matrix<SReal,3,4>   -> The matrix in wich it will write the data
//     *          SReal               -> Average length between the current segment and the next one.
//     *
//     * Output : The matrix J is initialized as said in the paper. (J = grad(C))
//     */
//    void computedOmega(const Quaternionr& q, Eigen::Matrix<SReal,3,4>& J, SReal averageSegmentLength);
//    vec3 m_stretchingAndShearingKs;
//    vec3 m_bendingAndTwistingKs;
//    sofa::core::objectmodel::Data<SReal> m_radius;
//    sofa::core::objectmodel::Data<SReal> m_maxError;
//};

#endif
