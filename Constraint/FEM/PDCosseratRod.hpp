#ifndef PDCOSSERATROD_HPP
#define PDCOSSERATROD_HPP

#include "PBDElasticRod.hpp"
/*
 *
 * This class implement the following paper
 * https://igl.ethz.ch/projects/cosserat-rods/CosseratRods-SCA2018.pdf
 *
 * Cosserat Rods with Projective Dynamics
 */
class PDCosseratRod : public PBDElasticRod
{
    typedef sofa::defaulttype::Vec3 vec3;
public:
    PDCosseratRod() : PBDElasticRod()
    {}

    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(PBDObject<sofa::defaulttype::Rigid3Types>& object, WriteCoord& p) override;

    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    /*
     * Inputs : ElasticRodData          -> Data structure representing the rod informations
     *          vector<Quaternion>     -> Free orientation
     *          vector<Vector3r>        -> Displacement
     *          vector<Quaternion>     -> Rotation correction
     *          PBDObject               -> Object on wich we will solve the constraint
     *          WriteCoord              -> Free positions on wich we apply the dispalcement
     *          uint                    -> index of the segment
     *
     * Output : Compute and apply the correction
     */
    static void solveLinearSystem( PDCosseratRodData& cRod,
                                   std::vector<Quaternion>& u,
                                   std::vector<vec3> dx,
                                   std::vector<Quaternion> dq,
                                   PBDObject<sofa::defaulttype::Rigid3Types>& object,
                                   WriteCoord& p,
                                   const uint e);

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

private:
    std::vector<vec3> m_dx;
    std::vector<Quaternion> m_dq;
};



#endif // PDCosseratRod_HPP
