#ifndef PBDVERTEXTOPOLOGY_HPP
#define PBDVERTEXTOPOLOGY_HPP

#include "../PBDBaseConstraintData.hpp"
typedef std::vector<std::vector<std::pair<uint,SReal>>> VertexTopologyData;

template <class T>
class PBDVertexTopology : public PBDBaseConstraintData<T>
{
    typedef sofa::component::container::MechanicalObject< T > Mech;
    typedef sofa::core::topology::BaseMeshTopology  Topo;
public:
    PBDVertexTopology(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  VertexTopologyData& data() {return m_data;}
    void initTopology(std::vector<std::vector<uint>>& t);

private:
    VertexTopologyData m_data;
};

template class PBDVertexTopology<sofa::defaulttype::Vec3Types>;
template class PBDVertexTopology<sofa::defaulttype::RigidTypes>;

#endif // PBDVERTEXTOPOLOGY_HPP
