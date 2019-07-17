#ifndef PBDVERTEXMass_HPP
#define PBDVERTEXMass_HPP

#include "PBDBaseConstraintData.hpp"
typedef std::vector<SReal> VertexMassData;

template < class T >
class PBDVertexMass : public PBDBaseConstraintData<T>
{
    typedef sofa::component::container::MechanicalObject< T > Mech;
    typedef sofa::core::topology::BaseMeshTopology  Topo;
public:
    PBDVertexMass(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  VertexMassData& m() {return m_mass;}
    inline  VertexMassData& w() {return m_weight;}
    inline  SReal m(uint i) {return m_mass[i];}
    inline  SReal w(uint i) {return m_weight[i];}

private:
    VertexMassData m_mass;
    VertexMassData m_weight;


};

template class PBDVertexMass<sofa::defaulttype::Vec3Types>;
template class PBDVertexMass<sofa::defaulttype::RigidTypes>;

#endif // PBDVERTEXMass_HPP
