#ifndef PBDBENDINGTOPOLOGY_HPP
#define PBDBENDINGTOPOLOGY_HPP

#include "../PBDBaseConstraintData.hpp"
#include "PBDTriangleAreaTopology.hpp"

typedef std::vector<std::vector<std::pair<uint[2],Eigen::Matrix4d>>> BendingTopologyData;

class PBDBendingTopology : public PBDBaseConstraintData
{
public:
    PBDBendingTopology(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  BendingTopologyData& bendingData()      {return m_bending_topology;}
    inline  TriangleAreaData&    triangleAreaData() {return m_triangle_rest_area.data() ;}

private:

    inline  void computeQ(const sofa::defaulttype::Vec3 *x[], Eigen::Matrix4d &Q,uint i, uint n);

    BendingTopologyData m_bending_topology;
    TriangleAreaTopology m_triangle_rest_area;
};

typedef PBDBendingTopology BendingTopology;

#endif // PBDBendingTOPOLOGY_HPP
