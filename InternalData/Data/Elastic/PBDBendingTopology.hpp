#ifndef PBDBENDINGTOPOLOGY_HPP
#define PBDBENDINGTOPOLOGY_HPP

#include "../PBDBaseConstraintData.hpp"
#include "PBDTriangleAreaTopology.hpp"

typedef std::vector<std::vector<std::pair<uint[2],Eigen::Matrix4d>>> BendingTopologyData;

/*
 * This class compute and store the data needed to solve the bending constraint described int this paper
 * https://cims.nyu.edu/gcl/papers/bergou2006qbm.pdf
 */
class PBDBendingTopology : public PBDBaseConstraintData<sofa::defaulttype::Vec3Types>
{

public:
    PBDBendingTopology(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  BendingTopologyData& bendingData()      {return m_bending_topology;}
    inline  TriangleAreaData&    triangleAreaData() {return m_triangle_rest_area.data() ;}

private:

    /*
     * Inputs : Vec3 *      -> Four vertices describing 2 adjacent triangles
     *          Matrix4d    -> The matrix defining the the bending (it has to bas full of zeros)
     *          uint        -> Index of the edge
     *          uint        -> Position in the neighborhood
     *
     * Output : Compute the matrix Q as described in the paper
     */
    inline  void computeQ(const sofa::defaulttype::Vec3 *x[], Eigen::Matrix4d &Q,uint i, uint n);

    BendingTopologyData m_bending_topology;
    PBDTriangleAreaTopology m_triangle_rest_area;
};


#endif // PBDBendingTOPOLOGY_HPP
