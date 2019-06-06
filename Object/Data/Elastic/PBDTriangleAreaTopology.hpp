#ifndef PBDTRIANGLEAREATOPOLOGY_HPP
#define PBDTRIANGLEAREATOPOLOGY_HPP

#include "../PBDBaseConstraintData.hpp"
typedef std::vector<std::vector<std::pair<float,float>>> TriangleAreaData;

/*
 * This class compute and store the data needed to solve (for the moment) the bending constraint described int this paper
 * https://cims.nyu.edu/gcl/papers/bergou2006qbm.pdf
 */
class PBDTriangleAreaTopology : public PBDBaseConstraintData
{
public:
    PBDTriangleAreaTopology(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    /*
     * Inputs : vector<vector<uint>> -> A void vector of vector of unsigned int
     *
     * Output : Create a "tree" of neighbors indices
     */
            void initTopology(std::vector<std::vector<uint>>& t);

    inline  TriangleAreaData&   data()                              {return m_data;}

private:

    TriangleAreaData m_data;
};

typedef PBDTriangleAreaTopology TriangleAreaTopology;

#endif // PBDTriangleAreaTOPOLOGY_HPP
