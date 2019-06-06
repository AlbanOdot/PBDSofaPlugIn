#ifndef PBDVERTEXTOPOLOGY_HPP
#define PBDVERTEXTOPOLOGY_HPP

#include "../PBDBaseConstraintData.hpp"
typedef std::vector<std::vector<std::pair<uint,SReal>>> VertexTopologyData;

class PBDVertexTopology : public PBDBaseConstraintData
{
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

private:
    VertexTopologyData m_data;
};

typedef PBDVertexTopology VertexTopology;

#endif // PBDVERTEXTOPOLOGY_HPP
