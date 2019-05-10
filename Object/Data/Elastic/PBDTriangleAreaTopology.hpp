#ifndef PBDTriangleAreaTOPOLOGY_HPP
#define PBDTriangleAreaTOPOLOGY_HPP

#include "../PBDBaseConstraintData.hpp"
typedef std::vector<std::vector<std::pair<float,float>>> TriangleAreaData;

class PBDTriangleAreaTopology : public PBDBaseConstraintData
{
public:
    PBDTriangleAreaTopology(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  TriangleAreaData&   data()                                              {return m_data;}
            void                initTopology(std::vector<std::vector<uint>>& t);

private:

    TriangleAreaData m_data;
};

typedef PBDTriangleAreaTopology TriangleAreaTopology;

#endif // PBDTriangleAreaTOPOLOGY_HPP
