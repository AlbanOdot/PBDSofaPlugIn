#ifndef PBDBEAMELEMENT_HPP
#define PBDBEAMELEMENT_HPP

#include "../PBDBaseConstraintData.hpp"
typedef struct s_infos{uint info[3];} infos;
class PBDBeamElement : public PBDBaseConstraintData
{
public:
    PBDBeamElement(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  uint    beginIdx(uint segment)      {return m_indicies[segment].info[0];}
    inline  uint    endIdx(uint segment)        {return m_indicies[segment].info[1];}
    inline  uint    segmentIdx(uint segment)    {return m_indicies[segment].info[2];}
    inline  SReal   length(uint segment)        {return m_averageLength[segment];}

private:
    std::vector<infos>    m_indicies;
    std::vector<SReal>    m_averageLength;
};

#endif // PBDBeamElement_HPP
