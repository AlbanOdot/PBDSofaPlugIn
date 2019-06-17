#ifndef PBDBEAMELEMENT_HPP
#define PBDBEAMELEMENT_HPP

#include "../PBDBaseConstraintData.hpp"
typedef struct s_infos{uint info[3];} infos;
class PBDBeamElement : public PBDBaseConstraintData
{
public:
    PBDBeamElement(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  uint               beginIdx(uint segment)       {return m_indices[segment].info[0];}
    inline  uint               endIdx(uint segment)         {return m_indices[segment].info[1];}
    inline  uint               segmentIdx(uint segment)     {return m_indices[segment].info[2];}
    inline  SReal              length(uint segment)         {return m_averageLength[segment];}
    inline  std::vector<SReal>& length()                    {return m_averageLength;}

protected:
    std::vector<infos>    m_indices;
    std::vector<SReal>    m_averageLength;
};

#endif // PBDBeamElement_HPP
