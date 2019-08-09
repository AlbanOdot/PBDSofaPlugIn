#ifndef PBDBEAMELEMENT_HPP
#define PBDBEAMELEMENT_HPP

#include "../PBDBaseConstraintData.hpp"
typedef struct s_infos{uint info[3];} infos;

class PBDBeamElement : public PBDBaseConstraintData<sofa::defaulttype::RigidTypes>
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

    enum COLOR{BLACK=true,RED=false};


    inline  SReal               length(uint segment)              {return m_averageLength[segment];}
    inline  std::vector<SReal>& length()                          {return m_averageLength;}
    inline  bool                color(uint segment)      const    {return m_color[segment];}
    inline  SReal               length(uint segment)     const    {return m_averageLength[segment];}
    inline  const std::vector<SReal>& length()                 const    {return m_averageLength;}

protected:
    std::vector<SReal>    m_averageLength;
    std::vector<bool>     m_color;
};

#endif // PBDBeamElement_HPP
