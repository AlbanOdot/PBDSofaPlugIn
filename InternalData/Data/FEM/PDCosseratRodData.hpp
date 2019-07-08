#ifndef PDCOSSERATRODDATA_HPP
#define PDCOSSERATRODDATA_HPP

#include "./PBDElasticRodData.hpp"

class PDCosseratRodData : public PBDElasticRodData
{


public:
    PDCosseratRodData(Mech * m = nullptr, Topo* t = nullptr) : PBDElasticRodData(m,t){}

    /*
     * Inputs : SReal   -> Young's modulus
     *          SReal   -> Poisson's ratio
     *          SReal   -> Radius of the rod
     *
     * Output : Setup ws and wbt
     */
            void setupW(const SReal E, const SReal nu, const SReal r);

    virtual void update() override;

    inline  std::vector<SReal>  ws()                                        {return m_ws;}
    inline  SReal               ws(uint i)                                  {return m_ws[i];}
    inline  std::vector<SReal>  wbt()                                       {return m_wbt;}
    inline  SReal               wbt(uint i)                                 {return m_wbt[i];}

private:
    std::vector<SReal>      m_ws;
    std::vector<SReal>      m_wbt;
};


#endif // PBDElasticRod_HPP
