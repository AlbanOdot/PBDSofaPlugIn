#ifndef PBDELASTICRODDATA_HPP
#define PBDELASTICRODATA_HPP

#include "./PBDBeamElement.hpp"

class PBDElasticRodData : public PBDBeamElement
{


public:
    PBDElasticRodData(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  std::vector<SReal>  wq()                                        {return m_wq;}
    inline  SReal               wq(uint i)                                  {return m_wq[i];}
            void                applyFixedPoint(const std::vector<uint>& k) {for(auto& i : k){if(i < m_wq.size())m_wq[i] = 0;}}



private:
    std::vector<SReal>      m_wq;

};

typedef PBDElasticRodData ElasticRodData;

#endif // PBDElasticRod_HPP
