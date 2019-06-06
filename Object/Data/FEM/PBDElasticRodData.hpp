#ifndef PBDELASTICRODDATA_HPP
#define PBDELASTICRODATA_HPP

#include "./PBDBeamElement.hpp"

class PBDElasticRodData : public PBDBeamElement
{


public:
    PBDElasticRodData(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  std::vector<SReal>  wq()                                        {return m_wq;}
    inline  SReal               wq(uint i)                                  {return m_wq[i];}
    /*
     * Inputs : vector<uint>        -> Set of index corresponding to the fixed point of PBDFixedPoint constraint
     *
     * Output : Will set to zero the rotation weightso that the quaternion always correspond to the actual edge orientation
     */
            void                applyFixedPoint(const std::vector<uint>& k) {for(auto& i : k){if(i < m_wq.size())m_wq[i] = 0;}}



private:
    std::vector<SReal>      m_wq;

};

typedef PBDElasticRodData ElasticRodData;

#endif // PBDElasticRod_HPP
