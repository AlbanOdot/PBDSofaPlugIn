#ifndef PBDSTIFFRODDATA_HPP
#define PBDSTIFFRODDATA_HPP

#include "./PBDBeamElement.hpp"

class PBDStiffRodData : public PBDBeamElement
{
    typedef Eigen::Matrix<SReal,6,1> Vec6;
public:
    PBDStiffRodData(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  std::vector<Vec6>&      invMass()                                           {return m_invMass;}
    inline  Vec6&                   invMass(uint i)                                     {return m_invMass[i];}
    inline  std::vector<Vec6>&      lambda()                                            {return m_lambda;}
    inline  Vec6&                   lambda(uint i)                                      {return m_lambda[i];}
    inline  std::vector<Vec6>&      alpha()                                             {return m_invAlpha;}
    inline  Vec6&                   alpha(uint i)                                       {return m_invAlpha[i];}
            void                    setInvMass(const std::vector<Vec6>& invMassDiag);
            void                    setInvAlpha(const std::vector<Vec6>& inv);
            void                    applyFixedPoint(const std::vector<uint>& k)         {for(auto& i : k){if(i < m_invMass.size())  m_invMass[i].setZero ();}}

private:
    std::vector<Vec6>       m_invMass;
    std::vector<Vec6>       m_lambda;
    std::vector<Vec6>       m_invAlpha;
};

typedef PBDStiffRodData StiffRodData;

#endif // PBDStiffRodData_HPP
