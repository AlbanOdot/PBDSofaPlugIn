#ifndef PBDSTIFFRODDATA_HPP
#define PBDSTIFFRODDATA_HPP

#include "./PBDBeamElement.hpp"

class PBDStiffRodData : public PBDBeamElement
{
    typedef Eigen::Matrix<SReal,6,6> Mat6;
    typedef Eigen::Matrix<SReal,6,1> Vec6;
public:
    PBDStiffRodData(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  std::vector<Mat6>&      invMass()                                   {return m_invMass;}
    inline  Mat6&                   invMass(uint i)                             {return m_invMass[i];}
    inline  std::vector<Vec6>&      lambda()                                    {return m_lambda;}
    inline  Vec6&                   lambda(uint i)                              {return m_lambda[i];}
    inline  void                    setInvMass(std::vector<Vec6> invMassDiag);

private:
    std::vector<Mat6>       m_invMass;
    std::vector<Vec6>       m_lambda;
};

typedef PBDStiffRodData StiffRodData;

#endif // PBDStiffRodData_HPP
