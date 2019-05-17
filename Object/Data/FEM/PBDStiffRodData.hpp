#ifndef PBDSTIFFRODDATA_HPP
#define PBDSTIFFRODDATA_HPP

#include "./PBDBeamElement.hpp"
#include "../../../Common/Common.hpp"

class PBDStiffRodData : public PBDBeamElement
{
    typedef Eigen::Matrix<SReal,6,1> Vec6;
public:
    PBDStiffRodData(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  std::vector<Vector3r>&      invInertia()                                        {return m_invInertia;}
    inline  Vector3r&                   invInertia(uint i)                                  {return m_invInertia[i];}
    inline  std::vector<Vector3r>&      lambdabt()                                          {return m_lambdabt;}
    inline  Vector3r&                   lambdabt(uint i)                                    {return m_lambdabt[i];}
    inline  std::vector<Vector3r>&      lambdar()                                           {return m_lambdar;}
    inline  Vector3r&                   lambdar(uint i)                                     {return m_lambdar[i];}
    inline  std::vector<Vector3r>&      alpha()                                             {return m_invAlpha;}
    inline  Vector3r&                   alpha(uint i)                                       {return m_invAlpha[i];}
    inline  SReal                       stiffness(uint i)                                   {return m_stiffness[i];}
    inline  std::vector<SReal>&         stiffness()                                         {return m_stiffness;}
            void                    setInvInertia(const std::vector<Vector3r>& invInertiaDiag);
            void                    setInvAlpha(const std::vector<Vector3r>& inv);
            void                    setStiffness(const std::vector<SReal>& stiffness);
            void                    applyFixedPoint(const std::vector<uint>& k)         {for(auto& i : k){if(i < m_invInertia.size())  m_invInertia[i].setZero ();}}

private:
    std::vector<Vector3r>       m_invInertia;//This contain only the bending part of the Inertia matrix
    std::vector<Vector3r>       m_lambdabt;//Bending and twisting part of lambda
    std::vector<Vector3r>       m_lambdar;//Zero stretching part ( R for real)
    std::vector<Vector3r>       m_invAlpha;
    std::vector<SReal>          m_stiffness;
};

typedef PBDStiffRodData StiffRodData;

#endif // PBDStiffRodData_HPP
