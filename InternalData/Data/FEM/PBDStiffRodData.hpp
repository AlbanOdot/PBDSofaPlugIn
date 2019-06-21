#ifndef PBDSTIFFRODDATA_HPP
#define PBDSTIFFRODDATA_HPP

#include "./PBDBeamElement.hpp"
#include "../../../Common/Common.hpp"

class PBDStiffRodData : public PBDBeamElement
{
    typedef Eigen::Matrix<SReal,6,1> Vec6;
public:
    PBDStiffRodData(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  std::vector<Vector6r>&      massMatrix()                                                {return m_massMatrix;}
    inline  Vector6r&                   massMatrix(uint i)                                          {return m_massMatrix[i];}
    inline  std::vector<Vector6r>&      lambda()                                                    {return m_lambda;}
    inline  Vector6r&                   lambda(uint i)                                              {return m_lambda[i];}
    inline  std::vector<Vector6r>&      alpha()                                                     {return m_alpha;}
    inline  Vector6r&                   alpha(uint i)                                               {return m_alpha[i];}
            void                        setMassMatrix(const std::vector<Vector6r>& massMatrix);
            void                        setAlpha(const std::vector<Vector6r>& inv);
            void                        applyFixedPoint(const std::vector<uint>& k)                 {for(auto& i : k){if(i < m_massMatrix.size())  m_massMatrix[i].setZero ();}}

private:
    std::vector<Vector6r>       m_massMatrix;
    std::vector<Vector6r>       m_lambda;
    std::vector<Vector6r>       m_alpha;
};

typedef PBDStiffRodData StiffRodData;

#endif // PBDStiffRodData_HPP
