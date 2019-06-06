#ifndef PBDTETRAHEDRONBASIS_HPP
#define PBDTETRAHEDRONBASIS_HPP

#include "../PBDBaseConstraintData.hpp"

typedef std::vector<std::pair<float,Eigen::Matrix3d>> TetrahedronBasisData;

class PBDTetrahedronBasis : public PBDBaseConstraintData
{
public:
    PBDTetrahedronBasis(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  TetrahedronBasisData&    data() {return m_data ;}

private:
    TetrahedronBasisData m_data;
};

typedef PBDTetrahedronBasis TetrahedronBasis;

#endif // PBDTETRAHEDRONBASIS_HPP
