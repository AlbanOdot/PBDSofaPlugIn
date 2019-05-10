#ifndef PBDTETRAHEDRONBASIS_HPP
#define PBDTETRAHEDRONBASIS_HPP

#include "../PBDBaseConstraintData.hpp"
#include "../Elastic/PBDTriangleAreaTopology.hpp"

typedef std::vector<std::pair<float,Eigen::Matrix3d>> TetrahedronBasisData;

class PBDTetrahedronBasis : public PBDBaseConstraintData
{
public:
    PBDTetrahedronBasis(Mech * m = nullptr, Topo* t = nullptr);
    virtual void init() override;
    virtual void update() override;

    inline  TetrahedronBasisData&    data() {return m_data ;}

private:
    TetrahedronBasisData m_data;
};

typedef PBDTetrahedronBasis TetrahedronBasis;

#endif // PBDTETRAHEDRONBASIS_HPP
