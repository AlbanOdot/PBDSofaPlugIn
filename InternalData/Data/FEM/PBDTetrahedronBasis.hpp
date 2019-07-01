#ifndef PBDTETRAHEDRONBASIS_HPP
#define PBDTETRAHEDRONBASIS_HPP

#include "../PBDBaseConstraintData.hpp"



class PBDTetrahedronBasis : public PBDBaseConstraintData<sofa::defaulttype::Vec3Types>
{
    typedef sofa::defaulttype::Vec3 Vec3;
    typedef sofa::defaulttype::Matrix3 Matrix3;
    typedef std::vector<std::pair<SReal,Matrix3>> TetrahedronBasisData;
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

#endif // PBDTETRAHEDRONBASIS_HPP
