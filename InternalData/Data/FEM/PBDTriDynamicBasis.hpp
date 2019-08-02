#ifndef PBDTRIBASIS_HPP
#define PBDTRIBASIS_HPP

#include "../PBDBaseConstraintData.hpp"



class PBDTriDynamicBasis : public PBDBaseConstraintData<sofa::defaulttype::Vec3Types>
{
    typedef sofa::defaulttype::Vec3 Vec3;
    typedef sofa::defaulttype::Vec2 Vec2;
    typedef sofa::defaulttype::Matrix2 Matrix2;
    typedef std::vector<std::pair<SReal,Matrix2>> TriangleBasisData;
public:
    PBDTriDynamicBasis(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  TriangleBasisData&    data() {return m_data ;}

private:
    TriangleBasisData m_data;
};

#endif // PBDTETRAHEDRONBASIS_HPP
