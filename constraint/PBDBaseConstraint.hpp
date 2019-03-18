#ifndef PBDANIMATIONLOOP_HPP
#define PBDANIMATIONLOOP_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/simulation/Node.h>

class PBDBaseConstraint : public sofa::core::behavior::BaseObject
{
    //Eigen
    typedef Eigen::Vector3f Vec;
    typedef Eigen::Vector4f Vec4;
    typedef std::vector<Vec> VecList;
    typedef std::vector<Vec4> Vec4List;
    typedef Eigen::Matrix3f Mat;
    typedef Eigen::Matrix4f Mat4;

public:
    SOFA_ABSTRACT_CLASS(BaseAnimationLoop, objectmodel::BaseObject);
    SOFA_BASE_CAST_IMPLEMENTATION(BaseAnimationLoop)

protected:
    BaseAnimationLoop();

    virtual ~BaseAnimationLoop();

    /// Stores starting time of the simulation
    SReal m_resetTime;

    /// Save the initial state for later uses in reset()
    virtual void storeResetState() override;


private:
    BaseAnimationLoop(const BaseAnimationLoop& n) ;
    BaseAnimationLoop& operator=(const BaseAnimationLoop& n) ;

public:
    /// Main computation method.
    ///
    /// Specify and execute all computations for computing a timestep, such
    /// as one or more collisions and integrations stages.
    virtual void solve() = 0;

    /// Returns starting time of the simulation
    SReal getResetTime() const;

    virtual bool insertInNode( objectmodel::BaseNode* node ) override;
    virtual bool removeInNode( objectmodel::BaseNode* node ) override;

};
};
#endif
