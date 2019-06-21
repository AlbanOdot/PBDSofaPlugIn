#ifndef PBDORIENTATION_HPP
#define PBDORIENTATION_HPP

#include "./PBDBaseConstraintData.hpp"
#include "../../Common/Common.hpp"


class PBDOrientation : public PBDBaseConstraintData<sofa::defaulttype::RigidTypes>
{
    typedef sofa::defaulttype::Vec3 Vec3;

public:
    PBDOrientation(Mech * m = nullptr, Topo* t = nullptr);
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() override;

    inline  std::vector<Quaternion>&   freeOrientation()                                   {return m_freeOrientation;}
    inline  Quaternion&                freeOrientation(uint i)                             {return m_freeOrientation[i];}
    inline  std::vector<Quaternion>&   restDarboux()                                       {return m_restDarboux;}
    inline  Quaternion&                restDarboux(uint i)                                 {return m_restDarboux[i];}
    inline  std::vector<Vec3>&      angularSpeed()                                      {return m_angularSpeed;}
    inline  Vec3&                   angularSpeed(uint i)                                {return m_angularSpeed[i];}
    inline  std::vector<Vec3>&      torque()                                            {return m_torque;}
    inline  Vec3&                   torque(uint i)                                      {return m_torque[i];}
    inline  std::vector<Vec3>&      inertia()                                           {return m_inertia;}
    inline  Vec3&                   inertia(uint i)                                     {return m_inertia[i];}

    void                        setAngularVelocity(const std::vector<Vec3> &as);
    void                        setTorque(const std::vector<Vec3>&);
    void                        setInertia(const std::vector<Vec3>&);


private:
    std::vector<Quaternion>  m_freeOrientation;
    std::vector<Quaternion>  m_restDarboux;
    std::vector<Vec3>    m_angularSpeed;
    std::vector<Vec3>    m_torque;
    std::vector<Vec3>    m_inertia;

};

typedef PBDOrientation Orientation;

#endif // PBDOrientation_HPP
