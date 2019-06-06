#ifndef PBDORIENTATION_HPP
#define PBDORIENTATION_HPP

#include "./PBDBaseConstraintData.hpp"
#include "../../Common/Common.hpp"

class PBDOrientation : public PBDBaseConstraintData
{
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

    inline  std::vector<Quaternionr>&   orientation()                                       {return m_orientation;}
    inline  Quaternionr&                orientation(uint i)                                 {return m_orientation[i];}
    inline  std::vector<Quaternionr>&   freeOrientation()                                   {return m_freeOrientation;}
    inline  Quaternionr&                freeOrientation(uint i)                             {return m_freeOrientation[i];}
    inline  std::vector<Quaternionr>&   restDarboux()                                       {return m_restDarboux;}
    inline  Quaternionr&                restDarboux(uint i)                                 {return m_restDarboux[i];}
    inline  std::vector<Vector3r>&      angularSpeed()                                      {return m_angularSpeed;}
    inline  Vector3r&                   angularSpeed(uint i)                                {return m_angularSpeed[i];}
    inline  std::vector<Vector3r>&      torque()                                            {return m_torque;}
    inline  Vector3r&                   torque(uint i)                                      {return m_torque[i];}
    inline  std::vector<Matrix3r>&      inertia()                                           {return m_inertia;}
    inline  Matrix3r&                   inertia(uint i)                                     {return m_inertia[i];}

            void                        setAngularVelocity(const std::vector<Vector3r> &as);
            void                        setTorque(const std::vector<Vector3r>&);
            void                        setInertia(const std::vector<Vector3r>&);


private:
    std::vector<Quaternionr> m_orientation;
    std::vector<Quaternionr> m_freeOrientation;
    std::vector<Quaternionr> m_restDarboux;
    std::vector<Vector3r>    m_angularSpeed;
    std::vector<Vector3r>    m_torque;
    std::vector<Matrix3r>    m_inertia;

};

typedef PBDOrientation Orientation;

#endif // PBDOrientation_HPP
