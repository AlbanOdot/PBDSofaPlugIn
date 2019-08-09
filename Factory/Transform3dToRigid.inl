#ifndef SOFA_COMPONENT_CONTROLLER_Transform3dToRigid_INL
#define SOFA_COMPONENT_CONTROLLER_Transform3dToRigid_INL

#include "Transform3dToRigid.hpp"
#include <sofa/core/visual/VisualParams.h>
#include <SofaOpenglVisual/OglModel.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/Quater.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaOpenglVisual/OglModel.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaConstraint/BilateralInteractionConstraint.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <sofa/helper/AdvancedTimer.h>
#include <iomanip>      // std::setprecision
#include "../Common/Common.hpp"

namespace sofa
{

    namespace core
    {

        namespace behavior
        {

            using namespace core::behavior;
            using namespace sofa::defaulttype;


            Transform3dToRigid::Transform3dToRigid()

                : d_out_pos(initData(&d_out_pos, "out_position", "out position"))
            {
                this->f_listening.setValue(true);
            }

            void Transform3dToRigid::init() {

                this->getContext()->get(m_container);
                if (m_container == NULL) {
                    serr << " Error cannot find the edgeContainer" << sendl;
                    return;
                }


                helper::vector<Rigid> res;
                uint m = m_container->getNbPoints();

                Vector3 from(0,0,1);
                for  (int i=0; i<m-1; i++) {

                    Rigid R;
                    Vector3 P0(m_container->getPX(i),m_container->getPY(i),m_container->getPZ(i));
                    Vector3 P1(m_container->getPX(i+1),m_container->getPY(i+1),m_container->getPZ(i+1));
                    Vector3 to = P1-P0;
                    to.normalize ();
                    Quaternion dq;
                    dq.setFromUnitVectors (from,to);
                    R.getOrientation () = i == 0 ? dq : dq * res[i-1].getOrientation ();
                    R.getCenter () = P0;
                    res.emplace_back(R);
                    from = to;
                }
                Rigid R;
                R.getOrientation () = res[m-2].getOrientation ();
                R.getCenter () = Vector3(m_container->getPX(m-1),m_container->getPY(m-1),m_container->getPZ(m-1));
                res.emplace_back(R);
                d_out_pos.setValue(res);
            }


        } //constraintset

    } //component

}//Sofa

#endif
