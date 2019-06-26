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
                int m = m_container->getNbPoints();
                for  (int i=0; i<m-1; i++) {

                    Rigid R;
                    Vector3 P0(m_container->getPX(i),m_container->getPY(i),m_container->getPZ(i));
                    Vector3 P1(m_container->getPX(i+1),m_container->getPY(i+1),m_container->getPZ(i+1));
                    Vector3 X = P1-P0;
                    X.normalize();
                    R.getCenter()=P0;
                    Quaternionr q; q.setFromTwoVectors(Vector3r(0,0,1),Vector3r(X[0],X[1],X[2]));
                    q.normalize ();//R.getOrientation().setFromUnitVectors (Vector3(0,0,1),X);
                    R.getOrientation () = Quaternion(q.x (),q.y (),q.z (),q.w ());
                    res.emplace_back(R);

                }

                int i = m-1;
                {
                    Rigid R;
                    R.getCenter()=Vector3(m_container->getPX(i),m_container->getPY(i),m_container->getPZ(i));
                    R.getOrientation() = res[res.size ()-1].getOrientation ();
                    res.emplace_back(R);
                }
                d_out_pos.setValue(res);
            }


        } //constraintset

    } //component

}//Sofa

#endif
