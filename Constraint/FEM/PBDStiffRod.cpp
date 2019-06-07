#include "PBDStiffRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDStiffRodClass = sofa::core::RegisterObject("Constraint that correct stiff rod.")
                       .add< PBDStiffRod >();

typedef sofa::defaulttype::Vec3 V3;

void PBDStiffRod::bwdInit ()
{

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    //This constant is here to make the unit of the scene in GPa -> 193 GPA = stainless steel
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*(2.0*mu));
    m_bendingAndTwistingKs /= m_nbIter.getValue ();

}


void PBDStiffRod::solve(PBDObject &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject::STIFFROD) || !object.hasDataType(PBDObject::ORIENTED))
    {
        initObject(object);
    }
    auto& rod = object.stiffRod ();
    auto& u   = object.orientation ().freeOrientation ();
    Matrix6r J1,J2;J1.setIdentity();J2.setIdentity();J2 = -J2;
    Eigen::Matrix<SReal,12,1> rhs; rhs.setZero();

    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint element = 0; element < rod.length().size () - 1; ++element)
        {
            uint begin = rod.beginIdx (element);
            uint end = rod.endIdx (element);
            auto x1 = Vector3r(p[begin][0],p[begin][1],p[begin][2]);
            auto x2 = Vector3r(p[end][0],p[end][1],p[end][2]);

            //Compute the right hand side
            Vector6r C;
            //C(x)
            C.block<3,1>(0,0) = x1-x2;
            C.block<3,1>(3,0) = (2.0/rod.length(begin)) * ((u[begin].conjugate()*u[end]).vec() - object.orientation().restDarboux(begin).vec());
            //C(x) + alpha * lambda
            rhs.block<6,1>(6,0) = C + rod.alpha (begin).asDiagonal ()* rod.lambda (begin);

            //Compute the left hand side
            //Compute J
            //Compute J1
            Eigen::Matrix<SReal, 3, 4> dO;
            computedOmega (u[begin],dO,rod.length (begin));
            Eigen::Matrix<SReal, 4, 3> G;
            computeG (u[begin],G);
            J1.block<3,3>(3,3) = dO*G;

            //Compute the matrix H of Equation 22
            Eigen::Matrix<SReal,12,12> lhs;
            lhs.block<6,6>(0,0) = rod.massMatrix (begin).asDiagonal ();
            lhs.block<6,6>(0,6) = -J1.transpose();
            lhs.block<6,6>(6,0) = -J1;
            lhs.block<6,6>(6,6) = -1.0 * rod.alpha(begin).asDiagonal ();

            //LDLT decomposition for fast solve
            Eigen::Matrix<SReal,12,1> solution(lhs.ldlt ().solve(rhs));

            // x = x + dx
            Vector6r dx (solution.block<6,1>(0,0));
            p[begin] += sofa::defaulttype::Vec3(dx(0),dx(1),dx(2));
            //Quaternionr q(0.0,dx(3),dx(4),dx(5));
            u[begin].coeffs () += G * dx.block<3,1>(0,0);//q.coeffs ();
            u[begin].normalize ();

            //lambda = lambda + dLambda
            rod.lambda (begin) += solution.block<6,1>(6,0); //dLambda


            //Compute J2
            computedOmega (u[end],dO,rod.length (end));
            computeG (u[end],G);
            J2.block<3,3>(3,3) = dO*G;

            //Compute the matrix H of Equation 22
            lhs.block<6,6>(0,0) = rod.massMatrix (end).asDiagonal ();
            lhs.block<6,6>(0,6) = J2.transpose();
            lhs.block<6,6>(6,0) = J2;
            lhs.block<6,6>(6,6) = -1.0 * rod.alpha(end).asDiagonal ();

            //LDLT decomposition for fast solve
            solution =lhs.ldlt ().solve(rhs);

            // x = x + dx
            dx = solution.block<6,1>(0,0);
            p[end] -= sofa::defaulttype::Vec3(dx(0),dx(1),dx(2));
            //Quaternionr q(0.0,dx(3),dx(4),dx(5));
            u[end].coeffs () +=  G * dx.block<3,1>(0,0);//q.coeffs ();
            u[end].normalize ();

            //lambda = lambda + dLambda
            rod.lambda (end) += solution.block<6,1>(6,0); //dLambda

        }
        u[u.size ()-1] = u[u.size()-2];
    }
}

void PBDStiffRod::computeG(const Quaternionr &q, Eigen::Matrix<SReal, 4, 3> &G)
{
    // w component at index 3
    G << static_cast<Real>(0.5)*q.w(),  static_cast<Real>(0.5)*q.z(), -static_cast<Real>(0.5)*q.y(),
            -static_cast<Real>(0.5)*q.z(),  static_cast<Real>(0.5)*q.w(),  static_cast<Real>(0.5)*q.x(),
            static_cast<Real>(0.5)*q.y(), -static_cast<Real>(0.5)*q.x(),  static_cast<Real>(0.5)*q.w(),
            -static_cast<Real>(0.5)*q.x(), -static_cast<Real>(0.5)*q.y(), -static_cast<Real>(0.5)*q.z();
}


void PBDStiffRod::computedOmega(const Quaternionr &q, Eigen::Matrix<SReal, 3, 4> &J, SReal averageSegmentLength)
{
    J <<  q.w(), q.z(),-q.y(),-q.x(),
            -q.z(), q.w(), q.x(),-q.y(),
            q.y(),-q.x(), q.w(),-q.z();
    J *= static_cast<Real>(2.0) / averageSegmentLength;
}

void PBDStiffRod::initObject(PBDObject& object)
{
    if(!object.hasDataType(PBDObject::ORIENTED))
        object.computeOrientation ();
    object.computeStiffRod();

    const auto& m = object.mass ();
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    std::vector<Vector6r> invMass;
    std::vector<Vector6r> alpha;
    SReal stiffness = 1e-8;
    const auto& l = object.stiffRod ().length ();

    for(uint i = 0; i < m.size(); ++i)
    {
        alpha.emplace_back(Vector6r());
        alpha[i] << stiffness,stiffness,stiffness,1.0/(E*I),1.0/(E*I),1.0/(mu*2.0*I);
        alpha[i] /= dt2;
        invMass.emplace_back(Vector6r());
        invMass[i] << object.mass (i), object.mass (i), object.mass (i),l[i]*I,l[i]*I,l[i]*2*I;
    }
    object.stiffRod().setAlpha (alpha);
    object.stiffRod().setMassMatrix(invMass);
    object.stiffRod().applyFixedPoint(m_indices.getValue());
}

/*
uint begin = rod.beginIdx (element);
uint end = rod.endIdx (element);
auto& x1 = p[begin];
auto& x2 = p[end];

//Compute C

//Compute the zero stretch part of the constraint
const auto& stretchViolation = x1 - x2;
Vector3r Cr = Vector3r(stretchViolation[0],stretchViolation[1],stretchViolation[2]);

//Compute the bending and twisting violation part of the constraint
Vector3r Cbt = (2.0/rod.length(element)) * ((u[element].conjugate ()*u[element+1]).vec() - object.orientation().restDarboux(element).vec()).eval();

// fill right hand side of the linear equation system (Equation (19))
Vector6r rhs;
rhs.block<3, 1>(0, 0) = - Cr - rod.stiffness(element) * rod.lambdar(element);
rhs.block<3, 1>(3, 0) = -Cbt - rod.alpha (element).cwiseProduct(rod.lambdabt(element));

// compute matrix of the linear equation system (using Equations (25), (26), and (28) in Equation (19))
Matrix6r JMJT(Matrix6r::Zero());

// compute stretch block
JMJT.block<3, 3>(0, 0) = (object.invMass (begin) + object.invMass (end))*Matrix3r::Identity ();

//Compute J1 (gradC) and J2
Matrix3r J1,J2;
Eigen::Matrix<SReal,4,3> G1,G2;
Eigen::Matrix<SReal,3,4> dOmega1,dOmega2;

computeG (u[begin],G1);
computeJacobian (u[begin],dOmega1,rod.length (element));
J1 = dOmega1 * G1;

computeG (u[end],G2);
computeJacobian (u[end],dOmega2,rod.length (element));
J2 = dOmega2 * G2;

// compute bending and torsion block
//Maybe need to change rod.invInertia to the actual weight (0 or 1)
//            SReal weight1 = rod.invInertia(begin)(0,0) == 0.0 ? 0.0 : 1.0;
//            SReal weight2 = rod.invInertia(end)(0,0) == 0.0 ? 0.0 : 1.0;
Matrix3r MInvJT0(rod.invInertia(begin).asDiagonal()* J1.transpose());
Matrix3r MInvJT1(rod.invInertia(end).asDiagonal()* J2.transpose());

Matrix3r JMJTOmega(Matrix3r::Zero());
if (rod.invInertia(begin)(0,0) != 0.0)
{
    JMJTOmega = J1*MInvJT0;
}

if (rod.invInertia(end)(0,0) != 0.0)
{
    JMJTOmega += J2*MInvJT1;
}
JMJT.block<3, 3>(3, 3) = JMJTOmega;

// add compliance
JMJT(0, 0) += rod.stiffness(element);
JMJT(1, 1) += rod.stiffness(element);
JMJT(2, 2) += rod.stiffness(element);
JMJT(3, 3) += rod.alpha (element)(0);
JMJT(4, 4) += rod.alpha (element)(1);
JMJT(5, 5) += rod.alpha (element)(2);

// solve linear equation system (Equation 19)
auto decomposition(JMJT.ldlt());
Vector6r deltaLambda(decomposition.solve(rhs));
Vector3r dlr(deltaLambda.block<3,1>(0,0));
Vector3r dlbt(deltaLambda.block<3,1>(3,0));

//update Lagrange multiplier
rod.lambdar(element) += dlr;
rod.lambdabt(element) += dlbt;

if (rod.invInertia(begin)(0,0) != 0.0)
{

    Vector3r dx = object.invMass(begin) * dlr;
    if( begin == 95 )
    {
        std::cout << "rhs : "<<(u[element].conjugate ()*u[element+1]).vec() - object.orientation().restDarboux(element).vec()<<std::endl;
    }
    p[begin] += sofa::defaulttype::Vec3(dx(0),dx(1),dx(2));
    u[begin].coeffs () += G1 * (MInvJT0 * dlbt);
    u[begin].normalize ();

}
if (rod.invInertia(end)(0,0) != 0.0)
{
    Vector3r dx = object.invMass(end) * dlr;
    p[end] -= sofa::defaulttype::Vec3(dx(0),dx(1),dx(2));
    u[end].coeffs () += G1 * (MInvJT1 * dlbt);
    u[end].normalize ();
}*/


/*******************************
 *
            //Compute the left hand side and solve
            if(rod.massMatrix (begin)(0,0) > 0.0)
            {
                Matrix6r lhs1(rod.alpha (begin).asDiagonal ());
                Vector6r dlambda1; dlambda1.setZero ();

                //Compute J1
                Eigen::Matrix<SReal, 3, 4> dO;
                computedOmega (u[begin],dO,rod.length (begin));
                Eigen::Matrix<SReal, 4, 3> G;
                computeG (u[begin],G);
                J1.block<3,3>(3,3) = dO*G;

                // JMJT + alpha (alpha is added at the start);
                lhs1 += J1*rod.massMatrix(begin).asDiagonal()*J1.transpose();

                //Solve Eq 19 for dLambda
                dlambda1 = lhs1.ldlt ().solve(rhs);

                //Solve Eq 20
                Vector6r dx = rod.massMatrix (begin).asDiagonal ().inverse () * J1 * dlambda1;

                //Integration
                rod.lambda (begin) += dlambda1;
                p[begin] += sofa::defaulttype::Vec3(dx(0),dx(1),dx(2));
                u[begin].coeffs() += G * dx.block<3,1>(3,0);
                u[begin].normalize ();
            }

            if(rod.massMatrix(end)(0,0) > 0.0)
            {
                Matrix6r lhs2(rod.alpha (begin).asDiagonal ());
                Vector6r dlambda2; dlambda2.setZero ();

                //Compute J2
                Eigen::Matrix<SReal, 3, 4> dO;
                computedOmega (u[end],dO,rod.length (end));
                Eigen::Matrix<SReal, 4, 3> G;
                computeG (u[end],G);
                J1.block<3,3>(3,3) = dO*G;

                // JMJT + alpha (alpha is added at the start);
                lhs2 += J2*rod.massMatrix (end).asDiagonal ()*J2.transpose ();

                //Solve Eq 19 for dLambda
                dlambda2 = lhs2.ldlt ().solve(rhs);

                //Solve Eq 20
                Vector6r dx = rod.massMatrix (end).asDiagonal ().inverse () * J2 * dlambda2;

                //Integration
                rod.lambda (end) += dlambda2;
                p[end] -= sofa::defaulttype::Vec3(dx(0),dx(1),dx(2));
                u[end].coeffs() += G * dx.block<3,1>(3,0);
                u[end].normalize ();
            }



*/
