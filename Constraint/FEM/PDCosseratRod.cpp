#include "PDCosseratRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PDCosseratRodClass = sofa::core::RegisterObject("Constraint that correct elastic rod.")
                         .add< PDCosseratRod >();

void PDCosseratRod::bwdInit ()
{

    m_mass = PBDVertexMass<sofa::defaulttype::RigidTypes>(m_mechanicalObject.getValue (), m_topology.getValue ());
    m_orientation = PBDOrientation(m_mechanicalObject.getValue(), m_topology.getValue());
    m_cosserat_rod = PDCosseratRodData(m_mechanicalObject.getValue (), m_topology.getValue());

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*2.0*mu);
    m_bendingAndTwistingKs /= m_nbIter.getValue ();

    m_cosserat_rod.applyFixedPoint(m_indices.getValue ());
    m_cosserat_rod.setupW(m_young_modulus.getValue (),m_poisson_ratio.getValue (), m_radius.getValue ());
    m_orientation.setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});

    Quaternionr q; q.coeffs ().setZero ();
    vec3 x(0,0,0);
    for(uint e = 0; e < m_cosserat_rod.wq().size (); ++e )
    {
        m_dq.emplace_back(q);
        m_dx.emplace_back(x);
    }

}

void PDCosseratRod::solve(sofa::simulation::Node * node)
{

    WriteCoordR p = m_pbdObject->getFreePosition ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        int lastBlack;
        //Compute first
        uint a = m_cosserat_rod.beginIdx (0);
        uint z = m_cosserat_rod.endIdx(0);
        m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
        m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
        //Ligne 9
        correction(m_cosserat_rod,
                   m_orientation.freeOrientation(),
                   m_orientation.restDarboux(m_cosserat_rod.beginIdx (0)),
                   m_mass.w (m_cosserat_rod.beginIdx (0)),
                   m_mass.w (m_cosserat_rod.endIdx (0)),
                   p,
                   m_orientation.inertia (0),
                   0);
        //Ligne 10
        solveLinearSystem(m_cosserat_rod,
                          m_orientation.freeOrientation(),
                          m_mass.w (m_cosserat_rod.beginIdx (0)),
                          m_mass.w (m_cosserat_rod.endIdx (0)),
                          p,
                          0);

        for(uint e = 1 ; e < m_cosserat_rod.wq().size (); ++e )
        {
            if(m_cosserat_rod.color (e) == PBDBeamElement::BLACK)
            {
                uint a = m_cosserat_rod.beginIdx (e);
                uint z = m_cosserat_rod.endIdx(e);
                m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
                m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
                //Ligne 9
                correction(m_cosserat_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(m_cosserat_rod.beginIdx (e)),
                           m_mass.w (m_cosserat_rod.beginIdx (e)),
                           m_mass.w (m_cosserat_rod.endIdx (e)),
                           p,
                           m_orientation.inertia (e),
                           e);
                //Ligne 10
                solveLinearSystem(m_cosserat_rod,
                                  m_orientation.freeOrientation(),
                                  m_mass.w (m_cosserat_rod.beginIdx (e)),
                                  m_mass.w (m_cosserat_rod.endIdx (e)),
                                  p,
                                  e);
                lastBlack = e;
            }

        }
        int previousRed = -1;
        for(int e = static_cast<int>(m_cosserat_rod.wq().size ()-1) ; e >= 0; --e )
        {
            if(m_cosserat_rod.color (e) == PBDBeamElement::RED)
            {
                uint a = m_cosserat_rod.beginIdx (e);
                uint z = m_cosserat_rod.endIdx(e);
                m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
                m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
                //Ligne 9
                correction(m_cosserat_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(m_cosserat_rod.beginIdx (e)),
                           m_mass.w (m_cosserat_rod.beginIdx (e)),
                           m_mass.w (m_cosserat_rod.endIdx (e)),
                           p,
                           m_orientation.inertia (e),
                           e);
                //Ligne 10
                solveLinearSystem(m_cosserat_rod,
                                  m_orientation.freeOrientation(),
                                  m_mass.w (m_cosserat_rod.beginIdx (e)),
                                  m_mass.w (m_cosserat_rod.endIdx (e)),
                                  p,
                                  e);

                previousRed = e;
                if(lastBlack >= 0){
                    p[lastBlack].getCenter () += m_dx[lastBlack];
                    m_orientation.freeOrientation()[lastBlack].coeffs() += m_dq[lastBlack].coeffs ();
                    lastBlack -= 2;
                }
                if(previousRed >= 0)
                {
                    p[e].getCenter () += m_dx[previousRed];
                    m_orientation.freeOrientation()[e].coeffs () += m_dq[previousRed].coeffs ();
                }
            }
        }
    }
}


