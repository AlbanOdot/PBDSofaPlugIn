#ifndef PBDOBJECT_HPP
#define PBDOBJECT_HPP

#include <SofaBaseMechanics/MechanicalObject.h>

#include "./Data/PBDVertexMass.hpp"
#include "./Data/PBDOrientation.hpp"

#include "./Data/Elastic/PBDVertexTopology.hpp"
#include "./Data/Elastic/PBDBendingTopology.hpp"

#include "./Data/FEM/PBDTetrahedronBasis.hpp"
#include "./Data/FEM/PBDStiffRodData.hpp"
#include "./Data/FEM/PBDElasticRodData.hpp"

class PBDObject
{
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;


public:
    enum PBDDataType        {STRETCH = 1, BENDING = 2, TETRAHEDRON = 4, MASS = 8, ORIENTED = 16, ELASTICROD = 32, STIFFROD = 64};
    enum PBDIntegrationType {NORMAL = 1,  ANGULAR = 2};

    PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
              sofa::core::topology::BaseMeshTopology * topo = nullptr);

    //Sofa data accessors
    inline const ReadCoord & rest() const                                                                           {return m_rest[0];}
    inline       WriteCoord position()                                                                              {return m_mechanicalObject->writePositions ();}
    inline       WriteDeriv velocity()                                                                              {return m_mechanicalObject->writeVelocities ();}
    inline       unsigned long int dataType()                                                               const   {return m_dataType;}
    inline       unsigned long int integrationType()                                                        const   {return m_integration_type;}
    inline       bool hasDataType(PBDDataType t)                                                            const   {return m_dataType & t;}
    inline       bool integrate(PBDIntegrationType t)                                                       const   {return m_integration_type & t;}

    //PBDDatatype setters
    inline       void setTopology(sofa::core::topology::BaseMeshTopology * topology);
                 void computeStretchTopology();
                 void computeBendingTopology();
                 void computeTetrahedraBasis();
                 void computeElasticRod();
                 void computeOrientation();
                 void computeQ(const sofa::defaulttype::Vec3 *x[4],Eigen::Matrix4d& Q, std::pair<float,float>& area);
                 void computeStiffRod();
                 void setupAngularVelocity(const std::vector<Vector3r>&);
                 void applyFixedPoint(const std::vector<uint>& idx);

     //PBDDatatype accessors
     inline       VertexMassData& mass()                                                                             {return m_mass.m();}
     inline       VertexMassData& invMass()                                                                          {return m_mass.w();}
     inline       SReal& mass(uint i)                                                                                {return m_mass.m(i);}
     inline       SReal& invMass(uint i)                                                                             {return m_mass.w(i);}
     inline       VertexTopology& topology()                                                                         {return m_stretch_topology;}
     inline       BendingTopology& bendTopology()                                                                    {return m_bending_topology;}
     inline       sofa::core::topology::BaseMeshTopology * sofaTopology()                                            {return m_sofa_topology;}
     inline       sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * object()            {return m_mechanicalObject;}
     inline       TetrahedronBasis& tetrahedraBases()                                                                {return m_tetra_bases;}
     inline       PBDOrientation& orientation()                                                                      {return m_orientation;}
     inline       ElasticRodData& elasticRod()                                                                       {return m_elasticRod;}
     inline       StiffRodData& stiffRod()                                                                           {return m_stiffRod;}
protected:

    //General
    PBDOrientation  m_orientation;
    VertexMass      m_mass;

    //Elastic
    VertexTopology  m_stretch_topology;//<<Topology optimized to apply constraints on vertex
    BendingTopology m_bending_topology;

    //FEM
    TetrahedronBasis    m_tetra_bases;
    ElasticRodData          m_elasticRod;
    StiffRodData            m_stiffRod;

    //SOFA
    sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * m_mechanicalObject;
    sofa::core::topology::BaseMeshTopology *                                       m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing
    std::vector<ReadCoord>                                                         m_rest;

    //Truth tellers
    unsigned long int   m_dataType;
    unsigned long int   m_integration_type;

};
#endif
