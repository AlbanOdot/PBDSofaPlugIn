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
#include "./Data/FEM/PDCosseratRodData.hpp"

template < class T >
class PBDObject
{
    typedef sofa::defaulttype::Vec3 Vec3;
    typedef typename T::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef typename T::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;


public:
    enum PBDDataType        {STRETCH = 1, BENDING = 2, TETRAHEDRON = 4, MASS = 8, ORIENTED = 16, ELASTICROD = 32, STIFFROD = 64, COSSERATROD = 128};
    enum PBDIntegrationType {NORMAL = 1,  ANGULAR = 2};

    PBDObject(sofa::component::container::MechanicalObject< T > * mobj,
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
    void computeStiffRod();
    void computeCosseratRod();


    //PBDDatatype accessors
    inline       VertexMassData& mass()                                                                             {return m_mass.m();}
    inline       VertexMassData& invMass()                                                                          {return m_mass.w();}
    inline       SReal& mass(uint i)                                                                                {return m_mass.m(i);}
    inline       SReal& invMass(uint i)                                                                             {return m_mass.w(i);}
    inline       PBDVertexTopology<T>& topology()                                                                         {return m_stretch_topology;}
    inline       PBDBendingTopology& bendTopology()                                                                    {return m_bending_topology;}
    inline       sofa::core::topology::BaseMeshTopology * sofaTopology()                                            {return m_sofa_topology;}
    inline       sofa::component::container::MechanicalObject< T > * object()            {return m_mechanicalObject;}
    inline       PBDTetrahedronBasis& tetrahedraBases()                                                                {return m_tetra_bases;}
    inline       PBDOrientation& orientation()                                                                      {return m_orientation;}
    inline       ElasticRodData& elasticRod()                                                                       {return m_elasticRod;}
    inline       StiffRodData& stiffRod()                                                                           {return m_stiffRod;}
    inline       PDCosseratRodData& cosseratRod()                                                                   {return m_PD_CosseratRod;}
protected:

    //General
    PBDOrientation  m_orientation;
    PBDVertexMass<T>      m_mass;

    //Elastic
    PBDVertexTopology<T>  m_stretch_topology;//<<Topology optimized to apply constraints on vertex
    PBDBendingTopology m_bending_topology;

    //FEM
    PBDTetrahedronBasis   m_tetra_bases;
    ElasticRodData      m_elasticRod;
    StiffRodData        m_stiffRod;
    PDCosseratRodData   m_PD_CosseratRod;

    //SOFA
    sofa::component::container::MechanicalObject< T > * m_mechanicalObject;
    sofa::core::topology::BaseMeshTopology *                                       m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing
    std::vector<ReadCoord>                                                         m_rest;

    //Truth tellers
    unsigned long int   m_dataType;
    unsigned long int   m_integration_type;

};

template class PBDObject<sofa::defaulttype::Vec3Types>;
template class PBDObject<sofa::defaulttype::RigidTypes>;
#endif
