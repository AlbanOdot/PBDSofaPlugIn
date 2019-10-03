#include "./PBDCollisionConstraint.hpp"
#include <sofa/core/topology/Topology.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <SofaMeshCollision/TriangleModel.h>
typedef sofa::defaulttype::Vec3Types v3;
typedef sofa::defaulttype::Rigid3Types r3;
using namespace sofa;
using namespace core;
using namespace simulation;
using namespace collision;
using namespace component;
using namespace core::topology;
using namespace component::topology;
using namespace visual;
using namespace defaulttype;
template<class T, class U>
bool PBDCollisionConstraint<T,U>::solve(sofa::simulation::Node *)
{
    if(m_detection == nullptr)
        return false;
    const helper::vector<DetectionOutput>* contacts = dynamic_cast<helper::vector<DetectionOutput>*>(m_detection);
    for( const DetectionOutput& collision : *contacts)
    {
        int type1 = collision.elem.first.getCollisionModel()->getEnumType ();
        int type2 = collision.elem.second.getCollisionModel()->getEnumType ();
        if(type1==type2)
        {
            switch(type1)
            {
            case CollisionModel::POINT_TYPE:
            {
                solvePP (collision);
            break;
            }
            case CollisionModel::LINE_TYPE:
            {
                solveEE(collision);
            break;
            }
            case CollisionModel::TRIANGLE_TYPE:
            {
                solveTT (collision);
            break;
            }
            default:
            break;
            }
        }else{
            switch(type1)
            {
            case CollisionModel::POINT_TYPE:
            {
                std::cout << "Point <-> ";
                switch(type2)
                {
                case CollisionModel::LINE_TYPE:
                    std::cout << "Edge"<<std::endl;
                    solvePE(collision);
                break;
                case CollisionModel::TRIANGLE_TYPE:
                    std::cout << "Triangle"<<std::endl;
                    solvePT (collision);
                break;
                default:
                break;
                }
                break;
            }
            case CollisionModel::LINE_TYPE:
            {
                std::cout << "Edge <-> ";
                switch(type2)
                {
                case CollisionModel::POINT_TYPE:
                    std::cout << "Point"<<std::endl;
                    solvePE (collision,true);
                break;
                case CollisionModel::TRIANGLE_TYPE:
                    std::cout << "Triangle"<<std::endl;
                    solveET (collision);
                break;
                default:
                break;
                }
            break;
            }
            case CollisionModel::TRIANGLE_TYPE:
            {
                std::cout << "Triangle <-> ";
                switch(type2)
                {
                case CollisionModel::POINT_TYPE:
                    std::cout << "Point"<<std::endl;
                    solvePT (collision,true);
                break;
                case CollisionModel::LINE_TYPE:
                    std::cout << "Edge"<<std::endl;
                    solveET(collision,true);
                break;
                default:
                break;
                }
            break;
            }
            default:
            break;
            }
        }
    }
    return true;
}


// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(const Vec3& p, Vec3& a, Vec3& b, Vec3& c, Vec3& bary)
{
    Vec3 v0 = b - a, v1 = c - a, v2 = p - a;
    SReal d00 = dot(v0, v0);
    SReal d01 = dot(v0, v1);
    SReal d11 = dot(v1, v1);
    SReal d20 = dot(v2, v0);
    SReal d21 = dot(v2, v1);
    SReal invDenom = 1.0 / (d00 * d11 - d01 * d01);
    bary[1] = (d11 * d20 - d01 * d21) * invDenom;
    bary[2] = (d00 * d21 - d01 * d20) * invDenom;
    bary[0] = 1.0 - bary[1] - bary[2];
}


/**
 *
 *
 * Point to point collision
 *
 *
 */
template<>
void PBDCollisionConstraint<v3,v3>::solvePP(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
        p1[static_cast<uint>(data.elem.first.getIndex ())] -= displacement;
    if( data.elem.second.getCollisionModel()->isSimulated ())
        p2[static_cast<uint>(data.elem.second.getIndex ())] += displacement;
}

template<>
void PBDCollisionConstraint<r3,r3>::solvePP(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
        p1[static_cast<uint>(data.elem.first.getIndex ())].getCenter () -= displacement;
    if( data.elem.second.getCollisionModel()->isSimulated ())
        p2[static_cast<uint>(data.elem.second.getIndex ())].getCenter () += displacement;
}

template<>
void PBDCollisionConstraint<v3,r3>::solvePP(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
        p1[static_cast<uint>(data.elem.first.getIndex ())] -= displacement;
    if( data.elem.second.getCollisionModel()->isSimulated ())
        p2[static_cast<uint>(data.elem.second.getIndex ())].getCenter () += displacement;
}
template<>
void PBDCollisionConstraint<r3,v3>::solvePP(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
        p1[static_cast<uint>(data.elem.first.getIndex ())].getCenter () -= displacement;
    if( data.elem.second.getCollisionModel()->isSimulated ())
        p2[static_cast<uint>(data.elem.second.getIndex ())] += displacement;

}


/**
 *
 *
 * Edge to edge collision
 *
 *
 */
template<>
void PBDCollisionConstraint<v3,v3>::solveEE(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]] -= displacement;
        p1[e1[1]] -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]] += displacement;
        p2[e2[1]] += displacement;
    }
}

template<>
void PBDCollisionConstraint<r3,r3>::solveEE(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]].getCenter () -= displacement;
        p1[e1[1]].getCenter () -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]].getCenter () += displacement;
        p2[e2[1]].getCenter () += displacement;
    }
}

template<>
void PBDCollisionConstraint<v3,r3>::solveEE(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]] -= displacement;
        p1[e1[1]] -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]].getCenter () += displacement;
        p2[e2[1]].getCenter () += displacement;
    }
}

template<>
void PBDCollisionConstraint<r3,v3>::solveEE(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]].getCenter () -= displacement;
        p1[e1[1]].getCenter () -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]] += displacement;
        p2[e2[1]] += displacement;
    }
}


/**
 *
 *
 * Triangle to triangle collision
 *
 *
 */
template<>
void PBDCollisionConstraint<v3,v3>::solveTT(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getTriangles()[static_cast<uint>(data.elem.first.getIndex ())];
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getTriangles()[static_cast<uint>(data.elem.second.getIndex ())];
    sofa::defaulttype::Vec3 displacement =  data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]] -= displacement;
        p1[e1[1]] -= displacement;
        p1[e1[2]] -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]] += displacement;
        p2[e2[1]] += displacement;
        p2[e2[2]] += displacement;
    }
}

template<>
void PBDCollisionConstraint<r3,r3>::solveTT(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();

    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]].getCenter () -= displacement;
        p1[e1[1]].getCenter () -= displacement;
        p1[e1[2]].getCenter () -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]].getCenter () += displacement;
        p2[e2[1]].getCenter () += displacement;
        p2[e2[2]].getCenter () += displacement;
    }
}

template<>
void PBDCollisionConstraint<v3,r3>::solveTT(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]] -= displacement;
        p1[e1[1]] -= displacement;
        p1[e1[2]] -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]].getCenter () += displacement;
        p2[e2[1]].getCenter () += displacement;
        p2[e2[2]].getCenter () += displacement;
    }
}

template<>
void PBDCollisionConstraint<r3,v3>::solveTT(const sofa::core::collision::DetectionOutput & data)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
    const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if( data.elem.first.getCollisionModel()->isSimulated ())
    {
        p1[e1[0]].getCenter () -= displacement;
        p1[e1[1]].getCenter () -= displacement;
        p1[e1[2]].getCenter () -= displacement;
    }
    if( data.elem.second.getCollisionModel()->isSimulated ())
    {
        p2[e2[0]] += displacement;
        p2[e2[1]] += displacement;
        p2[e2[2]] += displacement;
    }
}

/**
 *
 *
 * Point to edge collision
 *
 *
 */
template<>
void PBDCollisionConstraint<v3,v3>::solvePE(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]] += displacement;
            p2[e2[1]] += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]] -= displacement;
            p1[e1[1]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())] += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<r3,r3>::solvePE(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]].getCenter () += displacement;
            p2[e2[1]].getCenter () += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]].getCenter () -= displacement;
            p1[e1[1]].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())].getCenter () += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<v3,r3>::solvePE(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]].getCenter () += displacement;
            p2[e2[1]].getCenter () += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]] -= displacement;
            p1[e1[1]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())].getCenter () += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<r3,v3>::solvePE(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]] += displacement;
            p2[e2[1]] += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]].getCenter () -= displacement;
            p1[e1[1]].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())] += displacement;
        }
    }
}


/**
 *
 *
 * Point to edge collision
 *
 *
 */
template<>
void PBDCollisionConstraint<v3,v3>::solvePT(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteDeriv v1 = m_obj1->getFreeVelocity ();
    WriteCoord p2 = m_obj2->getFreePosition ();
    WriteDeriv v2 = m_obj2->getFreeVelocity ();

    if(!swap)
    {
////        const SReal bary0 = 1.0 - data.point[0][0] - - data.point[0][0];
////        const sofa::defaulttype::Vec3& v1 = bary0*data.
//        if( data.elem.first.getCollisionModel()->isSimulated ())
//        {
//            p1[static_cast<uint>(data.elem.first.getIndex ())] -= displacement;
//        }
//        if( data.elem.second.getCollisionModel()->isSimulated ())
//        {
//            const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
//            p2[e2[0]] += displacement;
//            p2[e2[1]] += displacement;
//            p2[e2[2]] += displacement;
        //}
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        uint i = static_cast<uint>(data.elem.second.getIndex ());
        Vec3 bary;
        Barycentric(data.point[0],p1[e1[0]],p1[e1[1]],p1[e1[2]],bary);
        SReal w0 = m_obj1->getMass ()->w (i);
        SReal w1 = m_obj1->getMass ()->w (e1[0]);
        SReal w2 = m_obj1->getMass ()->w (e1[1]);
        SReal w3 = m_obj1->getMass ()->w (e1[2]);
        SReal lambda = - data.value / (w0 + (bary[0] * bary[0]) * w1 + (bary[1] * bary[1]) *w2 + (bary[2] * bary[2]) *w3);

        if( data.elem.first.getCollisionModel()->isSimulated ())
        {

            p1[e1[0]] -= (w1*lambda * bary[0]) * data.normal;
            p1[e1[1]] -= (w2*lambda * bary[1]) * data.normal;
            p1[e1[2]] -= (w3*lambda * bary[2]) * data.normal;
            v1[e1[0]] = lambda * bary[0]*0.9*(v1[e1[0]] - data.normal * 2.0*(dot(v1[e1[0]],data.normal)));
            v1[e1[1]] = lambda * bary[1]*0.9*(v1[e1[1]] - data.normal * 2.0*(dot(v1[e1[1]],data.normal)));
            v1[e1[2]] = lambda * bary[2]*0.9*(v1[e1[2]] - data.normal * 2.0*(dot(v1[e1[2]],data.normal)));
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {

            p2[i] += (0.5*w0*lambda) * data.normal;
            v2[i] = 0.9*(v2[i] - data.normal * 2.0*(dot(v2[i],data.normal)));
        }
    }
}

template<>
void PBDCollisionConstraint<r3,r3>::solvePT(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]].getCenter () += displacement;
            p2[e2[1]].getCenter () += displacement;
            p2[e2[2]].getCenter () += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]].getCenter () -= displacement;
            p1[e1[1]].getCenter () -= displacement;
            p1[e1[1]].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())].getCenter () += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<v3,r3>::solvePT(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]].getCenter () += displacement;
            p2[e2[1]].getCenter () += displacement;
            p2[e2[2]].getCenter () += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]] -= displacement;
            p1[e1[1]] -= displacement;
            p1[e1[2]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())].getCenter () += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<r3,v3>::solvePT(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    sofa::defaulttype::Vec3 displacement = (0.5 * (data.value / m_nbIter.getValue())) * data.normal;
    if(!swap)
    {
        const auto& e2 =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[static_cast<uint>(data.elem.first.getIndex ())].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[e2[0]] += displacement;
            p2[e2[1]] += displacement;
            p2[e2[2]] += displacement;
        }
    }else{
        const auto& e1 =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            p1[e1[0]].getCenter () -= displacement;
            p1[e1[1]].getCenter () -= displacement;
            p1[e1[2]].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            p2[static_cast<uint>(data.elem.second.getIndex ())] += displacement;
        }
    }
}

/**
 *
 *
 * Edge to edge collision
 *
 *
 */
template<>
void PBDCollisionConstraint<v3,v3>::solveET(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    SReal w1pw2 = w1 + w2 + 1e-6;
    SReal d = data.value / m_nbIter.getValue();
    SReal W1 = w1 * d /w1pw2;
    SReal W2 = d-W1;
    if(!swap)
    {
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            const auto& e =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[e[0]] -= displacement;
            p1[e[1]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            const auto& t =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[t[0]] += displacement;
            p2[t[1]] += displacement;
            p2[t[2]] += displacement;
        }
    }else{
        std::cout << data.elem.first.getCollisionModel ()->name<<std::endl;
        std::cout << data.elem.second.getCollisionModel ()->name<<std::endl;
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            const auto& t =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[t[0]] -= displacement;
            p1[t[1]] -= displacement;
            p1[t[2]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            const auto& e =  data.elem.second.getCollisionModel()->getTopology()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[e[0]] += displacement;
            p2[e[1]] += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<r3,r3>::solveET(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    SReal w1pw2 = w1 + w2 + 1e-6;
    SReal d = data.value / m_nbIter.getValue();
    SReal W1 = w1 * d /w1pw2;
    SReal W2 = d-W1;
    if(!swap)
    {
        const auto& e =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        const auto& t =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[e[0]].getCenter () -= displacement;
            p1[e[1]].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[t[0]].getCenter () += displacement;
            p2[t[1]].getCenter () += displacement;
            p2[t[2]].getCenter () += displacement;
        }
    }else{
        const auto& e =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        const auto& t =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[t[0]].getCenter() -= displacement;
            p1[t[1]].getCenter() -= displacement;
            p1[t[2]].getCenter() -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[e[0]].getCenter() += displacement;
            p2[e[1]].getCenter() += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<v3,r3>::solveET(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoord p1 = m_obj1->getFreePosition();
    WriteCoordR p2 = m_obj2->getFreePosition ();
    SReal w1pw2 = w1 + w2 + 1e-6;
    SReal d = data.value / m_nbIter.getValue();
    SReal W1 = w1 * d /w1pw2;
    SReal W2 = d-W1;
    if(!swap)
    {
        const auto& e =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        const auto& t =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[e[0]] -= displacement;
            p1[e[1]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[t[0]].getCenter () += displacement;
            p2[t[1]].getCenter () += displacement;
            p2[t[2]].getCenter () += displacement;
        }
    }else{
        const auto& e =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        const auto& t =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[t[0]] -= displacement;
            p1[t[1]] -= displacement;
            p1[t[2]] -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[e[0]].getCenter() += displacement;
            p2[e[1]].getCenter() += displacement;
        }
    }
}

template<>
void PBDCollisionConstraint<r3,v3>::solveET(const sofa::core::collision::DetectionOutput & data,bool swap)
{
    WriteCoordR p1 = m_obj1->getFreePosition();
    WriteCoord p2 = m_obj2->getFreePosition ();
    SReal w1pw2 = w1 + w2 + 1e-6;
    SReal d = data.value / m_nbIter.getValue();
    SReal W1 = w1 * d /w1pw2;
    SReal W2 = d-W1;
    if(!swap)
    {
        const auto& e =  data.elem.first.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.first.getIndex ()));
        const auto& t =  data.elem.second.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.second.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {

            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[e[0]].getCenter () -= displacement;
            p1[e[1]].getCenter () -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[t[0]] += displacement;
            p2[t[1]] += displacement;
            p2[t[2]] += displacement;
        }
    }else{
        const auto& e =  data.elem.second.getCollisionModel()->getTopology ()->getEdge (static_cast<uint>(data.elem.second.getIndex ()));
        const auto& t =  data.elem.first.getCollisionModel()->getTopology()->getTriangle(static_cast<uint>(data.elem.first.getIndex ()));
        if( data.elem.first.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W1 * data.normal;
            p1[t[0]].getCenter() -= displacement;
            p1[t[1]].getCenter() -= displacement;
            p1[t[2]].getCenter() -= displacement;
        }
        if( data.elem.second.getCollisionModel()->isSimulated ())
        {
            sofa::defaulttype::Vec3 displacement =   W2 * data.normal;
            p2[e[0]] += displacement;
            p2[e[1]] += displacement;
        }
    }
}

