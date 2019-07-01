#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "./Common.hpp"

// ------------------------------------------------------------------------------------

/**
 *
 *
 * This file is entirely copy pasted from (with a few changes)
 *
 * https://github.com/InteractiveComputerGraphics/PositionBasedDynamics
 *
 *
 *
 *
 *
 **/
class MathFunctions
{
    typedef sofa::defaulttype::Matrix3 Matrix3;
    typedef sofa::defaulttype::Vec3 Vec3;
private:
    static void jacobiRotate(Matrix3 &A,
                             Matrix3 &R,
                             int p,
                             int q);

public:
    static Real infNorm(const Matrix3 &A);
    static Real oneNorm(const Matrix3 &A);

    static void eigenDecomposition(const Matrix3 &A,
                                   Matrix3 &eigenVecs,
                                   Vec3 &eigenVals);

    static void polarDecomposition(const Matrix3 &A,
                                   Matrix3 &R,
                                   Matrix3 &U,
                                   Matrix3 &D);

    static void polarDecompositionStable(const Matrix3 &M,
                                         const Real tolerance,
                                         Matrix3 &R);

    static void svdWithInversionHandling(const Matrix3 &A,
                                         Vec3 &sigma,
                                         Matrix3 &U,
                                         Matrix3 &VT);

    static Real cotTheta(const Vec3 &v, const Vec3 &w);

    /** Implementation of the paper: \n
         * Matthias Müller, Jan Bender, Nuttapong Chentanez and Miles Macklin,
         * "A Robust Method to Extract the Rotational Part of Deformations",
         * ACM SIGGRAPH Motion in Games, 2016
         */
    static void extractRotation(const Matrix3 &A, Quaternion &q, const unsigned int maxIter);
};


#endif
