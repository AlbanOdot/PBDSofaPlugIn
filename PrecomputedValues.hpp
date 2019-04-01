/*
 *
 *
 *
 *
 * This file contains a lot a precomputed values to avoid redundant computational cost
 *
 *
 *
 *
 *
 */
#include <cmath>


//T R I G O N O M E T R Y

namespace trigo{

static const unsigned int size0 = 10000;
static const double pi_over_size = M_PI / (double)size0;
static const double size_over_pi = (double)size0 / M_PI;

double sinV[size0];
double cosV[size0];
double tanV[size0];
double sininV[size0];
double cosinV[size0];
double taninV[size0];

//constexpr void computevalues()
//{
//    for(int i = 0; i < size0; ++i)
//    {
//        double angle = pi_over_size * (double)i;
//        sinV[i] = std::sin(angle);
//        cosV[i] = std::cos(angle);
//        tanV[i] = std::tan(angle);
//        sininV[i] = i > 0 ? 1.f/sinV[i] : 1e100;
//        cosinV[i] = angle != M_PI_2 && angle != 3*M_PI_2 ? 1.0/cosV[i] : 1e100;
//        taninV[i] = i != 0 ? 1.0/tanV[i] : 1e100;
//    }
//}

}//trigo
