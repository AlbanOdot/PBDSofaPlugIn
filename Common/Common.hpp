#ifndef COMMON_H
#define COMMON_H

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

#include <Eigen/Dense>
#include <float.h>
#include <sofa/simulation/Node.h>
//#define USE_DOUBLE
#define MIN_PARALLEL_SIZE 64

typedef SReal Real;

#define REAL_MAX DBL_MAX
#define REAL_MIN DBL_MIN
#define RealParameter DoubleParameter
#define RealParameterType ParameterBase::DOUBLE
#define RealVectorParameterType ParameterBase::VEC_DOUBLE

using Vector2r = Eigen::Matrix<SReal, 2, 1>;
using Vector3r = Eigen::Matrix<SReal, 3, 1>;
using Vector4r = Eigen::Matrix<SReal, 4, 1>;
using Vector5r = Eigen::Matrix<SReal, 5, 1>;
using Vector6r = Eigen::Matrix<SReal, 6, 1>;
using Matrix2r = Eigen::Matrix<SReal, 2, 2>;
using Matrix3r = Eigen::Matrix<SReal, 3, 3>;
using Matrix4r = Eigen::Matrix<SReal, 4, 4>;
using Matrix6r = Eigen::Matrix<SReal, 4, 4>;
using AlignedBox2r = Eigen::AlignedBox<SReal, 2>;
using AlignedBox3r = Eigen::AlignedBox<SReal, 3>;
using AngleAxisr = Eigen::AngleAxis<SReal>;
using Quaternionr = Eigen::Quaternion<SReal>;

//allocators to be used in STL collections containing Eigen structures
using Alloc_Vector2r = Eigen::aligned_allocator<Vector2r>;
using Alloc_Vector3r = Eigen::aligned_allocator<Vector3r>;
using Alloc_Vector4r = Eigen::aligned_allocator<Vector4r>;
using Alloc_Matrix2r = Eigen::aligned_allocator<Matrix2r>;
using Alloc_Matrix3r = Eigen::aligned_allocator<Matrix3r>;
using Alloc_Matrix4r = Eigen::aligned_allocator<Matrix4r>;
using Alloc_AlignedBox2r = Eigen::aligned_allocator<AlignedBox2r>;
using Alloc_AlignedBox3r = Eigen::aligned_allocator<AlignedBox3r>;
using Alloc_AngleAxisr = Eigen::aligned_allocator<AngleAxisr>;
using Alloc_Quaternionr = Eigen::aligned_allocator<Quaternionr>;

#if EIGEN_ALIGN
	#define PDB_MAKE_ALIGNED_OPERATOR_NEW EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	#define REPORT_MEMORY_LEAKS

#if defined(WIN32) || defined(_WIN32) || defined(WIN64)
#ifdef _DEBUG
	// Enable memory leak detection for Eigen new
	#undef PDB_MAKE_ALIGNED_OPERATOR_NEW
	#define PDB_MAKE_ALIGNED_OPERATOR_NEW	EIGEN_MAKE_ALIGNED_OPERATOR_NEW \
		void *operator new(size_t size, int const block_use, char const*  file_name, int const line_number) { \
		\
			return _aligned_malloc_dbg(size, 16, file_name, line_number); \
		} \
		void *operator new[](size_t size, int const block_use, char const*  file_name, int const line_number) { \
			return operator new(size, block_use, file_name, line_number); \
		}\
		void operator delete(void* block, int const block_use, char const*  file_name, int const line_number) { \
		\
			return _aligned_free_dbg(block); \
		} \
		void operator delete[](void* block, int const block_use, char const*  file_name, int const line_number) { \
			return operator delete(block, block_use, file_name, line_number); \
		}	
#endif
#endif
#else
	#define PDB_MAKE_ALIGNED_OPERATOR_NEW

#if defined(WIN32) || defined(_WIN32) || defined(WIN64)
	// Enable memory leak detection
#ifdef _DEBUG
	#define _CRTDBG_MAP_ALLOC 
	#include <stdlib.h>
	#include <crtdbg.h>
	#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__) 	
	#define REPORT_MEMORY_LEAKS _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#else
	#define REPORT_MEMORY_LEAKS
#endif
#else
	#define REPORT_MEMORY_LEAKS
	#define DEBUG_NEW new
#endif

#endif


#endif

#if defined(WIN32) || defined(_WIN32) || defined(WIN64)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE __attribute__((always_inline))
#endif
