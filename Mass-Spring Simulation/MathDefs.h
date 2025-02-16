

#if !defined(MATHDEFS_H__INCLUDED_)
#define MATHDEFS_H__INCLUDED_

#define M_PI        3.14159265358979323846f
#define HALF_PI	    1.57079632679489661923f

/// Trig Macros ///////////////////////////////////////////////////////////////
#define DEGTORAD(A)	((A * M_PI) / 180.0f)
#define RADTODEG(A)	((A * 180.0f) / M_PI)
///////////////////////////////////////////////////////////////////////////////

typedef struct
{
	union {
		float x;
		float u;
		float r;
	};
	union {
		float y;
		float v;
		float g;
	};
	union {
		float z;
		float w;
		float b;
	};
} tVector;

// NOT DECLARED AS float[4][4] BECAUSE OPENGL ACCESSES THIS STRANGLY
typedef struct
{
	float m[16];
} tMatrix;

// SOME STRUCTURES TO HELP ME ACCESS VERTEX DATA IN AN ARRAY
typedef struct
{
	float r,g,b;
	float x,y,z;
} tColoredVertex;

typedef struct
{
	float u,v;
	float x,y,z;
} tTexturedVertex;

typedef struct
{
	float u,v;
	float r,g,b;
	float x,y,z;
} tTexturedColoredVertex;

typedef struct
{
	float nx,ny,nz;
	float x,y,z;
} tNormalVertex;

typedef struct
{
	float u,v;
	float nx,ny,nz;
	float x,y,z;
} tTexturedNormalVertex;


/// Quaternion Definitions ////////////////////////////////////////////////////
typedef struct
{
	float x,y,z,w;
} tQuaternion;
///////////////////////////////////////////////////////////////////////////////

#define MAKEVECTOR(a,vx,vy,vz)	a.x = vx; a.y = vy; a.z = vz;

void	MultVectorByMatrix(tMatrix *mat, tVector *v,tVector *result);
double	VectorSquaredLength(tVector *v); 
double	VectorLength(tVector *v); 
void	NormalizeVector(tVector *v); 
double	DotProduct(tVector *v1, tVector *v2);
void	CrossProduct(tVector *v1, tVector *v2, tVector *result);
double	VectorSquaredDistance(tVector *v1, tVector *v2);
void	ScaleVector(tVector *v, float scale, tVector *result);
void	VectorSum(tVector *v1, tVector const *v2, tVector *result);
void	VectorDifference(tVector *v1, tVector const *v2, tVector *result);

#endif // !defined(MATH_H__INCLUDED_)

