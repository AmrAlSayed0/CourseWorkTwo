
#include "stdafx.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <assert.h>
#include <math.h>

#include <cmath>
#include <tuple>
#include <iomanip>
#include "Clothy.h"
#include "PhysEnv.h"
#include "SimProps.h"
#include "SetVert.h"
#include "AddSpher.h"

#include "System.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#pragma warning (disable:4244)      // I NEED TO CONVERT FROM DOUBLE TO FLOAT

/////////////////////////////////////////////////////////////////////////////
// CPhysEnv

// INITIALIZE THE SIMULATION WORLD
CPhysEnv::CPhysEnv()
{
	m_IntegratorType = EULER_INTEGRATOR;

	m_Pick[0] = -1;
	m_Pick[1] = -1;
	m_ParticleSys[0] = NULL;
	m_ParticleSys[1] = NULL;
	m_ParticleSys[2] = NULL;	// RESET BUFFER
	// THESE TEMP PARTICLE BUFFERS ARE NEEDED FOR THE MIDPOINT AND RK4 INTEGRATOR
	for (int i = 0; i < 5; i++)
		m_TempSys[i] = NULL;
	m_ParticleCnt = 0;
	m_Contact = NULL;
	m_Spring = NULL;
	m_SpringCnt = 0;
	m_MouseForceActive = FALSE;

	m_UseGravity = TRUE;
	m_DrawSprings = TRUE;
	m_DrawStructural = TRUE;	// By default only draw structural springs
	m_DrawBend = FALSE;
	m_DrawShear = FALSE;
	m_DrawVertices = TRUE;
	m_CollisionActive = TRUE;
	m_CollisionRootFinding = FALSE;		// I AM NOT LOOKING FOR A COLLISION RIGHT AWAY

	MAKEVECTOR(m_Gravity, 0.0f, -0.2f, 0.0f)
		m_UserForceMag = 100.0;
	m_UserForceActive = FALSE;
	m_MouseForceKs = 2.0f;	// MOUSE SPRING CONSTANT
	m_Kd = 0.04f;	// DAMPING FACTOR
	m_Kr = 0.1f;		// 1.0 = SUPERBALL BOUNCE 0.0 = DEAD WEIGHT
	m_Ksh = 5.0f;		// HOOK'S SPRING CONSTANT
	m_Ksd = 0.1f;		// SPRING DAMPING CONSTANT

	// CREATE THE SIZE FOR THE SIMULATION WORLD
	m_WorldSizeX = 15.0f;
	m_WorldSizeY = 15.0f;
	m_WorldSizeZ = 15.0f;

	m_CollisionPlane = (tCollisionPlane*)malloc(sizeof(tCollisionPlane) * 6);
	m_CollisionPlaneCnt = 6;

	// MAKE THE TOP PLANE (CEILING)
	MAKEVECTOR(m_CollisionPlane[0].normal, 0.0f, -1.0f, 0.0f)
		m_CollisionPlane[0].d = m_WorldSizeY / 2.0f;

	// MAKE THE BOTTOM PLANE (FLOOR)
	MAKEVECTOR(m_CollisionPlane[1].normal, 0.0f, 1.0f, 0.0f)
		m_CollisionPlane[1].d = m_WorldSizeY / 2.0f;

	// MAKE THE LEFT PLANE
	MAKEVECTOR(m_CollisionPlane[2].normal, -1.0f, 0.0f, 0.0f)
		m_CollisionPlane[2].d = m_WorldSizeX / 2.0f;

	// MAKE THE RIGHT PLANE
	MAKEVECTOR(m_CollisionPlane[3].normal, 1.0f, 0.0f, 0.0f)
		m_CollisionPlane[3].d = m_WorldSizeX / 2.0f;

	// MAKE THE FRONT PLANE
	MAKEVECTOR(m_CollisionPlane[4].normal, 0.0f, 0.0f, -1.0f)
		m_CollisionPlane[4].d = m_WorldSizeZ / 2.0f;

	// MAKE THE BACK PLANE
	MAKEVECTOR(m_CollisionPlane[5].normal, 0.0f, 0.0f, 1.0f)
		m_CollisionPlane[5].d = m_WorldSizeZ / 2.0f;

	m_SphereCnt = 0;
    testFile.open ( testFileName , ios_base::out );    

}

CPhysEnv::~CPhysEnv()
{
	if (m_ParticleSys[0])
		free(m_ParticleSys[0]);
	if (m_ParticleSys[1])
		free(m_ParticleSys[1]);
	if (m_ParticleSys[2])
		free(m_ParticleSys[2]);
	for (int i = 0; i < 5; i++)
	{
		if (m_TempSys[i])
			free(m_TempSys[i]);
	}
	if (m_Contact)
		free(m_Contact);
	free(m_CollisionPlane);
	free(m_Spring);

	free(m_Sphere);

    if ( testFile.is_open () )
    {
        testFile.flush ();
        testFile.close ();
    }
}

void CPhysEnv::RenderWorld()
{
	tParticle* tempParticle;
	tSpring* tempSpring;

	// FIRST DRAW THE WORLD CONTAINER  
	glColor3f(1.0f, 1.0f, 1.0f);
	// do a big linestrip to get most of edges
	glBegin(GL_LINE_STRIP);
	glVertex3f(-m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glVertex3f(m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glVertex3f(m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glVertex3f(-m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glVertex3f(-m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glVertex3f(-m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glEnd();
	// fill in the stragglers
	glBegin(GL_LINES);
	glVertex3f(m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glVertex3f(m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);

	glVertex3f(m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glVertex3f(m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);

	glVertex3f(-m_WorldSizeX / 2.0f, m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glVertex3f(-m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glEnd();

	// draw floor
	glDisable(GL_CULL_FACE);
	glBegin(GL_QUADS);
	glColor3f(0.0f, 0.0f, 0.5f);
	glVertex3f(-m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glVertex3f(m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, -m_WorldSizeZ / 2.0f);
	glVertex3f(m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glVertex3f(-m_WorldSizeX / 2.0f, -m_WorldSizeY / 2.0f, m_WorldSizeZ / 2.0f);
	glEnd();
	glEnable(GL_CULL_FACE);


	if (m_ParticleSys)
	{
		if (m_Spring && m_DrawSprings)
		{
			glBegin(GL_LINES);
			glColor3f(0.0f, 0.8f, 0.8f);
			tempSpring = m_Spring;
			for (int loop = 0; loop < m_SpringCnt; loop++)
			{
				// Only draw normal springs or the cloth "structural" ones
				if ((tempSpring->type == MANUAL_SPRING) ||
					(tempSpring->type == STRUCTURAL_SPRING && m_DrawStructural) ||
					(tempSpring->type == SHEAR_SPRING && m_DrawShear) ||
					(tempSpring->type == BEND_SPRING && m_DrawBend))
				{
					glVertex3fv((float*)&m_CurrentSys[tempSpring->p1].pos);
					glVertex3fv((float*)&m_CurrentSys[tempSpring->p2].pos);
				}
				tempSpring++;
			}
			if (m_MouseForceActive)	// DRAW MOUSESPRING FORCE
			{
				if (m_Pick[0] > -1)
				{
					glColor3f(0.8f, 0.0f, 0.8f);
					glVertex3fv((float*)&m_CurrentSys[m_Pick[0]].pos);
					glVertex3fv((float*)&m_MouseDragPos[0]);
				}
				if (m_Pick[1] > -1)
				{
					glColor3f(0.8f, 0.0f, 0.8f);
					glVertex3fv((float*)&m_CurrentSys[m_Pick[1]].pos);
					glVertex3fv((float*)&m_MouseDragPos[1]);
				}
			}
			glEnd();
		}
		if (m_DrawVertices)
		{
			glBegin(GL_POINTS);
			tempParticle = m_CurrentSys;
			for (int loop = 0; loop < m_ParticleCnt; loop++)
			{
				if (loop == m_Pick[0])
					glColor3f(0.0f, 0.8f, 0.0f);
				else if (loop == m_Pick[1])
					glColor3f(0.8f, 0.0f, 0.0f);
				else
					glColor3f(0.8f, 0.8f, 0.0f);
				glVertex3fv((float*)&tempParticle->pos);
				tempParticle++;
			}
			glEnd();
		}
	}

	if (m_SphereCnt > 0 && m_CollisionActive)
	{
		glColor3f(0.5f, 0.0f, 0.0f);
		for (int loop = 0; loop < m_SphereCnt; loop++)
		{
			glPushMatrix();
			glTranslatef(m_Sphere[loop].pos.x, m_Sphere[loop].pos.y, m_Sphere[loop].pos.z);
			glScalef(m_Sphere[loop].radius, m_Sphere[loop].radius, m_Sphere[loop].radius);
			glCallList(OGL_AXIS_DLIST);
			glPopMatrix();
		}
	}
}

void CPhysEnv::GetNearestPoint(int x, int y)
{
	/// Local Variables ///////////////////////////////////////////////////////////
	float* feedBuffer;
	int hitCount;
	tParticle* tempParticle;
	int loop;
	///////////////////////////////////////////////////////////////////////////////
		// INITIALIZE A PLACE TO PUT ALL THE FEEDBACK INFO (3 DATA, 1 TAG, 2 TOKENS)
	feedBuffer = (float*)malloc(sizeof(GLfloat) * m_ParticleCnt * 6);
	// TELL OPENGL ABOUT THE BUFFER
	glFeedbackBuffer(m_ParticleCnt * 6, GL_3D, feedBuffer);
	(void)glRenderMode(GL_FEEDBACK);	// SET IT IN FEEDBACK MODE

	tempParticle = m_CurrentSys;
	for (loop = 0; loop < m_ParticleCnt; loop++)
	{
		// PASS THROUGH A MARKET LETTING ME KNOW WHAT VERTEX IT WAS
		glPassThrough((float)loop);
		// SEND THE VERTEX
		glBegin(GL_POINTS);
		glVertex3fv((float*)&tempParticle->pos);
		glEnd();
		tempParticle++;
	}
	hitCount = glRenderMode(GL_RENDER); // HOW MANY HITS DID I GET
	//fout<<"hit count "<<hitCount<,endl;

	CompareBuffer(hitCount, feedBuffer, (float)x, (float)y);		// CHECK THE HIT 
	free(feedBuffer);		// GET RID OF THE MEMORY
}

///////////////////////////////////////////////////////////////////////////////
// Function:	CompareBuffer
// Purpose:		Check the feedback buffer to see if anything is hit
// Arguments:	Number of hits, pointer to buffer, point to test
///////////////////////////////////////////////////////////////////////////////
void CPhysEnv::CompareBuffer(int size, float* buffer, float x, float y)
{
	/// Local Variables ///////////////////////////////////////////////////////////
	GLint count;
	GLfloat token, point[3];
	int loop, currentVertex, result = -1;
	long nearest = -1, dist;
	///////////////////////////////////////////////////////////////////////////////
	count = size;
	while (count)
	{
		token = buffer[size - count];	// CHECK THE TOKEN
		count--;
		if (token == GL_PASS_THROUGH_TOKEN)	// VERTEX MARKER
		{
			currentVertex = (int)buffer[size - count]; // WHAT VERTEX
			count--;
		}
		else if (token == GL_POINT_TOKEN)
		{
			// THERE ARE THREE ELEMENTS TO A POINT TOKEN
			for (loop = 0; loop < 3; loop++)
			{
				point[loop] = buffer[size - count];
				count--;
			}
			dist = ((x - point[0]) * (x - point[0])) + ((y - point[1]) * (y - point[1]));
			if (result == -1 || dist < nearest)
			{
				nearest = dist;
				result = currentVertex;
			}
		}
	}

	if (nearest < 50.0f)
	{
		if (m_Pick[0] == -1)
			m_Pick[0] = result;
		else if (m_Pick[1] == -1)
			m_Pick[1] = result;
		else
		{
			m_Pick[0] = result;
			m_Pick[1] = -1;
		}
	}
}
////// CompareBuffer //////////////////////////////////////////////////////////
std::string CPhysEnv::ParticleCsvLine ( tParticle * particle )
{
    std::stringstream ss;
    ss << std::showpos << std::setprecision ( 10 );
    ss << particle->pos.x;
    ss << ",";
    ss << particle->pos.y;
    ss << ",";
    ss << particle->pos.z;
    ss << ",";
    ss << particle->v.x;
    ss << ",";
    ss << particle->v.y;
    ss << ",";
    ss << particle->v.z;
    ss << ",";
    ss << particle->f.x;
    ss << ",";
    ss << particle->f.y;
    ss << ",";
    ss << particle->f.z;
    return ss.str ();
}
//std::tuple < float , float , float > CPhysEnv::CalculateError ( bool reverse ) const
//{
//    if ( ! OUTPUT_TO_FILE )
//    {
//        return {
//            0 ,
//            0 ,
//            0
//        };
//    }
//    tParticle* target;
//    tParticle* current;
//    if ( reverse )
//    {
//        current = m_TargetSys;
//        target = m_CurrentSys;
//    }
//    else
//    {
//        current = m_CurrentSys;
//        target = m_TargetSys;
//    }
//    /* Calculate average error for positions over all the particles and dimensions */
//    float positionApproximateError = 0.0f;
//    int numOfPositionSums = 0;
//    for ( int i = 0; i < this->m_ParticleCnt; ++i )
//    {
//        if ( target [ i ].pos.x != 0 )
//        {
//            positionApproximateError += std::fabs ( ( target [ i ].pos.x - current [ i ].pos.x ) / target [ i ].pos.x );
//            ++numOfPositionSums;
//        }
//        if ( target [ i ].pos.y != 0 )
//        {
//            positionApproximateError += std::fabs ( ( target [ i ].pos.y - current [ i ].pos.y ) / target [ i ].pos.y );
//        }
//        if ( target [ i ].pos.z != 0 )
//        {
//            positionApproximateError += std::fabs ( ( target [ i ].pos.z - current [ i ].pos.z ) / target [ i ].pos.z );
//            ++numOfPositionSums;
//        }
//    }
//    float positionRelativeApproximateError;
//    if ( numOfPositionSums == 0 )
//    {
//        positionRelativeApproximateError = 0;
//    }
//    else
//    {
//        positionRelativeApproximateError = positionApproximateError / numOfPositionSums;
//    }
//    /* Calculate average error for velocities over all the particles and dimensions */
//    float velocityApproximateError = 0.0f;
//    int numOfVelocitySums = 0;
//    for ( int i = 0; i < this->m_ParticleCnt; ++i )
//    {
//        if ( target [ i ].v.x != 0 )
//        {
//            velocityApproximateError += std::fabs ( ( target [ i ].v.x - current [ i ].v.x ) / target [ i ].v.x );
//            ++numOfVelocitySums;
//        }
//        if ( target [ i ].v.y != 0 )
//        {
//            velocityApproximateError += std::fabs ( ( target [ i ].v.y - current [ i ].v.y ) / target [ i ].v.y );
//        }
//        if ( target [ i ].v.z != 0 )
//        {
//            velocityApproximateError += std::fabs ( ( target [ i ].v.z - current [ i ].v.z ) / target [ i ].v.z );
//            ++numOfVelocitySums;
//        }
//    }
//    float velocityRelativeApproximateError;
//    if ( numOfVelocitySums == 0 )
//    {
//        velocityRelativeApproximateError = 0;
//    }
//    else
//    {
//        velocityRelativeApproximateError = velocityApproximateError / numOfVelocitySums;
//    }
//    /* Calculate average error for forces over all the particles and dimensions */
//    float forceApproximateError = 0.0f;
//    int numOfForceSums = 0;
//    for ( int i = 0; i < this->m_ParticleCnt; ++i )
//    {
//        if ( target [ i ].f.x != 0 )
//        {
//            forceApproximateError += std::fabs ( ( target [ i ].f.x - current [ i ].f.x ) / target [ i ].f.x );
//            ++numOfForceSums;
//        }
//        if ( target [ i ].f.y != 0 )
//        {
//            forceApproximateError += std::fabs ( ( target [ i ].f.y - current [ i ].f.y ) / target [ i ].f.y );
//        }
//        if ( target [ i ].f.z != 0 )
//        {
//            forceApproximateError += std::fabs ( ( target [ i ].f.z - current [ i ].f.z ) / target [ i ].f.z );
//            ++numOfForceSums;
//        }
//    }
//    float forceRelativeApproximateError;
//    if ( numOfForceSums == 0 )
//    {
//        forceRelativeApproximateError = 0;
//    }
//    else
//    {
//        forceRelativeApproximateError = forceApproximateError / numOfForceSums;
//    }
//    /* Return the result as a tuple */
//    return {
//        positionRelativeApproximateError ,
//        velocityRelativeApproximateError ,
//        forceRelativeApproximateError
//    };
//}
std::tuple < float , float , float > CPhysEnv::CalculateError ( bool reverse ) const
{
    if ( ! OUTPUT_TO_FILE )
    {
        return {
            0 ,
            0 ,
            0
        };
    }
    tParticle* target;
    tParticle* current;
    if ( reverse )
    {
        current = m_TargetSys;
        target = m_CurrentSys;
    }
    else
    {
        current = m_CurrentSys;
        target = m_TargetSys;
    }
    float posX = 0.0f;
    float posY = 0.0f;
    float posZ = 0.0f;
    for ( int i = 0; i < this->m_ParticleCnt; ++i )
    {
        posX += current [ i ].pos.x;
        posY += current [ i ].pos.y;
        posZ += current [ i ].pos.z;
    }
    float averagePosition = std::sqrt ( std::pow ( posX / m_ParticleCnt , 2 ) + std::pow ( posY / m_ParticleCnt , 2 ) + std::pow ( posZ / m_ParticleCnt , 2 ) );
    /* Calculate average error for velocities over all the particles and dimensions */
    float vX = 0.0f;
    float vY = 0.0f;
    float vZ = 0.0f;
    for ( int i = 0; i < this->m_ParticleCnt; ++i )
    {
        vX += current [ i ].v.x;
        vY += current [ i ].v.y;
        vZ += current [ i ].v.z;
    }
    float averageVelocity = std::sqrt ( std::pow ( vX / m_ParticleCnt , 2 ) + std::pow ( vY / m_ParticleCnt , 2 ) + std::pow ( vZ / m_ParticleCnt , 2 ) );
    float fX = 0.0f;
    float fY = 0.0f;
    float fZ = 0.0f;
    for ( int i = 0; i < this->m_ParticleCnt; ++i )
    {
        fX += current [ i ].f.x;
        fY += current [ i ].f.y;
        fZ += current [ i ].f.z;
    }
    float averageForce = std::sqrt ( std::pow ( fX / m_ParticleCnt , 2 ) + std::pow ( fY / m_ParticleCnt , 2 ) + std::pow ( fZ / m_ParticleCnt , 2 ) );
    /* Return the result as a tuple */
    return {
        averagePosition ,
        averageVelocity ,
        averageForce
    };
}
void CPhysEnv::OutputErrorToCsV ( std::tuple < float , float , float > error , float time )
{
    std::stringstream ss;
    ss << std::showpos << std::setprecision ( std::numeric_limits < float >::digits + 1 );
    ss << std::get < 0 > ( error );
    ss << ",";
    ss << std::get < 1 > ( error );
    ss << ",";
    ss << std::get < 2 > ( error );
    ss << ",";
    ss << time;
    ss << ",";
    switch ( m_IntegratorType )
    {
        case EULER_INTEGRATOR: ss << "EULER";
            break;
        case MIDPOINT_INTEGRATOR: ss << "MIDPOINT";
            break;
        case HEUN_INTEGRATOR: ss << "HEUN";
            break;
        case RK4_INTEGRATOR: ss << "RK4";
            break;
        case RK5_INTEGRATOR: ss << "RK5";
            break;
        case RK4_ADAPTIVE_INTEGRATOR: ss << "RK4_ADAPTIVE";
            break;
        default: ss << "DEFAULT";
    }
    ss << '\n';
    testFile << ss.rdbuf ();

}
void CPhysEnv::SetWorldParticles(tTexturedVertex* coords, int particleCnt)
{
	tParticle* tempParticle;

	if (m_ParticleSys[0])
		free(m_ParticleSys[0]);
	if (m_ParticleSys[1])
		free(m_ParticleSys[1]);
	if (m_ParticleSys[2])
		free(m_ParticleSys[2]);
	for (int i = 0; i < 5; i++)
	{
		if (m_TempSys[i])
			free(m_TempSys[i]);
	}
	if (m_Contact)
		free(m_Contact);
	// THE SYSTEM IS DOUBLE BUFFERED TO MAKE THINGS EASIER
	m_CurrentSys = (tParticle*)malloc(sizeof(tParticle) * particleCnt);
	m_TargetSys = (tParticle*)malloc(sizeof(tParticle) * particleCnt);

	m_ParticleSys[2] = (tParticle*)malloc(sizeof(tParticle) * particleCnt);
	for (int i = 0; i < 5; i++)
	{
		m_TempSys[i] = (tParticle*)malloc(sizeof(tParticle) * particleCnt);
	}
	m_ParticleCnt = particleCnt;

	// MULTIPLIED PARTICLE COUNT * 2 SINCE THEY CAN COLLIDE WITH MULTIPLE WALLS
	m_Contact = (tContact*)malloc(sizeof(tContact) * particleCnt * 2);
	m_ContactCnt = 0;

	tempParticle = m_CurrentSys;
	for (int loop = 0; loop < particleCnt; loop++)
	{
		MAKEVECTOR(tempParticle->pos, coords->x, coords->y, coords->z)
			MAKEVECTOR(tempParticle->v, 0.0f, 0.0f, 0.0f)
			MAKEVECTOR(tempParticle->f, 0.0f, 0.0f, 0.0f)
			tempParticle->oneOverM = 1.0f;							// MASS OF 1
		tempParticle++;
		coords++;
	}

	// COPY THE SYSTEM TO THE SECOND ONE ALSO
	memcpy(m_TargetSys, m_CurrentSys, sizeof(tParticle) * particleCnt);
	// COPY THE SYSTEM TO THE RESET BUFFER ALSO
	memcpy(m_ParticleSys[2], m_CurrentSys, sizeof(tParticle) * particleCnt);

	m_ParticleSys[0] = m_CurrentSys;
	m_ParticleSys[1] = m_TargetSys;
}

///////////////////////////////////////////////////////////////////////////////
// Function:	FreeSystem
// Purpose:		Remove all particles and clear it out
///////////////////////////////////////////////////////////////////////////////
void CPhysEnv::FreeSystem()
{
	m_Pick[0] = -1;
	m_Pick[1] = -1;
	if (m_ParticleSys[0])
	{
		m_ParticleSys[0] = NULL;
		free(m_ParticleSys[0]);
	}
	if (m_ParticleSys[1])
	{
		free(m_ParticleSys[1]);
		m_ParticleSys[1] = NULL;
	}
	if (m_ParticleSys[2])
	{
		free(m_ParticleSys[2]);
		m_ParticleSys[2] = NULL;	// RESET BUFFER
	}
	for (int i = 0; i < 5; i++)
	{
		if (m_TempSys[i])
		{
			free(m_TempSys[i]);
			m_TempSys[i] = NULL;	// RESET BUFFER
		}
	}
	if (m_Contact)
	{
		free(m_Contact);
		m_Contact = NULL;
	}
	if (m_Spring)
	{
		free(m_Spring);
		m_Spring = NULL;
	}
	m_SpringCnt = 0;
	m_ParticleCnt = 0;
}
////// FreeSystem //////////////////////////////////////////////////////////////

void CPhysEnv::LoadData(FILE* fp)
{
	fread(&m_UseGravity, sizeof(BOOL), 1, fp);
	fread(&m_UseDamping, sizeof(BOOL), 1, fp);
	fread(&m_UserForceActive, sizeof(BOOL), 1, fp);
	fread(&m_Gravity, sizeof(tVector), 1, fp);
	fread(&m_UserForce, sizeof(tVector), 1, fp);
	fread(&m_UserForceMag, sizeof(float), 1, fp);
	fread(&m_Kd, sizeof(float), 1, fp);
	fread(&m_Kr, sizeof(float), 1, fp);
	fread(&m_Ksh, sizeof(float), 1, fp);
	fread(&m_Ksd, sizeof(float), 1, fp);
	fread(&m_ParticleCnt, sizeof(int), 1, fp);
	m_CurrentSys = (tParticle*)malloc(sizeof(tParticle) * m_ParticleCnt);
	m_TargetSys = (tParticle*)malloc(sizeof(tParticle) * m_ParticleCnt);
	m_ParticleSys[2] = (tParticle*)malloc(sizeof(tParticle) * m_ParticleCnt);
	for (int i = 0; i < 5; i++)
	{
		m_TempSys[i] = (tParticle*)malloc(sizeof(tParticle) * m_ParticleCnt);
	}
	m_ParticleSys[0] = m_CurrentSys;
	m_ParticleSys[1] = m_TargetSys;
	m_Contact = (tContact*)malloc(sizeof(tContact) * m_ParticleCnt);
	fread(m_ParticleSys[0], sizeof(tParticle), m_ParticleCnt, fp);
	fread(m_ParticleSys[1], sizeof(tParticle), m_ParticleCnt, fp);
	fread(m_ParticleSys[2], sizeof(tParticle), m_ParticleCnt, fp);
	fread(&m_SpringCnt, sizeof(int), 1, fp);
	m_Spring = (tSpring*)malloc(sizeof(tSpring) * (m_SpringCnt));
	fread(m_Spring, sizeof(tSpring), m_SpringCnt, fp);
	fread(m_Pick, sizeof(int), 2, fp);
	fread(&m_SphereCnt, sizeof(int), 1, fp);
	m_Sphere = (tCollisionSphere*)malloc(sizeof(tCollisionSphere) * (m_SphereCnt));
	fread(m_Sphere, sizeof(tCollisionSphere), m_SphereCnt, fp);
}

void CPhysEnv::SaveData(FILE* fp)
{
	fwrite(&m_UseGravity, sizeof(BOOL), 1, fp);
	fwrite(&m_UseDamping, sizeof(BOOL), 1, fp);
	fwrite(&m_UserForceActive, sizeof(BOOL), 1, fp);
	fwrite(&m_Gravity, sizeof(tVector), 1, fp);
	fwrite(&m_UserForce, sizeof(tVector), 1, fp);
	fwrite(&m_UserForceMag, sizeof(float), 1, fp);
	fwrite(&m_Kd, sizeof(float), 1, fp);
	fwrite(&m_Kr, sizeof(float), 1, fp);
	fwrite(&m_Ksh, sizeof(float), 1, fp);
	fwrite(&m_Ksd, sizeof(float), 1, fp);
	fwrite(&m_ParticleCnt, sizeof(int), 1, fp);
	fwrite(m_ParticleSys[0], sizeof(tParticle), m_ParticleCnt, fp);
	fwrite(m_ParticleSys[1], sizeof(tParticle), m_ParticleCnt, fp);
	fwrite(m_ParticleSys[2], sizeof(tParticle), m_ParticleCnt, fp);
	fwrite(&m_SpringCnt, sizeof(int), 1, fp);
	fwrite(m_Spring, sizeof(tSpring), m_SpringCnt, fp);
	fwrite(m_Pick, sizeof(int), 2, fp);
	fwrite(&m_SphereCnt, sizeof(int), 1, fp);
	fwrite(m_Sphere, sizeof(tCollisionSphere), m_SphereCnt, fp);
}

// RESET THE SIM TO INITIAL VALUES
void CPhysEnv::ResetWorld()
{
	memcpy(m_CurrentSys, m_ParticleSys[2], sizeof(tParticle) * m_ParticleCnt);
	memcpy(m_TargetSys, m_ParticleSys[2], sizeof(tParticle) * m_ParticleCnt);
}

void CPhysEnv::SetWorldProperties()
{
	CSimProps	dialog;
	dialog.m_CoefRest = m_Kr;
	dialog.m_Damping = m_Kd;
	dialog.m_GravX = m_Gravity.x;
	dialog.m_GravY = m_Gravity.y;
	dialog.m_GravZ = m_Gravity.z;
	dialog.m_SpringConst = m_Ksh;
	dialog.m_SpringDamp = m_Ksd;
	dialog.m_UserForceMag = m_UserForceMag;
	if (dialog.DoModal() == IDOK)
	{
		m_Kr = dialog.m_CoefRest;
		m_Kd = dialog.m_Damping;
		m_Gravity.x = dialog.m_GravX;
		m_Gravity.y = dialog.m_GravY;
		m_Gravity.z = dialog.m_GravZ;
		m_UserForceMag = dialog.m_UserForceMag;
		m_Ksh = dialog.m_SpringConst;
		m_Ksd = dialog.m_SpringDamp;
		for (int loop = 0; loop < m_SpringCnt; loop++)
		{
			m_Spring[loop].Ks = m_Ksh;
			m_Spring[loop].Kd = m_Ksd;
		}
	}
}

void CPhysEnv::SetVertexProperties()
{
	CSetVert	dialog;
	dialog.m_VertexMass = m_CurrentSys[m_Pick[0]].oneOverM;
	if (dialog.DoModal() == IDOK)
	{
		m_ParticleSys[0][m_Pick[0]].oneOverM = dialog.m_VertexMass;
		m_ParticleSys[0][m_Pick[1]].oneOverM = dialog.m_VertexMass;
		m_ParticleSys[1][m_Pick[0]].oneOverM = dialog.m_VertexMass;
		m_ParticleSys[1][m_Pick[1]].oneOverM = dialog.m_VertexMass;
		m_ParticleSys[2][m_Pick[0]].oneOverM = dialog.m_VertexMass;
		m_ParticleSys[2][m_Pick[1]].oneOverM = dialog.m_VertexMass;
	}
}

void CPhysEnv::ApplyUserForce(tVector* force)
{
	ScaleVector(force, m_UserForceMag, &m_UserForce);
	m_UserForceActive = TRUE;
}

///////////////////////////////////////////////////////////////////////////////
// Function:	SetMouseForce 
// Purpose:		Allows the user to interact with selected points by dragging
// Arguments:	Delta distance from clicked point, local x and y axes
///////////////////////////////////////////////////////////////////////////////
void CPhysEnv::SetMouseForce(int deltaX, int deltaY, tVector* localX, tVector* localY)
{
	/// Local Variables ///////////////////////////////////////////////////////////
	tVector tempX, tempY;
	///////////////////////////////////////////////////////////////////////////////
	ScaleVector(localX, (float)deltaX * 0.03f, &tempX);
	ScaleVector(localY, -(float)deltaY * 0.03f, &tempY);
	if (m_Pick[0] > -1)
	{
		VectorSum(&m_CurrentSys[m_Pick[0]].pos, &tempX, &m_MouseDragPos[0]);
		VectorSum(&m_MouseDragPos[0], &tempY, &m_MouseDragPos[0]);
	}
	if (m_Pick[1] > -1)
	{
		VectorSum(&m_CurrentSys[m_Pick[1]].pos, &tempX, &m_MouseDragPos[1]);
		VectorSum(&m_MouseDragPos[1], &tempY, &m_MouseDragPos[1]);
	}
}
/// SetMouseForce /////////////////////////////////////////////////////////////


void CPhysEnv::AddSpring()
{
	tSpring* spring;
	// MAKE SURE TWO PARTICLES ARE PICKED
	if (m_Pick[0] > -1 && m_Pick[1] > -1)
	{
		spring = (tSpring*)malloc(sizeof(tSpring) * (m_SpringCnt + 1));
		if (m_Spring != NULL)
			memcpy(spring, m_Spring, sizeof(tSpring) * m_SpringCnt);
		m_Spring = spring;
		spring = &m_Spring[m_SpringCnt++];
		spring->Ks = m_Ksh;
		spring->Kd = m_Ksd;
		spring->p1 = m_Pick[0];
		spring->p2 = m_Pick[1];
		spring->restLen =
			sqrt(VectorSquaredDistance(&m_CurrentSys[m_Pick[0]].pos,
				&m_CurrentSys[m_Pick[1]].pos));
		spring->type = MANUAL_SPRING;
	}
}

void CPhysEnv::AddSpring(int v1, int v2, float Ksh, float Ksd, int type)
{
	tSpring* spring;
	// MAKE SURE TWO PARTICLES ARE PICKED
	if (v1 > -1 && v2 > -1)
	{
		spring = (tSpring*)malloc(sizeof(tSpring) * (m_SpringCnt + 1));
		if (m_Spring != NULL)
		{
			memcpy(spring, m_Spring, sizeof(tSpring) * m_SpringCnt);
			free(m_Spring);
		}
		m_Spring = spring;
		spring = &m_Spring[m_SpringCnt++];
		spring->type = type;
		spring->Ks = Ksh;
		spring->Kd = Ksd;
		spring->p1 = v1;
		spring->p2 = v2;
		spring->restLen =
			sqrt(VectorSquaredDistance(&m_CurrentSys[v1].pos,
				&m_CurrentSys[v2].pos));
	}
}

void CPhysEnv::ComputeForces(tParticle* system)
{
	int loop;
	tParticle* curParticle, * p1, * p2;
	tSpring* spring;
	float		dist, Hterm, Dterm;
	tVector		springForce, deltaV, deltaP;

	curParticle = system;
	for (loop = 0; loop < m_ParticleCnt; loop++)
	{
		MAKEVECTOR(curParticle->f, 0.0f, 0.0f, 0.0f)		// CLEAR FORCE VECTOR

			if (m_UseGravity && curParticle->oneOverM != 0)
			{
				curParticle->f.x += (m_Gravity.x / curParticle->oneOverM);
				curParticle->f.y += (m_Gravity.y / curParticle->oneOverM);
				curParticle->f.z += (m_Gravity.z / curParticle->oneOverM);
			}

		if (m_UseDamping)
		{
			curParticle->f.x += (-m_Kd * curParticle->v.x);
			curParticle->f.y += (-m_Kd * curParticle->v.y);
			curParticle->f.z += (-m_Kd * curParticle->v.z);
		}
		else
		{
			curParticle->f.x += (-DEFAULT_DAMPING * curParticle->v.x);
			curParticle->f.y += (-DEFAULT_DAMPING * curParticle->v.y);
			curParticle->f.z += (-DEFAULT_DAMPING * curParticle->v.z);
		}
		curParticle++;
	}

	// CHECK IF THERE IS A USER FORCE BEING APPLIED
	if (m_UserForceActive)
	{
		if (m_Pick[0] != -1)
		{
			VectorSum(&system[m_Pick[0]].f, &m_UserForce, &system[m_Pick[0]].f);
		}
		if (m_Pick[1] != -1)
		{
			VectorSum(&system[m_Pick[1]].f, &m_UserForce, &system[m_Pick[1]].f);
		}
		MAKEVECTOR(m_UserForce, 0.0f, 0.0f, 0.0f);	// CLEAR USER FORCE
	}

	// NOW DO ALL THE SPRINGS
	spring = m_Spring;
	for (loop = 0; loop < m_SpringCnt; loop++)
	{
		p1 = &system[spring->p1];
		p2 = &system[spring->p2];
		VectorDifference(&p1->pos, &p2->pos, &deltaP);	// Vector distance 
		dist = VectorLength(&deltaP);					// Magnitude of deltaP

		Hterm = (dist - spring->restLen) * spring->Ks;	// Ks * (dist - rest)

		VectorDifference(&p1->v, &p2->v, &deltaV);		// Delta Velocity Vector
		Dterm = (DotProduct(&deltaV, &deltaP) * spring->Kd) / dist; // Damping Term

		ScaleVector(&deltaP, 1.0f / dist, &springForce);	// Normalize Distance Vector
		ScaleVector(&springForce, -(Hterm + Dterm), &springForce);	// Calc Force
		VectorSum(&p1->f, &springForce, &p1->f);			// Apply to Particle 1
		VectorDifference(&p2->f, &springForce, &p2->f);	// - Force on Particle 2
		spring++;					// DO THE NEXT SPRING
	}

	// APPLY THE MOUSE DRAG FORCES IF THEY ARE ACTIVE
	if (m_MouseForceActive)
	{
		// APPLY TO EACH PICKED PARTICLE
		if (m_Pick[0] > -1)
		{
			p1 = &system[m_Pick[0]];
			VectorDifference(&p1->pos, &m_MouseDragPos[0], &deltaP);	// Vector distance 
			dist = VectorLength(&deltaP);					// Magnitude of deltaP

			if (dist != 0.0f)
			{
				Hterm = (dist)*m_MouseForceKs;					// Ks * dist

				ScaleVector(&deltaP, 1.0f / dist, &springForce);	// Normalize Distance Vector
				ScaleVector(&springForce, -(Hterm), &springForce);	// Calc Force
				VectorSum(&p1->f, &springForce, &p1->f);			// Apply to Particle 1
			}
		}
		if (m_Pick[1] > -1)
		{
			p1 = &system[m_Pick[1]];
			VectorDifference(&p1->pos, &m_MouseDragPos[1], &deltaP);	// Vector distance 
			dist = VectorLength(&deltaP);					// Magnitude of deltaP

			if (dist != 0.0f)
			{
				Hterm = (dist)*m_MouseForceKs;					// Ks * dist

				ScaleVector(&deltaP, 1.0f / dist, &springForce);	// Normalize Distance Vector
				ScaleVector(&springForce, -(Hterm), &springForce);	// Calc Force
				VectorSum(&p1->f, &springForce, &p1->f);			// Apply to Particle 1
			}
		}
	}
}
/**
 * \brief Does the Integration for all the points in a system.
 * \param initial The initial positions, velocities and forces of the system.
 * \param source The source positions, velocities and forces of the system. This is where the slopes are calculated/taken from. Acceleration or 𝑑𝑣 / 𝑑𝑡 is Force / Mass. Force / Mass is the slope for calculating Velocity. Velocity is the slope for calculating Distance/Position
 * \param target The target positions, velocities and forces of the system. The result of the integration will be calculated and place here.
 * \param deltaTime The time step of the integration.
 */
void CPhysEnv::IntegrateSysOverTime(tParticle* initial, tParticle* source, tParticle* target, float deltaTime)
{
	///////////////////////////////////////////////////////////////////////////////
    for (int loop = 0; loop < m_ParticleCnt; loop++)
    {
		const float deltaTimeMass = deltaTime * initial->oneOverM;
        // DETERMINE THE NEW VELOCITY FOR THE PARTICLE
		// 𝑦ᵢ₊₁ = 𝑦ᵢ + 𝑓( 𝑥ᵢ , 𝑦ᵢ ) * 𝒉
		// 𝑓( 𝑥ᵢ , 𝑦ᵢ ) = 𝑑𝑣 / 𝑑𝑡  = 𝑎 ( 𝑡 ) = 𝐹 / 𝑚
		// ∂𝑣 / ∂𝑥 * ∂𝑥 / ∂𝑡
        target->v.x = initial->v.x + (source->f.x * deltaTimeMass);
		// ∂𝑣 / ∂𝑦 * ∂𝑦 / ∂𝑡
        target->v.y = initial->v.y + (source->f.y * deltaTimeMass);
		// ∂𝑣 / ∂𝑧 * ∂𝑧 / ∂𝑡
        target->v.z = initial->v.z + (source->f.z * deltaTimeMass);

		// The mass doesn't change
        target->oneOverM = initial->oneOverM;

        // SET THE NEW POSITION
		// This is a Time vs Velocity graph. If we integrate it we get distance which is used to calculate the new position.
		// 𝑦ᵢ₊₁ = 𝑦ᵢ + 𝑓( 𝑥ᵢ , 𝑦ᵢ ) * 𝒉
		// 𝑓( 𝑥ᵢ , 𝑦ᵢ ) = 𝑑𝑥 / 𝑑𝑡  = 𝑣 ( 𝑡 )
		// ∂𝑥 / ∂𝑡
		target->pos.x = initial->pos.x + (deltaTime * source->v.x);
		// ∂𝑦 / ∂𝑡
        target->pos.y = initial->pos.y + (deltaTime * source->v.y);
		// ∂𝑧 / ∂𝑡
        target->pos.z = initial->pos.z + (deltaTime * source->v.z);

        initial++;
        source++;
        target++;
    }
}
/**
 * \brief Uses the Forward Euler method to integrate the system.
 * \param DeltaTime The amount of time to integrate the system over.
 */
void CPhysEnv::EulerIntegrate(float DeltaTime)
{
    // Use the array of particles "m_CurrentSys" to fill the system of particles "cur"
    //tParticle* yn = m_CurrentSys;
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ )
    //tParticle* k1 = m_CurrentSys;
    //tParticle* ynp1 = m_TargetSys;
    // Read through the implementation of the function (implemented in PhysEnv.cpp) and make sure you understand how it works
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + 𝑘₁ * 𝒉
    //IntegrateSysOverTime ( yn , k1 , ynp1 , DeltaTime );
    IntegrateSysOverTime ( m_CurrentSys , m_CurrentSys , m_TargetSys , DeltaTime );
}
/**
 * \brief Uses the Explicit Midpoint method to integrate the system.
 * \param DeltaTime The amount of time to integrate the system over.
 */
void CPhysEnv::MidPointIntegrate ( float DeltaTime )
{
    const float halfDeltaT = DeltaTime / 2.0f;
    //tParticle * yn = m_CurrentSys;
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ )
    //tParticle * k1 = m_CurrentSys;
    //System k2 ( m_ParticleCnt );
    // Compute the state of the system at the half of the interval
    // 𝑘₂ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₁ * 𝒉 )
    //IntegrateSysOverTime ( yn , k1 , k2 , halfDeltaT );
    IntegrateSysOverTime ( m_CurrentSys , m_CurrentSys , m_TempSys [ 0 ] , halfDeltaT );
    // Evaluate derivatives at the half of the interval
    // The function ComputeForces will update the forces on each particle in the System "k2"
    //ComputeForces ( k2 );
    ComputeForces ( m_TempSys [ 0 ] );
    //tParticle * ynp1 = m_TargetSys;
    // Use these derivatives to compute the state at the end of the interval
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + 𝑘₂ * 𝒉
    //IntegrateSysOverTime ( yn , k2 , ynp1 , DeltaTime );
    IntegrateSysOverTime ( m_CurrentSys , m_TempSys [ 0 ] , m_TargetSys , DeltaTime );
}
float CPhysEnv::CalculateTwoSystemError ( tParticle * systemOne , tParticle * systemTwo , int particleCount ) const
{
    float error = 0.0f;
    for ( std::size_t i = 0; i < particleCount; ++i )
    {
        error += std::abs ( ( systemOne [ i ].pos.x - systemTwo [ i ].pos.x ) / systemOne [ i ].pos.x ) * 100;
        error += std::abs ( ( systemOne [ i ].pos.y - systemTwo [ i ].pos.y ) / systemOne [ i ].pos.y ) * 100;
        error += std::abs ( ( systemOne [ i ].pos.z - systemTwo [ i ].pos.z ) / systemOne [ i ].pos.z ) * 100;
    }
    error /= ( 3 * particleCount );
    return error;
}
/**
 * \brief Uses Heun's method to integrate the system.
 * \param DeltaTime The amount of time to integrate the system over.
 */
//void CPhysEnv::HeunIntegrate ( float DeltaTime )
//{
//    //The velocities and positions at the intial t
//    System y_0 ( m_CurrentSys , m_ParticleCnt ); //2,Velocity,Position
//    //The slope at the intial t
//    System y_prime_0 ( m_CurrentSys , m_ParticleCnt );//3,Slope
//    //The velocities and positions at the end t
//    System y_0_1 ( m_ParticleCnt );
//    IntegrateSysOverTime ( y_0 , y_prime_0 , y_0_1 , DeltaTime );//5,Velocity,Position
//    while ( true )
//    {
//        //The slope at the end t
//        System y_prime_1 ( y_0_1 );
//        ComputeForces ( y_prime_1 ); //6.402164,Slope
//        //Average slope at start and end
//        System y_bar_prime = 1.0f / 2.0f * ( y_prime_0 + y_prime_1 ); //4.701082,Slope
//        //The velocities and positions at the end t after correction
//        System y_i_1 ( m_ParticleCnt );
//        IntegrateSysOverTime ( y_0 , y_bar_prime , y_i_1 , DeltaTime ); //6.701082,Velocity,Position
//        const float error = CalculateTwoSystemError ( y_i_1 , y_0_1 , m_ParticleCnt );
//        if ( error > 0.01 )
//        {
//            y_0_1 = y_i_1;
//            continue;
//        }
//        break;
//    }
//    y_0_1.fillOut ( m_TargetSys );
//}
void CPhysEnv::HeunIntegrate ( float DeltaTime )
{
    const float stoppingCriterion = 0.00001f; //Percentage 
    //The velocities and positions at the intial t
    System y_0 ( m_CurrentSys , m_ParticleCnt ); //Position
    //The slope at the intial t
    System y_prime_0 ( m_CurrentSys , m_ParticleCnt );//Slope
    //The velocities and positions at the end t
    System y_0_1 ( m_ParticleCnt );
    IntegrateSysOverTime ( y_0 , y_prime_0 , y_0_1 , DeltaTime );//Position
    while ( true )
    {
        //The slope at the end t
        System y_prime_1 ( y_0_1 );
        ComputeForces ( y_prime_1 ); //Slope
        //Average slope at start and end
        System y_bar_prime = 1.0f / 2.0f * ( y_prime_0 + y_prime_1 ); //Slope
        //The velocities and positions at the end t after correction
        System y_i_1 ( m_ParticleCnt );
        IntegrateSysOverTime ( y_0 , y_bar_prime , y_i_1 , DeltaTime ); //Position
        const float error = CalculateTwoSystemError ( y_i_1 , y_0_1 , m_ParticleCnt );
        if ( error > stoppingCriterion )
        {
            y_0_1 = y_i_1;
            continue;
        }
        break;
    }
    y_0_1.fillOut ( m_TargetSys );
}
void CPhysEnv::Swap ( tParticle * source , tParticle * target ) const
{
    for ( int loop = 0; loop < m_ParticleCnt; loop++ )
    {
        target->v.x = source->v.x;
        target->v.y = source->v.y;
        target->v.z = source->v.z;
        target->oneOverM = source->oneOverM;
        target->pos.x = source->pos.x;
        target->pos.y = source->pos.y;
        target->pos.z = source->pos.z;
        source++;
        target++;
    }
}
/**
 * \brief Uses Fourth-Order Runge–Kutta (4th order) method to integrate the system.
 * \param DeltaTime The amount of time to integrate the system over.
 */
void CPhysEnv::RK4Integrate ( float DeltaTime )
{
    const float halfDeltaT = DeltaTime / 2.0f;
    //System yn ( m_CurrentSys , m_ParticleCnt );
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ ) I know I'm creating an extra copy that I don't really need but this is done for clarity not efficiency.
    System k1 ( m_CurrentSys , m_ParticleCnt );
    System k2 ( m_ParticleCnt );
    // 𝑘₂ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₁ * 𝒉 )
    //IntegrateSysOverTime ( yn , k1 , k2 , halfDeltaT );
    IntegrateSysOverTime ( m_CurrentSys , k1 , k2 , halfDeltaT );
    ComputeForces ( k2 );
    System k3 ( m_ParticleCnt );
    // 𝑘₃ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₂ * 𝒉 )
    //IntegrateSysOverTime ( yn , k2 , k3 , halfDeltaT );
    IntegrateSysOverTime ( m_CurrentSys , k2 , k3 , halfDeltaT );
    ComputeForces ( k3 );
    System k4 ( m_ParticleCnt );
    // 𝑘₄ = 𝑓( 𝑥ᵢ + 𝒉 , 𝑦ᵢ + 𝑘₃ * 𝒉 )
    //IntegrateSysOverTime ( yn , k3 , k4 , DeltaTime );
    IntegrateSysOverTime ( m_CurrentSys , k3 , k4 , DeltaTime );
    ComputeForces ( k4 );
    //System ynp1 ( m_ParticleCnt );
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + ( 1 / 6 ) * ( 𝑘₁ + 2𝑘₂ + 2𝑘₃ + 𝑘₄ ) * 𝒉
    //IntegrateSysOverTime ( yn , 1.0f / 6.0f * ( k1 + 2 * k2 + 2 * k3 + k4 ) , ynp1 , DeltaTime );
    //IntegrateSysOverTime ( m_CurrentSys , 1.0f / 6.0f * ( k1 + 2 * k2 + 2 * k3 + k4 ) , m_TargetSys , DeltaTime );
    IntegrateSysOverTime ( m_CurrentSys , 1.0f / 6.0f * ( k1 + 2 * ( k2 + k3 ) + k4 ) , m_TargetSys , DeltaTime );
    //ynp1.fillOut ( m_TargetSys );
}
void CPhysEnv::RK5Integrate ( float DeltaTime )
{
    const float quarterDelta = DeltaTime / 4;
    const float halfDelta = DeltaTime / 2;
    const float threeQuartersDelta = 3 * ( DeltaTime / 4 );
    //System yn ( m_CurrentSys , m_ParticleCnt );
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ ) I know I'm creating an extra copy that I don't really need but this is done for clarity not efficiency.
    System k1 ( m_CurrentSys , m_ParticleCnt );
    System k2 ( m_ParticleCnt );
    // 𝑘₂ = 𝑓( 𝑥ᵢ + ( 1 / 4 ) * 𝒉 , 𝑦ᵢ + ( 1 / 4 ) * 𝑘₁ * 𝒉 )
    //IntegrateSysOverTime ( yn , k1 , k2 , quarterDelta );
    IntegrateSysOverTime ( m_CurrentSys , k1 , k2 , quarterDelta );
    ComputeForces ( k2 );
    System k3 ( m_ParticleCnt );
    // 𝑘₃ = 𝑓( 𝑥ᵢ + ( 1 / 4 ) * 𝒉 , 𝑦ᵢ + ( 1 / 8 ) * ( 𝑘₁ +  𝑘₂) * 𝒉  )
    //    = 𝑓( 𝑥ᵢ + ( 1 / 4 ) * 𝒉 , 𝑦ᵢ + ( 1 / 4 ) * ( ( 1 / 2 ) * ( 𝑘₁ + 𝑘₂ ) ) * 𝒉 )
    //IntegrateSysOverTime ( yn , ( 1.0f / 2.0f ) * ( k1 + k2 ) , k3 , quarterDelta );
    IntegrateSysOverTime ( m_CurrentSys , ( 1.0f / 2.0f ) * ( k1 + k2 ) , k3 , quarterDelta );
    ComputeForces ( k3 );
    System k4 ( m_ParticleCnt );
    // 𝑘₄ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( ( -1 / 2 ) * 𝑘₂ + 𝑘₃ ) * 𝒉 )
    //    = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * ( −𝑘₂ + 2 * 𝑘₃ ) * 𝒉 )
    //IntegrateSysOverTime ( yn , -1 * k2 + 2 * k3 , k4 , halfDelta );
    IntegrateSysOverTime ( m_CurrentSys , -1 * k2 + 2 * k3 , k4 , halfDelta );
    ComputeForces ( k4 );
    System k5 ( m_ParticleCnt );
    // 𝑘₅ = 𝑓( 𝑥ᵢ + ( 3 / 4 ) * 𝒉 , 𝑦ᵢ + ( ( 3 / 16 ) * 𝑘₁ + ( 9 / 16 ) * 𝑘₄ ) * 𝒉 )
    //    = 𝑓( 𝑥ᵢ + ( 3 / 4 ) * 𝒉 , 𝑦ᵢ + ( 3 / 4 ) * ( ( 1 / 4 ) * 𝑘₁ + ( 3 / 4 ) * 𝑘₄ ) * 𝒉 )
    //IntegrateSysOverTime ( yn , ( 1.0f / 4.0f ) * k1 + ( 3.0f / 4.0f ) * k4 , k5 , threeQuartersDelta );
    //IntegrateSysOverTime ( m_CurrentSys , ( 1.0f / 4.0f ) * k1 + ( 3.0f / 4.0f ) * k4 , k5 , threeQuartersDelta );
    IntegrateSysOverTime ( m_CurrentSys , ( k1 + 3 * k4 ) / 4.0f , k5 , threeQuartersDelta );
    ComputeForces ( k5 );
    System k6 ( m_ParticleCnt );
    // 𝑘₆ = 𝑓( 𝑥ᵢ + 𝒉 , 𝑦ᵢ + ( ( -3 / 7 ) * 𝑘₁ + ( 2 / 7 ) * 𝑘₂ + ( 12 / 7 ) * 𝑘₃ + ( -12 / 7 ) * 𝑘₄ + ( 8 / 7 ) * 𝑘₅ ) * 𝒉 )
    //IntegrateSysOverTime ( yn , ( -3.0f / 7.0f ) * k1 + ( 2.0f / 7.0f ) * k2 + ( 12.0f / 7.0f ) * k3 + ( -12.0f / 7.0f ) * k4 + ( 8.0f / 7.0f ) * k5 , k6 , DeltaTime );
    //IntegrateSysOverTime ( m_CurrentSys , ( -3.0f / 7.0f ) * k1 + ( 2.0f / 7.0f ) * k2 + ( 12.0f / 7.0f ) * k3 + ( -12.0f / 7.0f ) * k4 + ( 8.0f / 7.0f ) * k5 , k6 , DeltaTime );
    IntegrateSysOverTime ( m_CurrentSys , ( -3.0f * k1 + 2.0f * k2 + 12.0f * k3 - 12.0f * k4 + 8.0f * k5 ) / 7.0f , k6 , DeltaTime );
    ComputeForces ( k6 );
    //System ynp1 ( m_ParticleCnt );
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + ( 1 / 90 ) * ( 7𝑘₁ + 32𝑘₃ + 12𝑘₄ + 32𝑘₅ + 7𝑘₆ ) * 𝒉
    //IntegrateSysOverTime ( yn , ( 1.0f / 90.0f ) * ( 7 * k1 + 32 * k3 + 12 * k4 + 32 * k5 + 7 * k6 ) , ynp1 , DeltaTime );
    //IntegrateSysOverTime ( m_CurrentSys , ( 1.0f / 90.0f ) * ( 7 * k1 + 32 * k3 + 12 * k4 + 32 * k5 + 7 * k6 ) , m_TargetSys , DeltaTime );
    IntegrateSysOverTime ( m_CurrentSys , ( 7 * k1 + 32 * k3 + 12 * k4 + 32 * k5 + 7 * k6 ) / 90.0f , m_TargetSys , DeltaTime );
    //ynp1.fillOut ( m_TargetSys );
}
void CPhysEnv::RK4AdaptiveIntegrate ( float DeltaTime )
{
    const float h1 = DeltaTime;
    const float halfH1 = h1 / 2.0f;
    const float h2 = halfH1;
    const float halfH2 = h2 / 2.0f;


    //The initial state of the system.
    //System yn ( m_CurrentSys , m_ParticleCnt );
    // THE SINGLE STEP
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ ) I know I'm creating an extra copy that I don't really need but this is done for clarity not efficiency.
    System k1_1 ( m_CurrentSys , m_ParticleCnt );
    System k2_1 ( m_ParticleCnt );
    // 𝑘₂ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₁ * 𝒉 )
    //IntegrateSysOverTime ( yn , k1_1 , k2_1 , halfH1 );
    IntegrateSysOverTime ( m_CurrentSys , k1_1 , k2_1 , halfH1 );
    ComputeForces ( k2_1 );
    System k3_1 ( m_ParticleCnt );
    // 𝑘₃ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₂ * 𝒉 )
    IntegrateSysOverTime ( m_CurrentSys , k2_1 , k3_1 , halfH1 );
    ComputeForces ( k3_1 );
    System k4_1 ( m_ParticleCnt );
    // 𝑘₄ = 𝑓( 𝑥ᵢ + 𝒉 , 𝑦ᵢ + 𝑘₃ * 𝒉 )
    IntegrateSysOverTime ( m_CurrentSys , k3_1 , k4_1 , h1 );
    ComputeForces ( k4_1 );
    System y1 ( m_ParticleCnt );
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + ( 1 / 6 ) * ( 𝑘₁ + 2𝑘₂ + 2𝑘₃ + 𝑘₄ ) * 𝒉
    IntegrateSysOverTime ( m_CurrentSys , 1.0f / 6.0f * ( k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1 ) , y1 , h1 );


    // THE FIRST HALF STEP
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ ) Extra copy. Clarity > Efficiency
    System k1Half_2 ( m_CurrentSys , m_ParticleCnt );
    // 𝑘₂ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₁ * 𝒉 )
    System k2Half_2 ( m_ParticleCnt );
    //IntegrateSysOverTime ( yn , k1Half_2 , k2Half_2 , halfH2 );
    IntegrateSysOverTime ( m_CurrentSys , k1Half_2 , k2Half_2 , halfH2 );
    ComputeForces ( k2Half_2 );
    System k3Half_2 ( m_ParticleCnt );
    // 𝑘₃ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₂ * 𝒉 )
    //IntegrateSysOverTime ( yn , k2Half_2 , k3Half_2 , halfH2 );
    IntegrateSysOverTime ( m_CurrentSys , k2Half_2 , k3Half_2 , halfH2 );
    ComputeForces ( k3Half_2 );
    System k4Half_2 ( m_ParticleCnt );
    // 𝑘₄ = 𝑓( 𝑥ᵢ + 𝒉 , 𝑦ᵢ + 𝑘₃ * 𝒉 )
    //IntegrateSysOverTime ( yn , k3Half_2 , k4Half_2 , h2 );
    IntegrateSysOverTime ( m_CurrentSys , k3Half_2 , k4Half_2 , h2 );
    ComputeForces ( k4Half_2 );
    System y2Half ( m_ParticleCnt );
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + ( 1 / 6 ) * ( 𝑘₁ + 2𝑘₂ + 2𝑘₃ + 𝑘₄ ) * 𝒉
    //IntegrateSysOverTime ( yn , 1.0f / 6.0f * ( k1Half_2 + 2 * k2Half_2 + 2 * k3Half_2 + k4Half_2 ) , y2Half , h2 );
    IntegrateSysOverTime ( m_CurrentSys , 1.0f / 6.0f * ( k1Half_2 + 2 * k2Half_2 + 2 * k3Half_2 + k4Half_2 ) , y2Half , h2 );
    ComputeForces ( y2Half );


    // THE SECOND HALF STEP
    // 𝑘₁ = 𝑓( 𝑥ᵢ , 𝑦ᵢ ) Extra copy. Clarity > Efficiency
    //System k1_2 = y2Half;
    System k2_2 ( m_ParticleCnt );
    // 𝑘₂ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₁ * 𝒉 )
    //IntegrateSysOverTime ( y2Half , k1_2 , k2_2 , halfH2 );
    IntegrateSysOverTime ( y2Half , y2Half , k2_2 , halfH2 );
    ComputeForces ( k2_2 );
    System k3_2 ( m_ParticleCnt );
    // 𝑘₃ = 𝑓( 𝑥ᵢ + ( 1 / 2 ) * 𝒉 , 𝑦ᵢ + ( 1 / 2 ) * 𝑘₂ * 𝒉 )
    IntegrateSysOverTime ( y2Half , k2_2 , k3_2 , halfH2 );
    ComputeForces ( k2_2 );
    System k4_2 ( m_ParticleCnt );
    // 𝑘₄ = 𝑓( 𝑥ᵢ + 𝒉 , 𝑦ᵢ + 𝑘₃ * 𝒉 )
    IntegrateSysOverTime ( y2Half , k3_2 , k4_2 , h2 );
    ComputeForces ( k4_2 );
    System y2 ( m_ParticleCnt );
    // 𝑦ᵢ₊₁ = 𝑦ᵢ + ( 1 / 6 ) * ( 𝑘₁ + 2𝑘₂ + 2𝑘₃ + 𝑘₄ ) * 𝒉
    //IntegrateSysOverTime ( y2Half , 1.0f / 6.0f * ( k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2 ) , y2 , h2 );
    IntegrateSysOverTime ( y2Half , 1.0f / 6.0f * ( y2Half + 2 * k2_2 + 2 * k3_2 + k4_2 ) , y2 , h2 );
    System errorDelta = y2 - y1;
    y2 = y2 + errorDelta / 15;
    //tParticle * ynp1 = m_TargetSys;
    //y2.fillOut ( ynp1 );
    y2.fillOut ( m_TargetSys );
}

///////////////////////////////////////////////////////////////////////////////

int CPhysEnv::CheckForCollisions(tParticle* system)
{
	// be optimistic!
	int collisionState = NOT_COLLIDING;
	float const depthEpsilon = 0.001f;

	int loop;
	tParticle* curParticle;

	m_ContactCnt = 0;		// THERE ARE CURRENTLY NO CONTACTS

	curParticle = system;
	for (loop = 0; (loop < m_ParticleCnt) && (collisionState != PENETRATING);
		loop++, curParticle++)
	{
		// CHECK THE MAIN BOUNDARY PLANES FIRST
		for (int planeIndex = 0;(planeIndex < m_CollisionPlaneCnt) &&
			(collisionState != PENETRATING);planeIndex++)
		{
			tCollisionPlane* plane = &m_CollisionPlane[planeIndex];

			float axbyczd = DotProduct(&curParticle->pos, &plane->normal) + plane->d;

			if (axbyczd < -depthEpsilon)
			{
				// ONCE ANY PARTICLE PENETRATES, QUIT THE LOOP
				collisionState = PENETRATING;
			}
			else
				if (axbyczd < depthEpsilon)
				{
					float relativeVelocity = DotProduct(&plane->normal, &curParticle->v);

					if (relativeVelocity < 0.0f)
					{
						collisionState = COLLIDING;
						m_Contact[m_ContactCnt].particle = loop;
						memcpy(&m_Contact[m_ContactCnt].normal, &plane->normal, sizeof(tVector));
						m_ContactCnt++;
					}
				}
		}
		if (m_CollisionActive)
		{
			// THIS IS NEW FROM MARCH - ADDED SPHERE BOUNDARIES
			// CHECK ANY SPHERE BOUNDARIES
			for (int sphereIndex = 0;(sphereIndex < m_SphereCnt) &&
				(collisionState != PENETRATING);sphereIndex++)
			{
				tCollisionSphere* sphere = &m_Sphere[sphereIndex];

				tVector	distVect;

				VectorDifference(&curParticle->pos, &sphere->pos, &distVect);

				float radius = VectorSquaredLength(&distVect);
				// SINCE IT IS TESTING THE SQUARED DISTANCE, SQUARE THE RADIUS ALSO
				float dist = radius - (sphere->radius * sphere->radius);

				if (dist < -depthEpsilon)
				{
					// ONCE ANY PARTICLE PENETRATES, QUIT THE LOOP
					collisionState = PENETRATING;
				}
				else
					if (dist < depthEpsilon)
					{
						// NORMALIZE THE VECTOR
						NormalizeVector(&distVect);

						float relativeVelocity = DotProduct(&distVect, &curParticle->v);

						if (relativeVelocity < 0.0f)
						{
							collisionState = COLLIDING;
							m_Contact[m_ContactCnt].particle = loop;
							memcpy(&m_Contact[m_ContactCnt].normal, &distVect, sizeof(tVector));
							m_ContactCnt++;
						}
					}
			}
		}

	}

	return collisionState;
}

void CPhysEnv::ResolveCollisions(tParticle* system)
{
	tContact* contact;
	tParticle* particle;		// THE PARTICLE COLLIDING
	float		VdotN;
	tVector		Vn, Vt;				// CONTACT RESOLUTION IMPULSE
	contact = m_Contact;
	for (int loop = 0; loop < m_ContactCnt; loop++, contact++)
	{
		particle = &system[contact->particle];
		// CALCULATE Vn
		VdotN = DotProduct(&contact->normal, &particle->v);
		ScaleVector(&contact->normal, VdotN, &Vn);
		// CALCULATE Vt
		VectorDifference(&particle->v, &Vn, &Vt);
		// SCALE Vn BY COEFFICIENT OF RESTITUTION
		ScaleVector(&Vn, m_Kr, &Vn);
		// SET THE VELOCITY TO BE THE NEW IMPULSE
		VectorDifference(&Vt, &Vn, &particle->v);
	}
}

void CPhysEnv::Simulate(float DeltaTime, BOOL running)
{
	float		CurrentTime = 0.0f;
	float		TargetTime = DeltaTime;
	tParticle* tempSys;
	int			collisionState;

	while (CurrentTime < DeltaTime)
	{
        if (running)
		{
			ComputeForces(m_CurrentSys);
			// IN ORDER TO MAKE THINGS RUN FASTER, I HAVE THIS LITTLE TRICK
			// IF THE SYSTEM IS DOING A BINARY SEARCH FOR THE COLLISION POINT,
			// I FORCE EULER'S METHOD ON IT. OTHERWISE, LET THE USER CHOOSE.
			// THIS DOESN'T SEEM TO EFFECT STABILITY EITHER WAY
			if ( m_CollisionRootFinding )
			{
				EulerIntegrate( TargetTime - CurrentTime );
            }
			else
			{
				switch (m_IntegratorType)
				{
				case EULER_INTEGRATOR:
					EulerIntegrate( TargetTime - CurrentTime );
					break;
				case MIDPOINT_INTEGRATOR:
					MidPointIntegrate( TargetTime - CurrentTime );
					break;
				case HEUN_INTEGRATOR:
					HeunIntegrate( TargetTime - CurrentTime );
					break;
				case RK4_INTEGRATOR:
					RK4Integrate( TargetTime - CurrentTime );
					break;
				 case RK5_INTEGRATOR:
					RK5Integrate( TargetTime - CurrentTime );
					 break;
				case RK4_ADAPTIVE_INTEGRATOR:
					RK4AdaptiveIntegrate( TargetTime - CurrentTime );
					break;
				}
			}
		}
		collisionState = CheckForCollisions(m_TargetSys);

		if (collisionState == PENETRATING)
		{
			// TELL THE SYSTEM I AM LOOKING FOR A COLLISION SO IT WILL USE EULER
			m_CollisionRootFinding = TRUE;
			// we simulated too far, so subdivide time and try again
			TargetTime = (CurrentTime + TargetTime) / 2.0f;

			// blow up if we aren't moving forward each step, which is
			// probably caused by interpenetration at the frame start
			assert(fabs(TargetTime - CurrentTime) > EPSILON);
		}
		else
		{
			// either colliding or clear
			if (collisionState == COLLIDING)
			{
				int Counter = 0;
				do
				{
					ResolveCollisions(m_TargetSys);
					Counter++;
				} while ((CheckForCollisions(m_TargetSys) ==
					COLLIDING) && (Counter < 100));

				assert(Counter < 100);
				m_CollisionRootFinding = FALSE;  // FOUND THE COLLISION POINT
			}

			// we made a successful step, so swap configurations
			// to "save" the data for the next step
            if ( OUTPUT_TO_FILE )
            {
                auto errors = CalculateError ();

            }
			CurrentTime = TargetTime;
			TargetTime = DeltaTime;

			// SWAP MY TWO PARTICLE SYSTEM BUFFERS SO I CAN DO IT AGAIN
			tempSys = m_CurrentSys;
			m_CurrentSys = m_TargetSys;
			m_TargetSys = tempSys;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Function:	AddCollisionSphere 
// Purpose:		Add a collision sphere to the system
///////////////////////////////////////////////////////////////////////////////
void CPhysEnv::AddCollisionSphere()
{
	/// Local Variables ///////////////////////////////////////////////////////////
	tCollisionSphere* temparray;
	CAddSpher		dialog;
	///////////////////////////////////////////////////////////////////////////////
	dialog.m_Radius = 2.0f;
	dialog.m_XPos = 0.0f;
	dialog.m_YPos = -3.0f;
	dialog.m_ZPos = 0.0f;
	if (dialog.DoModal())
	{
		temparray = (tCollisionSphere*)malloc(sizeof(tCollisionSphere) * (m_SphereCnt + 1));
		if (m_SphereCnt > 0)
		{
			memcpy(temparray, m_Sphere, sizeof(tCollisionSphere) * m_SphereCnt);
			free(m_Sphere);
		}

		MAKEVECTOR(temparray[m_SphereCnt].pos, dialog.m_XPos, dialog.m_YPos, dialog.m_ZPos)
			temparray[m_SphereCnt].radius = dialog.m_Radius;

		m_Sphere = temparray;
		m_SphereCnt++;
	}
}