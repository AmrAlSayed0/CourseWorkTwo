#pragma once
#include <utility>
#include <vector>
#include "MathDefs.h"
/* A particle is defined inside PhysEnv.h as
struct tParticle
{
	tVector pos;		// Position of Particle
    tVector v;			// Velocity of Particle
	tVector f;			// Force Acting on Particle
	float	oneOverM;	// 1 / Mass of Particle
};
*/
#include "PhysEnv.h"
/**
 * System of particles (masses that are connected together using springs) In this class, we encapsulate a vector of particles.
 *
 * When implementing different integrators, we often need to perform sum/average/... of the derivatives of the system at different time steps.
 * In our case, the derivatives are the velocity (first derivative) and acceleration/force (second derivative) of each particle in the system.
 *
 * Below are some helper functions that can make your life easier when doing such arithmetic operations
 */
class System
{
    std::vector < tParticle > particles_;
public:
    /**
     * \brief Create a system of n particles. The system is initially empty and acts as a place holder for later operations.
     * \param n The number of particles the system should hold.
     */
    explicit System ( const int n ): particles_ ( n )
    {}
    /**
     * \brief Creates a system of particles and fills it using the input array.
     * \param sys The array of particles to fill the system with.
     * \param n The number of particles in the array.
     */
    System ( tParticle * sys , const int n ): System ( n )
    {
        std::memcpy ( static_cast < void * > ( this->particles_.data () ) , static_cast < const void * > ( sys ) , static_cast < std::size_t > ( n  * sizeof ( tParticle ) ) );
    }
    /**
     * \brief Creates a system of particles from a vector of particles. The vector is copied not taken.
     * \param particles The particles the system should contain
     */
    explicit System ( std::vector < tParticle > & particles ) : particles_ ( particles )
    {}
    /**
     * \brief Copy constructor.
     * \param other The System to copy.
     */
    System ( const System & other ) = default;
    /**
     * \brief Move constructor.
     * \param other The System to move.
     */
    System ( System && other ) noexcept : particles_ ( std::exchange ( static_cast < std::vector < tParticle > & > ( other.particles_ ) , static_cast < std::vector < tParticle > && > ( std::vector < tParticle > ( 0 ) ) ) )
    {}
    ~System () = default;
    /**
     * \brief Copy assignment operator.
     * \param other The System to copy.
     * \return A reference to the left (this) system after the right (other) system is copied into it.
     */
    System & operator= ( const System & other )
    {
        if ( this == &other )
        {
            return static_cast < System & > ( *this );
        }
        particles_ = other.particles_;
        return static_cast < System & > ( *this );
    }
    /**
     * \brief Move assignment operator.
     * \param other The System to move.
     * \return A reference to the left (this) system after the right (other) system is moved into it.
     */
    System & operator= ( System && other ) noexcept
    {
        if ( this == &other )
        {
            return static_cast < System & > ( *this );
        }
        particles_ = std::move ( other.particles_ );
        return static_cast < System & > ( *this );
    }
    /**
     * \brief Overload of the addition operator. Adds the derivatives of the left and right systems and returns a new system.
     * \param other The other (right) system to add to the left one.
     * \return A new system that is the result of the addition.
     */
    System operator + ( const System & other ) const
    {
        System temp ( *this );
        return temp += other;
    }
    /**
     * \brief Overload of the addition assignment operator. Adds the derivatives of the right (other) system to the left (this) system.
     * \param other The other system to add to this system.
     * \return A reference to this system after the addition operation.
     */
    System & operator+= ( const System & other )
    {
        const std::size_t this_n = this->particles_.size ();
        if ( this_n > other.particles_.size () )
        {
            throw std::invalid_argument ( "The 2 systems must have the same number of particles" );
        }
        for ( std::size_t i = 0; i < this_n; ++i )
        {
            VectorSum ( &this->particles_ [ i ].pos , &other.particles_ [ i ].pos , &this->particles_ [ i ].pos );
            VectorSum ( &this->particles_ [ i ].v , &other.particles_ [ i ].v , &this->particles_ [ i ].v );
            VectorSum ( &this->particles_ [ i ].f , &other.particles_ [ i ].f , &this->particles_ [ i ].f );
            //this->particles_ [ i ].oneOverM += other.particles_ [ i ].oneOverM;
        }
        return static_cast < System & > ( *this );
    }
    /**
     * \brief Overload of the subtraction operator. Subtracts the derivatives of the right from the left systems and returns a new system.
     * \param other The other (right) system to subtract from the left one.
     * \return A new system that is the result of the subtraction.
     */
    System operator - ( const System & other ) const
    {
        System temp ( *this );
        return temp -= other;
    }
    /**
     * \brief Overload of the subtraction assignment operator. Subtracts the derivatives of the right (other) from the left (this) system.
     * \param other The other system to subtract from this system.
     * \return A reference to this system after the subtraction operation.
     */
    System & operator -= ( const System & other )
    {
        const std::size_t this_n = this->particles_.size ();
        if ( this_n > other.particles_.size () )
        {
            throw std::invalid_argument ( "The 2 systems must have the same number of particles" );
        }
        for ( std::size_t i = 0; i < this_n; ++i )
        {
            VectorDifference ( &this->particles_ [ i ].pos , &other.particles_ [ i ].pos , &this->particles_ [ i ].pos );
            VectorDifference ( &this->particles_ [ i ].v , &other.particles_ [ i ].v , &this->particles_ [ i ].v );
            VectorDifference ( &this->particles_ [ i ].f , &other.particles_ [ i ].f , &this->particles_ [ i ].f );
            //this->particles_ [ i ].oneOverM += other.particles_ [ i ].oneOverM;
        }
        return static_cast < System & > ( *this );
    }
    /**
     * \brief Overload of the multiplication operator. Scales all the vectors of the left system by a given factor.
     * \param k The factor to scale by.
     * \return A new system that is the result of the scaling.
     */
    System operator * ( const float k ) const
    {
        System temp ( *this );
        return temp *= k;
    }
    /**
     * \brief Overload of the multiplication assignment operator. Scales all the vectors of the left (this) system by a given factor.
     * \param k The factor to scale by.
     * \return A reference to this system after the scaling operation.
     */
    System & operator *= ( const float k )
    {
        const std::size_t this_n = this->particles_.size ();
        for ( std::size_t i = 0; i < this_n; ++i )
        {
            ScaleVector ( &this->particles_ [ i ].pos , k , &this->particles_ [ i ].pos );
            ScaleVector ( &this->particles_ [ i ].v , k , &this->particles_ [ i ].v );
            ScaleVector ( &this->particles_ [ i ].f , k , &this->particles_ [ i ].f );
        }
        return static_cast < System & > ( *this );
    }
    /**
     * \brief Overload of the division operator. De-scales all the vectors of the left system by a given factor.
     * \param k The factor to de-scales by.
     * \return A new system that is the result of the de-scaling.
     */
    System operator / ( const float k ) const
    {
        System temp ( *this );
        return temp /= k;
    }
    /**
     * \brief Overload of the multiplication assignment operator. De-scales all the vectors of the left (this) system by a given factor.
     * \param k The factor to de-scales by.
     * \return A reference to this system after the de-scaling operation.
     */
    System & operator /= ( const float k )
    {
        const std::size_t this_n = this->particles_.size ();
        const float scaleFactor = 1 / k;
        for ( std::size_t i = 0; i < this_n; ++i )
        {
            ScaleVector ( &this->particles_ [ i ].pos , scaleFactor , &this->particles_ [ i ].pos );
            ScaleVector ( &this->particles_ [ i ].v , scaleFactor , &this->particles_ [ i ].v );
            ScaleVector ( &this->particles_ [ i ].f , scaleFactor , &this->particles_ [ i ].f );
        }
        return static_cast < System & > ( *this );
    }
    /**
     * \brief Functions like IntegrateSysOverTime and ComputeForces inside PhysEnv.cpp take their input as tParticle*.
     *
     * For compatibility, we need to be able to cast our system of particles to tParticle* to be able to pass our System as an argument to these functions
     */
    operator tParticle * ()
    {
        return this->particles_.data ();
    }
    /**
     * \brief Copies the particles of this system to the provided array. You might find this method useful. You might need to call s.fillOut(m_TargetSys) somewhere, which will copy the particles in the system "s" to the array of particles "m_TargetSys"
     * \param sys The array to copy the particles of this system to.
     */
    void fillOut ( tParticle * sys )
    {
        const std::size_t n = this->particles_.size ();
        std::memcpy ( static_cast < void * > ( sys ) , static_cast < void * > ( this->particles_.data () ) , n * sizeof ( tParticle ) );
    }
};
/**
 * \brief Overload of the multiplication operator. Scales all the vectors of the right system by a given factor.
 * \param left The factor to scale by.
 * \param right The System to scale.
 * \return A new system that is the result of the scaling operation.
 */
inline System operator * ( const float & left , const System & right )
{
    return right * left;
}
/**
 * \brief Overload of the division operator. De-scales all the vectors of the right system by a given factor.
 * \param left The factor to de-scales by.
 * \param right The System to de-scale.
 * \return A new system that is the result of the de-scaling.
 */
inline System operator / ( const float & left , const System & right )
{
    return right * ( 1 / left );
}
