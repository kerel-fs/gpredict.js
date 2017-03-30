/*
 * Unit SGP_Math
 *       Author:  Dr TS Kelso
 * Original Version:  1991 Oct 30
 * Current Revision:  1998 Mar 17
 *          Version:  3.00
 *        Copyright:  1991-1998, All Rights Reserved
 *
 *   ported to C by:  Neoklis Kyriazis  April 9 2001
 */

/* Returns sign of a double */
function Sign(arg)
{
  if( arg > 0 )
    return( 1 );
  else if( arg < 0 )
    return( -1 );
  else
    return( 0 );
} /* Function Sign*/

/*------------------------------------------------------------------*/

/* Returns square of a double */
function Sqr(arg)
{
  return( arg*arg );
} /* Function Sqr */

/*------------------------------------------------------------------*/

/* Returns cube of a double */
function Cube(arg)
{
  return( arg*arg*arg );
} /*Function Cube*/

/*------------------------------------------------------------------*/

/* Returns angle in radians from arg id degrees */
function Radians(arg)
{
  return( arg*de2ra );
} /*Function Radians*/

/*------------------------------------------------------------------*/

/* Returns angle in degrees from arg in rads */
function Degrees(arg)
{
  return( arg/de2ra );
} /*Function Degrees*/

/*------------------------------------------------------------------*/

/* Returns the arcsine of the argument */
var ArcSin = Math.asin;

/*------------------------------------------------------------------*/

/* Returns orccosine of rgument */
var ArcCos = Math.acos;

/** \brief Arccosine implementation. 
*
* Returns a value between zero and two pi.
* Borrowed from gsat 0.9 by Xavier Crehueras, EB3CZS.
* Optimized by Alexandru Csete.
*/
function arccos (x, y)
{
   if (x && y) {
       if (y > 0.0)
           return Math.acos (x/y);
       else if (y < 0.0)
           return pi + Math.acos (x/y);
   }

   return 0.0;
}

/*------------------------------------------------------------------*/

/* Calculates scalar magnitude of a vector_t argument */
function Magnitude(v)
{
  v.w = sqrt(Sqr(v.x) + Sqr(v.y) + Sqr(v.z));
} /*Procedure Magnitude*/

/*------------------------------------------------------------------*/

/* Adds vectors v1 and v2 together to produce v3 */
function Vec_Add(v1, v2, v3)
{
  v3.x = v1.x + v2.x;
  v3.y = v1.y + v2.y;
  v3.z = v1.z + v2.z;

  Magnitude(v3);
} /*Procedure Vec_Add*/

/*------------------------------------------------------------------*/

/* Subtracts vector v2 from v1 to produce v3 */
function Vec_Sub(v1, v2, v3)
{
  v3.x = v1.x - v2.x;
  v3.y = v1.y - v2.y;
  v3.z = v1.z - v2.z;

  Magnitude(v3);
} /*Procedure Vec_Sub*/

/*------------------------------------------------------------------*/

/* Multiplies the vector v1 by the scalar k to produce the vector v2 */
function Scalar_Multiply(k, v1, v2)
{
  v2.x = k * v1.x;
  v2.y = k * v1.y;
  v2.z = k * v1.z;
  v2.w = fabs(k) * v1.w;
} /*Procedure Scalar_Multiply*/

/*------------------------------------------------------------------*/

/* Multiplies the vector v1 by the scalar k */
function Scale_Vector(k, v)
{ 
  v.x *= k;
  v.y *= k;
  v.z *= k;
  Magnitude(v);
} /* Procedure Scale_Vector */

/*------------------------------------------------------------------*/

/* Returns the dot product of two vectors */
function Dot(v1, v2)
{
  return( v1.x*v2.x + v1.y*v2.y + v1.z*v2.z );
}  /*Function Dot*/

/*------------------------------------------------------------------*/

/* Calculates the angle between vectors v1 and v2 */
function Angle(v1, v2)
{
  Magnitude(v1);
  Magnitude(v2);
  return( ArcCos(Dot(v1,v2)/(v1.w*v2.w)) );
} /*Function Angle*/

/*------------------------------------------------------------------*/

/* Produces cross product of v1 and v2, and returns in v3 */
function Cross(v1, v2, v3)
{
  v3.x = v1.y*v2.z - v1.z*v2.y;
  v3.y = v1.z*v2.x - v1.x*v2.z;
  v3.z = v1.x*v2.y - v1.y*v2.x;
  Magnitude(v3);
} /*Procedure Cross*/

/*------------------------------------------------------------------*/

/* Normalizes a vector */
function Normalize( v )
{
  v.x /= v.w;
  v.y /= v.w;
  v.z /= v.w;
} /*Procedure Normalize*/

/*------------------------------------------------------------------*/

/* Four-quadrant arctan function */
var AcTan = Math.atan2;

/*------------------------------------------------------------------*/

/* Returns mod 2pi of argument */
function FMod2p(x)
{
  return (x % twopi);
} /* fmod2p */

/*------------------------------------------------------------------*/

/* Returns arg1 mod arg2 */
function Modulus(arg1, arg2)
{
  return (arg1 % arg2);
} /* modulus */

/*------------------------------------------------------------------*/

/* Returns fractional part of double argument */
function Frac( arg )
{
  return( arg - floor(arg) );
} /* Frac */

/*------------------------------------------------------------------*/

/* Returns argument rounded up to nearest integer */
var Round = Math.round;

/*------------------------------------------------------------------*/

/* Returns the floor integer of a double arguement, as double */
var Int = Math.floor;

/*------------------------------------------------------------------*/

/* Converts the satellite's position and velocity  */
/* vectors from normalised values to km and km/sec */ 
function Convert_Sat_State( pos, vel )
{
      Scale_Vector( xkmper, pos );
      Scale_Vector( xkmper*xmnpda/secday, vel );

} /* Procedure Convert_Sat_State */

/*------------------------------------------------------------------*/
