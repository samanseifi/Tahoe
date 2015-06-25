// All data use SI: dimension--m; density--Kg/m^3; pressure--Pa; time--second                  

#ifndef PARAMETER_H
#define PARAMETER_H

namespace dem { 

////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1. time integration method 

int NUM_STEP              = 250000;
int NUM_PRINT             = 100;

//  --- dynamic
long double TIMESTEP      = 5.0e-07; // time step
long double MASS_SCL      = 1;       // mass scaling
long double MNT_SCL       = 1;       // moment of inertia scaling
long double GRVT_SCL      = 1;       // gravity scaling
long double DMP_F         = 0;       // background viscous damping on mass   
long double DMP_M         = 0;       // background viscous damping on moment of inertial

/*
// --- dynamic relaxation and scaling
long double TIMESTEP      = 5.0e-06;
long double MASS_SCL      = 1.0e+01;
long double MNT_SCL       = 1.0e+01;
long double GRVT_SCL      = 0;//1.0e+03;
long double DMP_F         = 2.0/TIMESTEP;
long double DMP_M         = 2.0/TIMESTEP;
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2. normal damping and tangential friction
long double DMP_CNT       = 0.05;    // damping ratio of viscous damping for normal contact force, for both particle-particle and particle-boundary contact
long double FRICTION      = 0.5;     // constant coefficient of static friction between particles
long double BDRYFRIC      = 0.5;     // constant coefficient of static friction between particle and rigid wall
long double COHESION      = 0;       // cohesion between particles (5.0e+8)
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3. boundary displacement rate
long double COMPRESS_RATE = 7.0e-03; // 7.0e-03 for triaxial; 1.0e-03 for isotropic and odometer.
long double RELEASE_RATE  = 7.0e-03;
long double PILE_RATE     = 2.5e-01; // pile penetration velocity
long double STRESS_ERROR  = 2.0e-02; // tolerance of stress equilibrium on rigid walls
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// 4. other global variables
long     idum             = -1;      // random number seed
std::ofstream g_exceptioninf;        // record debugging information
int      g_iteration;                // iteration number 
////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace dem ends

#endif
