//////////////////////////////////////////////////////////////////////////////////////////
// The MEX part of the Word of Particles Simulation Engine to be called in the main loop
// of the WoP simulation. Implements MATLAB function WoP_Solver.
//
//  Usage:
//
//      result = WoP_Solver( isOddStep, isIgnorePP. isKillInactive, isKillOutOfBox )
//
//  where the 'result' is a structure with computed variables as fields: 
//
//      t, X, V, W, R, P, E, ccStat, N_p
//
//  The solver uses the following variables from the caller's scope:
//
//      M, R, X, V, W, NSFC, XSFC, HF_z, HF_ends,
//      g, k_D, k_M, k_W, v_fluid, gamma, mu_part, mu_sfc,
//      imodel, k_e, k_p, k_s, k_d
//      stepper, h
//
//  See the *WoP_Solver* class bellow for the description of all input, output
//  and referenced variables from the caller's scope.
//
//  CONVENTIONS:
//  Particles' variables are kept in rows while space dimensions are kept in columns.
//  Capital letters are used for M*N matrices or M*1 column vectors
//  and small letters are used for 1*N row vectors and scalars.
//
//  PERFORMANCE:
//  Example of performance comparison (P4/3GHz,Mem/2GB) between a MATLAB m-file and 
//  an equivalent implementation using the Matrix class.
//
//      N_p    t_f    h         NT    NT * N_p     m-file    mex/Matrix.h
//      10     10     0.001   10000    1e5          68 s       < 4 s
//     100     10     0.001   10000    1e6         580 s        20 s
//    1000     0.1    0.001     100    1e5          66 s        11 s
//   10000     0.001  0.001       1    1e4          36 s        20 s
//
//  Filename: WoP_Solver.cpp
//  Revision: 0.8
//  Date:     2012-04-12
//  Author:   Mikica B Kocic 

#include "mex.h"
#include "WoP_Matrix.h"

//////////////////////////////////////////////////////////////////////////////////////////
// The WoP_Solver class
// Implements an object-oriented wrapper the WoP_Solver mexFunction.
//
class WoP_Solver
{
    //////////////////////////////////////////////////////////////////////////////////////
    // Constants

    enum // Supported interaction models
    { 
        SpringModel  = 1,   // Linear-spring interaction model
        ImpulseModel = 2    // Impulse collision model
    };

    enum // Supported integrator methods
    { 
        ForwardEuler      = 1,   // Forward Euler integrator
        SemiImplicitEuler = 2,   // Semi-implicit Euler integrator
        Leapfrog          = 3    // Leapfrog integrator
    };

    //////////////////////////////////////////////////////////////////////////////////////

    MatrixHeap mResults;    // Heap for allocating Matrix objects used as results

    //////////////////////////////////////////////////////////////////////////////////////
    // Input variables
    //
    Logical mIgnorePP;      // Ignore particle-particle collisions
    Logical mKillInactive;  // Make inactive non-moving particles near ground surface
    Logical mKillOutOfBox;  // Make inactive particles going out of bounds

    //////////////////////////////////////////////////////////////////////////////////////
    // State variables from the caller's scope
    //
    Matrix M;               // Particle masses
    Matrix R;               // Particle radii; NaN for dead particles 
    Matrix X;               // Positions
    Matrix V;               // Velocities
    Matrix W;               // Angular velocities
    Matrix NSFC;            // Normal unit vectors to surfaces
    Matrix XSFC;            // Points on surfaces
    Matrix HF;              // Height field, z-component
    Matrix HF_ends;         // Height field x limits in 1st row and y limits in 2nd row
    Matrix X_sink;          // Sink position
    Matrix R_sink;          // Sink radius
    Matrix g;               // Acceleration of gravity
    Matrix k_D;             // Drag coefficient
    Matrix k_M;             // Magnus force coefficient
    Matrix k_W;             // Angular velocity damping coefficient
    Matrix mu_part;         // Friction coefficient, in contact between particles
    Matrix mu_sfc;          // Friction coefficient, in contact with surface
    Matrix v_fluid;         // Velocity of the surrounding fluid
    Matrix gamma;           // Viscous force gamma factor
    Matrix imodel;          // Interaction method
    Matrix k_e;             // Impulse collision model: restitution
    Matrix k_p;             // Impulse collision model: projections factor
    Matrix k_s;             // Linear-spring model: spring coefficient
    Matrix k_d;             // Linear-spring model: damping coefficient
    Matrix stepper;         // Integrator method
    Matrix h;               // Integrator time step
    Matrix t_0;             // Initial time
    Matrix t;               // Current time
    Matrix nstep;           // Current integrator step (variable 'n' in caller's scope)

    //////////////////////////////////////////////////////////////////////////////////////
    // Calculated and temporary variables
    //
    Matrix F_tot;           // Total force
    Matrix VJ_tot;          // Velocity jolts
    Matrix XJ_tot;          // Position projections
    Matrix W_dot;           // Angular acceleration
    Matrix zero;            // A 'zero' spatial vector

    //////////////////////////////////////////////////////////////////////////////////////
    // Returned results (in addition to solved state variables X and V)
    //
    Logical valid;          // True if time step is complete and ODE solved
    Matrix  P;              // Linear momentum
    Matrix  E;              // Energy (columns: kinetic, potential and total energy)
    Matrix  ccStat;         // Collision statistics; columns: part-part and part-surface

    //////////////////////////////////////////////////////////////////////////////////////
    // System parameters
    //
    int N_p;                // Number of particles in the system
	int N_dim;              // Number of spatial dimensions
    int N_sfc;              // Number of surfaces
    int N_p_alive;          // Number of active particles in the system
    bool mOddStep;          // Indicates if it is an odd half-step (used in leapfrog)
    double mLowestSpeed;    // Particles gets killed if their speed drops bellow the limit

public:

    //////////////////////////////////////////////////////////////////////////////////////
    // Returns true if equations are solved and state variables are valid.
    //
    bool IsSolutionValid () const 
    { 
        return valid;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Constructs an instance of the WoP solver from the mexFunction input arguments
    // and required variables from the caller's scope.
    //
    WoP_Solver(
        int nArgIn, const mxArray* argIn[] 
    )
        : mResults( "Results Heap" )
        //
        // Input variables ---------------------------------------------------------------
        //
        , mIgnorePP     ( nArgIn, argIn, 0, "IgnorePP"     )  // 1nd argument
        , mKillInactive ( nArgIn, argIn, 1, "KillInactive" )  // 2nd argument
        , mKillOutOfBox ( nArgIn, argIn, 2, "KillOutOfBox" )  // 3nd argument
        //
        // Variables from caller's workspace ---------------------------------------------
        //
        , M       ( "M"        )   // Particle masses
        , R       ( "R"        )   // Radii
        , X       ( "X"        )   // Positions
        , V       ( "V"        )   // Velocities
        , W       ( "W"        )   // Angular velocities, always 3D
        , NSFC    ( "NSFC"     )   // Normal unit vectors to surfaces
        , XSFC    ( "XSFC"     )   // Points on surfaces
        , HF      ( "HF_z"     )   // Height field z component
        , HF_ends ( "HF_ends"  )   // Height field x/y limits
        , X_sink  ( "X_sink"   )   // Sink position
        , R_sink  ( "R_sink"   )   // Sink radius
        , g       ( "g"        )   // Acceleration of gravity
        , k_D     ( "k_D"      )   // Drag coefficient
        , k_M     ( "k_M"      )   // Magnus force coefficient
        , k_W     ( "k_W"      )   // Magnus force coefficient
        , mu_part ( "mu_part"  )   // Friction coefficient, for particle-particle contact
        , mu_sfc  ( "mu_sfc"   )   // Friction coefficient, for particle-surface contact
        , gamma   ( "gamma"    )   // Viscous force gamma factor
        , v_fluid ( "v_fluid"  )   // Velocity of the surrounding fluid 
        , imodel  ( "imodel"   )   // Interaction model
        , k_e     ( "k_e"      )   // Impulse restitution
        , k_p     ( "k_p"      )   // Impulse projections factor
        , k_s     ( "k_s"      )   // Spring coefficient
        , k_d     ( "k_d"      )   // Damping coefficient
        , stepper ( "stepper"  )   // Integrator method
        , h       ( "h"        )   // Integrator time-step
        , nstep   ( "n"        )   // Current integrator step
        , t_0     ( "t_0"      )   // Initial time
        , t       ( "t"        )   // Elapsed time
    {
        dbgf( "\n%p : Entering WoP_Solver constructor\n\n", this );

        // Use mResults for storing Matrix objects
        //
        HeapStack::UseMxMalloc ();
        HeapStack results( mResults ); 

        // Get dimensions of the system
        //
        N_p       = X.GetM ();     // Number of particles
        N_dim     = X.GetN ();     // Number of spatial dimensions
        N_sfc     = NSFC.GetM ();  // Number of surfaces
        N_p_alive = 0;             // Number of active particles

        zero = Matrix( 1, N_dim, 0.0 ); // A 'zero' spatial vector

        // Verify that dimensions of input variables do agree
        //
        M       . VerifyDims( N_p,   N_dim,  "M dimension mismatch"               );
        R       . VerifyDims( N_p,   1,      "R dimension mismatch"               );
        X       . VerifyDims( N_p,   N_dim,  "X dimension mismatch"               );
        V       . VerifyDims( N_p,   N_dim,  "V dimension mismatch"               );
        W       . VerifyDims( N_p,   3,      "W dimension mismatch"               );
        NSFC    . VerifyDims( N_sfc, N_dim,  "NSFC dimension mismatch"            );
        XSFC    . VerifyDims( N_sfc, N_dim,  "XSFC dimension mismatch"            );
        k_D     . VerifyDims( N_p,   1,      "k_D dimension mismatch"             );
        k_M     . VerifyDims( N_p,   1,      "k_M dimension mismatch"             );
        k_W     . VerifyDims( N_p,   1,      "k_W dimension mismatch"             );
        g       . VerifyDims( 1,     N_dim,  "g should be a spatial vector"       );

        if ( ! v_fluid.IsEmpty () ) {
            v_fluid.VerifyDims( 1, N_dim, "v_fluid should be a spatial vector" );
        } else {
            v_fluid = zero;
        }

        if ( ! X_sink.IsEmpty () ) {
            X_sink.VerifyDims( 1, N_dim, "X_sink dimension mismatch"  );
            R_sink.VerifyDims( 1, 1,     "R_sink should be a scalar"  );
        }

        if ( N_p > 0 && N_dim > 0 ) 
        {
            k_e      . VerifyDims( 1, 1, "k_e should be a scalar" );
            k_p      . VerifyDims( 1, 1, "k_p should be a scalar" );
            k_s      . VerifyDims( 1, 1, "k_s should be a scalar" );
            k_d      . VerifyDims( 1, 1, "k_d should be a scalar" );
            gamma    . VerifyDims( 1, 1, "gamma should be a scalar" );
            mu_part  . VerifyDims( 1, 1, "mu_part should be a scalar" );
            mu_sfc   . VerifyDims( 1, 1, "mu_sfc should be a scalar" );
        }

        // Verify 'n' and 't', optionally
        //
        if ( N_p > 0 && N_dim > 0 && nArgIn > 0 )
        {
            t     . VerifyDims( 1, 1, "t should be a scalar" );
            nstep . VerifyDims( 1, 1, "n should be a scalar" );
        }

        // Initialize output and temporary variables on results stack
        //
        mOddStep = true;                  // Indicates if it is an odd half-step
        valid    = false;                 // ODE is not solved yet
        E        = Matrix( 1, 3, 0.0 );   // Accumulated energy for this time step
        ccStat   = Matrix( 1, 3, 0.0 );   // Accumulated collision statistics
        F_tot    . Resize( N_p, N_dim );  // Total force per particle
        VJ_tot   . Resize( N_p, N_dim );  // Velocity jolts
        XJ_tot   . Resize( N_p, N_dim );  // Position projections
        W_dot    . Resize( N_p, 3 );      // Angular acceleration, always 3D (!)
        mResults << R;                    // Move R to results stack

        // The lowest speed allowed for active particles depends on gravity
        // and integration time-step
        //
        mLowestSpeed = g.Norm () * h(0);

        // When called without arguments skip solving ODE and just calculate 
        // derived quantities from the state variables.
        //
        if ( nArgIn == 0 ) {
            valid = true;
        }

        dbgf( "%p : Exiting WoP_Solver constructor\n", this );
    }

    ~WoP_Solver ()
    {
        // mResults.ShowStatistics ();
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Displays current values of all variables (for debugging purposes)
    //
    void Display ()
    {
        disp( M       );  disp( R      ); 
        disp( X       );  disp( V      );  disp( W      );
        disp( NSFC    );  disp( XSFC   ); 
        disp( k_D     );  disp( k_M    );  disp( k_W    );
        disp( g       ); 
        disp( imodel  ); 
        disp( k_e     );  disp( k_p    ); 
        disp( k_s     );  disp( k_d    );
        disp( mu_part );  disp( mu_sfc );
        disp( v_fluid );  disp( gamma  );
        disp( stepper );  disp( h      );
        disp( F_tot   );  disp( VJ_tot );  disp( XJ_tot );  disp( W_dot );
        disp( P       );  disp( E      );  disp( ccStat );

        mexPrintf( "mOddStep:      %d\n", bool( mOddStep      ) );
        mexPrintf( "mIgnorePP:     %d\n", bool( mIgnorePP     ) );
        mexPrintf( "mKillInactive: %d\n", bool( mKillInactive ) );
        mexPrintf( "mKillOutOfBox: %d\n", bool( mKillOutOfBox ) );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // The equation solver methods
    //
    void CollisionDetectionAndResponse ();
    void HeightFieldCollisions ();
    void CalculateForces ();
    void SolveEquations ();
    void DerivedVariables ();
    void ReturnResults( int nArgOut, mxArray* argOut[] );
    
    //////////////////////////////////////////////////////////////////////////////////////
    // Private methods
    //
private:

    void CollisionResponse_Particles( 
        /* output */ Matrix& f_tot, Matrix& vj_tot, Matrix& xj_tot, double& E_p,
        /* input  */ bool isSpringModel, const Matrix& dX, const Matrix& dV, 
                     const Matrix& nX, const Matrix& m_ratio
    );

    void CollisionResponse_Surfaces( 
        /* output */ Matrix& f_tot, Matrix& vj_tot, Matrix& xj_tot, double& E_p,
        /* input  */ bool isSpringModel, const Matrix& dX, const Matrix& dV, 
                     const Matrix& nX 
    );
    
    void DisplayParticleInfo( int i, const char* title )
    {
        mexPrintf( "-------------------------------------------\n" );
        mexPrintf( "Particle %d  --> %s\n", i+1, title );
        mexPrintf( "-------------------------------------------\n" );
        mexPrintf( "t =  % 7.4e\n", t(0) );
        mexPrintf( "X = " ); X.Row( i ).Display ();
        mexPrintf( "V = " ); V.Row( i ).Display ();
        mexPrintf( "W = " ); W.Row( i ).Display ();
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
// Collision detection and collision response
// Search for colliding geometries and responds to collisions.
//
// Two contact models are implemented: 
// * the penalty method with linear-spring forces
// * the collision response impulse transfer method with position projections that 
//   resolve overlapping geometries, and 
// 
// See also: CollisionResponse_Particles and CollisionResponse_Particles methods
//
// NOTE: To avoid jittery behaviour in multiple-particle hard-core interaction, we should
// iterate the loop until velocity deflections (jolts) and position projections settle 
// down. Alternative is to use position projections factor k_p (in range between 0 and 1) 
// which softens overlappings over time (see usage of k_p bellow).
//
void WoP_Solver::CollisionDetectionAndResponse ()
{
    if ( valid || N_p <= 0 ) {
        return; // Refuse to proceede if ODE are already solved
    }

    dbgf( "\n%p : Entering WoP_Solver::CollisionDetectionAndResponse\n\n", this );

    //////////////////////////////////////////////////////////////////////////////////////
    // Get references to output variables' contents
    //
    double& pp_colc = ccStat(0);  // Particle-particle collision counter
    double& ps_colc = ccStat(1);  // Particle-surface collision counter

    //////////////////////////////////////////////////////////////////////////////////////
    // Verify interaction model
    //
    if ( int(imodel(0)) != SpringModel && int(imodel(0)) != ImpulseModel )
    {
        Matrix::SevereError( "WoP_Solver:invimodel",
            "Invalid interaction model %g", imodel(0)
        );
    }

    bool isSpringModel = int( imodel(0) ) == SpringModel;

    // Do always in linear-spring model but only on odd half-steps in impulse model 
    // with leapfrog integrator, i.e. don't check collisions for impulse model 
    // with leapfrog integrator on even half-steps. Note that A || B || C is 
    // equivalent to: A || ( ! A && B ) || ( !A && !B && C )
    //
    bool check_collisions = isSpringModel 
        || int(stepper(0)) != Leapfrog
        || mOddStep;

    if ( ! check_collisions )
    {
        // Nothing to do; just set forces and jolts to 0
        F_tot   = 0;
        VJ_tot  = 0;
        XJ_tot  = 0;
        return;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Compute forces and jolts for each particle (a)
    //
    for ( int i = 0; i < N_p; ++i ) //////////////////////////////////////////////////////           
    {
        if ( mxIsNaN( R(i) ) ) {
            continue; // Ignore dead particles
        }

        HeapStack stack; // Matrix instances go on a temporary heap

        //////////////////////////////////////////////////////////////////////////////////////
        // Kill all particles near the ground surface that either have speed bellow 
        // a certain limit or are in the sink circular area
        //
        if ( mKillInactive && X(i,N_dim-1) <= R(i,0) )
        {
            // Check if particle speed became too low
            //
            double speed = V.Row( i ).Norm ();

            if ( speed <= mLowestSpeed )
            {
                DisplayParticleInfo( i, "Stopped at z = 0" );
                R.SetRow( i, mxGetNaN () ); // Dead particle have NaN radius
                continue;
            }

            // Find if particle felt into the sink calculating norm of the displacement
            // and comparing it to the sink radius.
            //
            if ( ! X_sink.IsEmpty () )
            {
                Matrix displ = X.Row(i) - X_sink;
                displ(N_dim-1) = 0; // Calculate only horizontal distance

                if ( displ.Norm () <= R_sink(0) )
                {
                    ++ccStat(2); // Signal user that we have particle in the sink
                    DisplayParticleInfo( i, "In sink at z = 0" );
                    R.SetRow( i, mxGetNaN () ); // Dead particle have NaN radius
                    continue;
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////

        Matrix f_tot  = zero;  // Total force for for current particle (a)
        Matrix vj_tot = zero;  // Total velocity jolt for current particle
        Matrix xj_tot = zero;  // Total position projection for current particle

        //////////////////////////////////////////////////////////////////////////////////
        // Particle to other particles interaction

        // Do not check particle-particle overlaps in single particle system
        // or if particle-particle interactions are ignored
        // 
        if ( N_p > 1 && ! mIgnorePP )
        {
            // Displacement between current particle (a) and other particles (b)
            //
            //   x_ab = x_a - x_b
            //
            Matrix dXp = X.Row(i) - X;

            // Distance between current particle (a) and other particles
            //
            //   |x_ab| = | x_a - x_b |
            //
            Matrix norm_dXp = dXp.Norm_ByRow ();

            // Overlap i.e. penetration depth from other particles
            //
            //   dx_ab = r_a + r_b - | x_a - x_b |
            //
            Matrix dX = R.Row(i) + R - norm_dXp;

            // Find all overlaps, excluding self-overlap.
            // Particle (a) overlaps with (b) if dx_ab > 0
            //
            int count = 0;
            unsigned char* overlaps = dX.FindPositive( count ); 
            // Exclude the self-overlap. Note that overlaps is a column vector.
            count -= overlaps[i] != 0;
            overlaps[i] = 0;

            if ( count > 0 ) // if there was overlaps ====================================
            {
                // Update number of collisions between particles in this time-step.
                // Since the particle interacts always with some other particle, we
                // should count only half of the interactions per one particle.
                //
                pp_colc = pp_colc + double( count ) / 2;

                // From now on, consider only overlapping particles...
                //
                dX  = dX( overlaps, count );
                dXp = dXp( overlaps, count );

                // Unit direction from b to a:
                //
                //   n_ab = ( x_a - x_b ) / | x_a - x_b |
                //
                Matrix nX = dXp.UnitVector_ByRow ();

                // Relative velocities from other overlapping particles:
                //
                //   v_ab = v_a - v_b
                //
                Matrix dV = V.Row(i) - V( overlaps, count );

                Matrix m_ratio;
                if ( ! isSpringModel )
                {
                    // Ratio of the reduced mass m_ab to particle (a) mass m_a:
                    //
                    //   m_ratio = m_red_ab * m_b^-1
                    //           = m_a^-1 / ( m_a^-1 + m_b^-1 )
                    //           = m_b / ( m_a + m_b )
                    //
                    double m_i = M(i,0); // Current particle mass m_a
                    Matrix m_overlaps = M.Column(0)( overlaps, count ); // m_b
                    m_ratio = m_overlaps / ( m_i + m_overlaps );
                }
                
                // Finally, calculate the response
                //
                CollisionResponse_Particles(
                    /* output */ f_tot, vj_tot, xj_tot, /* E_p = */ E(1),
                    /* input  */ isSpringModel, dX, dV, nX, m_ratio
                );

            } // end if count > 0  =======================================================
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Particle to surfaces interaction
        // 
        {
            // Distance between current particle (a) and surfaces (s):
            //
            //   x_as = < x_a - r_s, n_s >
            //
            Matrix dXp = ( X.Row(i) - XSFC ).DotProd_ByRow( NSFC );

            // Penetration depth (overlap) of particle (a) to surface (s):
            //
            //   dx_as = r_a - x_as
            //
            Matrix dX = R.Row(i) - dXp;

            // Find all overlaping surfaces (s).
            // Particle (a) overlaps with surface (s) if dx_as > 0
            //
            int count = 0;
            unsigned char* overlaps = dX.FindPositive( count ); 

            if ( count > 0 ) // if there was overlaps ====================================
            {
                // Update number of collisions with surfaces in this time-step
                //
                ps_colc = ps_colc + double( count );

                // Inactivate particles touching box surfaces
                //
                if ( mKillOutOfBox )
                {
                    DisplayParticleInfo( i, "*** Out of Bounds ***" );
                    R.SetRow( i, mxGetNaN () ); // Inactive particle have NaN radius
                    continue;
                }

                // From now on, consider only overlapping surfaces...
                //
                dX = dX( overlaps, count );
                
                // Unit direction from the overlapping surface
                //
                //   n_ab = NSFC( overlapping surfaces )
                //
                Matrix nX = NSFC( overlaps, count );

                // Relative velocity of particle (a) to surface (s).
                // Surfaces are stationary i.e. v_s = 0.
                //
                //   v_ab = v_a - v_s = v_a - 0 = v_a
                //
                Matrix dV = V.Row(i) - Matrix( dX.GetM (), N_dim, 0.0 );

                // Finally, calculate the response
                //
                CollisionResponse_Surfaces(
                    /* output */ f_tot, vj_tot, xj_tot, /* E_p = */ E(1),
                    /* input  */ isSpringModel, dX, dV, nX
                );

            } // end if count > 0  =======================================================
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Save the total force and jolts for the particle (i)
        //
        F_tot .SetRow( i, f_tot  );
        VJ_tot.SetRow( i, vj_tot );
        XJ_tot.SetRow( i, xj_tot );

    }  // for each particle (i) //////////////////////////////////////////////////////////

    dbgf( "\n%p : Entering WoP_Solver::CollisionDetectionAndResponse\n\n", this );
}

//////////////////////////////////////////////////////////////////////////////////////////
// Detects and solves collisions with height field.
//
void WoP_Solver::HeightFieldCollisions ()
{
    if ( valid || N_p <= 0 ) {
        return; // Refuse to proceede if ODE are already solved
    }

    if ( HF.IsEmpty () || HF_ends.IsEmpty () || N_dim == 1 ) {
        return; // Cowardly refuse to deal with empty HF matrices or when in 1D
    }

    // Do always in linear-spring model but only on odd half-steps in impulse model 
    // with leapfrog integrator, i.e. don't check collisions for impulse model 
    // with leapfrog integrator on even half-steps. Note that A || B || C is 
    // equivalent to: A || ( ! A && B ) || ( !A && !B && C )
    //
    bool isSpringModel = int( imodel(0) ) == SpringModel;
    bool check_collisions = isSpringModel 
        || int(stepper(0)) != Leapfrog
        || mOddStep;

    if ( ! check_collisions ) {
        return; // Nothing to do
    }

    dbgf( "\n%p : Entering WoP_Solver::HeightFieldCollisions\n\n", this );

    //////////////////////////////////////////////////////////////////////////////////////
    // Compute jolts for each particle (i)
    //
    for ( int i = 0; i < N_p; ++i ) //////////////////////////////////////////////////////           
    {
        if ( mxIsNaN( R(i) ) ) {
            continue; // Ignore dead particles
        }

        HeapStack stack; // Matrix instances go on a temporary heap

        // Split particle position into x, y and z components
        //
        Matrix xp = X.Row(i);
        double x = xp(0);
        double y = N_dim == 3 ? xp(1) : 0;
        double z = N_dim == 3 ? xp(2) : xp(1);
        
        // Find z of the height field at particle (x,y) position
        //
        double hf = HF.Interpolate( x, y, HF_ends );
        if ( mxIsInf( hf ) ) {
            // If there is no HF at x,y assume ground level height to be 0
            hf = 0;
        }

        // Distance between current particle (a) and the height field
        //
        //   x_as = x_a(z) - hf( x_a.x, x_a.y )
        //
        double dXp = z - hf;

        // Penetration depth (overlap) of particle (a) to height field surface
        // Works only if particle radius is much smaller than surface roughness.
        //
        //   dx_as = r_a - x_as
        //
        Matrix dX( 1, 1, R(i) - dXp );

        if ( dX(0) <= 0 ) {
            continue; // No overlap
        }
        
        dbgf( "Particle %d overlaps hf = %g at (%g,%g) z = %g\n", i, hf, x, y, z );

        // Relative velocity of particle (a) to surface (s).
        // Surfaces are stationary i.e. v_s = 0.
        //
        //   v_ab = v_a - v_s = v_a - 0 = v_a
        //
        Matrix dV = V.Row(i);

        // Kill all particles near the ground surface that either have speed bellow 
        // a certain limit or are in the sink circular area
        //
        if ( mKillInactive )
        {
            // If particle velocity drops bellow the limit in contct with the surface
            // make it inactive (dead) by setting radius to NaN.
            //
            if ( dV.Norm () <= mLowestSpeed )
            {
                DisplayParticleInfo( i, "Stopped on HF" );
                R.SetRow( i, mxGetNaN () ); // Dead particle have NaN radius
                continue;
            }

            // Find if particle felt into the sink calculating norm of the displacement
            // and comparing it to the sink radius.
            //
            if ( ! X_sink.IsEmpty () )
            {
                Matrix displ =  X.Row(i) - X_sink;
                displ(N_dim-1) = 0; // Calculate only horizontal distance

                if ( displ.Norm () <= R_sink(0) )
                {
                    ++ccStat(2); // Signal user we have particle in the sink
                    DisplayParticleInfo( i, "In sink on HF" );
                    R.SetRow( i, mxGetNaN () ); // Dead particle have NaN radius
                    continue;
                }
            }
        }

        // Find differentials at point (x,y) approximated for delta x or y as particle 
        // radius centered around (x,y). This assumes that the particle radius
        // is smaller than height field surface roughness!
        //
        double delta = R(i,0)/3;
        double z_prim_x = HF.Interpolate( x + delta/2, y, HF_ends )
                        - HF.Interpolate( x - delta/2, y, HF_ends );
        double z_prim_y = HF.Interpolate( x, y + delta/2, HF_ends )
                        - HF.Interpolate( x, y - delta/2, HF_ends );

        // Now, calculate tangent vectors along x and y.
        // We should handle 2-D and 3-D cases separately. Note that we should not
        // suppose to be here in case of 1-D.
        //
        Matrix t_x = Matrix( 1, N_dim, 0.0 );
        Matrix t_y = Matrix( 1, N_dim, 0.0 );
        switch( N_dim )
        {
            case 2:
                // Tangent t_x: Derivative z_prim_x for constant y
                t_x(0) = delta;    // delta X
                t_x(1) = z_prim_x; // delta Z
                break;
            case 3:
                // Tangent t_x: Derivative z_prim_x for constant y
                t_x(0) = delta;    // delta X
                t_x(1) = 0;        // constant Y
                t_x(2) = z_prim_x; // delta Z
                // Tangent t_y: Derivative z_prim_y for constant x
                t_y(0) = 0;        // contant X
                t_y(1) = delta;    // delta Y
                t_y(2) = z_prim_y; // delta Z
                break;
        }

        // Normal vector at the n-1 dimensional surface is outer product of the 
        // tangents t_1 ... t_m to the surface, where m = n - 1 in n-dimensional case:
        //
        //    n_ab = [ t_1, ..., t_m ] / | [ t_1, ..., t_m ] |
        //
        // Note that in 1-dimensional space there are no 0-dimensional 'hypersurfaces',
        // which agrees with the initial sanity-check for this method that N_dim != 1.
        // Ref: http://en.wikipedia.org/wiki/Normal_%28geometry%29
        //
        Matrix nX;
        switch ( N_dim )
        {
            case 2:
                // Calculate n = ( -y, x ) i.e. as 90 degree rotation counter-clockwise
                // (We can ue t_y as temporary variable for this.)
                //
                t_y(0) = - t_x(1);
                t_y(1) =   t_x(0);
                nX = t_y.UnitVector_ByRow ();
                break;

            case 3:
                // Crossproduct n_ab = ( t_x × t_y ) / | t_x × t_y |
                //
                nX = ( t_x.Cross( t_y ) ).UnitVector_ByRow ();
                break;
        }

        // Finally, calculate the response (total force and jolts)
        //
        Matrix f_tot  = F_tot .Row( i ); // Total force for for current particle (a)
        Matrix vj_tot = VJ_tot.Row( i ); // Total velocity jolt for current particle
        Matrix xj_tot = XJ_tot.Row( i ); // Total position projection for current particle
        double& E_p   = E(1);            // References to potential energy

        CollisionResponse_Surfaces(
            /* output */ f_tot, vj_tot, xj_tot, E_p,
            /* input  */ isSpringModel, dX, dV, nX
        );

        // and save the total force and jolts for the particle (i)
        //
        F_tot .SetRow( i, f_tot  );
        VJ_tot.SetRow( i, vj_tot );
        XJ_tot.SetRow( i, xj_tot );

    }  // for each particle (i) //////////////////////////////////////////////////////////

    dbgf( "\n%p : Entering WoP_Solver::HeightFieldCollisions\n\n", this );
}

//////////////////////////////////////////////////////////////////////////////////////////
// Collision response for particle-particle collisions (private method)

void WoP_Solver::CollisionResponse_Particles
( 
    /* output */ Matrix& f_tot, Matrix& vj_tot, Matrix& xj_tot, double& E_p,
    /* input  */ bool isSpringModel, const Matrix& dX, const Matrix& dV, 
                 const Matrix& nX, const Matrix& m_ratio
    )
{
    // Get the projection of v_ab on n_ab (scalar product)
    //
    //   v_n_ab = < v_ab, n_ab >
    //
    Matrix dV_n = dV.DotProd_ByRow( nX );

    // Get the unit of tangential component of v_ab
    //
    //   v_t_ab = n_ab × ( v_ab × n_ab ) 
    //          = v_ab - n_ab * < v_ab, n_ab >
    //          = v_ab - n_ab * v_n_ab
    //
    //   t_ab = v_t_ab / |v_t_ab|
    //
    Matrix dV_t = dV - nX * dV_n;
    Matrix tX = dV_t.UnitVector_ByRow ();

    if ( isSpringModel ) // --------------------------------------------------
    {
        // Total penalty forces from overlapping particles (sum per column)
        //
        //   f_ab = k_s * < x_ab, n_ab > - k_d * < v_n_ab, n_ab >
        //
        f_tot += ( ( k_s * dX - k_d * dV_n ) * nX ).Sum_ByColumn ();

        // Accumulate the total potential energy.
        // Since we are summing spring energy twice (for each side of the
        // spring i.e. for each particle) we should sum only 1/2 of the
        // usual E_p = 1/2 * k_s * |dX|^2.
        //
        E_p += 0.25 * k_s(0) * dX.SquaredNorm ();
    }
    else // ImpulseModel  ----------------------------------------------------
    {
        // Suppress impulse jolt for separating particles
        // and calculate base velocity jolt: vj = j/m
        //
        //   v_n_ab( v_n_ab > 0 ) = 0    i.e. remove separating contacts
        //   v_j = -( 1 + k_e ) * m_ratio * v_ab
        //
        Matrix vj = -( 1.0 + k_e ) * m_ratio * dV_n.ReplacePositive( 0 );

        // Total impulse jolt from overlapping particles (sum per column)
        //
        //   vj_ab = v_j * n_ab - CLAMP( mu * v_j * t_ab, v_t_ab )
        //
        vj_tot += ( vj * nX - ( vj * mu_part * tX ).Clamp( dV_t )
                  ).Sum_ByColumn ();

        // Position projection from overlapping particles (multiplied by a
        // position projectsion factor k_p where k_p = 0 disables
        // projections, and k_p = 1 fully enables projections while value
        // inbetween softens overlappings)
        //
        //   xj_ab = k_p * m_ratio < x_ab, n_ab >
        //
        if ( k_p(0) != 0 ) {
            xj_tot += ( ( m_ratio * dX ) * nX ).Sum_ByColumn () * k_p;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////
// Collision response for particle-surface collisions (private method)

void WoP_Solver::CollisionResponse_Surfaces
( 
    /* output */ Matrix& f_tot, Matrix& vj_tot, Matrix& xj_tot, double& E_p,
    /* input  */ bool isSpringModel, const Matrix& dX, const Matrix& dV, const Matrix& nX 
    )
{
    // Get the projection of v_ab on n_ab (scalar product)
    //
    //   v_n_ab = < v_ab, n_ab >
    //
    Matrix dV_n = dV.DotProd_ByRow( nX );

    // Get the unit of tangential component of v_ab
    //
    //   v_t_ab = n_ab × ( v_ab × n_ab ) 
    //          = v_ab - n_ab * < v_ab, n_ab >
    //          = v_ab - n_ab * v_n_ab
    //
    //   t_ab = v_t_ab / |v_t_ab|
    //
    Matrix dV_t = dV - nX * dV_n;
    Matrix tX = dV_t.UnitVector_ByRow ();

    if ( isSpringModel ) // --------------------------------------------------
    {
        // Total penalty forces from overlapping surfaces (sum per column)
        //
        //   f_ab = k_s * < x_ab, n> - k_d * < v_n_ab, n >
        //
        f_tot += ( ( k_s * dX - k_d * dV_n ) * nX ).Sum_ByColumn ();

        // Accumulate the total potential energy E_p = 0.5 * k_s * |dX|^2
        //
        E_p += 0.5 * k_s(0) * dX.SquaredNorm ();
    }
    else // ImpulseModel -----------------------------------------------------
    {
        // Suppress impulse jolt when separating from the surface
        // and calculate base velocity jolt: vj = j/m
        //
        //   v_n_ab( v_n_ab > 0 ) = 0       i.e. remove separating contacts
        //   v_j = -( 1 + k_e ) * v_n_ab
        //
        Matrix vj = -( 1.0 + k_e ) * dV_n.ReplacePositive( 0 );

        // Total impulse jolt from overlapping surfaces (sum per column)
        //
        //   vj_ab = v_j * n_ab - CLAMP( mu * v_j * t_ab, v_t_ab )
        //
        vj_tot += ( vj * nX - ( vj * mu_sfc * tX ).Clamp( dV_t )
                  ).Sum_ByColumn ();

        // Position projection from overlappings surfaces
        //
        //   xj_ab = k_p * < x_ab, n_ab >
        //
        if ( k_p(0) != 0 ) {
            xj_tot += k_p * ( dX * nX ).Sum_ByColumn ();
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////
// Compute Forces and Constraints
// Calculate external forces on the objects and internal forces (except penalty
// forces which were computed earlier).
//
void WoP_Solver::CalculateForces ()
{
    if ( valid || N_p <= 0 ) {
        return; // Refuse to proceede if ODE are already solved
    }

    dbgf( "\n%p : Entering WoP_Solver::CalculateForces\n\n", this );

    N_p_alive = 0; // Number of active particles

    //////////////////////////////////////////////////////////////////////////////////////
    // Compute forces for each particle (i)
    //
    for ( int i = 0; i < N_p; ++i ) //////////////////////////////////////////////////////           
    {
        if ( mxIsNaN( R(i) ) )
        {
            // Clear all forces, acceleration, jolts and velocities for dead particles
            F_tot  .SetRow( i, 0 );
            VJ_tot .SetRow( i, 0 );
            XJ_tot .SetRow( i, 0 );
            W_dot  .SetRow( i, 0 );
            V      .SetRow( i, 0 );
            W      .SetRow( i, 0 );
            continue;
        }

        ++N_p_alive;       // Count active particles
        
        HeapStack stack;     // Matrix instances go on a temporary heap

        Matrix f_tot = F_tot.Row( i );  // Load total force

        //////////////////////////////////////////////////////////////////////////////////
        // Calculate the gravitational force:
        //
        //   f_g = M * G
        //
        f_tot += M.Row(i) * g;

        //////////////////////////////////////////////////////////////////////////////////
        // For viscous forces and Magnus force we need to calculate velocity relative 
        // to surrounding fluid
        //
        Matrix v    = V.Row(i) - v_fluid;  // Relative velocity of the object to the fluid
        double v_norm = v.Norm ();         // |V_rel|

        //////////////////////////////////////////////////////////////////////////////////
        // Calculate the viscous force, where V is velocity relative to the fluid:
        //
        //   f_visc = - k_D * |V| ^ gamma * V/|V|
        //          = - k_D * |V| ^ (gamma-1) * V
        //
        f_tot -= k_D.Row(i) * pow( v_norm, gamma(0) - 1 ) * v;

        //////////////////////////////////////////////////////////////////////////////////
        // Calculate the Magnus force, where V is velocity relative to the fluid:
        //
        //   f_M = - k_M * |V| * ( W × V )
        //
        Matrix w = W.Row(i);
        Matrix w_cross_v;
        switch( N_dim )
        {
            case 3:
                // Cross product
                w_cross_v = w.Cross( v );
                break;
            case 2:
                // Assuming rotation in XZ plane around Y axis (in N_dim = 2 kept vector 
                // components are X and Z) we have W(i,0) = W(i,1) = 0 so:
                w_cross_v    = Matrix( 1, 2 );
                w_cross_v(0) = - w(2) * v(1);
                w_cross_v(1) =   w(2) * v(0);
                break;
        }
        f_tot += k_M.Row(i) * v_norm * w_cross_v;

        //////////////////////////////////////////////////////////////////////////////////
        // Save the total force
        //
        F_tot.SetRow( i, f_tot ); 

        //////////////////////////////////////////////////////////////////////////////////
        // Calculate the angular velocity damping acceleration, 
        // where V is velocity relative to the fluid:
        //
        //    w_dot = - k_W * v_w_term * W
        //    where v_w_term = ( R*|W| + |V| )^(gamma-1)
        //
        double v_w_term = pow( R(i) * w.Norm () + v_norm, gamma(0) - 1 );
        W_dot.SetRow( i, - k_W.Row(i) * v_w_term * w );

    }  // for each particle (i) //////////////////////////////////////////////////////////

    dbgf( "%p : Exiting WoP_Solver::CalculateForces\n", this );
}

//////////////////////////////////////////////////////////////////////////////////////////
// Steps Forward and Updates State Variables
// Advances the simulation data from the current point of time to the next and
// also computes the new state variables (solving a system of equations).

void WoP_Solver::SolveEquations ()
{
    if ( valid ) {
        return; // Refuse to proceede if ODE are already solved.
    }

    dbgf( "\n%p : Entering WoP_Solver::SolveEquations\n\n", this );

    HeapStack stack; // Matrix instances go on a temporary heap

    Matrix A = F_tot / M; // Calculate the acceleration

    // NOTE: at this time we should also calculate angular acceleration W_dot from 
    // total torque. Torque is already divided by the moment of inertia i.e. it is 
    // included in k_W coefficient.

    switch( int(stepper(0)) )
    {
        case ForwardEuler: ///////////////////////////////////////////////////////////////

            X = X + V * h + VJ_tot;     // Solve position
            V = V + A * h + XJ_tot;     // Solve velocity
            W = W + W_dot * h;          // Solve angular velocity
            t = t_0 + nstep * h;        // Solve time
            valid = true;               // ODEs are solved and state variables are valid
            break;

        case SemiImplicitEuler: //////////////////////////////////////////////////////////

            V = V + A * h + VJ_tot;     // Solve velocity
            X = X + V * h + XJ_tot;     // Solve position
            W = W + W_dot * h;          // Solve angular velocity
            t = t_0 + nstep * h;        // Solve time
            valid = true;               // ODEs are solved and state variables are valid
            break;

        case Leapfrog: ///////////////////////////////////////////////////////////////////

            // In the leapfrog , one call of SolveEquations() corresponds to one
            // half-step (h/2).

            V = V + A * h/2;            // Solve velocity each half-step
            W = W + W_dot * h/2;        // Solve angular velocity

            if ( mOddStep ) {           // If it is the first (odd) half-step:
                mOddStep = false;       // - mark next step as the final (even) half-step
                V = V + VJ_tot;         // - apply velocity jolts
                X = X + V * h + XJ_tot; // - solve position and do position projections
                t = t_0 + ( nstep - 0.5 ) * h; // - do half-step of elapsed time
                E(1) = 0;               // - reset accumulated E_p_n
                valid = false;          // - ODE are still not solved (it's not full-step)
            }
            else {                      // If it is the second (even) half-step:
                t = t_0 + nstep * h;    // - do full-step of elapsed time
                valid = true;           // - ODEs are solved and state variables are valid
            }

            break;

        default: /////////////////////////////////////////////////////////////////////

            Matrix::SevereError( "WoP_Solver:invstepper",
                "Invalid integrator method %g", stepper(0)
            );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Save the results (move the results to the results heap)

    mResults << t << X << V << W;

    dbgf( "%p : Exiting WoP_Solver::SolveEquations\n", this );
}

//////////////////////////////////////////////////////////////////////////////////////
// Computes Derived Quantities
// After the state variables are updated, we can update all derived quantities.

void WoP_Solver::DerivedVariables ()
{
    if ( ! valid ) {
        return; // Don't update derived quantities if ODE are not solved
    }

    dbgf( "\n%p : Entering WoP_Solver::DerivedVariables\n\n", this );

    HeapStack stack; // Matrix instances go on a temporary heap

    //////////////////////////////////////////////////////////////////////////////////////
    // Get references to output variables' contents

    double& E_k   = E(0);     // Kinetic energy
    double& E_p   = E(1);     // Potential energy
    double& E_tot = E(2);     // Total energy

    //////////////////////////////////////////////////////////////////////////////////////
    // Calculate the linear momentum

    P = M * V;  // Note: The total linear momentum of the system is: sum( P, 1 )

    //////////////////////////////////////////////////////////////////////////////////////
    // Compute the kinetic and the total energy of the system
    // Note that the potential energy of the penalty forces is accumulated earlier.

    E_k   = 0.5 * ( P * V ).Sum ();      // Kinetic energy
    E_p   = E_p - ( M * g * X ).Sum ();  // Potential energy
    E_tot = E_k + E_p;                   // Total energy

    //////////////////////////////////////////////////////////////////////////////////////
    // Save the results

    mResults << P;

    dbgf( "%p : Exiting WoP_Solver::DerivedVariables\n", this );
}

//////////////////////////////////////////////////////////////////////////////////////////
// Returns results to the caller of the mexFunction a structure with computed variables.
//
void WoP_Solver::ReturnResults( int nArgOut, mxArray* argOut[] )
{
    dbgf( "\n%p : Entering WoP_Solver::ReturnResults\n\n", this );

    if ( nArgOut == 0 ) {
        mResults.ShowStatistics ();
        HeapStack::ShowStatistics ();
    }

    if ( nArgOut >= 1 ) 
    {
        // Return a structure with computed variables
        //
        static const mwSize dims = { 1 };
        static const char* fields[] = 
        { 
            "t", "X", "V", "W", "R", "P", "E", "ccStat", "N_p"
        };
        mxArray* result = mxCreateStructArray( 1, &dims, 
                sizeof( fields ) / sizeof( fields[0] ), fields );

        // "valid", mxSetField( result, 0,  "valid",  valid   );
        mxSetField( result, 0,  "t",      t       );
        mxSetField( result, 0,  "X",      X       );
        mxSetField( result, 0,  "V",      V       );
        mxSetField( result, 0,  "W",      W       );
        mxSetField( result, 0,  "R",      R       );
        mxSetField( result, 0,  "P",      P       );
        mxSetField( result, 0,  "E",      E       );
        mxSetField( result, 0,  "ccStat", ccStat  );
        mxSetField( result, 0,  "N_p",    mxCreateDoubleScalar( N_p_alive ) );

        argOut[0] = result;
    }
    
    dbgf( "%p : Exiting WoP_Solver::ReturnResults\n", this );
}

//////////////////////////////////////////////////////////////////////////////////////////
// The main entry routine for the MATLAB function WoP_Solver.
// See the *WoP_Solver* class above for the description of input and output variables.
//
void mexFunction
(
    int nArgOut, mxArray* argOut[],       // left-hand side (varargout) of the function
    int nArgIn, const mxArray* argIn[]    // right-hand side (varargin) of the function
    )
{
    WoP_Solver solver( nArgIn, argIn );      // Setup solver from input variables
    
    while ( ! solver.IsSolutionValid () )    // Cycle half-steps until ODE are solved
    {
        solver.CollisionDetectionAndResponse (); // Particle-particle and particle-surface
        solver.HeightFieldCollisions ();         // Hight field collision det. & response
        solver.CalculateForces ();               // Calculate internal and external forces
        solver.SolveEquations ();                // Solve ODE of motion
    }

    solver.DerivedVariables ();              // Compute derived quantities
    solver.ReturnResults( nArgOut, argOut ); // Return the results
}
