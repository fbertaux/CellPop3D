

#include "contactdetector.hpp"

struct IntegratorDopr5
{
    // error tolerances
    Doub atol , rtol ;
    // numerical limit
    Doub EPS ;
    
    // storing the init state of a step (updated when step successful)
    VecDoub y_init_step ;
    VecDoub dydt_init_step ;

	// result of a step TRIAL
	VecDoub y_end_step , y_error_step ;
	VecDoub dydt_end_step ;
	Doub h_step , h_next_step ;
    Doub error_step , error_previous_step ;
    Bool reject_previous_step ;

	// computed and used during a step TRIAL
	VecDoub dydt_substep_2 , dydt_substep_3 , dydt_substep_4 , dydt_substep_5 , dydt_substep_6 ;
 	VecDoub y_substep ;

    // for contact detection
    ContactDetector contactdetector ;

 	// construction and re-initialization
    IntegratorDopr5 ( Doub atol , Doub rtol ) ;
    void reInitialize ( Int sys_size ) ;

    // higher-level routines : advance by a fixed time
    void integrateForDuration ( Doub duration , State* state ) ;

 	// high-level routine : realize one step (but can 'try' many steps)
	void doStep ( Doub h_init_step , State* state ) ;

	// mid-level routine : try a step
 	bool doStepTrial ( State* state ) ;

 	// low-level routines : within a step trial
 	void computeStep ( State* state ) ;	
    void computeError () ;
 	bool checkSuccess () ;
 
 	// the derivative function
 	void computeDerivatives ( State* state , VecDoub* y , VecDoub* dydt ) ;
} ;
