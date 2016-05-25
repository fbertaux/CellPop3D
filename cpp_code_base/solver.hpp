
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "integrator.hpp"
#include "random.hpp"


struct Solver
{
	// construction
	Solver ( State* s ) ;

	// the state that has to evolve
	State* state ;

    // the total simulated time
    Doub simulated_time ;

	// main routine to call for simulation
	void advanceSimulation ( Doub duration ) ;

	// sub-routine realizing a simulation step
	Doub doFullStep ( Doub h_max ) ;

	// error tolerance parameters
	Doub event_tol ;
	Doub stoch_tol ;
	Doub abs_ode_tol ;
	Doub rel_ode_tol ;

    // the vector storing propensities of stochastic events
    VecDoub propensities ;

    // the random generator for events
    RanGen ran ;

	// ode integration
	void integrate ( Doub h ) ;
	IntegratorDopr5 integrator ;

	// the specific events methods
    #include "stoch_events.hpp"
    #include "det_events.hpp"
} ;


#endif

