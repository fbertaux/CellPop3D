
#include "solver.hpp"


Solver::Solver ( State* s ) : 
    state(s) , simulated_time (0.) ,
	event_tol (0.1) , stoch_tol (0.1) ,
	abs_ode_tol (1e-6) , rel_ode_tol (1e-6) ,
    ran ( RanGen(0) ) ,
	integrator ( IntegratorDopr5 (abs_ode_tol,rel_ode_tol) )
	{}

void Solver::advanceSimulation ( Doub duration )
{
	Doub duration_done = 0. ;
	while ( duration_done < duration )
	{
		duration_done += doFullStep ( duration - duration_done ) ;
	}
    simulated_time += duration ;
}

Doub Solver::doFullStep ( Doub h_max )
{
	applyDeterministicEvents () ;
	Doub max_propensity = computePropensities () ;
    Doub h_stoch = stoch_tol / max_propensity ;
	Doub h = MIN ( h_max , MIN ( h_stoch , event_tol ) ) ;
    // cout << "\t\t timestep chosen = " << h << endl ;
	applyStochasticEvents ( h ) ;
	integrate ( h ) ;
    return h ;
}

void Solver::integrate ( Doub h ) { integrator.integrateForDuration ( h , state ) ; }

