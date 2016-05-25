
#include "solver.hpp"

int main ()
{

	// create initial state
	State state ;

	// create solver
	Solver solver ( &state ) ;

	// simulate
    Doub total_duration = 100 ;
    Doub interval = 10. ;
    for (Uint i = 0 ; i < total_duration / interval ; i++ )
	{
        solver.advanceSimulation ( interval ) ;
        cout << endl << "time = " << solver.simulated_time << endl ;
	}

	return 0 ;
}
