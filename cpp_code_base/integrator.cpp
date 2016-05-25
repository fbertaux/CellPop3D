
#include "integrator.hpp"


// parameters
static const Doub a21=0.2,a31=3.0/40.0,
a32=9.0/40.0,a41=44.0/45.0,a42=-56.0/15.0,a43=32.0/9.0,a51=19372.0/6561.0,
a52=-25360.0/2187.0,a53=64448.0/6561.0,a54=-212.0/729.0,a61=9017.0/3168.0,
a62=-355.0/33.0,a63=46732.0/5247.0,a64=49.0/176.0,a65=-5103.0/18656.0,
a71=35.0/384.0,a73=500.0/1113.0,a74=125.0/192.0,a75=-2187.0/6784.0,
a76=11.0/84.0,e1=71.0/57600.0,e3=-71.0/16695.0,e4=71.0/1920.0,
e5=-17253.0/339200.0,e6=22.0/525.0,e7=-1.0/40.0 ;
static const Doub beta=0.0,alpha=0.2-beta*0.75,safe=0.9,minscale=0.2,
maxscale=10.0 ;


IntegratorDopr5::IntegratorDopr5 ( Doub abs_tol , Doub rel_tol ) : 
    atol(abs_tol) , rtol(rel_tol) , EPS(numeric_limits<Doub>::epsilon()) ,
    error_previous_step(1e-4) , reject_previous_step(false) ,
    contactdetector(Params::length_x,Params::length_y,2.*Params::R_max)
    {}

void IntegratorDopr5::reInitialize (Int sys_size)
{
    y_init_step.resize (sys_size) ;
    dydt_init_step.resize (sys_size) ;
    y_substep.resize (sys_size) ;
    y_end_step.resize (sys_size) ;
    dydt_end_step.resize (sys_size) ;
    dydt_substep_2.resize (sys_size) ;
    dydt_substep_3.resize (sys_size) ;
    dydt_substep_4.resize (sys_size) ;
    dydt_substep_5.resize (sys_size) ;
    dydt_substep_6.resize (sys_size) ;
    y_error_step.resize (sys_size) ;
}

void IntegratorDopr5::integrateForDuration ( Doub duration , State* state )
{
    Doub duration_done = 0. ;
    state->countVarsAndAssignIndexes () ;
    reInitialize ( state->num_vars ) ;
    state->exportVarsInVector ( &y_init_step ) ;
    computeDerivatives ( state , &y_init_step , &dydt_init_step ) ;
    Doub h = duration ;
    while ( duration_done < duration )
    {
        if ( h > duration - duration_done )
            h = duration - duration_done ;
        doStep ( h , state ) ;
        duration_done += h_step ;
        h = h_next_step ;
    }
    // cout << "\t\t age first cell before importing = " << y_init_step[ state->cells[0]->I_age ] << endl ;
    state->importVarsFromVector ( &y_init_step ) ;
}

void IntegratorDopr5::doStep ( Doub h_init_step , State* state )
{
	h_step = h_init_step ;
	bool success = false ;
    int num_trials = 0 ;
	while ( ! success )
	{
		success = doStepTrial (state) ;
        num_trials ++ ;
	}
    // copy the result in the 'initial' vectors
    for ( Int i = 0 ; i < y_init_step.size () ; i++ )
    {
        y_init_step[i] = y_end_step[i] ;
        dydt_init_step[i] = dydt_end_step[i] ;
    }
    // cout << "\t\t num trials = " << num_trials << endl ;
}

bool IntegratorDopr5::doStepTrial ( State* state )
{
	computeStep (state) ;
    computeError () ;
	return checkSuccess () ;
}

void IntegratorDopr5::computeStep ( State* state )
{
    // cout << "\t\t\t computing step. h = " << h_step << endl ;
	Int i ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_substep[i] = y_init_step[i] + h_step * a21 * dydt_init_step[i] ;
   // cout << "\t\t\t y_substep_1 [ first cell age ] = " << y_substep[state->cells[0]->I_age] << endl ;

	computeDerivatives ( state , &y_substep , &dydt_substep_2 ) ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_substep[i] = y_init_step[i] + h_step * ( a31 * dydt_init_step[i] + a32 * dydt_substep_2[i] ) ;

	computeDerivatives ( state , &y_substep , &dydt_substep_3 ) ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_substep[i] = y_init_step[i] + h_step * ( a41 * dydt_init_step[i] + a42 * dydt_substep_2[i] + a43 * dydt_substep_3[i] ) ;

	computeDerivatives ( state , &y_substep , &dydt_substep_4 ) ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_substep[i] = y_init_step[i] + h_step * ( a51 * dydt_init_step[i] + a52 * dydt_substep_2[i] + a53 * dydt_substep_3[i] + a54 * dydt_substep_4[i] ) ;

	computeDerivatives ( state , &y_substep , &dydt_substep_5 ) ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_substep[i] = y_init_step[i] + h_step * ( a61 * dydt_init_step[i] + a62 * dydt_substep_2[i] + a63 * dydt_substep_3[i] + a64 * dydt_substep_4[i] + a65 * dydt_substep_5[i] ) ;
   // cout << "\t\t\t y_substep_5 [ first cell age ] = " << y_substep[state->cells[0]->I_age] << endl ;

	computeDerivatives ( state , &y_substep , &dydt_substep_6 ) ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_end_step[i] = y_init_step[i] + h_step * ( a71 * dydt_init_step[i] + a73 * dydt_substep_3[i] + a74 * dydt_substep_4[i] + a75 * dydt_substep_5[i] + a76 * dydt_substep_6[i] ) ;
   // cout << "\t\t\t y_end_step [ first cell age ] = " << y_end_step[state->cells[0]->I_age] << endl ;

	computeDerivatives ( state , &y_end_step , &dydt_end_step ) ;
    for ( i = 0 ; i < y_init_step.size () ; i++ )
		y_error_step[i] = h_step * ( e1 * dydt_init_step[i] + e3 * dydt_substep_3[i] + e4 * dydt_substep_4[i] + e5 * dydt_substep_5[i] + e6 * dydt_substep_6[i] + e7 * dydt_end_step[i] ) ;

   // cout << "\t\t\t computed step. y_step [ first cell age ] = " << y_end_step[state->cells[0]->I_age] << endl ;
}

void IntegratorDopr5::computeError ()
{
    Doub err = 0.0 , sk ;
    for ( Int i=0 ; i < y_init_step.size () ; i++ )
    {
        sk = atol + rtol * MAX ( abs(y_init_step[i]) , abs(y_end_step[i]) ) ;
//        cout << "\t\t sk = " << sk << endl ;
        err += SQR (y_error_step[i]/sk) ;
    }
    error_step = sqrt (err/y_init_step.size()) ;
    // cout << "\t computing error = " << error_step << endl ;
}

bool IntegratorDopr5::checkSuccess ()
{
    Doub scale;
    if (error_step <= 1.0)
    {
        if (error_step == 0.0)
            scale=maxscale;
        else
        {
            scale=safe*pow(error_step,-alpha)*pow(error_previous_step,beta) ;
            if (scale<minscale) scale=minscale ;
            if (scale>maxscale) scale=maxscale ;
        }
        if (reject_previous_step)
            h_next_step=h_step*MIN(scale,1.0) ;
        else
            h_next_step=h_step*scale;
        error_previous_step=MAX(error_step,1.0e-4) ;
        reject_previous_step = false ;
        return true ;
    }
    else
    {
        scale=MAX(safe*pow(error_step,-alpha),minscale) ;
        h_step *= scale ;
        reject_previous_step = true ;
        return false ;
    }
}

#include "derivatives.hpp"


