
#include "contactdetector.hpp"


ContactDetector::ContactDetector (Doub lx , Doub ly , Doub lv) :
	nx(lx/lv) , ny(ly/lv) , l_voxel(lv) , 
	length_x (lx) ,
	length_y (ly)
{
	cout << "ContactDetection grid : nx = " << nx << " _ ny = " << ny << endl ;
	x_first_vox = - l_voxel * (nx/2. - 0.5) ;
	y_first_vox = - l_voxel * (ny/2. - 0.5) ;
	cout << "ContactDetection grid : x_first_vox = " << x_first_vox << " _ y_first_vox = " << y_first_vox << endl ;
	// voxels construction
	for ( Int I = 0 ; I < nx*ny ; I++ ) voxels.push_back ( new ContVoxel ) ;
	for ( Int ix = 0 ; ix < nx ; ix++ )
	{
		for ( Int iy = 0 ; iy < ny ; iy++ )
		{
			Int I = ny*ix + iy ;
			// left-right relationship
			if ( ix+1 < nx )
			{
				voxels[I]->neighbors.push_back ( voxels[ny*(ix+1) + iy] ) ;
				voxels[ny*(ix+1) + iy]->neighbors.push_back ( voxels[I] ) ;
			}
			// down-up relationship
			if ( iy+1 < ny )
			{
				voxels[I]->neighbors.push_back ( voxels[ny*ix + iy+1] ) ;
				voxels[ny*ix + iy+1]->neighbors.push_back ( voxels[I] ) ;
			}
			// diag-left relationship
			if ( ix+1 < nx && iy+1 < ny )
			{
				voxels[I]->neighbors.push_back ( voxels[ny*(ix+1) + iy+1] ) ;
				voxels[ny*(ix+1) + iy+1]->neighbors.push_back ( voxels[I] ) ;				
			}
			// diag-right relationship
			if ( ix > 0 && iy+1 < ny )
			{
				voxels[I]->neighbors.push_back ( voxels[ny*(ix-1) + iy+1] ) ;
				voxels[ny*(ix-1) + iy+1]->neighbors.push_back ( voxels[I] ) ;				
			}
		}
	}
}

void ContactDetector::resetCellMapping ()
{
	for ( Uint v = 0 ; v < voxels.size () ; v++ ) voxels[v]->I_cells.clear () ;
}


void ContactDetector::locateCellsInVoxels (  State* state , VecDoub* y )
{
	// Doub X , Y ;
	Int ix , iy ;
	for ( Uint c = 0 ; c < state->cells.size () ; c++ )
	{
		// de-activate grid detection for testing
		// voxels[0]->I_cells.push_back (c) ; continue ;

		// cout << "\t X = " << (*y)[ state->cells[c]->I_X ] << " _ Y = " << (*y)[ state->cells[c]->I_Y ] << endl ; 
		ix = ( (*y)[ state->cells[c]->I_X ] - x_first_vox ) / l_voxel ;
		iy = ( (*y)[ state->cells[c]->I_Y ] - y_first_vox ) / l_voxel ;
		// if outside, put in closest boundary box
		if ( ix > nx-1 ) { ix = nx-1 ; } else if ( ix < 0 ) { ix = 0 ; }
		if ( iy > ny-1 ) { iy = ny-1 ; } else if ( iy < 0 ) { iy = 0 ; }		
		// cout << "\t ix = " << ix << " - iy = " << iy << endl ;
		// cout << "\t I = " << ny*ix + iy  << endl << endl ;
		voxels[ ny*ix + iy ]->I_cells.push_back (c) ;
	}
}


void ContactDetector::computeVelocities (  State* state , VecDoub* y , VecDoub* dydt )
{
	// iterate on voxels
	for ( Uint v = 0 ; v < voxels.size () ; v++ )
	{
		// iterate on cells in this voxels
		for (list <Int>::iterator it = voxels[v]->I_cells.begin() ; it != voxels[v]->I_cells.end() ; ++it )
		{
			// iterate on cells in same voxel
			for (list <Int>::iterator itsv = voxels[v]->I_cells.begin() ; itsv != voxels[v]->I_cells.end() ; ++itsv )
			{
				computeVelocity ( state , y , dydt , *it , *itsv ) ;
			}
			// iterate neighbor voxels
			for (list <ContVoxel*>::iterator itnv = voxels[v]->neighbors.begin() ; itnv!= voxels[v]->neighbors.end() ; ++itnv )
			{
				// iterate on cells in neighbor voxel
				for (list <Int>::iterator itnc = (*itnv)->I_cells.begin() ; itnc != (*itnv)->I_cells.end() ; ++itnc )
				{
					computeVelocity ( state , y , dydt , *it , *itnc ) ;
				}
			}
		}
	}
}


Doub
extendedHertzForce ( const Doub &A_rad , const Doub &B_rad , const Doub &AB_dist )
{
	// FB , 2012/11/13 : GOOD
	if (AB_dist>A_rad+B_rad) {return 0;} // added 2012/11/18
    return 4. / 3. / ( (1. - Params::poisson_ratio * Params::poisson_ratio)/Params::young_modulus + (1. - Params::poisson_ratio * Params::poisson_ratio)/Params::young_modulus) * sqrt( A_rad * B_rad / (A_rad + B_rad) ) * pow(A_rad + B_rad - AB_dist, 1.5 ) - M_PI * Params::adhesivity * A_rad * B_rad / (A_rad + B_rad) ;
}

void ContactDetector::computeVelocity ( State* state , VecDoub* y , VecDoub* dydt , Int c1 , Int c2 )
{
	if (c1 >= c2) return ;
	// cout << "\t c1 = " << c1 <<  " __ c2 = " << c2 << endl ;
	Doub contdist2 = SQR ( (*y)[ state->cells[c1]->I_R ] + (*y)[ state->cells[c2]->I_R ] ) ;
	// cout << "\t contact dist ^ 2 = " << contdist2 << endl ;
	Doub dx = (*y)[ state->cells[c1]->I_X ] - (*y)[ state->cells[c2]->I_X ] ;
	// periodicity: modulation of dx
	// modulationX (dx) ;
	// cout << "\t dx = " << dx << endl ;
    if (dx*dx>contdist2) { return; }

    Doub dy = (*y)[ state->cells[c1]->I_Y ] - (*y)[ state->cells[c2]->I_Y ] ;
	// periodicity: modulation of dy
	// modulationY (dy) ;
	// cout << "\t dy = " << dy << endl ;
    if (dy*dy>contdist2) { return; }

    if ( dx*dx+dy*dy  < contdist2 )
    {
    	Doub dist = sqrt ( dx*dx+dy*dy ) ;
    	Doub force = extendedHertzForce ( (*y)[ state->cells[c1]->I_R ] , (*y)[ state->cells[c2]->I_R ] , dist ) ;
    	Doub fric1 = Params::friction * 4. * M_PI * pow ( (*y)[ state->cells[c1]->I_R ] , 2.) ;
    	Doub fric2 = Params::friction * 4. * M_PI * pow ( (*y)[ state->cells[c2]->I_R ] , 2.) ;
    	// cout << "\t\t contact ! dist = " << dist << " , force = " << force << endl ;

    	(*dydt)[ state->cells[c1]->I_X ] += dx / dist * force / fric1 ;
    	(*dydt)[ state->cells[c1]->I_Y ] += dy / dist * force / fric1 ;

    	(*dydt)[ state->cells[c2]->I_X ] -= dx / dist * force / fric2 ;
    	(*dydt)[ state->cells[c2]->I_Y ] -= dy / dist * force / fric2 ;
    	// cout << "\t\t\t dx/dt = " << dx / dist * force / fric1 << endl ; 
    	// cout << "\t\t\t dy/dt = " << dy / dist * force / fric1 << endl ; 
        }
    // cout << endl ;
}





