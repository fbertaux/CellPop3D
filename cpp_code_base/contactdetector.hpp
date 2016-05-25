
#include "state.hpp"
#include "parameters.hpp"

struct ContVoxel
{
	list <ContVoxel*> neighbors ;
	list <Int> I_cells ;
} ;


struct ContactDetector
{
	ContactDetector (Doub lx , Doub ly , Doub lv) ;
	const Int nx , ny ;
	const Doub l_voxel ;
	const Doub length_x , length_y ;
	Doub x_first_vox , y_first_vox ;
	vector <ContVoxel*> voxels ;

	void locateCellsInVoxels ( State* state , VecDoub* y ) ;
	void resetCellMapping () ;
	void computeVelocities ( State* state , VecDoub *y , VecDoub *dydt ) ;
	void computeVelocity ( State* state , VecDoub* y , VecDoub* dydt , Int c1 , Int c2 ) ;
} ;