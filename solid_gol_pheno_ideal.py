#!/usr/bin/python

import modeling
import simulating


### model construction

model = modeling.Model ("solid_gol_pheno_ideal")

# agents
model.addAgent ( name="Cell" , unique=False , properties=[ "R" , "X" , "Y" , "CC_length" , "is_alive" , "density_sensing" , "death_type" ] )

# event: cell division
model.addDeterministicEvent ( name="cell_division" , 
							  agent="Cell" ,
							  kind="creation" , 
							  parameters={"CC_avg":3.0,"CC_std":0.25,"R_birth":0.5**(1./3.)} ,
							  trigger_expression="R * R * R > 2. * R_birth * R_birth * R_birth" ,
							  realization_function=
							  		["new.R = R_birth" ,
							  		 "R = new.R" ,
							  		 "new.CC_length = ran.norm ( CC_avg , CC_std )" ,
							  		 "CC_length = ran.norm ( CC_avg , CC_std )" ,
							  		 "new.X = X" , "new.Y = Y" , "new.is_alive = is_alive" ,
							  		 "new.density_sensing = density_sensing" ] )


# event: cell death (no destruct because dead cells remain there physically !)
model.addDeterministicEvent ( name="cell_death" ,
						   agent="Cell" ,  
						   parameters={ "low_threshold":1.5 , "high_threshold":10. , "comm_distance":5. } ,
						   trigger_expression="density_sensing < low_threshold || density_sensing > high_threshold" ,
						   realization_function=["is_alive = 0."] )


# continuous dynamics : cell growth
model.addContinuousChange ( name="cell_growth" , agent="Cell" , property="R" , rate_expression="is_alive * Params::R_birth * ( pow(2.,1./3.) - 1. ) / CC_length" )






### model simulation (i.e. generation of cpp code for simulation)

simulating.writeAllCode ( model=model , target_folder="solid_gol_pheno_ideal" )



