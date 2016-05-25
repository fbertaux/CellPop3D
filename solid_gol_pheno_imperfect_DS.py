#!/usr/bin/python

import modeling
import simulating


### model construction

model = modeling.Model ("solid_gol_pheno_imperfectDS")

# agents
model.addAgent ( name="Cell" , unique=False , properties= [ "CC_length" , "is_alive" , "death_type" ] )
model.addAgent ( name="Medium" , unique=True , properties= [ "messenger" ] )

# resolve spatially
model.resolveAsGrid ( agent="Medium" , two_dim=True , length_x=50.0 , length_y=50.0 , step_length=0.5 )

# currently being implemented
model.resolveAsSphere ( agent="Cell" , two_dim=True , length_x=50.0 , length_y=50.0 , max_radius=1.0 )


# continuous dynamics : messenger degradation 
model.addContinuousChange ( name="mess_degrad" , agent="Medium" , property="messenger" , 
							parameters={ "mess_degrad_rate":0. } ,
							rate_expression="- mess_degrad_rate * messenger" )

# continuous dynamics : messenger production
model.addContinuousChange ( name="mess_prod" , agent_source="Cell" , agent="Medium" , property="messenger" ,
							parameters={ "mess_prod_rate_per_cell":1.0 } ,
							rate_expression="mess_prod_rate_per_cell" )

# continuous dynamics : messenger diffusion
model.addContinuousChange ( name="mess_diff" , agent="Medium" , property="messenger" ,
							parameters={ "mess_diff_coeff":10. } ,
							rate_expression="mess_diff_coeff * ( neighbor.messenger - messenger ) / state->medium_l / state->medium_l" ) 

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
							  		 "new.X = X" , "new.Y = Y" , "new.is_alive = is_alive" ] )


# event: cell death (no destruct because dead cells remain there physically !)
model.addDeterministicEvent ( name="cell_death_LT" ,
						   agent="Cell" ,  
						   parameters={ "low_threshold":1.5 } ,
						   trigger_expression="state->giveMedium ( X , Y )->messenger < low_threshold" ,
						   realization_function=["is_alive = 0.","death_type = -1."] )
model.addDeterministicEvent ( name="cell_death_HT" ,
						   agent="Cell" ,  
						   parameters={ "high_threshold":10. } ,
						   trigger_expression="state->giveMedium ( X , Y )->messenger > high_threshold" ,
						   realization_function=["is_alive = 0.","death_type = 1."] )



# continuous dynamics : cell growth
model.addContinuousChange ( name="cell_growth" , agent="Cell" , property="R" , rate_expression="is_alive * Params::R_birth * ( pow(2.,1./3.) - 1. ) / CC_length" )



# AND CELL PHYSICS ?
# model.addContinuousChange ( name="cell_cell_physics" , agent="Cell" , property="X" , rate_expression="" )


### model simulation (i.e. generation of cpp code for simulation)

simulating.writeAllCode ( model=model , target_folder="solid_gol_pheno_imperfectDS" )



