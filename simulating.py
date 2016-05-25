

import modeling
import os.path
import shutil


import simulating_state as state
import simulating_ode as ode
import simulating_det_events as det_events
import simulating_stoch_events as stoch_events
import simulating_contacts as contacts
import simulating_parameters as params

#####################################


def writeAllCode ( model , target_folder ) :

	if not os.path.exists ( target_folder ) : os.mkdir ( target_folder )
	if not os.path.exists ( target_folder + "/source" ) : os.mkdir ( target_folder + "/source" )
	if not os.path.exists ( target_folder + "/build" ) : os.mkdir ( target_folder + "/build" )
	
	# write the model specific code
	state.writeStateCodeHeader ( model , target_folder + "/source" )
	state.writeStateCodeSource ( model , target_folder + "/source" )
	ode.writeOdeCode ( model , target_folder + "/source" )
	det_events.writeDetEventCode ( model , target_folder + "/source" )
	stoch_events.writeStochEventCode ( model , target_folder + "/source" )
	contacts.writeContactsCode ( model , target_folder + "/source" )
	params.writeParameterCode ( model , target_folder + "/source" )

	# copy the generic code
	shutil.copy ( "cpp_code_base/nr3.hpp" , target_folder + "/source" )
	shutil.copy ( "cpp_code_base/random.hpp" , target_folder + "/source" )
	shutil.copy ( "cpp_code_base/integrator.hpp" , target_folder + "/source" )
	shutil.copy ( "cpp_code_base/integrator.cpp" , target_folder + "/source" )
	shutil.copy ( "cpp_code_base/solver.cpp" , target_folder + "/source" )
	shutil.copy ( "cpp_code_base/solver.hpp" , target_folder + "/source" )
	shutil.copy ( "cpp_code_base/compile.sh" , target_folder )


	# shutil.copy ( "cpp_code_base/contactdetector.cpp" , target_folder + "/source" )
	# shutil.copy ( "cpp_code_base/contactdetector.hpp" , target_folder + "/source" )

	# shutil.copy ( "cpp_code_base/main.cpp" , target_folder + "/source" )



