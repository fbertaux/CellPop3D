############# PARAMETERS ###################

def giveParametersDeclaration ( model ) :
	param_names = []
	params_def = "struct Params\n{\n"
	event_types = [ "det_events" , "stoch_events" , "ode_dynamics" ]
	for event_type in event_types :
		for event in eval ( "model." + event_type ) :
			for param in event.parameters.keys () :
				if param not in param_names :
					param_names.append ( param )
					params_def += "\tstatic Doub " + param + " ; \n"
	# add params for sphere contact detection
	# if model.giveSphereAgents () : for now no clear separation aspatial / spatial
	params_def += "\n"
	params_def += "\tstatic Doub length_x ;\n"
	params_def += "\tstatic Doub length_y ;\n"
	params_def += "\tstatic Doub R_max ;\n"
	params_def += "} ;\n\n"
	return params_def

def giveParametersDefinition ( model ) :
	param_names = []
	params_decl = ""
	event_types = [ "det_events" , "stoch_events" , "ode_dynamics" ]
	for event_type in event_types :
		for event in eval ( "model." + event_type ) :
			for param in event.parameters.keys () :
				if param not in param_names :
					param_names.append ( param )
					params_decl += "Doub Params::" + param + " = " + str ( event.parameters[param] ) + " ; \n"
	# params for sphere contact detc (TO DO: non default, depend on the "resolved as sphere")
	params_decl += "\n"
	lx = [ ag.length_x for ag in model.giveSphereAgents () ]
	ly = [ ag.length_y for ag in model.giveSphereAgents () ]
	r_max = [ ag.max_radius for ag in model.giveSphereAgents () ]
	# do latter case 3Ds
	if lx :
		params_decl += "Doub Params::length_x = " + str(max(lx)) + " ;\n"
		params_decl += "Doub Params::length_y = " + str(max(ly)) + " ;\n"
		params_decl += "Doub Params::R_max = " + str(max(r_max)) + " ;\n"
	params_decl += "\n\n"					
	return params_decl

def giveUnitsDeclaration ( model ) :
	units_decl = "\n#ifndef PARAMETERS_HPP\n#define PARAMETERS_HPP\n\n"
	units_decl += "\n#include \"nr3.hpp\"\n\n"
	units_decl += "struct Units\n{\n"
	units_decl += "\tstatic Doub DistanceUnit ; // in meters (SI)\n"
	units_decl += "\tstatic Doub EnergyUnit ; // in J = N.m (SI)\n"
	units_decl += "\tstatic Doub TimeUnit ; // in seconds (SI)\n"
	units_decl += "\tstatic Doub ForceUnit ; // in N (SI)\n"
	units_decl += "\tstatic Doub PressureUnit ; // in Pa (SI)\n"
	units_decl += "} ;\n\n"
	return units_decl

def giveUnitsDefinition ( model ) :
	# TO DO : distance units = R_max ?
	units_def = "\n#include \"parameters.hpp\"\n\n"
	units_def += "Doub Units::DistanceUnit = 5.e-6 ; // in meters (SI)\n"
	units_def += "Doub Units::EnergyUnit = 1e-16 ; // in J = N.m (SI)\n"
	units_def += "Doub Units::TimeUnit =  3600. ; // in seconds (SI)\n"
	units_def += "Doub Units::ForceUnit = Units::EnergyUnit/Units::DistanceUnit ; // in N (SI)\n"
	units_def += "Doub Units::PressureUnit = Units::ForceUnit/(Units::DistanceUnit*Units::DistanceUnit) ; // in Pa (SI)\n"
	units_def += "\n\n"
	return units_def


def giveSimParametersDeclaration ( model ) :
	sim_params_decl = "struct SimParams\n{\n"
	sim_params_decl += "\tstatic Doub total_duration ;\n"
	sim_params_decl += "\tstatic Doub output_period ;\n"
	sim_params_decl += "} ;\n\n"
	return sim_params_decl

def giveSimParametersDefinition ( model ) :
	sims_params_def = "Doub SimParams::total_duration = 0. ;\n"
	sims_params_def += "Doub SimParams::output_period = 0. ;\n\n"
	return sims_params_def

def writeParameterCode ( model , target_folder ) :
	writer = open ( target_folder + "/parameters.hpp" , "w" )
	writer.write ( giveUnitsDeclaration ( model ) )
	writer.write ( giveParametersDeclaration ( model ) )
	writer.write ( giveSimParametersDeclaration ( model) )	
	writer.write ( "#endif\n" )
	writer = open ( target_folder + "/parameters.cpp" , "w" )
	writer.write ( giveUnitsDefinition ( model ) )
	writer.write ( giveParametersDefinition ( model ) )	
	writer.write ( giveSimParametersDefinition ( model) )


