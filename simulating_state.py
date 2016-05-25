############## STATE CODE ###############

import simulating_helpers as helpers


def giveAgentDeclarations ( model ) :
	agents_decl = ""
	for agent in model.agents :
		agents_decl += "struct " + agent.name + "\n{\n"
		agents_decl += "\t" + agent.name + " () ;\n"
		for prop in agent.properties :
			agents_decl += "\tDoub " + prop + " ; Uint I_" + prop + " ;\n"
		# if is a grid, add a list of neighbors
		if hasattr (agent,"grid") :
			agents_decl += "\tlist <" + agent.name +"*> neighbors ; \n" 
		agents_decl += "} ;\n\n"
	return agents_decl

def giveStateDeclaration ( model ) :
	state_decl = "struct State\n{\n"
	state_decl += "\tState () ;\n"
	for agent in model.agents :
		if agent.unique :
			state_decl += "\t" + agent.name + " " + agent.name.lower () + " ;\n"
		else :
			state_decl += "\tvector <" + agent.name + "*> " + agent.name.lower () + "s ;\n"
			# if is a grid, add the grid info
			if hasattr (agent,"grid") :
				state_decl += "\tDoub " + agent.name.lower() + "_l ;\n"
				state_decl += "\tDoub " + agent.name.lower() + "_lx ;\n"
				state_decl += "\tDoub " + agent.name.lower() + "_ly ;\n"
				state_decl += "\tInt " + agent.name.lower() + "_nx ;\n"
				state_decl += "\tInt " + agent.name.lower() + "_ny ;\n"
				state_decl += "\tDoub " + agent.name.lower() + "_xf , " + agent.name.lower() + "_yf ;\n"
				state_decl += "\t" + agent.name + "* give" + agent.name + " ( Doub x , Doub y ) ;\n"
	state_decl += "\n\tUint num_vars ;\n"
	state_decl += "\tvoid countVarsAndAssignIndexes () ;\n"
	state_decl += "\tvoid exportVarsInVector ( VecDoub* vars_vec ) ;\n"
	state_decl += "\tvoid importVarsFromVector ( VecDoub* vars_vec ) ;\n"
	state_decl += "} ;\n\n"
	return state_decl

def writeStateCodeHeader ( model , target_folder ) :
	writer = open ( target_folder + "/state.hpp" , "w" )
	writer.write ( "#ifndef STATE_HPP\n#define STATE_HPP\n\n")
	writer.write ( "\n#include \"nr3.hpp\"\n\n" )
	writer.write ( giveAgentDeclarations ( model ) )
	writer.write ( giveStateDeclaration ( model ) )
	writer.write ( "\n#endif\n\n")

def giveStateDefinitions ( model ) :
	## constructor
	definitions = "State::State ()\n{\n"
	# iterate on agents in search for grid
	for agent in model.agents :
		if hasattr (agent,"grid") :
			# infos about the grid
			definitions += "\t" + agent.name.lower() + "_l = " + str(agent.step_length)  + " ;\n"
			definitions += "\t" + agent.name.lower() + "_lx = " + str(agent.length_x)  + " ;\n"
			definitions += "\t" + agent.name.lower() + "_ly = " + str(agent.length_y)  + " ;\n"
			definitions += "\t" + agent.name.lower() + "_nx = " + agent.name.lower() + "_lx / " +   agent.name.lower() + "_l ;\n"
			definitions += "\t" + agent.name.lower() + "_ny = " + agent.name.lower() + "_ly / " +   agent.name.lower() + "_l ;\n"
			definitions += "\t" + agent.name.lower() + "_xf = - " + agent.name.lower() + "_l * ( " + agent.name.lower() + "_nx/2. - 0.5) ;\n"
			definitions += "\t" + agent.name.lower() + "_yf = - " + agent.name.lower() + "_l * ( " + agent.name.lower() + "_ny/2. - 0.5) ;\n"
			# creation of voxels
			definitions += "\tfor ( Int I = 0 ; I < " + agent.name.lower() + "_nx*" + agent.name.lower() + "_ny ; I++ ) "
			definitions += agent.name.lower() + "s.push_back ( new " + agent.name + " ) ;\n"
			# creation of topological relationships between voxels
			definitions += "\tfor (Int ix = 0 ; ix < " + agent.name.lower() + "_nx ; ix++ )\n\t{\n"
			definitions += "\t\tfor (Int iy = 0 ; iy < " + agent.name.lower() + "_ny ; iy++ )\n\t\t{\n"
			definitions += "\t\t\tInt I = " + agent.name.lower() + "_ny*ix + iy ;\n"
			definitions += "\t\t\tif ( ix+1 < " + agent.name.lower() + "_nx )\n\t\t\t{\n"
			definitions += "\t\t\t\t" + agent.name.lower() + "s[I]->neighbors.push_back ( " + agent.name.lower()
			definitions += "s[" + agent.name.lower() + "_ny*(ix+1) + iy] ) ;\n"
			definitions += "\t\t\t\t" + agent.name.lower() + "s[" + agent.name.lower() + "_ny*(ix+1) + iy]->neighbors.push_back ( "
			definitions += agent.name.lower() + "s[I] ) ;\n"
			definitions += "\t\t\t}\n"
			definitions += "\t\t\tif ( iy+1 < " + agent.name.lower() + "_ny )\n\t\t\t{\n"
			definitions += "\t\t\t\t" + agent.name.lower() + "s[I]->neighbors.push_back ( " + agent.name.lower()
			definitions += "s[" + agent.name.lower() + "_ny*ix + iy+1] ) ;\n"
			definitions += "\t\t\t\t" + agent.name.lower() + "s[" + agent.name.lower() + "_ny*ix + iy+1]->neighbors.push_back ( "
			definitions += agent.name.lower() + "s[I] ) ;\n"
			definitions += "\t\t\t}\n"
			definitions += "\t\t}\n\t}\n"			
	definitions += "}\n\n"
	## methods for grid location
	for agent in model.agents :
		if hasattr (agent,"grid") :
			definitions += agent.name + "* State::give" + agent.name + " ( Doub x , Doub y )\n{\n"
			definitions += "\tInt ix , iy ;\n"
			definitions += "\tix = ( x - " + agent.name.lower() + "_xf ) / " + agent.name.lower() + "_l ;\n"
			definitions += "\tiy = ( y - " + agent.name.lower() + "_yf ) / " + agent.name.lower() + "_l ;\n"
			definitions += "\tif ( ix > " + agent.name.lower() + "_nx-1 ) { cout << \"asked location out of grid.\" << endl ; exit(1) ; } "
			definitions += "else if ( ix < 0 ) { cout << \"asked location out of grid.\" << endl ; exit(1) ; }\n"
			definitions += "\tif ( iy > " + agent.name.lower() + "_ny-1 ) { cout << \"asked location out of grid.\" << endl ; exit(1) ; } "
			definitions += "else if ( iy < 0 ) { cout << \"asked location out of grid.\" << endl ; exit(1) ; }\n"
			definitions += "\treturn " + agent.name.lower() + "s[" + agent.name.lower() + "_ny*ix+iy] ;\n"
			definitions += "}\n\n"
	return definitions

def giveAgentDefinitions ( model ) :
	definitions = ""
	for agent in model.agents :
		definitions += agent.name + "::" + agent.name + " ()\n{\n"
		for prop in agent.properties :
			definitions += "\t" + prop + " = 0. ;\n"
		definitions += "}\n\n"
	return definitions

def giveCountVarsDefinition ( model ) :
	definition = "void State::countVarsAndAssignIndexes ()\n{\n\tnum_vars = 0 ;\n"
	for agent in model.agents :
		if agent.unique :
			for prop in agent.properties :
				definition += "\t" + agent.name.lower () + ".I_" + prop + " = num_vars ; num_vars++ ;\n"
		else :
			ag_iter = helpers.giveAgentIterationCode ( agent=agent , indent=1 )
			definition += ag_iter["line_for"] + "\t{\n"
			for prop in agent.properties :
				definition += "\t\t" + ag_iter["ag_vector"] + "[" + ag_iter["ag_idx"] + "]->I_" + prop + " = num_vars ; num_vars++ ;\n"
			definition += "\t}\n"
	definition += "}\n\n"
	return definition

def giveExportVarsDefinition ( model ) :
	definition = "void State::exportVarsInVector ( VecDoub* vars_vec )\n{\n"
	for agent in model.agents :
		if agent.unique :
			for prop in agent.properties :
				definition += "\t(*vars_vec) [ " + agent.name.lower () + ".I_" + prop + " ] = "
				definition += agent.name.lower () + "." + prop + " ;\n"
		else :
			ag_iter = helpers.giveAgentIterationCode ( agent=agent , indent=1 )
			definition += ag_iter["line_for"] + "\t{\n"
			for prop in agent.properties :
				definition += "\t\t(*vars_vec) [ " + ag_iter["ag_vector"] + "[" + ag_iter["ag_idx"] + "]->I_" + prop + " ] = "
				definition += ag_iter["ag_vector"] + "[" + ag_iter["ag_idx"] + "]->" + prop + " ;\n"
			definition += "\t}\n"
	definition += "}\n\n"
	return definition

def giveImportVarsDefinition ( model ) :
	definition = "void State::importVarsFromVector ( VecDoub* vars_vec )\n{\n"
	for agent in model.agents :
		if agent.unique :
			for prop in agent.properties :
				definition += "\t" + agent.name.lower () + "." + prop + " = (*vars_vec) [ "
				definition += agent.name.lower () + ".I_" + prop + " ] ;\n"
		else :
			ag_iter = helpers.giveAgentIterationCode ( agent=agent , indent=1 )
			definition += ag_iter["line_for"] + "\t{\n"
			for prop in agent.properties :
				definition += "\t\t" + ag_iter["ag_vector"] + "[" + ag_iter["ag_idx"] + "]->" + prop + " = (*vars_vec) [ "
				definition += ag_iter["ag_vector"] + "[" + ag_iter["ag_idx"] + "]->I_" + prop + " ] ;\n"
			definition += "\t}\n"
	definition += "}\n\n"
	return definition

def writeStateCodeSource ( model , target_folder ) :
	writer = open ( target_folder + "/state.cpp" , "w" )
	writer.write ( "\n#include \"state.hpp\"\n\n" )
	writer.write ( giveStateDefinitions (model) )
	writer.write ( giveAgentDefinitions (model) )
	writer.write ( giveCountVarsDefinition ( model ) )
	writer.write ( giveExportVarsDefinition ( model ) )
	writer.write ( giveImportVarsDefinition ( model ) )
