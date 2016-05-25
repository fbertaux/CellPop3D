############## STOCH EVENTS CODE ###############

import simulating_helpers as helpers

def givePropensityFunc ( stoch_event ) :
	prop_func = "Doub " + stoch_event.name + "_propensity ( "
	prop_func += stoch_event.agent + "* " + stoch_event.agent.lower () + " )\n{\n"
	prop_func += "\tDoub prop = "
	if hasattr ( stoch_event , "propensity_expression" ) :
		prop_func += stoch_event.propensity_expression
	else :
		for param in stoch_event.parameters :
			prop_func += "Params::" + param + " "
	prop_func += " ;\n\treturn prop ;\n}\n"
	return prop_func

def givePropensityFuncs ( model ) :
	prop_funcs = ""
	for stoch_event in model.stoch_events :
		prop_funcs += givePropensityFunc ( stoch_event ) + "\n"
	return prop_funcs

# TO-DO : if no stoch rule, nothing...
def giveComputePropensitiesDefinition ( model ) :
	# if no stoch events, simply return a value such that h_stoch = event_tol
	if not [stoch_event for stoch_event in model.stoch_events ] :
		return "Doub computePropensities () { return stoch_tol / event_tol ; }\n\n"
	definition = "Doub computePropensities ()\n{\n"
	definition += "\tDoub max_prop = 0. ; Uint i_prop = 0 ;\n"
	propensities_num_expr = []
	ag_propensities_exprs = []
	# iterate on agents that are targeted by stochastic events
	for agent in model.agents :
		events = [ stoch_event for stoch_event in model.stoch_events if stoch_event.agent == agent.name ]
		if len (events) != 0 :
			# computation of number of propensities and of those propensities
			if agent.unique :
				propensities_num_expr.append ( str ( len (events) ) )
				ag_propensities_expr = ""
				for event in events :
					ag_propensities_expr += "\tpropensities [i_prop] = " + event.name + "_propensity ( state->" + agent.name.lower () + " ) ;\n"
					ag_propensities_expr += "\tif ( propensities [i_prop] > max_prop ) max_prop = propensities [i_prop] ;\n\ti_prop ++ ;\n" 
			else :
				propensities_num_expr.append ( "state->" + agent.name.lower () + "s.size () * " + str ( len (events) ) )
				ag_propensities_expr = helpers.giveAgentIterationCode ( agent , "state->" , 1 )["line_for_rev_it"]
				ag_propensities_expr += "\t{\n\t\t" + agent.name + "* " + agent.name.lower () + " = *rit ; \n"
				for event in events :
					ag_propensities_expr += "\t\tpropensities [i_prop] = " + event.name + "_propensity ( " + agent.name.lower () + " ) ;\n"
					ag_propensities_expr += "\t\tif ( propensities [i_prop] > max_prop ) max_prop = propensities [i_prop] ;\n\t\ti_prop ++ ;\n" 
				ag_propensities_expr += "\t}\n"
			ag_propensities_exprs.append ( ag_propensities_expr )
	if not propensities_num_expr :
		propensities_num_expr = ["0"]
	definition += "\tpropensities.resize ( " + " + ".join ( propensities_num_expr ) + " ) ;\n"
	definition += "".join ( ag_propensities_exprs )
	definition += "\treturn max_prop ;\n}\n\n"
	return definition

def giveApplyStochEventsDefinition ( model ) :
	# if no stoch events, empty function
	if not [stoch_event for stoch_event in model.stoch_events ] :
		return "void applyStochasticEvents ( Doub h ) {}\n\n"
	definition = "void applyStochasticEvents ( Doub h )\n{\n\tUint i_prop = 0 ;\n"
	# iterate on agents that are targeted by stochastic events
	for agent in model.agents :
		Ag = agent.name
		ag = Ag.lower ()
		events = [ stoch_event for stoch_event in model.stoch_events if stoch_event.agent == agent.name ]
		if len (events) != 0 :
			# if unique agent, no death or creation event
			if agent.unique :
				definition += "\tif ( ran.doub01 () < h * propensities[i_prop] ) "
				definition += event.name + "_realization ( state->" + agent.name.lower () + " ) ; \n\ti_prop++ ;\n" 
			else :
				ags = Ag.lower () + "s"
				ag_vector = "state->" + ags
				ag_idx = ags[0:2]
				# list to collect the agents to keep
				definition += "\tlist <" + Ag + "*> kept_" + ags + " ;\n"
				# iteration on agent instances
				definition += "\tUint num_" + ags + " = " + ag_vector + ".size () ;\n"
				definition += "\tfor ( Uint " + ag_idx + " = 0 ; " + ag_idx + " < num_" + ags + " ; " + ag_idx + "++ )\n\t{\n"
				definition += "\t\t" + Ag + "* " + ag + " = " + ag_vector + ".back () ; " + ag_vector + ".pop_back () ;\n"
				# iterate on events targeting this agent
				for event in events :
					# if normal event kind (no agent creation or destruction), trigger realization function
					if event.kind == "normal" :
						definition += "\t\tif ( ran.doub01 () < h * propensities[i_prop] ) "
						definition += event.name + "_realization ( " + agent.name.lower () + " ) ;\n"
						definition += "\t\ti_prop++ ;\n"
					# if agent destruction, simply delete
					if event.kind == "destruction" :
						definition += "\t\tif ( ran.doub01 () < h * propensities[i_prop] ) "
						definition += "{ delete " + agent.name.lower () + " ; continue ; }\n"
						definition += "\t\ti_prop++ ;\n"
					# if agent creation
					if event.kind == "creation" :
						definition += "\t\tif ( ran.doub01 () < h * propensities[i_prop] )\n\t\t{\n"
						definition += "\t\t\t" + agent.name +"* new_" + ag + " = new " + Ag + " ;\n"
						definition += "\t\t\t" + event.name + "_realization ( " + ag + " , new_" + ag + " ) ;\n"
						definition += "\t\t\tkept_" + ag + "s.push_back ( new_" +  ag + " ) ;\n\t\t}\n"
						definition += "\t\ti_prop++ ;\n"						
				# if still there, it means no death, the agent should be kept
				definition += "\t\tkept_" + agent.name.lower () + "s.push_back ( " +  agent.name.lower () + " ) ;\n"
				definition += "\t}\n"
				# add the agent kept
				definition += "\tfor ( list <" + Ag + "*>::iterator it = kept_" + ags + ".begin () ; it != kept_" + ags + ".end () ; ++it ) "
				definition += "state->" + ags + ".push_back ( *it ) ;\n"
	definition += "}\n\n"
	return definition

def giveStochEventFunc ( stoch_event ) :
	real_func = ""
	if stoch_event.kind != "destruction" :
		real_func = "void " + stoch_event.name + "_realization ( "
		if stoch_event.kind == "normal" :
			real_func += stoch_event.agent + "* " + stoch_event.agent.lower () + " )\n{\n"
		if stoch_event.kind == "creation" :
			real_func += stoch_event.agent + "* " + stoch_event.agent.lower () + " , "
			real_func += stoch_event.agent + "* new_" + stoch_event.agent.lower () + " )\n{\n"
		if hasattr ( stoch_event , "realization_function" ) :
			for expr in stoch_event.realization_function :
				real_func += "\t" + expr + " ;\n"
		else :
			for param in stoch_event.parameters :
				real_func += "Params::" + param + " " 
			real_func += "\n"
	real_func += "}\n\n"
	return real_func

def giveStochEventFuncs ( model ) :
	real_funcs = ""
	for stoch_event in model.stoch_events :
		real_funcs += giveStochEventFunc ( stoch_event )
	return real_funcs

def writeStochEventCode ( model , target_folder ) :
	writer = open ( target_folder + "/stoch_events.hpp" , "w" )
	writer.write ( givePropensityFuncs ( model ) )
	writer.write ( giveStochEventFuncs ( model ) )
	writer.write ( giveComputePropensitiesDefinition ( model ) )
	writer.write ( giveApplyStochEventsDefinition ( model ) )




