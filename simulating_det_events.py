
############## DET EVENTS CODE #################
def giveTriggerFunc ( det_event ) :
	trig_func = "Bool " + det_event.name + "_trigger ( "
	trig_func += det_event.agent + "* " + det_event.agent.lower () + " )\n{\n"
	trig_func += "\treturn "
	if hasattr ( det_event , "trigger_expression" ) :
		trig_func += det_event.trigger_expression
	else :
		for param in det_event.parameters :
			trig_func += "Params::" + param + " " 
	trig_func += " ;\n}\n\n"
	return trig_func 

def giveTriggerFuncs ( model ) :
	trig_funcs = ""
	for det_event in model.det_events :
		trig_funcs += giveTriggerFunc ( det_event )
	return trig_funcs

def giveDetEventFunc ( det_event ) :
	real_func = ""
	if det_event.kind != "destruction" :
		real_func = "void " + det_event.name + "_realization ( "
		if det_event.kind == "normal" :
			real_func += det_event.agent + "* " + det_event.agent.lower () + " )\n{\n"
		if det_event.kind == "creation" :
			real_func += det_event.agent + "* " + det_event.agent.lower () + " , "
			real_func += det_event.agent + "* new_" + det_event.agent.lower () + " )\n{\n"
		if hasattr ( det_event , "realization_function" ) :
			for expr in det_event.realization_function :
				real_func += "\t" + expr + " ;\n"
		else :
			for param in det_event.parameters :
				real_func += "Params::" + param + " " 
			real_func += "\n"
		real_func += "}\n\n"
	return real_func

def giveDetEventFuncs ( model ) :
	real_funcs = ""
	for det_event in model.det_events :
		real_funcs += giveDetEventFunc ( det_event )
	return real_funcs	

def giveApplyDetEventsDefinition ( model ) :
	definition = "void applyDeterministicEvents ()\n{\n"
	# iterate on agents that are targeted by det events
	for agent in model.agents :
		Ag = agent.name
		ag = Ag.lower ()
		events = [ det_event for det_event in model.det_events if det_event.agent == agent.name ]
		if len (events) != 0 :
			# if unique agent, no death or creation event
			if agent.unique :
				definition += "\tif ( " + event.name + "_trigger ( state->" + ag + " ) ) "
				definition += event.name + "_realization ( state->" + ag + " ) ; \n" 
			else :
				ags = ag + "s"
				ag_vector = "state->" + ag + "s"
				ag_idx = ag[0:2]
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
						definition += "\t\tif ( " + event.name + "_trigger ( " + ag + " ) ) "
						definition += event.name + "_realization ( " + ag + " ) ;\n"
					# if agent destruction, simply delete
					if event.kind == "destruction" :
						definition += "\t\tif ( " + event.name + "_trigger ( " + ag + " ) ) "
						definition += "{ delete " + ag + " ; continue ; }\n"
					# if agent creation
					if event.kind == "creation" :
						definition += "\t\tif ( " + event.name + "_trigger ( " + ag + " ) )\n\t\t{\n"
						definition += "\t\t\t" + agent.name +"* new_" + ag + " = new " + Ag + " ;\n"
						definition += "\t\t\t" + event.name + "_realization ( " + ag + " , new_" + ag + " ) ;\n"
						definition += "\t\t\tkept_" + ag + "s.push_back ( new_" +  ag + " ) ;\n\t\t}\n"
				# if not out of the for loop, the agent is alive
				definition += "\t\tkept_" + ag+ "s.push_back ( " +  ag + " ) ;\n"						
				definition += "\t}\n"
				# add the agent kept
				definition += "\tfor ( list <" + Ag + "*>::iterator it = kept_" + ags + ".begin () ; it != kept_" + ags + ".end () ; ++it ) "
				definition += "state->" + ags + ".push_back ( *it ) ;\n"
	definition += "}\n\n"
	return definition


def writeDetEventCode ( model , target_folder ) :
	writer = open ( target_folder + "/det_events.hpp" , "w" )
	writer.write ( giveTriggerFuncs ( model ) )
	writer.write ( giveDetEventFuncs ( model ) )
	writer.write ( giveApplyDetEventsDefinition ( model ) )
