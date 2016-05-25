############## ODE CODE ####################

import simulating_helpers as helpers

# derivatives for rules with a source being a shared(unique) agent
def giveDerivativesSourceIsUnique ( model , change ) :
	src_ag_unq = [ ag for ag in model.agents if ag.name == change.agent_source ][0]
	trg_ag_unq = [ ag for ag in model.agents if ag.name == change.agent ][0]
	if not trg_ag_unq.unique : raise Exception ("If source agent unique, target also should be unique.")
	derivs = "\t(*dydt)[ " + "state->" + trg_ag_unq.name.lower () + ".I_" + change.property + " ] += "
	if hasattr ( change , "rate_expression" ) :
		derivs += change.rate_expression
	else :
		for param in change.parameters :
			derivs += "Params::" + param + " "
	derivs += " ;\n"
	# add before a comment line indicating the rule in question
	derivs = "\t// " + change.name + "\n" + derivs	
	return derivs

# derives for rules with a source being a "non-shared" ( for now englobes sphr/non-sphr/etc...)
def giveDerivativesSourceNonUnique ( model , change ) :
	src_ag = [ ag for ag in model.agents if ag.name == change.agent_source ][0]
	trg_ag = [ ag for ag in model.agents if ag.name == change.agent ][0]
	src_ag_istc = src_ag.name.lower ()
	# reference to target agent ?
	# case 1 : target is self
	if trg_ag is src_ag : 
		trg_ag_istc = src_ag_istc
	# case 2 : target is unique
	elif trg_ag.unique : 
		trg_ag_istc = "state->" + trg_ag.name.lower ()
	# case 3 : target is grid
	elif hasattr ( trg_ag , "grid" ) : 
		# if source is not self, then access by X,Y
		trg_ag_istc = "state->give" + trg_ag.name + " ( " + src_ag_istc + "->X , " + src_ag_istc + "->Y )"
	# expression giving the left hand side on derivative :
	lhs_deriv = "\t\t(*dydt)[ " + trg_ag_istc + "->I_" + change.property + " ] += "
	# interesting case : there is an expression for the rhs
	if hasattr ( change , "rate_expression" ) :
		# if source is a grid, we might have "neighbor iteration" in the expression
		if hasattr ( src_ag , "grid" ) :
			# check if we have neighbor appearing in the expression
			if "neighbor." in change.rate_expression :
				# iteration on neighbor
				derivs = "\t\tfor (list <"+src_ag.name+"*>::iterator neighbor = "+src_ag.name.lower()+"->neighbors.begin() ; neighbor != "+src_ag.name.lower()+"->neighbors.end() ; ++neighbor )\n\t\t{\n"
				# translate the ref to neighbor.property in the expression
				rate_expr = change.rate_expression
				for prop in src_ag.properties :
					rate_expr = rate_expr.replace ("neighbor."+prop,"(*y)[ (*neighbor)->I_"+prop+" ]")
				# write the deriv
				derivs += "\t" + lhs_deriv + rate_expr + " ;\n"
				# finish iteration
				derivs += "\t\t}\n"
			else :
				derivs = lhs_deriv + change.rate_expression + " ;\n"
		# if no grid, no iteration (sphere contact...)
		else :
			derivs = lhs_deriv + change.rate_expression + " ;\n"
	# if no expression for rhs, just write params
	else :
		for param in change.parameters :
			derivs += "Params::" + param + " ;\n"
	# add before a comment line indicating the rule in question
	derivs = "\t\t// " + change.name + "\n" + derivs
	return derivs


def giveDerivatives ( model ) :
	derivs = "void IntegratorDopr5::computeDerivatives ( State* state , VecDoub *y , VecDoub *dydt )\n{\n"
	# set the dydt vec to zero
	derivs += "\n\tfor ( Int i = 0 ; i < y->size () ; i++ ) (*dydt)[i] = 0. ;\n"
	# if there is sphere agents, detect contacts
	if model.giveSphereAgents () :
		derivs += "\t// contact detection\n"
		derivs += "\tcontactdetector.detectContacts ( state , y ) ;\n"
	# iterate on agents, potential sources of ODE change
	for source_agent in model.agents :
		# continuous changes for which current agent is the source
		changes = [ chge for chge in model.ode_dynamics if chge.agent_source == source_agent.name ]
		if changes :
			derivs += "\n\t// differential dynamics triggered by agent(s) of type : " + source_agent.name + "\n"
		# case 1 : source is unique agent
		if source_agent.unique :
			for change in changes :
				derivs += giveDerivativesSourceIsUnique ( model , change )
		# case 2 : source is normal (grid or sphere or aspatial)
		else :
			# common loop on this non unique source agent
			ag_iter = helpers.giveAgentIterationCode ( source_agent , "state->" , 1 )
			derivs += ag_iter["line_for"] + "\t{\n"
			derivs += "\t\t" + source_agent.name + "* " + source_agent.name.lower () + " = state->" + source_agent.name.lower () + "s[" + source_agent.name.lower()[0:2] + "] ;\n"
			# line for each change
			for change in changes :
				derivs += giveDerivativesSourceNonUnique ( model , change )
			# end of the common loop on this non unique source agent
			derivs += "\t}\n"
	derivs += "}\n\n"
	return derivs

def writeOdeCode ( model , target_folder ) :
	writer = open ( target_folder + "/derivatives.hpp" , "w" )
	writer.write ( giveDerivatives ( model ) )
