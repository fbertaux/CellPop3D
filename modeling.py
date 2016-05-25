

#### A model is first made of different types of agents ####
class Agent (object) :	
	def __init__ ( self , name , properties=[] , unique=False ) :
		self.name = name
		self.properties = properties
		self.unique = unique
	def display (self) :
		print "agent type name = %s , unique = %i" % ( self.name , self.unique )
		for prop in self.properties :
			print "\tproperty : %s" % prop

#### A model dynamics can have stochastic events, triggered via propensity
class StochasticEvent (object) :
	def __init__ ( self , name , agent , kind="normal" , parameters={} ) :
		self.name = name 
		self.agent = agent
		self.kind = kind
		self.parameters = parameters
	def display (self) :
		print "stochastic event \"%s\" acting on agent %s , of kind \"%s\", with parameters %s" % ( self.name , self.agent , self.kind , str (self.parameters) )

#### A model dynamics can have deterministic events, triggered by a boolean function of the state
class DeterministicEvent (object) :
	def __init__ ( self , name , agent , kind="normal" , parameters={} ) :
		self.name = name 
		self.agent = agent
		self.kind = kind
		self.parameters = parameters
	def display (self) :
		print "deterministic event \"%s\" acting on agent %s , of kind \"%s\", with parameters %s" % ( self.name , self.name , self.kind , str (self.parameters) )

#### Model state can evolve continuously via ODE rules (or computation)
class ContinuousChange (object) :
	def __init__ ( self , name , agent , property , agent_source=None , parameters={} ) :
		self.name = name
		self.agent = agent
		self.property = property
		self.agent_source = agent_source
		if agent_source == None :
			self.agent_source = agent
		self.parameters = parameters

#### A model is a set of agent types and dynamics affecting them ####
class Model (object) :

	def __init__ ( self , name="Model" ) :
		self.name = name 
		self.agents = []
		self.stoch_events = []
		self.det_events = []
		self.ode_dynamics = []

	def addAgent ( self , name , properties=[] , unique=False ) :
		self.agents.append ( Agent ( name=name , properties=properties , unique=unique ) )

	def resolveAsGrid ( self , agent , length_x , length_y , step_length , two_dim=False , length_z=0. ) :
		agent = [ ag for ag in self.agents if ag.name == agent ][0]
		if not agent.unique : raise Exception ("agent to resolve as grid should be unique")
		agent.unique = False
		agent.grid = True
		agent.length_x=length_x
		agent.length_y=length_y
		agent.two_dim=two_dim
		agent.length_z=length_z
		agent.step_length=step_length

	def resolveAsSphere ( self , agent , length_x , length_y , two_dim=False , max_radius=1. , length_z=0. ) : 
		agent = [ ag for ag in self.agents if ag.name == agent ][0]
		if agent.unique : raise Exception ("agent to resolve as sphere should not be unique")
		agent.sphere = True
		agent.length_x=length_x
		agent.length_y=length_y
		agent.two_dim=two_dim
		agent.length_z=length_z
		agent.max_radius=max_radius
		props = [ "R" , "X" , "Y" ]
		if not two_dim:
			props.append ( "Z" )
		for prop in props :
			if prop not in agent.properties :
				agent.properties.append (prop)
			else :
				raise Exception ("R,X,Y or Z already used as a property ?!")


	def addStochasticEvent ( self , name , agent , kind="normal" , parameters={} , propensity_expression=None , realization_function=None ) :
		stoch_event = StochasticEvent ( name=name , agent=agent , kind=kind , parameters=parameters )
		self.stoch_events.append ( stoch_event )
		if propensity_expression is not None :
			stoch_event.propensity_expression = self.translateExpression ( propensity_expression , stoch_event.agent ) 
		if realization_function is not None :
			stoch_event.realization_function = [ self.translateExpression ( expr , stoch_event.agent ) for expr in realization_function ]
	def addDeterministicEvent ( self , name , agent , kind="normal" , parameters={} , trigger_expression=None , realization_function=None ) :
		det_event = DeterministicEvent ( name=name , agent=agent , kind=kind , parameters=parameters )
		self.det_events.append ( det_event )
		if trigger_expression is not None :
			det_event.trigger_expression = self.translateExpression ( trigger_expression , det_event.agent ) 
		if realization_function is not None :
			det_event.realization_function = [ self.translateExpression ( expr , det_event.agent ) for expr in realization_function ]
	def addContinuousChange ( self , name , agent , property , agent_source=None , parameters={} , rate_expression=None ) :
		change = ContinuousChange ( name=name , agent=agent , property=property , agent_source=agent_source , parameters=parameters ) 
		self.ode_dynamics.append ( change )
		if rate_expression is not None :
			change.rate_expression = self.translateExpression ( rate_expression , change.agent , True )
	def translateExpression ( self , expr , target_agent , for_ode=False ) :
		trans_dict = {}
		# agent props
		for ag in self.agents :
			# search how to access to this agent 
			if ag.name == target_agent : # if it is the target
				if ag.unique : # if unique, access via state
					ag_ref = "state->" + ag.name.lower () + "."
				else : # not unique
					ag_ref = ag.name.lower () + "->"
			else : # non targeted agent: can be accessed only if unique !
				if ag.unique :
					ag_ref = "state->" + ag.name.lower () + "."
				else : 
					ag_ref = ag.name.lower () + "->"
			# iterate on properties
			for prop in ag.properties :
				if for_ode :
					trans_dict[prop] = "(*y)[ " + ag_ref + "I_" + prop + " ]"
				else :
					trans_dict[prop] = ag_ref + prop
					trans_dict["new."+prop] = "new_" + ag_ref + prop
		# parameters
		params = [ event.parameters for event in self.det_events ]
		params.extend ( [ event.parameters for event in self.stoch_events ] )
		params.extend ( [ event.parameters for event in self.ode_dynamics ] )
		params = [ param.keys () for param in params ]
		for param_set in params :
			for param in param_set :
				trans_dict[param] = "Params::" + param
		# do the translation
		trans_expr = []
		atoms = expr.split (" ")
		for atom in atoms :
			if trans_dict.has_key (atom) :
				trans_expr.append ( trans_dict[atom] )
			else :
				trans_expr.append ( atom )
		return " ".join ( trans_expr )

	def giveGridAgents ( self ) :
		return [ ag for ag in self.agents if hasattr (ag,"grid") ]

	def giveSphereAgents (self) :
		return [ ag for ag in self.agents if hasattr (ag,"sphere") ]

	def display (self) :
		for agent in self.agents : 
			agent.display ()
		for stoch_event in self.stoch_events :
			stoch_event.display ()
		for det_event in self.det_events :
			det_event.display ()

