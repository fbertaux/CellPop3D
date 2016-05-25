################################################
# CellPop, Version 1.0
# Francois Bertaux, Inria Paris-Rocquencourt
# francois.bertaux@inria.fr
# March 2015
################################################


#### imports
import modeling
import simulating
import math

#### some check functions
def checkName (name) :
	if name == "" : raise Exception ("Empty name forbidden.")



## function translating structure of CellPop model into a HMA model
def translateStructureCellPopModel (cpmodel,hmamodel) :
	env_properties = []
	# iterate on cell types
	for cell_type in cpmodel.cell_types.values () :
		# get properties of that cell type
		properties = [ np.name for np in cell_type.nativeProteins ] + [ mp.name for mp in cell_type.modifiedProteins ]
		# remove properties belonging to environment
		for prop in properties :
			if "Env." in prop :
				properties.remove (prop)
				env_properties.append ( prop.replace("Env.","") )
		# add properties for gene & mRNA of native proteins
		for np in cell_type.nativeProteins :
			properties += [ np.name + "_gene" , np.name + "_mRNA"  ]
		# add properties age and cell_cycle length if cell type has a division rule
		if cpmodel.divisions.has_key (cell_type.name) :
			properties += [ "age" , "cell_cycle_length" ]
		hmamodel.addAgent ( name=cell_type.name[0].upper()+cell_type.name[1:] , properties=properties )
	hmamodel.addAgent ("Medium",unique=True,properties=env_properties)

## function translating cellular reactions of CellPop model into a HMA model
def translateReactionsCellPopModel (cpmodel,hmamodel) :
	# iterate on cell types
	for cell_type in cpmodel.cell_types.values () :
		ct_name = cell_type.name[0].upper()+cell_type.name[1:]
		# reactions for native prot fluctuations
		for np in cell_type.nativeProteins :
			np_name = np.name
			# gene stoch switches
			hmamodel.addStochasticEvent (name=np_name+"_gene_on_switch",agent=ct_name,parameters={np_name+"_kon":np.kon},
				propensity_expression="( 1. - "+np_name+"_gene ) * "+np_name+"_kon",realization_function=[np_name+"_gene = 1."])
			hmamodel.addStochasticEvent (name=np_name+"_gene_off_switch",agent=ct_name,parameters={np_name+"_koff":np.koff},
				propensity_expression=np_name+"_gene * "+np_name+"_koff",realization_function=[np_name+"_gene = 0."])
			# mRNA prod and degrad
			hmamodel.addStochasticEvent (name=np_name+"_mrna_prod",agent=ct_name,parameters={np_name+"_ksm":np.ksm},
				propensity_expression=np_name+"_gene * "+np_name+"_ksm",realization_function=[np_name+"_mRNA += 1."])
			hmamodel.addStochasticEvent (name=np_name+"_mrna_deg",agent=ct_name,parameters={np_name+"_rm":np.rm},
				propensity_expression=np_name+"_mRNA * "+np_name+"_rm",realization_function=[np_name+"_mRNA -= 1."])
			# protein prod and degrad
			hmamodel.addContinuousChange (name=np_name+"_prod_and_degrad",agent=ct_name,property=np_name,parameters={np_name+"_ksp":np.ksp,np_name+"_rp":np.rp},
				rate_expression=np_name+"_mRNA * "+np_name+"_ksp - "+np_name+"_rp * "+np_name)
		# other cellular reactions
		for c_reac in cell_type.signalingReactions :
			# get the common part of rate expression
			param = c_reac.name+"_"+c_reac.rate[0]
			rate_expression = param
			for reactant in c_reac.reactants :
				if "Env." in reactant :
					rate_expression += " * " + reactant.replace ("Env.","")
				else :
					rate_expression += " * " + reactant
			# add continuous change for each participant in the reaction
			for reactant in c_reac.reactants :
				if "Env." in reactant :
					hmamodel.addContinuousChange (name=c_reac.name+"_change_"+reactant.replace("Env.",""),agent_source=ct_name,agent="Medium",property=reactant.replace("Env.",""),parameters={param:c_reac.rate[1]},
						rate_expression="- "+rate_expression)			
				else :
					hmamodel.addContinuousChange (name=c_reac.name+"_change_"+reactant,agent=ct_name,property=reactant,parameters={param:c_reac.rate[1]},
						rate_expression="- "+rate_expression)			
			for product in c_reac.products :
				if "Env." in product :
					hmamodel.addContinuousChange (name=c_reac.name+"_change_"+product.replace("Env.",""),agent_source=ct_name,agent="Medium",property=product.replace("Env.",""),parameters={param:c_reac.rate[1]},
						rate_expression=rate_expression)			
				else :
					hmamodel.addContinuousChange (name=c_reac.name+"_change_"+product,agent=ct_name,property=product,parameters={param:c_reac.rate[1]},
						rate_expression=rate_expression)
		# degradation of modified proteins
		for mp in cell_type.modifiedProteins :
			if not "Env." in mp.name :
				hmamodel.addContinuousChange (name=mp.name+"_degrad",agent=ct_name,property=mp.name,parameters={mp.name+"_degrad_rate":mp.degRate},
					rate_expression="- "+mp.name+" * "+mp.name+"_degrad_rate")						
	# degradation of environment species
	for env_species in cpmodel.environment_species_degrad.keys () :
		hmamodel.addContinuousChange (name=env_species+"_degrad",agent="Medium",property=env_species,parameters={env_species+"_degrad_rate":math.log(2.)/cpmodel.environment_species_degrad[env_species]},
			rate_expression="- "+env_species+" * "+env_species+"_degrad_rate")

## function translating death rules of CellPop model into HMA model
def translateDeathRulesCellPopModel ( cpmodel , hmamodel ) :
	for cell_type in cpmodel.death_rules.keys () :
		death_rule_bool_expr = cpmodel.death_rules[cell_type]["expr"]
		death_rule_params = cpmodel.death_rules[cell_type]["params"]
		hmamodel.addDeterministicEvent (name="death_rule_"+cell_type,agent=cell_type,kind="destruction",parameters=death_rule_params,trigger_expression=death_rule_bool_expr)

def translateDivisionsCellPopModel ( cpmodel , hmamodel ) :
	for cell_type in cpmodel.divisions.keys () :
		## do the division event
		avg = cpmodel.divisions[cell_type]["avg_length"]
		std = cpmodel.divisions[cell_type]["std_length"]
		real_func = ["new.age = age - cell_cycle_length" ,"age = new.age" ,
					 "new.cell_cycle_length = ran.norm ( "+cell_type+"_cell_cycle_length_avg"+" , "+cell_type+"_cell_cycle_length_std )" ,
					 "cell_cycle_length = ran.norm ( "+cell_type+"_cell_cycle_length_avg"+" , "+cell_type+"_cell_cycle_length_std )"]
		# add the state copy instructions
		properties = [ np.name for np in cpmodel.cell_types[cell_type].nativeProteins ] + [ mp.name for mp in cpmodel.cell_types[cell_type].modifiedProteins ]
		for prop in [ prop for prop in properties if not "Env." in prop ] :
			real_func.append ( "new."+prop+" = "+prop )
		for np in cpmodel.cell_types[cell_type].nativeProteins :
			real_func.append ( "new."+np.name+"_gene = "+np.name+"_gene")
			real_func.append ( "new."+np.name+"_mRNA = "+np.name+"_mRNA")
		hmamodel.addDeterministicEvent (name="division_"+cell_type,agent=cell_type,kind="creation",parameters={cell_type+"_cell_cycle_length_avg":avg,cell_type+"_cell_cycle_length_std":std},
			trigger_expression="age > cell_cycle_length",
			realization_function=real_func)
		## add the aging rule
		hmamodel.addContinuousChange (name=cell_type+"_aging",agent=cell_type,property="age",rate_expression="1.")	

### function translating a CellPop model into a HMA model
def translateCellPopModel ( model ) :
	hmamodel = modeling.Model ( model.name )
	translateStructureCellPopModel ( model , hmamodel )
	translateReactionsCellPopModel ( model , hmamodel )
	translateDeathRulesCellPopModel ( model , hmamodel )
	translateDivisionsCellPopModel ( model , hmamodel )
	return hmamodel


####  classes to describe a CellPop model
class Model (object) :
	def __init__ ( self , name="model" ) :
		self.name = name
		self.cell_types = {}
		self.environment_species_degrad = {}
		self.death_rules = {}
		self.divisions = {}
	def addCellType ( self , name ) :
		checkName (name)
		if name in [ ct.name for ct in self.cell_types.keys() ] :
			raise Exception ("Cell Type with that name exists already")
		self.cell_types[name] = CellType (name)
	def addFluctuatingProtein (self,standard=False,**kwargs) :
		if standard :
			self.cell_types[kwargs["cell_type"]].addNativeProteinStdFluct (name=kwargs["name"],EP=kwargs["EP"],dilutionHalfLife=kwargs["dilutionHalfLife"])
		else :
			self.cell_types[kwargs["cell_type"]].addNativeProtein (name=kwargs["name"],HLP=kwargs["HLP"],EP=kwargs["EP"],
																   HLM=kwargs["HLM"],EM=kwargs["EM"],Ton=kwargs["Ton"],Toff=kwargs["Toff"])
	def addCellularReaction (self,**kwargs) :
		self.cell_types[kwargs["cell_type"]].addReaction (name=kwargs["name"],reactants=kwargs["reactants"],products=kwargs["products"],rate=kwargs["rate"])
	def addCellularReversibleReaction (self,**kwargs) :
		self.cell_types[kwargs["cell_type"]].addReversibleReaction (name=kwargs["name"],reactants=kwargs["reactants"],products=kwargs["products"],rates=kwargs["rates"])
	def addCellularCatalyticReaction (self,**kwargs) :
		self.cell_types[kwargs["cell_type"]].addCatalyticReaction (name=kwargs["name"],substrate=kwargs["substrate"],catalyst=kwargs["catalyst"],product=kwargs["product"],rates=kwargs["rates"])
	def setAllModifiedProteinDegradation (self,**kwargs) :
		self.cell_types[kwargs["cell_type"]].setAllModifiedProteinDegradation (halfLife=kwargs["halfLife"])
	def setModifiedProteinDegradation (self,**kwargs) :
		self.cell_types[kwargs["cell_type"]].setModifiedProteinDegradation (name=kwargs["name"],halfLife=kwargs["halfLife"])
	def setEnvironmentSpeciesDegradation (self,name,halfLife) :
		self.environment_species_degrad[name]=halfLife
	def setCellDeathRule (self,cell_type,death_bool_expression,params={}) :
		self.death_rules[cell_type]={}
		self.death_rules[cell_type]["params"]=params
		self.death_rules[cell_type]["expr"]=death_bool_expression
	def defineDivision (self,cell_type,cell_cycle_length_avg,cell_cycle_length_std) :
		self.divisions[cell_type]={}
		self.divisions[cell_type]["avg_length"]=cell_cycle_length_avg
		self.divisions[cell_type]["std_length"]=cell_cycle_length_std

class NativeProtein (object) :
	def __init__ ( self , name , kon , koff , ksm , rm , ksp , rp ) :
		checkName (name)
		self.name = name
		self.kon = kon
		self.koff = koff
		self.ksm = ksm
		self.rm = rm
		self.ksp = ksp
		self.rp = rp

class ModifiedProtein (object) :
	def __init__ ( self , name , degRate ) :
		checkName (name)
		self.name = name
		self.degRate = degRate

class SignalingReaction (object) :
	def __init__ ( self , name , reactants , products , rate ) :
		self.name = name
		self.rate = rate
		self.reactants = reactants
		self.products = products

class CellType (object) :
	def __init__ ( self , name = "Cell" ) :
		self.name = name
		self.nativeProteins = []
		self.modifiedProteins = []
		self.signalingReactions = []
	def addNativeProtein ( self , name , Ton , Toff , EM , HLM , EP , HLP ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		for s in ["*",":"] : 
			if s in name : raise Exception ("Symbols * or , not allowed in names.")
		if name in [ prot.name for prot in self.nativeProteins ] : raise Exception ("Name already taken by a nativeProt.")
		if name in [ prot.name for prot in self.modifiedProteins ] : raise Exception ("Name already taken by a modifiedProt.")
		koff = 1./Ton
		kon = 1./Toff
		EG = kon / (kon+koff)
		rm = math.log(2.) / HLM 
		ksm = EM * rm / EG
		rp = math.log(2.) / HLP
		ksp = EP * rp / EM
		nativeProt = NativeProtein ( name=name , kon=kon , koff=koff , ksm=ksm , rm=rm , ksp=ksp , rp=rp )
		self.nativeProteins.append (nativeProt)
	def addNativeProteinStdFluct ( self , name , EP , dilutionHalfLife ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		if name in [ prot.name for prot in self.nativeProteins ] : raise Exception ("Name already taken by a nativeProt.")
		if name in [ prot.name for prot in self.modifiedProteins ] : raise Exception ("Name already taken by a modifiedProt.")
		self.addNativeProtein ( name=name , Ton=0.1 , Toff=2.58 , EM=17. , HLM=9. , EP=EP , HLP=dilutionHalfLife)		
	def addModifiedProtein ( self , name , degRate=0. ) :
		if name == "" : raise Exception ("Empty name forbidden.")
		if name == "" : raise Exception ("Empty name forbidden.")
		for s in ["*",":"] : 
			if s in name : raise Exception ("Symbols * or , not allowed in names.")
		if name in [ prot.name for prot in self.nativeProteins ] : raise Exception ("Name already taken by a nativeProt.")
		if name in [ prot.name for prot in self.modifiedProteins ] : raise Exception ("Name already taken by a modifiedProt.")
		mprot = ModifiedProtein ( name , degRate )
		self.modifiedProteins.append (mprot)
	def setModifiedProteinDegradation ( self , name , halfLife ) :
		mprot_matches = [ mprot for mprot in self.modifiedProteins if mprot.name == name ]
		if len(mprot_matches) == 0 : raise Exception ("Not existing modifiedProt: %s" % name)
		mprot_matches[0].degRate = math.log(2.) / halfLife
	def setAllModifiedProteinDegradation ( self , halfLife ) :
		for mprot in self.modifiedProteins : self.setModifiedProteinDegradation ( mprot.name , halfLife )
	def addReaction ( self , name , reactants , products , rate ) :
		if len (reactants) == 0 : raise Exception ("I don't accept 0-th order reaction (for now at least.)")
		for reactant in reactants :
			if reactant not in [ prot.name for prot in self.nativeProteins ] + [ mprot.name for mprot in self.modifiedProteins ] :
				self.addModifiedProtein ( name=reactant )
		for product in products :
			if product not in [ prot.name for prot in self.nativeProteins ] + [ mprot.name for mprot in self.modifiedProteins ] :
				self.addModifiedProtein ( name=product )
		reaction = SignalingReaction ( name=name , reactants=reactants , products=products  , rate=rate)
		self.signalingReactions.append (reaction)
	def addReversibleReaction ( self , name , reactants , products , rates ) :
		self.addReaction ( name=name+"_forward" , reactants=reactants , products=products , rate=rates[0] )
		self.addReaction ( name=name+"_backward" , reactants=products , products=reactants , rate=rates[1] )
	def addCatalyticReaction ( self , name , substrate , catalyst , product , rates ) :
		self.addReaction ( name=name+"_binding" , reactants=[substrate,catalyst] , products=[substrate+"_"+catalyst] , rate=rates[0] )
		self.addReaction ( name=name+"_unbinding" , reactants=[substrate+"_"+catalyst] , products=[substrate,catalyst] , rate=rates[1] )
		self.addReaction ( name=name+"_catalysis" , reactants=[substrate+"_"+catalyst] , products=[catalyst,product] , rate=rates[2] )
	def giveProteinIndexFromName ( self , name ) :
		if name not in [ nprot.name for nprot in self.nativeProteins ] + [ mprot.name for mprot in self.modifiedProteins ] : raise Exception ("Protein does not exist.")
		np_matches = [ nprot for nprot in self.nativeProteins if nprot.name==name ]
		if len(np_matches) > 0 : idxProt = self.nativeProteins.index(np_matches[0])
		mp_matches = [ mprot for mprot in self.modifiedProteins if mprot.name==name ]
		if len(mp_matches) > 0 : idxProt = len(self.nativeProteins)+self.modifiedProteins.index(mp_matches[0])
		return idxProt


