############## SOME HELPERS FUNCTIONS #########

def giveAgentIterationCode ( agent , state_prefix="" , indent=0 ) :
	ag_name = agent.name.lower ()
	ag_idx = ag_name[0:2]
	ag_vector = state_prefix + ag_name + "s"
	line_for = "\t"*indent + "for ( Uint " + ag_idx + " = 0 ; "
	line_for += ag_idx + " < " + ag_vector + ".size () ; " + ag_idx + "++ )\n"
	line_for_rev_it = "\t"*indent + "for ( vector <" + agent.name + "*>::reverse_iterator rit = "
	line_for_rev_it += ag_vector + ".rbegin () ; rit != " + ag_vector + ".rend () ; ++rit )\n"
	return { "line_for":line_for , "line_for_rev_it":line_for_rev_it , "ag_name":ag_name , "ag_idx":ag_idx , "ag_vector":ag_vector }

