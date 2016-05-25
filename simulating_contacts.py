

##### CONTACT DETECTION CODE ######


def writeContactsCodeHeader ( model , target_folder ) :
	header = "#ifndef CONTACTS_HPP\n#define CONTACTS_HPP\n\n"
	header += "\n#include \"state.hpp\"\n"
	header += "#include \"parameters.hpp\"\n\n"
	
	# the class describing a voxel of the contact grid
	header += "struct ContactVoxel\n{\n"
	header += "\t" + "list <ContactVoxel*> neighbors ;\n"
	# the indexes of the spheres agent in that voxel
	for agent in model.giveSphereAgents () :
		header += "\t" + "list <Int> " + "I_" + agent.name.lower() + "s ;\n"
	header += "} ;\n\n"

	# the classes for storing contacts
	for agent_a in model.giveSphereAgents () :
		for agent_b in model.giveSphereAgents () :
			# self-self contact
			if agent_a is agent_b :
				header += "struct Contact_" + agent.name + "_" + agent.name + "\n{\n"
				header += "\t" + agent.name + "* A ;\n"
				header += "\t" + agent.name + "* B ;\n"
				header += "\tDoub dist_sqr , dist ;\n"
				header += "\tDoub dx , dy"
				if not model.giveSphereAgents ()[0].two_dim :
					header += " , dz"
				header += " ;\n"
				header += "} ;\n\n"
			# self-other contact
			if model.giveSphereAgents ().index (agent_a) < model.giveSphereAgents ().index (agent_b) :
				raise Exception ("Not yet implemented, contacts btw dfft sphere types")

	# the class describing the contact grid
	header += "struct ContactDetector\n{\n"
	# constructor
	header += "\t// constructor\n"
	header += "\tContactDetector ( Doub lx , Doub ly , "
	if not model.giveSphereAgents ()[0].two_dim :
		header += "Doub lz , "
	header += "Doub lv ) ;\n\n"
	# the fields of the contact grid
	header += "\t// fields\n"
	header += "\tconst Int nx , ny "
	if not model.giveSphereAgents ()[0].two_dim :
		header += ", nz "
	header += ";\n\tconst Doub l_voxel ;\n"
	header += "\t const Doub length_x , length_y "
	if not model.giveSphereAgents ()[0].two_dim :
		header += ", length_z "
	header += ";\n"
	header += "\tDoub x_first_vox , y_first_vox "
	if not model.giveSphereAgents ()[0].two_dim :
		header += ", z_first_vox "
	header += ";\n"
	header += "\tvector <ContactVoxel*> voxels ;\n"
	for agent_a in model.giveSphereAgents () :
		for agent_b in model.giveSphereAgents () :
			if agent_a is agent_b :
				header += "\tvector <Contact_" + agent_a.name + "_" + agent_a.name + "> contacts_"
				header += agent_a.name.lower() + "_" + agent_a.name.lower() + " ;\n"
			if model.giveSphereAgents ().index (agent_a) < model.giveSphereAgents ().index (agent_b) :
				raise Exception ("Not yet implemented, contacts btw dfft sphere types")
	# methods of the contact grid
	header += "\n\t// methods\n"
	header += "\tvoid detectContacts ( State* state , VecDoub* y ) ;\n"
	header += "\tvoid localizeSpheresInVoxels ( State* state , VecDoub* y ) ;\n"
	for agent_a in model.giveSphereAgents () :
		for agent_b in model.giveSphereAgents () :
			if agent_a is agent_b :
				header += "\tvoid detect_" + agent_a.name + "_" + agent_a.name + "_Contact ( State* state , VecDoub* y , "
				header += "Int " + agent_a.name.lower()[0:2] + "_1 , Int " + agent_a.name.lower()[0:2] + "_2 ) ;\n"
			if model.giveSphereAgents ().index (agent_a) < model.giveSphereAgents ().index (agent_b) :
				raise Exception ("Not yet implemented, contacts btw dfft sphere types")
	header += "} ;\n\n"
	writer = open ( target_folder + "/contactdetector.hpp" , "w" )
	writer.write ( header )
	writer.write ( "\n#endif\n\n")


def assembleContactGridConstructor ( model , target_folder ) :
	# beginning
	declaration = "ContactDetector::ContactDetector ( Doub lx , Doub ly"
	if not model.giveSphereAgents ()[0].two_dim :
		declaration += " , Doub lz"
	declaration += " , Doub lv ) :\n\tnx(lx/lv) , ny(ly/lv)"
	if not model.giveSphereAgents ()[0].two_dim :
		declaration += " , nz(lz/lv)"
	declaration += " , l_voxel(lv) ,\n\tlength_x(lx) , length_y(ly)"
	if not model.giveSphereAgents ()[0].two_dim :
		declaration += " , length_z(lz)"
	declaration += "\n{\n"	
	# body: computation of x,y,(z) of first voxel
	declaration += "\tx_first_vox = - l_voxel * (nx/2. - 0.5) ;\n"
	declaration += "\ty_first_vox = - l_voxel * (ny/2. - 0.5) ;\n"
	if not model.giveSphereAgents ()[0].two_dim :
		declaration += 	"\tz_first_vox = - l_voxel * (nz/2. - 0.5) ;\n"
	# body: voxels construction
	declaration += "\n\t// voxels construction\n"
	if model.giveSphereAgents ()[0].two_dim :
		declaration += "\tfor ( Int I = 0 ; I < nx*ny ; I++ ) voxels.push_back ( new ContactVoxel ) ;\n"
	else :
		declaration += "\tfor ( Int I = 0 ; I < nx*ny*nz ; I++ ) voxels.push_back ( new ContactVoxel ) ;\n"
	# body: creation of the topology
	declaration += "\n\t// construction of the topology\n"
	# beginning on the dimensional loop on voxels
	declaration += "\tfor ( Int ix = 0 ; ix < nx ; ix++ )\n\t{\n"
	declaration += "\t\tfor ( Int iy = 0 ; iy < ny ; iy++ )\n\t\t{\n"
	if not model.giveSphereAgents ()[0].two_dim :
		declaration += "\t\t\tfor ( Int iz = 0 ; iz < nz ; iz++ )\n\t\t\t{\n"
	# the corresponding unidimensional index
	if model.giveSphereAgents ()[0].two_dim :
		declaration += "\t\t\tInt I = ny*ix + iy ; // unidimensional index of the voxel\n"
	else :
		declaration += "\t\t\t\tInt I = ny*nz*ix + nz*iy + iz ; // unidimensional index of the voxel\n"
	# relationships: 2D case
	if model.giveSphereAgents ()[0].two_dim :
		declaration+= "\t\t\t// left-right relationship\n"
		declaration+= "\t\t\tif ( ix+1 < nx )\n"
		declaration+= "\t\t\t{\n"
		declaration+= "\t\t\t	voxels[I]->neighbors.push_back ( voxels[ny*(ix+1) + iy] ) ;\n"
		declaration+= "\t\t\t	voxels[ny*(ix+1) + iy]->neighbors.push_back ( voxels[I] ) ;\n"
		declaration+= "\t\t\t}\n"
		declaration+= "\t\t\t// down-up relationship\n"
		declaration+= "\t\t\tif ( iy+1 < ny )\n"
		declaration+= "\t\t\t{\n"
		declaration+= "\t\t\t	voxels[I]->neighbors.push_back ( voxels[ny*ix + iy+1] ) ;\n"
		declaration+= "\t\t\t	voxels[ny*ix + iy+1]->neighbors.push_back ( voxels[I] ) ;\n"
		declaration+= "\t\t\t}\n"
		declaration+= "\t\t\t// diag-left relationship\n"
		declaration+= "\t\t\tif ( ix+1 < nx && iy+1 < ny )\n"
		declaration+= "\t\t\t{\n"
		declaration+= "\t\t\t	voxels[I]->neighbors.push_back ( voxels[ny*(ix+1) + iy+1] ) ;\n"
		declaration+= "\t\t\t	voxels[ny*(ix+1) + iy+1]->neighbors.push_back ( voxels[I] ) ;\n"				
		declaration+= "\t\t\t}\n"
		declaration+= "\t\t\t// diag-right relationship\n"
		declaration+= "\t\t\tif ( ix > 0 && iy+1 < ny )\n"
		declaration+= "\t\t\t{\n"
		declaration+= "\t\t\t	voxels[I]->neighbors.push_back ( voxels[ny*(ix-1) + iy+1] ) ;\n"
		declaration+= "\t\t\t	voxels[ny*(ix-1) + iy+1]->neighbors.push_back ( voxels[I] ) ;\n"				
		declaration+= "\t\t\t}\n"
	# relationships: 3D case
	if not model.giveSphereAgents ()[0].two_dim :
		raise Exception ("Not finished yet, topology..")
	# end of the dimensional loop on voxels
	if not model.giveSphereAgents ()[0].two_dim :
		declaration += "\t\t\t}\n"
	declaration += "\t\t}\n\t}\n"
	declaration += "}\n\n"
	return declaration



def assembleSpheresLocalizationMethod ( model , target_folder ) :
	definition = "void ContactDetector::localizeSpheresInVoxels ( State* state , VecDoub* y )\n{\n"
	definition += "\t// reset previous localization\n"
	definition += "\tfor ( Uint v = 0 ; v < voxels.size () ; v++ )\n\t{\n"
	for agent in  model.giveSphereAgents () :
		definition += "\t\tvoxels[v]->I_" + agent.name.lower() + "s.clear () ;\n"
	definition += "\t}\n"
	definition += "\t// localization, per type of sphere agents\n"	
	definition += "\tInt ix , iy"
	if not model.giveSphereAgents ()[0].two_dim :
		definition += " , iz"
	definition += " ;\n"
	# one iteration per type of sphere
	for agent in model.giveSphereAgents () :
		definition += "\tfor ( Uint "+agent.name.lower()[0:2]+" = 0 ; "+agent.name.lower()[0:2]
		definition += " < state->"+agent.name.lower()+"s.size () ; "+agent.name.lower()[0:2]+"++ )\n\t{\n"
		# spatial index from X,Y(,Z)
		definition += "\t\t// get spatial indexing from position\n"
		definition += "\t\tix = ( (*y)[ state->"+agent.name.lower()+"s["+agent.name.lower()[0:2]+"]->I_X ] - x_first_vox ) / l_voxel ;\n"
		definition += "\t\tiy = ( (*y)[ state->"+agent.name.lower()+"s["+agent.name.lower()[0:2]+"]->I_Y ] - y_first_vox ) / l_voxel ;\n"
		if not model.giveSphereAgents ()[0].two_dim :
			definition += "\t\tiz = ( (*y)[ state->"+agent.name.lower()+"s["+agent.name.lower()[0:2]+"]->I_Z ] - z_first_vox ) / l_voxel ;\n"
		# if outside, put in closest boundary box
		definition += "\t\t// if outside limits, use the closest voxel\n"
		definition += "\t\tif ( ix > nx-1 ) { ix = nx-1 ; } else if ( ix < 0 ) { ix = 0 ; }\n"
		definition += "\t\tif ( iy > ny-1 ) { iy = ny-1 ; } else if ( iy < 0 ) { iy = 0 ; }\n"	
		if not model.giveSphereAgents ()[0].two_dim :
			definition += "\t\tif ( iz > nz-1 ) { iz = nz-1 ; } else if ( iz < 0 ) { iz = 0 ; }\n"	
		# add sphere to its voxel
		if model.giveSphereAgents ()[0].two_dim :
			definition += "\t\tvoxels[ ny*ix + iy ]->I_"+agent.name.lower()+"s.push_back ("+agent.name.lower()[0:2]+") ;\n"
		else :
			definition += "\t\tvoxels[ ny*nz*ix + nz*iy + iz ]->I_"+agent.name.lower()+"s.push_back ("+agent.name.lower()[0:2]+") ;\n"
		definition += "\t}\n"
	definition += "}\n\n"
	return definition

def assembleContactDetectionMethod ( model , target_folder ) :
	definition = "void ContactDetector::detectContacts ( State* state , VecDoub* y )\n{\n"
	# clearing previous contacts
	definition += "\t// clear previously detected contacts\n"
	for agent_a in model.giveSphereAgents () :
		for agent_b in model.giveSphereAgents () :
			if agent_a is agent_b :
				definition += "\tcontacts_" + agent_a.name.lower() + "_" + agent_b.name.lower() + ".clear () ;\n"
			if model.giveSphereAgents ().index (agent_a) < model.giveSphereAgents ().index (agent_b) :
				raise Exception ("Not yet implemented, contacts btw dfft sphere types")	
	# localization step
	definition += "\t// localization of spheres on the contact grid\n"
	definition += "\tlocalizeSpheresInVoxels ( state , y ) ;\n"
	# detection of self-self contacts
	for agent in model.giveSphereAgents () :
		definition += assembleSelfContactDetectionLoop ( model , target_folder , agent )
	definition += "}\n\n"
	return definition


def assembleSelfContactDetectionLoop ( model , target_folder , agent ) :
	ag_N = agent.name
	ag_n = agent.name.lower ()
	ag_1 = agent.name.lower()[0:2] + "_1"
	ag_2 = agent.name.lower()[0:2] + "_2"
	# loop on voxels of the contact grid
	definition = "\t// detection of " + ag_n + "-" + ag_n + " contacts\n"
	definition += "\tfor ( Uint v = 0 ; v < voxels.size () ; v++ )\n\t{\n"
	# loop on agents in that voxel
	definition += "\t\t// loop on " + ag_n + "s in that voxel\n"
	definition += "\t\tfor (list <Int>::iterator " + ag_1 + " = voxels[v]->I_"+ag_n+"s.begin() ; "
	definition += ag_1 + " != voxels[v]->I_" + ag_n + "s.end() ; ++" + ag_1 +" )\n\t\t{\n"
	# potential contacts in the same voxel
	definition += "\t\t\t// potential " + ag_n + "s in the same voxel\n"
	definition += "\t\t\tfor (list <Int>::iterator " + ag_2 + " = voxels[v]->I_"+ag_n+"s.begin() ; "
	definition += ag_2 + " != voxels[v]->I_" + ag_n + "s.end() ; ++" + ag_2 +" )\n\t\t\t{\n"
	definition += "\t\t\t\tdetect_" + ag_N + "_" + ag_N + "_Contact ( state , y , *" + ag_1 + " , *" + ag_2 + " ) ;\n\t\t\t}\n"
	# iteration on neighbor voxels
	definition += "\t\t\t// search neighbor voxels\n"
	definition += "\t\t\tfor (list <ContactVoxel*>::iterator vn = voxels[v]->neighbors.begin() ; vn!= voxels[v]->neighbors.end() ; ++vn )\n\t\t\t{\n"
	# potential contacts with agent in neighbor voxel
	definition += "\t\t\t\t// potential contacts with cells in neighbor voxel\n"
	definition += "\t\t\t\tfor (list <Int>::iterator " + ag_2 + " = (*vn)->I_"+ag_n+"s.begin() ; "
	definition += ag_2 + " != (*vn)->I_" + ag_n + "s.end() ; ++" + ag_2 +" )\n\t\t\t\t{\n"
	definition += "\t\t\t\t\tdetect_" + ag_N + "_" + ag_N + "_Contact ( state , y , *" + ag_1 + " , *" + ag_2 + " ) ;\n\t\t\t\t}\n"
	definition += "\t\t\t}\n"
	definition += "\t\t}\n"
	definition += "\t}\n"
	return definition


def assembleSelfContactDetectionMethod ( model , target_folder , agent ) :
	ag_1 = agent.name.lower()[0:2] + "_1"
	ag_2 = agent.name.lower()[0:2] + "_2"
	ag_vec = "state->" + agent.name.lower() + "s"
	ag_1_from_vec = ag_vec + "[" + ag_1 + "]"
	ag_2_from_vec = ag_vec + "[" + ag_2 + "]"
	definition = "void ContactDetector::detect_" + agent.name + "_" + agent.name + "_Contact ( State* state , VecDoub* y , "
	definition += "Int " + ag_1 + " , Int " + ag_2 + " )\n{\n"
	# avoid double check
	definition += "\t// don't check twice the same potential contact\n"
	definition += "\tif ( " + ag_1 + " < " + ag_2 + " ) return ;\n"
	# dist2 below which contact
	definition += "\t// compute the distance^2 below which there is contact\n"
	definition += "\tDoub contdist2 = SQR ( (*y)[ " + ag_1_from_vec +"->I_R ] + (*y)[ " + ag_2_from_vec + "->I_R ] ) ;\n\n"
	# dx computation, fast check
	definition += "\t// dx computation, fast check for no contact\n"
	definition += "\tDoub dx = (*y)[ " + ag_1_from_vec +"->I_X ] - (*y)[ " + ag_2_from_vec + "->I_X ] ;\n"
	definition += "\tif (dx*dx>contdist2) { return; }\n\n"
	# dy computation, fast check
	definition += "\t// dy computation, fast check for no contact\n"
	definition += "\tDoub dy = (*y)[ " + ag_1_from_vec + "->I_Y ] - (*y)[ " + ag_2_from_vec + "->I_Y ] ;\n"
	definition += "\tif (dy*dy>contdist2) { return; }\n\n"
	# check and store contact
	definition += "\t// check contact and store it\n"
	definition += "\tif ( dx*dx+dy*dy  < contdist2 )\n"
	definition += "\t{\n"
	definition += "\t\t// create the contact\n"
	definition += "\t\tContact_" + agent.name + "_" + agent.name + " contact ;\n"
	definition += "\t\tcontact.A = " + ag_1_from_vec + " ;\n"
	definition += "\t\tcontact.B = " + ag_2_from_vec + " ;\n"
	definition += "\t\tcontact.dx = dx ;\n"
	definition += "\t\tcontact.dy = dy ;\n"
	definition += "\t\tcontact.dist_sqr = dx*dx+dy*dy ;\n"
	definition += "\t\tcontact.dist = sqrt ( contact.dist_sqr ) ;\n"
	definition += "\t\t// add the contact\n"
	definition += "\t\tcontacts_" + agent.name.lower() + "_" + agent.name.lower() + ".push_back ( contact ) ;\n"
	definition += "\t}\n"
	definition += "}\n\n"
	return definition


def writeContactsCodeSource ( model , target_folder ) :
	source = "#include \"contactdetector.hpp\"\n\n"
	source += assembleContactGridConstructor ( model , target_folder )
	source += assembleSpheresLocalizationMethod ( model , target_folder )
	source += assembleContactDetectionMethod ( model , target_folder )
	for agent in model.giveSphereAgents () :
		source += assembleSelfContactDetectionMethod ( model , target_folder , agent )
	writer = open ( target_folder + "/contactdetector.cpp" , "w" )
	writer.write ( source )


def writeContactsCode ( model , target_folder ) :
	# if no sphere agents, no need for that code
	if  not model.giveSphereAgents () :
		return
	# otherwise write the code
	writeContactsCodeHeader ( model , target_folder )
	writeContactsCodeSource ( model , target_folder )
