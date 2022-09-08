'''Program to read a pdb file and extract the coordinates of the C-alpha atoms of each residue of a chain.
   Calculate the adjacency, degree and laplacean matrices from the coordinates and write them into files.
   Command line arguments:
	1)pdb file name.
	2)chain id of the chain to be extracted.
   Output files:
	_adjacency.txt: each line represents the row of adjacency matrix.
	_distance.txt: each line represents the row of distance matrix.
	_laplacean.txt: each line represents the row of laplacean matrix.
	_edgelist.txt: each line contains the two residues of an edge.
	_xyz.txt: each line contains the residue number as in pdb file, x coordinate, y coordinate and the z coordinate separeted by tab.
'''
import math
def distance(p1,p2):
   	return math.sqrt(math.pow((float(p1[0])-float(p2[0])),2)+math.pow((float(p1[1])-float(p2[1])),2)+math.pow((float(p1[2])-float(p2[2])),2))

def readpdb(pdbfilename,central_atom,chains,edgelist_file,xyz_file,weight,upper_cutoff,lower_cutoff,res_sep):
	fp=open(pdbfilename,'r');
	fp_co=open(xyz_file,'w');
	coordinates,residue=[],[];
	#extracting the coordinates of CA atoms from the pdb file
	for i in fp.readlines():
		#if more than 1 model is present, it will take the first model and terminate the for loop
		if i[0:6]=='ENDMDL':
			break;
		if((i[0:4]=='ATOM') and (i[26]==' ') and (i[16]==' ' or i[16]=='A' or i[16]=='1')):	
			#and condition to check for alternate position of an atom and take only one position
			temp_central_atom=central_atom;
			if(i[17:20].lstrip().rstrip()=='GLY'):#Always consider C-alpha for GLY
				temp_central_atom='CA';
			if i[12:16].lstrip().rstrip()==temp_central_atom and i[21].lstrip().rstrip() in chains:
				xyz=[];
				xyz.append(float(i[30:38].lstrip().rstrip()));
				xyz.append(float(i[38:46].lstrip().rstrip()));
				xyz.append(float(i[46:54].lstrip().rstrip()));
				coordinates.append(xyz);
				residue.append(i[21].lstrip().rstrip()+i[22:26].lstrip().rstrip());
				fp_co.write(i[21].lstrip().rstrip()+i[22:26].lstrip().rstrip()+'\t'+str(xyz[0])+'\t'+str(xyz[1])+'\t'+str(xyz[2])+'\t'+i[17:20].lstrip().rstrip()+'\n');
	fp_co.close();
	seqLength=len(coordinates);
	fp.close();
	#calculating the distance matrix
	distanceMatrix=[];
	for i in coordinates:
		tempDist=[];
		for j in coordinates:
			tempDist.append(distance(i,j));
		distanceMatrix.append(tempDist);

	#calculating the adjacency matrix and the degree matrix
	adjacencyMatrix=[];
	for i in range(seqLength):
		tempAd=[];
		for j in range(seqLength):
			if i==j:
				tempAd.append(0);
			#discard if residue separatation less than threshold (for long range interaction network)
			elif((residue[j][0]==residue[i][0]) and (int(residue[j][1:])-int(residue[i][1:])<int(res_sep))):
				tempAd.append(0);				
			elif(float(lower_cutoff)<=distanceMatrix[i][j]<=float(upper_cutoff)):
				if weight=='unweighted':
					tempAd.append(1);
				elif weight=='weighted':
					tempAd.append(float(1/distanceMatrix[i][j]));
			else:
				tempAd.append(0);
		adjacencyMatrix.append(tempAd);
	#creating the edgelist
	fp_el=open(edgelist_file[:-4]+'.txt','w');
	if(weight=='weighted'):#weighted network is saved as unweighted for k-clique
		fp_el_uw=open(edgelist_file[:-4]+'_uw.txt','w');
	for i in range(seqLength):
		for j in range(i+1,seqLength):
			if adjacencyMatrix[i][j]==1 and weight=='unweighted':
				fp_el.write(str(i)+' '+str(j)+'\n');
			if adjacencyMatrix[i][j]>0 and weight=='weighted':
				fp_el.write(str(i)+' '+str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
				fp_el_uw.write(str(i)+' '+str(j)+'\n');
	fp_el.close();
	if(weight=='weighted'):
		fp_el_uw.close();

	#creating the edgelist in Long Graph Layout (LGL) format
	fp_el=open(edgelist_file,'w');
	for i in range(seqLength):
		fp_el.write('# '+str(i)+'\n');
		for j in range(i+1,seqLength):
			if adjacencyMatrix[i][j]==1 and weight=='unweighted':
				fp_el.write(str(j)+'\n');
			elif adjacencyMatrix[i][j]>0 and weight=='weighted':
				fp_el.write(str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
	fp_el.close();
	return(0);


def readpdb_anyatom(pdbfilename,chains,edgelist_file,xyz_file,weight,upper_cutoff,lower_cutoff,res_sep):
	#An edge is drawn if any pair of atoms of the two residues are within cutoff distance
	atom_list=['N','CA','C','O','CB','CG','CG1','CG2','CD','CD1','CD2','SD','CE','CE1','CE2','CE3','CZ','CZ2','CZ3','NE','NE1','NE2','CH2','OG','OG1','SG','OH','OD1','OD2','ND1','ND2','OE1','OE2','NZ','NH1','NH2'];
	fp=open(pdbfilename,'r');
	coordinates={};
	residue,residue_name=[],[];
	#extracting the coordinates of CA atoms from the pdb file
	for i in fp.readlines():
		#if more than 1 model is present, it will take the first model and terminate the for loop
		if i[0:6]=='ENDMDL':
			break;
		if((i[0:4]=='ATOM') and (i[26]==' ') and (i[16]==' ' or i[16]=='A' or i[16]=='1')):	
			#and condition to check for alternate position of an atom and take only one position
			if((i[12:16].lstrip().rstrip() in atom_list) and (i[21].lstrip().rstrip() in chains)):
				xyz=[float(i[30:38].lstrip().rstrip()),float(i[38:46].lstrip().rstrip()),float(i[46:54].lstrip().rstrip())];
				resid=i[21].lstrip().rstrip()+i[22:26].lstrip().rstrip();
				if(resid in residue):
					coordinates[resid].append(xyz);
				else:
					residue.append(resid);
					residue_name.append(i[17:20].lstrip().rstrip());
					coordinates[resid]=[xyz];
	fp.close();
	seqLength=len(residue);

	fp_co=open(xyz_file,'w');
	#calculating the adjacency matrix and writing the coordinates of the residue to file
	adjacencyMatrix=[];
	for i in range(seqLength):
		tempAd=[];
		for j in range(seqLength):
			count=0;
			if(i==j):
				count=0;
			#discard if residue separatation less than threshold (for long range interaction network
			elif((residue[j][0]==residue[i][0]) and (int(residue[j][1:])-int(residue[i][1:])<int(res_sep))):
				count=0;
			else:
				for r1 in coordinates[residue[i]]:
					for r2 in coordinates[residue[j]]:
						if(float(lower_cutoff)<=distance(r1,r2)<=float(upper_cutoff)):
							count=count+1;
			tempAd.append(count);
		adjacencyMatrix.append(tempAd);
		#coordinates of geometric centre of all the atoms represent the residue in the network
		xyz_sum=[0,0,0];
		for r in coordinates[residue[i]]:
			xyz_sum[0]=xyz_sum[0]+r[0];
			xyz_sum[1]=xyz_sum[1]+r[1];
			xyz_sum[2]=xyz_sum[2]+r[2];
		n_atoms=len(coordinates[residue[i]]);
		gc_xyz=[float(xyz_sum[0])/n_atoms,float(xyz_sum[1])/n_atoms,float(xyz_sum[2])/n_atoms];
		fp_co.write(residue[i]+'\t'+str(gc_xyz[0])+'\t'+str(gc_xyz[1])+'\t'+str(gc_xyz[2])+'\t'+residue_name[i]+'\n');
	fp_co.close();
			
	#creating the edgelist
	fp_el=open(edgelist_file[:-4]+'.txt','w');
	if(weight=='weighted'):#weighted network is saved as unweighted for k-clique
		fp_el_uw=open(edgelist_file[:-4]+'_uw.txt','w');
	for i in range(seqLength):
		for j in range(i+1,seqLength):
			if adjacencyMatrix[i][j]>0:
				if weight=='unweighted':
					fp_el.write(str(i)+' '+str(j)+'\n');
				elif weight=='weighted':
					fp_el.write(str(i)+' '+str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
					fp_el_uw.write(str(i)+' '+str(j)+'\n');
	fp_el.close();
	if(weight=='weighted'):
		fp_el_uw.close();

	#creating the edgelist in Long Graph Layout (LGL) format
	fp_el=open(edgelist_file,'w');
	for i in range(seqLength):
		fp_el.write('# '+str(i)+'\n');
		for j in range(i+1,seqLength):
			if adjacencyMatrix[i][j]>0:
				if weight=='unweighted':
					fp_el.write(str(j)+'\n');
				elif weight=='weighted':
					fp_el.write(str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
	fp_el.close();
	return(0);

def readpdb_centroid(pdbfilename,chains,edgelist_file,xyz_file,weight,upper_cutoff,lower_cutoff,res_sep):
	#An edge is drawn if centroids(center of mass) of side chains are within cutoff distance
	atom_list=['CB','CG','CG1','CG2','CD','CD1','CD2','SD','CE','CE1','CE2','CE3','CZ','CZ2','CZ3','NE','NE1','NE2','CH2','OG','OG1','SG','OH','OD1','OD2','ND1','ND2','OE1','OE2','NZ','NH1','NH2'];
	fp=open(pdbfilename,'r');
	atomic_mass={'C':12,'N':14,'O':16,'S':32};
	coordinates,res_atom_seq={},{};
	residue,residue_name=[],[];
	#extracting the coordinates of CA atoms from the pdb file
	for i in fp.readlines():
		#if more than 1 model is present, it will take the first model and terminate the for loop
		if i[0:6]=='ENDMDL':
			break;
		if((i[0:4]=='ATOM') and (i[26]==' ') and (i[16]==' ' or i[16]=='A' or i[16]=='1')):	
			#and condition to check for alternate position of an atom and take only one position
			current_atom=i[12:16].lstrip().rstrip();
			if(((current_atom in atom_list) and (i[21].lstrip().rstrip() in chains)) or ((i[17:20].lstrip().rstrip()=='GLY') and (current_atom=="CA") and (i[21].lstrip().rstrip() in chains))):
				#or condition for GLY
				xyz=[float(i[30:38].lstrip().rstrip()),float(i[38:46].lstrip().rstrip()),float(i[46:54].lstrip().rstrip())];
				resid=i[21].lstrip().rstrip()+i[22:26].lstrip().rstrip();
				if(resid in residue):
					res_atom_seq[resid].append(current_atom);
					coordinates[resid].append(xyz);
				else:
					residue.append(resid);
					res_atom_seq[resid]=[current_atom];
					residue_name.append(i[17:20].lstrip().rstrip());
					coordinates[resid]=[xyz];
	fp.close();

	fp_co=open(xyz_file,'w');
	seqLength=len(residue);
	centroids=[];
	#Calculating centroid and writing the coordinates to file
	for i in range(seqLength):
		sum_mx=[0,0,0]; #Sum of (mass*coordinate)
		sum_m=0; #Sum of mass
		for atm_itr in range(len(res_atom_seq[residue[i]])):
			sum_mx[0]=sum_mx[0]+(atomic_mass[res_atom_seq[residue[i]][atm_itr][0]]*coordinates[residue[i]][atm_itr][0]);
			sum_mx[1]=sum_mx[1]+(atomic_mass[res_atom_seq[residue[i]][atm_itr][0]]*coordinates[residue[i]][atm_itr][1]);
			#atomic_mass[res_atom_seq[residue[i]][atm_itr][0]] indicates the first character of the atom identifier
			sum_mx[2]=sum_mx[2]+(atomic_mass[res_atom_seq[residue[i]][atm_itr][0]]*coordinates[residue[i]][atm_itr][2]);
			sum_m=sum_m+atomic_mass[res_atom_seq[residue[i]][atm_itr][0]];
		gc_xyz=[float(sum_mx[0])/sum_m,float(sum_mx[1])/sum_m,float(sum_mx[2])/sum_m];
		centroids.append(gc_xyz);
		fp_co.write(residue[i]+'\t'+str(gc_xyz[0])+'\t'+str(gc_xyz[1])+'\t'+str(gc_xyz[2])+'\t'+residue_name[i]+'\n');
	fp_co.close();

	#calculating the distance matrix
	distanceMatrix=[];
	for i in centroids:
		tempDist=[];
		for j in centroids:
			tempDist.append(distance(i,j));
		distanceMatrix.append(tempDist);

	#calculating the adjacency matrix
	adjacencyMatrix=[];
	for i in range(seqLength):
		tempAd=[];
		for j in range(seqLength):
			if i==j:
				tempAd.append(0);
			#discard if residue separatation less than threshold (for long range interaction network
			elif((residue[j][0]==residue[i][0]) and (int(residue[j][1:])-int(residue[i][1:])<int(res_sep))):
				tempAd.append(0);
			elif(float(lower_cutoff)<=distanceMatrix[i][j]<=float(upper_cutoff)):
				if weight=='unweighted':
					tempAd.append(1);
				elif weight=='weighted':
					tempAd.append(float(1/distanceMatrix[i][j]));
			else:
				tempAd.append(0);
		adjacencyMatrix.append(tempAd);
	#creating the edgelist
	fp_el=open(edgelist_file[:-4]+'.txt','w');
	if(weight=='weighted'):#weighted network is saved as unweighted for k-clique
		fp_el_uw=open(edgelist_file[:-4]+'_uw.txt','w');
	for i in range(seqLength):
		for j in range(i+1,seqLength):
			if adjacencyMatrix[i][j]>0:
				if weight=='unweighted':
					fp_el.write(str(i)+' '+str(j)+'\n');
				elif weight=='weighted':
					fp_el.write(str(i)+' '+str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
					fp_el_uw.write(str(i)+' '+str(j)+'\n');
	fp_el.close();
	if(weight=='weighted'):
		fp_el_uw.close();

	#creating the edgelist in Long Graph Layout (LGL) format
	fp_el=open(edgelist_file,'w');
	for i in range(seqLength):
		fp_el.write('# '+str(i)+'\n');
		for j in range(i+1,seqLength):
			if adjacencyMatrix[i][j]>0:
				if weight=='unweighted':
					fp_el.write(str(j)+'\n');
				elif weight=='weighted':
					fp_el.write(str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
	fp_el.close();
	return(0);

def readpdb_sidechain(pdbfilename,chains,edgelist_file,xyz_file,weight,upper_cutoff,res_sep,distance_matrix_file):
	#An edge is drawn if interaction strength between two residues is within cutoff
	
	atom_list=['CB','CG','CG1','CG2','CD','CD1','CD2','SD','CE','CE1','CE2','CE3','CZ',
	'CZ2','CZ3','NE','NE1','NE2','CH2','OG','OG1','SG','OH','OD1','OD2','ND1','ND2','OE1',
	'OE2','NZ','NH1','NH2'];

	NORM={'ALA':55.7551, 'ARG':93.7891, 'ASN':73.4097, 'ASP':75.1507,'CYS':54.9528, 'GLN':78.1301,
	'GLU':78.8288, 'GLY': 47.3129, 'HIS': 83.7357, 'ILE':67.9452, 'LEU':72.2517, 'LYS': 69.6096, 'MET': 69.2569,
	'PHE': 93.3082, 'PRO': 51.3310, 'SER': 61.3946, 'THR': 63.7075, 'TRP': 106.703, 'TYR': 100.719, 'VAL': 62.3673}

	fp=open(pdbfilename,'r');
	coordinates={};
	residue,sequence=[],[];
	#extracting the coordinates of CA atoms from the pdb file
	for i in fp.readlines():
		#if more than 1 model is present, it will take the first model and terminate the for loop
		if i[0:6]=='ENDMDL':
			break;
		if((i[0:4]=='ATOM') and (i[26]==' ') and (i[16]==' ' or i[16]=='A' or i[16]=='1')):	
			#and condition to check for alternate position of an atom and take only one position
			current_atom=i[12:16].lstrip().rstrip();
			if(((current_atom in atom_list) and (i[21].lstrip().rstrip() in chains)) or ((i[17:20].lstrip().rstrip()=='GLY') and (current_atom=="CA") and (i[21].lstrip().rstrip() in chains))):
				#or condition for GLY
				xyz=[float(i[30:38].lstrip().rstrip()),float(i[38:46].lstrip().rstrip()),float(i[46:54].lstrip().rstrip())];
				resid=i[21].lstrip().rstrip()+i[22:26].lstrip().rstrip();
				if(resid in residue):
					coordinates[resid].append(xyz);
				else:
					residue.append(resid);
					coordinates[resid]=[xyz];
					sequence.append(i[17:20].lstrip().rstrip());#Tracking the sequence information
	fp.close();
	seqLength=len(residue);

	fp_co=open(xyz_file,'w');
	#calculating the adjacency matrix and writing the coordinates of the residue to file
	adjacencyMatrix=[];
	for i in range(seqLength):
		tempAd=[];
		for j in range(seqLength):
			count=0;
			if(i==j):
				count=0;
			#discard if residue separatation less than threshold (for long range interaction network
			elif((residue[j][0]==residue[i][0]) and (int(residue[j][1:])-int(residue[i][1:])<int(res_sep))):
				count=0;
			else:
				for r1 in coordinates[residue[i]]:
					for r2 in coordinates[residue[j]]:
						if distance(r1,r2)<=4.5:
							#Hard threshold of 4.5 is used as given by Vishveshwara(1999,2005)
							count=count+1;
			interaction_strength=float(count*100)/math.sqrt(NORM[sequence[i]]*NORM[sequence[j]]);
			tempAd.append(interaction_strength);
		adjacencyMatrix.append(tempAd);
		#coordinates of geometric centre of sidechain atoms represent the residue in the network
		xyz_sum=[0,0,0];
		for r in coordinates[residue[i]]:
			xyz_sum[0]=xyz_sum[0]+r[0];
			xyz_sum[1]=xyz_sum[1]+r[1];
			xyz_sum[2]=xyz_sum[2]+r[2];
		n_atoms=len(coordinates[residue[i]]);
		gc_xyz=[float(xyz_sum[0])/n_atoms,float(xyz_sum[1])/n_atoms,float(xyz_sum[2])/n_atoms];
		fp_co.write(residue[i]+'\t'+str(gc_xyz[0])+'\t'+str(gc_xyz[1])+'\t'+str(gc_xyz[2])+'\t'+sequence[i]+'\n');
	fp_co.close();

	#writing the distance matrix
	fp_dm=open(distance_matrix_file,'w');
	for i in range(seqLength):
		for j in range(seqLength):
			fp_dm.write(str(adjacencyMatrix[i][j])+'\t');
		fp_dm.write('\n');
	fp_dm.close();

	#creating the edgelist
	fp_el=open(edgelist_file[:-4]+'.txt','w');
	if(weight=='weighted'):#weighted network is saved as unweighted for k-clique
		fp_el_uw=open(edgelist_file[:-4]+'_uw.txt','w');
	for i in range(seqLength):
		for j in range(i+1,seqLength):
			if((i!=j) and (adjacencyMatrix[i][j]>=float(upper_cutoff))):
				if weight=='unweighted':
					fp_el.write(str(i)+' '+str(j)+'\n');
				elif weight=='weighted':
					fp_el.write(str(i)+' '+str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
					fp_el_uw.write(str(i)+' '+str(j)+'\n');
	fp_el.close();
	if(weight=='weighted'):
		fp_el_uw.close();

	#creating the edgelist in Long Graph Layout (LGL) format
	fp_el=open(edgelist_file,'w');
	for i in range(seqLength):
		fp_el.write('# '+str(i)+'\n');
		for j in range(i+1,seqLength):
			if((i!=j) and (adjacencyMatrix[i][j]>=float(upper_cutoff))):
				if weight=='unweighted':
					fp_el.write(str(j)+'\n');
				elif weight=='weighted':
					fp_el.write(str(j)+' '+str(adjacencyMatrix[i][j])+'\n');
	fp_el.close();
	return(0);
