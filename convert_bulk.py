#file to convert the internal files generated for network analysis to the downloadable files
#input: 1) file name (eg. 1234_edgelist.txt)
#2) jobid

def convert_bulkfiles(name,jobid,NAPS_path):
	file_type=name.split('.')[0].split('_')[-1];
	file_name=name.split('_'+file_type[0])[0];

	# NAPS_path='/var/www/html/NAPS/files_generated/bulkupload/'+jobid+'/';

	#Reading the residue names
	fp_xyz=open(NAPS_path+jobid+'/xyz/'+file_name+'_xyz.txt');
	residue=[];
	for line in fp_xyz.readlines():
		residue.append(line.split()[0]);
	fp_xyz.close();

	if(file_type=='edgelist'):#For edgelist file
		fp_in=open(NAPS_path+jobid+'/edgelist/'+name,'r');
		fp_out=open(NAPS_path+jobid+'/download/'+file_name+'_edgelist.txt','w');
		# fp_out=open(NAPS_path+jobid+'/download/'+file_name+'_edgelist.txt','wb');
		for line in fp_in.readlines():
			sp=line.split();
			if(len(sp)>2):#for weighted networks
				fp_out.write(str(residue[int(sp[0])])+'\t'+str(residue[int(sp[1])])+'\t'+sp[2]+'\r\n');
			else:
				fp_out.write(str(residue[int(sp[0])])+'\t'+str(residue[int(sp[1])])+'\r\n');
		fp_in.close();
		fp_out.close();
	elif(file_type=='centrality'):#For centrality file
		fp_in=open(NAPS_path+'centrality/'+name,'r');
		fp_out=open(NAPS_path+'download/'+file_name+'_centrality.txt','w');
		for line in fp_in.readlines():
			sp=line.split();
			if(sp[0]=='Node'):
				fp_out.write(line[:-1]+'\r\n');
			else:
				fp_out.write(str(residue[int(sp[0])-1]));
				for itr in range(1,len(sp)):
					fp_out.write('\t'+sp[itr]);
				fp_out.write('\r\n');
		fp_in.close();
		fp_out.close();
	elif(file_type=='clique'):#For clique file
		fp_in=open(NAPS_path+'clique/'+name,'r');
		fp_out=open(NAPS_path+'download/'+file_name+'_clique.txt','w');
		for line in fp_in.readlines():
			l=line[:-1];#Choppinf off the newline char at the end of line
			cl=l.split(':');
			sp=cl[1].split(',');
			fp_out.write(cl[0]+':');
			for itr in range(len(sp)-1):
				fp_out.write(str(residue[int(sp[itr])])+',');
			fp_out.write(str(residue[int(sp[-1])])+'\r\n');
		fp_in.close();
		fp_out.close();
	elif(file_type=='path'):#For shortest path file
		fp_in=open(NAPS_path+'shortest_path/'+name,'r');
		fp_out=open(NAPS_path+'download/'+file_name+'_path.txt','w');
		line=fp_in.readline();
		fp_out.write(line);
		for line in fp_in.readlines():
			l=line[:-1];#Choppinf off the newline char at the end of line
			sp=l.split('-');
			for itr in range(len(sp)-1):
				fp_out.write(residue[int(sp[itr])]+'-');
			fp_out.write(residue[int(sp[-1])]+'\r\n');
		fp_in.close();
		fp_out.close();
	elif(file_type=='lap'):
		fp_in=open(NAPS_path+'laplacian/'+name,'r');
		ssevc=[];
		for line in fp_in.readlines():
			sp=line.split();
			if(sp[0]!='Node'):
				ssevc.append(sp[1]);
		fp_in.close();
		fp_in=open(NAPS_path+'centrality/'+name[:10]+'_centrality.txt','r');
		alevc=[];
		for line in fp_in.readlines():
			sp=line.split('\t');
			if(sp[0]!='Node'):
				alevc.append(sp[5]);
		fp_in.close();
		fp_out=open(NAPS_path+'download/'+file_name+'_lap.txt','w');
		fp_out.write('Node\tEigenvector centrality\tSecond smallest eigenvector of laplacian');
		for itr in range(len(ssevc)):
			fp_out.write('\r\n'+str(residue[itr])+'\t'+str(alevc[itr])+'\t'+str(ssevc[itr]));
		fp_out.close();
	if(file_type=='edgecentrality'):#For edgecentrality
		fp_in=open(NAPS_path+'edgecentrality/'+name,'r');
		fp_out=open(NAPS_path+'download/'+file_name+'_edgecentrality.txt','w');
		line=fp_in.readline();
		fp_out.write(line);
		for line in fp_in.readlines():
			sp=line.split();
			res_ind=sp[0].split('_')
			fp_out.write(str(residue[int(res_ind[0])])+'_'+str(residue[int(res_ind[1])])+'\t'+str(sp[1])+'\r\n');
		fp_in.close();
		fp_out.close();
