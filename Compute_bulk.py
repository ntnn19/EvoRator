'''Construct network and compute network parameters and node centralities
command line arguments:
1. pdb list in text format
2. analysis type: basic/centrality/clique
3. graph type: ca/cb
4. weight: weighted/unweighted
5. lower cutoff (eg. 0)
6. upper cutoff (eg. 7)
7. residue separation (eg. 1)
'''
import os
import sys
NAPS_path='/var/www/html/NAPS/files_generated/bulkupload/';

# jobid=sys.argv[1].split('.')[0];
jobid='/Compute_bulk_input';
print(jobid)
jobid_filename=sys.argv[1];
print(jobid_filename)
analysis=sys.argv[2];
graph_type=sys.argv[3];
weight=sys.argv[4];
upper_cutoff=sys.argv[5];
lower_cutoff=sys.argv[6];
residue_separation=sys.argv[7];
NAPS_path=sys.argv[8];
print(NAPS_path)
print(NAPS_path+jobid)

if not os.path.exists(NAPS_path+jobid):
	print(31)
	os.makedirs(NAPS_path+jobid);
if not os.path.exists(NAPS_path+jobid+'/pdb'):
	os.makedirs(NAPS_path+jobid+'/pdb');
if not os.path.exists(NAPS_path+jobid+'/edgelist'):
	os.makedirs(NAPS_path+jobid+'/edgelist');
if not os.path.exists(NAPS_path+jobid+'/xyz'):
	os.makedirs(NAPS_path+jobid+'/xyz');
if not os.path.exists(NAPS_path+jobid+'/distance_matrix'):
	os.makedirs(NAPS_path+jobid+'/distance_matrix');
if not os.path.exists(NAPS_path+jobid+'/download'):
	os.makedirs(NAPS_path+jobid+'/download');

print(NAPS_path+jobid_filename)
# fp=open(NAPS_path+jobid_filename);
fp=open(jobid_filename);
for line in fp.readlines():
	sp=line.split();
	pdbid=sp[0];
	chain=sp[1:];
	chain_txt='';
	for ch in chain:
		chain_txt=chain_txt+ch;
	print('cp '+ NAPS_path+"/"+pdbid+'.pdb '+NAPS_path+jobid+'/pdb/')
	os.system('cp '+ NAPS_path+"/"+pdbid+'.pdb '+NAPS_path+jobid+'/pdb/');
	#The following is to generate the edgelist file
	import readpdb_chain;
	if graph_type=="ca":

		o=readpdb_chain.readpdb(NAPS_path+jobid+'/pdb/'+pdbid+'.pdb','CA',chain,NAPS_path+jobid+'/edgelist/'+pdbid+'_'+chain_txt+'_edgelist.lgl',NAPS_path+jobid+'/xyz/'+pdbid+'_'+chain_txt+'_xyz.txt',weight,upper_cutoff,lower_cutoff,residue_separation);
	elif graph_type=="cb":
		o=readpdb_chain.readpdb(NAPS_path+jobid+'/pdb/'+pdbid+'.pdb','CB',chain,NAPS_path+jobid+'/edgelist/'+pdbid+'_'+chain_txt+'_edgelist.lgl',NAPS_path+jobid+'/xyz/'+pdbid+'_'+chain_txt+'_xyz.txt',weight,upper_cutoff,lower_cutoff,residue_separation);
	elif graph_type=="any":
		o=readpdb_chain.readpdb_anyatom(NAPS_path+jobid+'/pdb/'+pdbid+'.pdb',chain,NAPS_path+jobid+'/edgelist/'+pdbid+'_'+chain_txt+'_edgelist.lgl',NAPS_path+jobid+'/xyz/'+pdbid+'_'+chain_txt+'_xyz.txt',weight,upper_cutoff,lower_cutoff,residue_separation);
	elif graph_type=="centroid":
		o=readpdb_chain.readpdb_centroid(NAPS_path+jobid+'/pdb/'+pdbid+'.pdb',chain,NAPS_path+jobid+'/edgelist/'+pdbid+'_'+chain_txt+'_edgelist.lgl',NAPS_path+jobid+'/xyz/'+pdbid+'_'+chain_txt+'_xyz.txt',weight,upper_cutoff,lower_cutoff,residue_separation);
	elif graph_type=="sidechain":
		o=readpdb_chain.readpdb_sidechain(NAPS_path+jobid+'/pdb/'+pdbid+'.pdb',chain,NAPS_path+jobid+'/edgelist/'+pdbid+'_'+chain_txt+'_edgelist.lgl',NAPS_path+jobid+'/xyz/'+pdbid+'_'+chain_txt+'_xyz.txt',weight,upper_cutoff,residue_separation,NAPS_path+jobid+'/distance_matrix/'+pdbid+'_'+chain_txt+'_dm.txt');

	if((analysis=='basic') or (analysis=='centrality')):#Compute the graph properties and centrality measures
		if not os.path.exists(NAPS_path+jobid+'/parameters'):
			os.makedirs(NAPS_path+jobid+'/parameters');
		if not os.path.exists(NAPS_path+jobid+'/centrality'):
			os.makedirs(NAPS_path+jobid+'/centrality');
		if weight=="unweighted":
			os.system('/var/www/html/NAPS/c_codes/executables/centrality.out '+NAPS_path+jobid+'/edgelist/'+pdbid+'_'+chain_txt+'_edgelist.txt '+NAPS_path+jobid+'/centrality/'+pdbid+'_'+chain_txt+'_centrality.txt '+NAPS_path+jobid+'/parameters/'+pdbid+'_'+chain_txt+'_parameters.txt');
		elif weight=="weighted":
			import networkx_centrality_weighted;
			networkx_centrality_weighted.centrality_weighted(pdbid+'_'+chain_txt,NAPS_path+jobid+'/');

	import convert_bulk;
	if(analysis=='edgelist'):
		convert_bulk.convert_bulkfiles(pdbid+'_'+chain_txt+'_edgelist.txt',jobid,NAPS_path);
	elif(analysis=='centrality'):
		convert_bulk.convert_bulkfiles(pdbid+'_'+chain_txt+'_centrality.txt',jobid);

import zipfile;
from os.path import basename

zf = zipfile.ZipFile(NAPS_path+jobid+".zip", "w");
folder_name='download';
if(analysis=='basic'):
	folder_name='parameters';
for dirname, subdirs, files in os.walk(NAPS_path+jobid+'/'+folder_name):
	#zf.write(dirname);
	for filename in files:
		zf.write(os.path.join(dirname, filename), basename(os.path.join(dirname, filename)));
zf.close();
