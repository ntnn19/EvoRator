
# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = "TAU BioSequence \<bioSequence\@tauex.tau.ac.il\>"
ADMIN_USER_NAME = "bioSequence"
ADMIN_PASSWORD = "elana"
# use constant SMTP_SERVER =                "ex1.tau.ac.il"
SMTP_SERVER = "mxout.tau.ac.il"

# the name of the list of all running processes
QUEUING_JOBS = "/bioseq/bioSequence_scripts_and_constants/queuing_jobs.list"
RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/running_jobs.list"
SUBMITTED_JOBS = "/bioseq/bioSequence_scripts_and_constants/submitted_jobs.list"
JOBS_ON_BIOSEQ_NODE = "/bioseq/bioSequence_scripts_and_constants/jobs_on_bioc.01_node.list"
JOBS_WAITING_BIOSEQ_NODE = "/bioseq/bioSequence_scripts_and_constants/jobs_waiting_bioc.01_node.list"
CONSURF_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/consurf_running_jobs.list"
SELECTON_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/selecton_running_jobs.list"
CONSEQ_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/conseq_running_jobs.list"
PEPITOPE_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/pepitope_running_jobs.list"

# Databases urls
PROTEOPEDIA = "http://proteopedia.org/wiki/index.php/"
PDB_DB = "http://www.rcsb.org/pdb/explore/explore.do?structureId="
RCSB_WGET = "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/"
RCSB = "http://www.rcsb.org/"
PISA_WGET = "wget http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/multimer.pdb?"

# CGIs paths
CONSURF_CGI_DIR = "/var/www/cgi-bin/ConSurf"

# general paths
SERVERS_RESULTS_DIR = "/bioseq/data/results/"

SERVERS_LOGS_DIR = "/bioseq/data/logs/"
# use constant SEND_EMAIL_DIR =             "/db1/Local/src/sendEmail" # path on biocluster
SEND_EMAIL_DIR = "/bioseq/bioSequence_scripts_and_constants/sendEmail"
SEND_EMAIL_DIR_IBIS = "/bioseq/bioSequence_scripts_and_constants/sendEmail"  # path on ibis
DAEMON_LOG_FILE = "/bioseq/bioSequence_scripts_and_constants/daemon.log"
UPDATE_RUN_TIME_LOG_FILE = "/bioseq/bioSequence_scripts_and_constants/update_runTime.log"
CONSURF_CGI = "/var/www/cgi-bin/ConSurf"  # on ibis
BIOSEQ_TEMP = "/bioseq/temp/"

# servers urls:
SELECTON_URL = "http://selecton.tau.ac.il"
CONSEQ_URL = "http://conseq.tau.ac.il/"
CONSURF_URL = "http://consurf.tau.ac.il/"
NEW_CONSURF_URL = "http://consurf.tau.ac.il/"  # "http://consurftest.tau.ac.il/"
EPITOPIA_URL = "http://epitopia.tau.ac.il/"
PEPITOPE_URL = "http://pepitope.tau.ac.il/"
QMF_URL = "http://quasimotifinder.tau.ac.il/"
PATCHFINDER_URL = "http://patchfinder.tau.ac.il/"
# use constant FASTML_URL =      "http://ibis.tau.ac.il/fastml/"
# use constant FASTML_URL =      "http://fastmlng.tau.ac.il/"
FASTML_URL = "http://fastml.tau.ac.il/"
CRISTA_URL = "http://crista.tau.ac.il/"
RECONST_URL = "http://fastml.tau.ac.il/reconst/"
GAIN_LOSS_URL = "http://gloome.tau.ac.il/"
CONSURF_DB_URL = "http://consurfdb.tau.ac.il/"
# use constant GILAD_SERVER_URL ="http://consurftest.tau.ac.il/Gilad/"
GILAD_SERVER_URL = "http://mud.tau.ac.il/"
MCPep_URL = "http://bental.tau.ac.il/MCPep/"
ConFind_URL = "http://contemplateng.tau.ac.il/"
CAPTCHER_URL = "http://bental.tau.ac.il/CAPTCHER/"
GUIDANCE_URL = "http://guidance.tau.ac.il/"
GUIDANCE_INDELS_URL = "http://guidance.tau.ac.il/indels/"
SPECBOOST_URL = "http://bental.tau.ac.il/specBoost/"
PROMAYA_URL = "http://bental.tau.ac.il/ProMaya/"
HOMOLOGY_SEARCH_URL = "http://fastml.tau.ac.il/HomologySearch/"
COPAP_URL = "http://copap.tau.ac.il/"

TraitRateProp_URL = "http://traitrate.tau.ac.il/prop/"
# use constant spartaabc_URL=>			"http://spartaabc.tau.ac.il/webserver/"
spartaabc_URL = "http://spartaabc.tau.ac.il/"

ASAP_URL = "http://asap.tau.ac.il/"

# servers logs:
CONSURF_LOG = "/bioseq/ConSurf_old/consurf.log"
CONSURF_NEW_LOG = "/bioseq/ConSurf/consurf.log"
SELECTON_LOG = "/bioseq/Selecton/selecton.log"
EPITOPIA_LOG = "/bioseq/epitopia/epitopia.log"
CONSEQ_LOG = "/bioseq/ConSeq/conseq.log"
PEPITOPE_LOG = "/bioseq/pepitope/pepitope.log"
RECONST_LOG = "/bioseq/ReConst_Server/reconst.log"
MCPep_LOG = "/bioseq/MCPep/mcpep.log"
ConFind_LOG = "/bioseq/ConFind/ConFind.log"
CAPTCHER_LOG = "/bioseq/ConFind/CAPTCHER.log"
Guidance_LOG = "/bioseq/Guidance/guidance.log"
Guidance_Indels_LOG = "/bioseq/GuidanceIndels/guidance_Indels.log"
MuD_LOG = "/bioseq/Gilad_Server/MuD.log"
FASTML_LOG = "/bioseq/FastML/fastml.log"
CRISTA_LOG = "/bioseq/crista/crista.log"
SPECBOOST_LOG = "/bioseq/specBoost/specBoost.log"
GAIN_LOSS_LOG = "/bioseq/GainLoss/GainLoss.log"
PROMAYA_LOG = "/bioseq/ProMaya/ProMaya.log"
COPAP_LOG = "/bioseq/CoPAP/CoPAP.log"
TraitRateProp_LOG = "/bioseq/TraitRate/Prop/server/TraitRateProp_runs.log"
spartaabc_LOG = "/bioseq/spartaabc/server/spartaabc_runs.log"
ASAP_LOG = "/bioseq/asap/ASAP_runs.log"

# use constant COPAP_LOG =     "/groups/pupko/haim/bioseq/CoPAP/logs"

# servers results urls:
# servers urls:
SELECTON_RESULTS_URL = SELECTON_URL + "/results/"

# external databases
# use constant PQS=                 "/bioseq/data/results/PQS/"
PQS = "/bioseq/PQS/"
PDB_DIVIDED = "/bioseq/PDB/data/structures/divided/pdb/"
SWISSPROT_DB = "/bioseq/BLAST/Proteins/swissprot"
UNIPROT_DB = "/bioseq/BLAST/Proteins/uniprot"
CLEAN_UNIPROT_DB = "/bioseq/BLAST/Proteins/clean_uniprot"
UNIREF90_DB = "/bioseq/BLAST/Proteins/uniref90"  # "/groups/bioseq.home/HAIM/UNIREF90/uniref90"

SWISSPROT_DB_FASTA = "/bioseq/FASTA/uniprot/swissprot.fa"
UNIPROT_DB_FASTA = "/bioseq/FASTA/uniprot/trembl.fa"
CLEAN_UNIPROT_DB_FASTA = "/bioseq/FASTA/uniprot/clean_uniprot.fa"
UNIREF90_DB_FASTA = "/bioseq/FASTA/uniprot/uniref90.fa"  # "/groups/bioseq.home/HAIM/UNIREF90/uniref90"
PDBAA_NCBI = "/bioseq/BLAST/Proteins/pdbaa"
CULLED_PDB = "/bioseq/BLAST/Proteins/pdbaaent_dun"  # "/groups/bioseq.home/HAIM/PDBAA/pdbaaent"  # TO CHANGE TO: /bioseq/BLAST/dunbrack.fccc.edu/Guoli/culledpdb/pdbaaent_dun
PDB_DUNBRACK = "/bioseq/BLAST/Proteins/pdbaa_dun"  # "/groups/bioseq.home/HAIM/PDBAA/pdbaa"     # TO CHANGE TO: /bioseq/BLAST/dunbrack.fccc.edu/Guoli/culledpdb/pdbaa_dun
NR_PROT_DB = "/bioseq/BLAST/Proteins/nr"
NR_PROT_DB_FASTA = "/bioseq/FASTA/nr/nr.fa"
NR_NUC_DB = "/bioseq/BLAST/Nucleotides/nt"
NR_NUC_DB_FASTA = "/bioseq/FASTA/nt/nt.fa"
UNIPROT_DAT_INDEX = "/bioseq/data/results/GB_CDS/uniprot.dat.bp_index"
PDB_TO_UNIPROT = "/bioseq/data/results/PDB_to_UNIPROT/idmapping_PDB_UNIPROTKB.dat"  # "/bioseq/idmapping_PDB_UNIPROTKB.dat"
PDB_TO_UNIPROT_test = "/bioseq/idmapping_PDB_UNIPROTKB.dat"
# internal databases
EPITOPIA_DATA = "/bioseq/epitopia/data"

# external programs
BLASTALL = "/opt/bio/ncbi/bin/blastall"  # "/opt/Bio/ncbi/bin/blastall" # on the lecs
BLASTPGP = "blastpgp"  # "/opt/Bio/ncbi/bin/blastpgp" # on the lecs
CS_BLAST = "/share/apps/csblast-2.1.0-linux64/csblast_static"  # on the lecs
# use constant MUSCLE_LECS =        "/share/apps/bin/muscle"  # on the lecs
MUSCLE_LECS = "module load muscle muscle"  # on power
MUSCLE = "/usr/local/bin/muscle"  # on the biocluster
MUSCLE_3_6 = "/bioseq/Programs/muscle_3.6_from_BIOCLUSTER/muscle3.6/muscle"  # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)
# use constant CLUSTALW_LECS =      "/share/apps/bin/clustalw" # on the lecs
CLUSTALW_LECS = "module load clustalw/2.1 clustalw"  # on power
CLUSTALW = "/usr/local/bin/clustalw"  # on the biocluster
CLUSTALW_1_82 = "/bioseq/Programs/ClustalW_1.82/clustalw1.82/clustalw"  # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)
CLUSTALW_1_81 = "/bioseq/Programs/ClustalW_1.81/clustalw1.81/clustalw"  # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)
CLUSTALW_2_0_10 = "/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2"  # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)

MAFFT_LINSI = "/usr/local/bin/mafft-linsi"  # on the biocluster
MAFFT = "/usr/local/bin/mafft"  # on the biocluster
# use constant MAFFT_GUIDANCE =     "/groups/pupko/privmane/bin/mafft" #v6.711b
# use constant MAFFT_LINSI_GUIDANCE =	    "/groups/pupko/privmane/bin/mafft --localpair --maxiterate 1000" #v6.711b
# use constant MAFFT_GUIDANCE =     "/bioseq/Programs/MAFFT_6.711b/mafft" #v6.711b
MAFFT_GUIDANCE = "/bioseq/Programs/MAFFT_6.833/bin/mafft"  # v6.833
# use constant MAFFT_GUIDANCE =      "/bioseq/Programs/MAFFT_6.857/bin/mafft" #v6.857 !!! make sure: 'setenv MAFFT_BINARIES /bioseq/Programs/MAFFT_6.857/mafft-6.857-with-extensions/binaries' BEFORE
# use constant MAFFT_LINSI_GUIDANCE =	    "/bioseq/Programs/MAFFT_6.711b/mafft --localpair --maxiterate 1000" #v6.711b
MAFFT_LINSI_GUIDANCE = "/bioseq/Programs/MAFFT_6.833/bin/mafft --localpair --maxiterate 1000"  # v6.833
# use constant MAFFT_LINSI_GUIDANCE ="/bioseq/Programs/MAFFT_6.857/bin/mafft --localpair --maxiterate 1000" #v6.857 !!! make sure: 'setenv MAFFT_BINARIES /bioseq/Programs/MAFFT_6.857/mafft-6.857-with-extensions/binaries' BEFORE
MAFFT_v7_222 = '/bioseq/Programs/MAFFT_7.222/installation/bin/mafft'  # v7.222 # IMPORTANT: one must run the command: setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}" ahead of this mafft command so all components will be found...
# use constant PRANK_LECS =         "/share/apps/bin/prank" # on the lecs
PRANK_LECS = "prank"  # on power
PRANK = "/usr/local/bin/prank"  # on the biocluster
T_COFFEE = "/share/apps/T-COFFEE-8.47/bin/binaries/linux/t_coffee"  # requiers setenv PATH /share/apps/T-COFFEE-8.47/bin/binaries/linux:$PATH
# use constant PAGAN_LECS =         "/share/apps/pagan-msa/bin/pagan" # requires:  "module load gcc/gcc461" before!!
PAGAN_LECS = "/share/apps/pagan-msa/pagan/bin/pagan"

RNA_FOLD = "/bioseq/Programs/ViennaRNA/ViennaRNA-2.2.0_Installation/bin/RNAfold"
COLOR_RNAFOLD_CONSURF = "perl /bioseq/Programs/ViennaRNA/ViennaRNA-2.2.0_Installation/ViennaRNA/bin/relplot_ConSurf.pl"
COLOR_CBS_RNAFOLD_CONSURF = "perl /bioseq/Programs/ViennaRNA/ViennaRNA-2.2.0_Installation/ViennaRNA/bin/relplot_ConSurf_CBS.pl"
HH_PRED = "perl /bioseq/Programs/HHSuite/hhpred/hhpred.pl"
TREE_VIEWER_DIR = "/bioseq/ConSurf_old/treeViewer/"
PACC_path = "/bioseq/ConSeq/external_scripts/PACC/"
RATE4SITE_BIOC_VER = "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4site.exe"
RATE4SITE_SLOW_BIOC_VER = "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4siteSlow.exe"
RATE4SITE = "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4site.exe"
RATE4SITE_SLOW = "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4siteSlow.exe"
RATE4SITE_SLOW_LECS = "/share/apps/bin/rate4site_slow"
RATE4SITE_LOCAL = "/bioseq/rate4site/rate4site"
RATE4SITE_SLOW_LOCAL = "/bioseq/rate4site/rate4site.doubleRep"
RATE4SITE_WITH_LG = "/bioseq/rate4site/With_LG/rate4site"
RATE4SITE_WITH_LG_SLOW = "/bioseq/rate4site/With_LG/rate4site.doubleRep"
RUBY = "/share/apps/bin/ruby"  # "/usr/bin/ruby"
# use constant CD_HIT_DIR =         "/db1/Local/src/cd-hit_redundency/"
CD_HIT_DIR = "/bioseq/cd_hit/"
PREDICT_PACC = "/bioseq/ConSeq/external_scripts/PACC/run.sh"
MSA_to_HSSP = "/bioseq/ConSeq/external_scripts/PACC/MSA2hssp.pl"
# use constant SEMPHY =             "/groups/pupko/privmane/alignment/run/semphy" #on Biocluster
SEMPHY = "/bioseq/Programs/Semphy/semphy.doubleRep"

# internal programs
EPITOPIA_EXECUTABLES = "/bioseq/epitopia/executables"

# constant values
BLAST_MAX_HOMOLOGUES_TO_DISPLAY = 500
BLAST_PDB_MAX_HOMOLOGUES_TO_DISPLAY = 25
CONSURF_PIPE_FORM = "/bioseq/ConSurf_old/consurf_pipe.form"
SELECTON_MAX_NUCLEOTIDE = 15000
MAX_WALLTIME = "96:00:00"

# Queue Details
BIOSEQ_NODE = "bioc01.tau.ac.il"  # Node on BioCluster dedicated to Bioseq runs (Not part of the queue)
# use constant MAX_QUEUE_RUNS =         60
MAX_QUEUE_RUNS = 999

# external links
RCSB_WEB = "http://www.rcsb.org/"
PYMOL_WEB = "http://pymol.sourceforge.net/"
CHIMERA_WEB = 'http://www.rbvi.ucsf.edu/chimera/'
CHIMERA_SAVING_FIGURE = 'http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/print.html'
CHIMERA_DOWNLOAD = CHIMERA_WEB + "download.html"
MSA_CONVERT = 'http://www.ebi.ac.uk/Tools/sfc/'
MSA_FORMATS = 'http://consurf.tau.ac.il/2016/quick_help.php#UserMSA'

# redirect pages
CONSURF_REDIRECT_PAGE = CONSURF_URL+ "too_many_runs.html"
SELECTON_REDIRECT_PAGE = SELECTON_URL + "/too_many_runs.html"
CONSEQ_REDIRECT_PAGE = CONSEQ_URL + "too_many_runs.html"
PEPITOPE_REDIRECT_PAGE = PEPITOPE_URL + "too_many_runs.html"

# faq pages
CONSURF_TREE_FAQ = CONSURF_URL + 'quick_help.html#note5'

# Files Name Conventions
TEMPLATES_LIST_FILE = "List_of_Templates"
PISA_ERRORS_FILE = "PISA_Errors"

SPECIFIC_PDB_CHAIN_EXE = "/bioseq/pepitope/Exec/specificPdbChain"
PEPSURF_EXE = "/bioseq/pepitope/Exec/pepSurf.exe"
MAPITOPE_EXE = "/bioseq/pepitope/Exec/mapitope"
MAPITOPE_EXE_FOR_TESTS = "/bioseq/pepitope/Exec/mapitope_test_version_09Aug11"
NACCESS = "/bioseq/Programs/naccess2.1.1/naccess"
SURFRACER_DIR = "/bioseq/pepitope/"
SURFRACER = "surfrace5_0_linux_64bit"

# ---------------------------------------------
# sub
# print_to_output
# {
#     my $OutHtmlFile = shift
# my $server_name = shift
# my $run_name = shift
# my $recipient = shift
#
# open
# OUTPUT, ">>$OutHtmlFile"
# # flock OUTPUT, 2
# print
# OUTPUT
# "\n<p><font size=+3 color='red'>ERROR! $server_name session has been terminated: </font>\n<br><b>A system error occured during the calculation. Please try to run $server_name again in a few minutes.</b>\n</p>\n"
# print
# OUTPUT
# "<H3><center>For assistance please <a href=\"mailto:".ADMIN_EMAIL.
# "?subject=".$server_name.
# "%20Run%20No:%20".$run_name.
# "\">contact us</aand mention this number: $run_name</H3>\n"
# # flock OUTPUT, 8
# close
# OUTPUT
# & send_mail($server_name, $recipient, $run_name, "error", "error") if ($recipient ne "NO")
# & stop_reload($OutHtmlFile)
# }
# # ---------------------------------------------
#
# # in case the desired mail report on error: the vars $email_subject and $email_message should be 'error'
# sub
# send_mail
# {  # to user
#     my $server_name = shift
# my $recipient = shift
# my $run_name = shift
# my $email_subject = shift
# my $email_message = shift
# my $email_attach = shift
# my $from_server = ""
# $from_server = shift
# my $OutputURL
# my $mail
#
# if ($server_name eq "Selecton")
# {$OutputURL = SELECTON_URL.
# "/results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "ConSeq") {$OutputURL = CONSEQ_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "Epitopia") {$OutputURL = EPITOPIA_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "pepitope") {$OutputURL = PEPITOPE_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "ConSurf") {$OutputURL = CONSURF_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "QuasiMotiFinder") {$OutputURL = QMF_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "fastml") {$OutputURL = FASTML_URL.
# "results/$run_name".
# "/output.html"}
#
# $email_subject = "Error in $server_name running" if $email_subject
# eq
# "error"
# $email_message = "Hello!\n\nUnfortunately there was an error while running the $server_name server.\nPlease click on the following link to see more details\nWe apologize for the inconvenience\n\n$OutputURL\n" if $email_message
# eq
# "error"
# chdir
# SEND_EMAIL_DIR
# chdir
# SEND_EMAIL_DIR_IBIS if ($from_server eq "ibis")
# $mail = 'perl sendEmail.pl -f \''.ADMIN_EMAIL.
# '\' -t \''.$recipient.
# '\' -u \''.$email_subject.
# '\' -s '.SMTP_SERVER.
# ' -m \''.$email_message.
# "\'"
# # $mail ='perl sendEmail.pl -f \''.ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u \''.$email_subject.'\' -xu '.ADMIN_USER_NAME.' -xp '.ADMIN_PASSWORD.' -s '.SMTP_SERVER.' -m \''.$email_message."\'"
# if ($email_attach ne '')
# {$mail. = " -a $email_attach"}
# `$mail
# `
# }
# # ---------------------------------------------
# sub
# stop_reload
# {
#     my $OutHtmlFile = shift
#
# sleep
# 10
# open
# OUTPUT, "<$OutHtmlFile"
# my @ output = < OUTPUT >
# close
# OUTPUT
# open
# OUTPUT, ">$OutHtmlFile"
# foreach
# my $line( @ output){  # we remove the refresh lines and the button which codes for Selecton cancelled job
#     unless($line = ~ / REFRESH / i or $line = ~ / NO - CACHE / i or $line = ~ / ACTION =\"\/cgi\/kill_process.cgi/ or
# $line = ~ / VALUE =\"Cancel Selecton Job\"/ or $line =~ /TYPE=hidden NAME=\"pid\"/ or
# $line = ~ / TYPE = hidden
# NAME =\"selecton_http\"/ or $line =~ /TYPE=hidden NAME=\"run_no\"/ or
# $line = ~ / <!--job_ /){
#     print
# OUTPUT $line
# }
# }
# close
# OUTPUT
# }
# ---------------------------------------------
# sub
# print_Q_status_in_html
# {
#     my $html_file = shift
# my $_status = shift
# my $_time = shift
# my $_estimated_run_time = shift
#
# my($line, $line1, $line2)
# my $out = "/bioseq/ELANA/from_GENERAL_CONST.txt"
#
# $_time = "" if ($_time eq "no")
# unless(open
# HTML, "+>>".$html_file) {
# return "print_Q_status_in_html : Could not open file $html_file to update the status. Status is: $_status  reason: $!\n"}
# else {
# # flock HTML, 2
# seek
# HTML, 0, 0  # rewind the pointer to the beginning
# my @ html_lines = < HTML >  # read the contents into the array
# truncate
# HTML, 0  # remove all the information, The 0 represents the size of the file that we want
# foreach( @ html_lines){
# if (/ < !--job_stat--<.+ Your job status is: < \ / a (.+) < br / ){
# if ($_status ne "")
# {
#     s /$1 /$_status /
# }
# }
# elsif( / <!--job_pass -->The
# time
# that
# passed
# since
# submitting
# the
# query is: (.+) < br / ){
# if ($_time ne "")
# {
#     s /$1 /$_time /
# }
# }
# elsif( / <!--(job_time - -)
# Estimated
# run
# time is: (-->) / and $_estimated_run_time
# ne
# "none"){
# $line = $_
# $line1 = $1
# $line2 = $2
# if ($_estimated_run_time =~ m / \d+:\d+:\d+:\d+ /)
# {
# $_estimated_run_time. = " days"
# }
# elsif($_estimated_run_time = ~ m /\d +:\d +:\d + / ) {
# $_estimated_run_time. = " hours"
# }
# elsif($_estimated_run_time = ~ m /\d +:\d + / ){
# $_estimated_run_time. = " minutes"
# }
# $_ = $line  # since we make another RE comparison, the original values of $_ and $1 are changing, therefore we must save them at the beginning and change them back here.
# s /$line2 /$_estimated_run_time < br /  # the reason we first substitue the second part, is that the first part creates an expression --which might be wrongly replaced with this value
# s /$line1 /$line1 /
# }
# }
# print
# HTML $_
# foreach( @ html_lines)
# # flock HTML, 8
# close
# HTML
# return "OK"
# }
# }
#
#
# # in case the desired mail report on error: the vars $email_subject and $email_message should be 'error'
# sub
# send_mail2
# {  # to user
#     my $server_name = shift
# my $recipient = shift
# my $run_name = shift
# my $email_subject = shift
# my $email_message = shift
# my $email_attach = shift
# my $from_server = shift
# my $OutputURL
# my $mail
#
# if ($server_name eq "Selecton")
# {$OutputURL = SELECTON_URL.
# "/results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "ConSeq") {$OutputURL = CONSEQ_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "Epitopia") {$OutputURL = EPITOPIA_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "pepitope") {$OutputURL = PEPITOPE_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "ConSurf") {$OutputURL = CONSURF_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "QuasiMotiFinder") {$OutputURL = QMF_URL.
# "results/$run_name".
# "/output.html"}
# elsif($server_name
# eq
# "fastml") {$OutputURL = FASTML_URL.
# "results/$run_name".
# "/output.html"}
#
# $email_subject = "Error in $server_name running" if $email_subject
# eq
# "error"
# $email_message = "Hello!\n\nUnfortunately there was an error while running the $server_name server.\nPlease click on the following link to see more details\nWe apologize for the inconvenience\n\n$OutputURL\n" if $email_message
# eq
# "error"
# chdir
# SEND_EMAIL_DIR
# chdir
# SEND_EMAIL_DIR_IBIS if ($from_server eq "ibis")
# $mail = 'perl sendEmail.pl -f \''.ADMIN_EMAIL.
# '\' -t \''.$recipient.
# '\' -u \''.$email_subject.
# '\' -s '.SMTP_SERVER.
# ' -m \''.$email_message.
# "\'"
# # $mail ='perl sendEmail.pl -f \''.ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u \''.$email_subject.'\' -xu '.ADMIN_USER_NAME.' -xp '.ADMIN_PASSWORD.' -s '.SMTP_SERVER.' -m \''.$email_message."\'"
# if ($email_attach ne '')
# {$mail. = " -a $email_attach"}
# $mail = 'sh -c \' $mail 2>/dev/null\''
# `$mail
# `
# }
