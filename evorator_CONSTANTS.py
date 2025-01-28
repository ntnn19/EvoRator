#!/groups/pupko/modules/python-anaconda3.6.5/bin/python

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'natannag@mail.tau.ac.il' #'orenavram@gmail.com' #'shiranos@gmail.com' 'evolseq@mail.tau.ac.il'
SMTP_SERVER = 'mxout.tau.ac.il'
OWNER_EMAIL = 'natannagar89@gmail.com'

QUEUE_NAME = 'pupkoweb'
#PIPELINE_NAME = 'ModelTeller'

# general paths
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'

WEBSERVER_NAME = 'evorator'
#WEBSERVER_NAME = 'EVORATOR'
WEBSERVER_URL = 'https://evorator.tau.ac.il'
#MODELTELLER_LOG = '/bioseq/modelteller/MODELTELLER_runs.log'
#APACHE_ERROR_LOG = '/var/log/httpd/modelteller.error_log'
CONSURFDB_DB = '/var/www/html/consurfdb/DB'
PDB_DIVIDED = "/bioseq/PDB/data/structures/divided/pdb/"
IDENTICAL_2_UNIQUE_DICT = '/bioseq/mutevorator/identical_to_unique_dict.txt'
RELOAD_INTERVAL = 30
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

HR_STYLE = 'style="height:1px;border:none;color:#333;background-color:#333;"'

EVORATOR_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'evorator')
EVORATOR_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'evorator')
EVORATOR_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')
EVORATOR_HTML_DIR = '/data/www/html/evorator'

# EVORATOR_EXEC = '/groups/pupko/natannag/natan_git/EvoRator'
EVORATOR_EXEC = os.getcwd()

CATALYTIC_RES_DB = os.path.join(EVORATOR_EXEC,'data', 'catalytic_sites.csv')
#TRAINED_SVR = os.path.join(EVORATOR_EXEC, 'regression_results','trained_classifier_keep','estimator_4_ORFan_analysis_keep.pkl')
#TRAINED_XTR = os.path.join(EVORATOR_EXEC, 'results_4_webserver','EvoRator_default_all_features_final.joblib')
TRAINED_SVR = os.path.join(EVORATOR_EXEC, 'results_4_webserver','EvoRator_default_all_features_SVR_final_110921.joblib')
TRAINED_SVR_perfo = os.path.join(EVORATOR_EXEC, 'evorator_model','EvoRator_SVR_final_model_PERFO_281121.joblib')
TRAINED_SVR_perfgr = os.path.join(EVORATOR_EXEC, 'evorator_model','EvoRator_SVR_final_model_PERFGR.joblib')
FITTED_PSSM_PREPROCESSOR = os.path.join(EVORATOR_EXEC, 'trained_pssm_model','preprocessor_pssm_280722.pkl')
TRAINED_PSSM_ANN = os.path.join(EVORATOR_EXEC, 'trained_pssm_model','trained_ANN_pssm_280722.pkl')
#MAIN_SCRIPT = os.path.join(EVORATOR_EXEC, 'evorator.py')
#MAIN_SCRIPT = os.path.join(EVORATOR_EXEC, 'evorator_final.py')
MAIN_SCRIPT = os.path.join(EVORATOR_EXEC, 'evorator_final_backup_270722.py')
CALC_GDV_SCRIPTS = os.path.join(EVORATOR_EXEC, 'calc_gdv')
CONVERT2LEDA_GDV_SCRIPT = os.path.join(CALC_GDV_SCRIPTS, 'convert_edgelist_2_LEDA_and_calc_GDV.py')
COUNT_GDV_SCRIPT = os.path.join(CALC_GDV_SCRIPTS, 'count.py')
ORCA_GDV_SCRIPT = os.path.join(CALC_GDV_SCRIPTS, 'orca','orca.exe')
PSSM_PRED_EXE = os.path.join(EVORATOR_EXEC, 'preprocessor_pssm_backup_290722.py')
HUANG_DATA_PATH = os.path.join(EVORATOR_EXEC, 'predict_huang' ,'huang_with_evorator_features_and_predictions.csv')
HUANG_PSSM_PATH = os.path.join(EVORATOR_EXEC, 'predict_huang', 'huang_aa_variance_pssm.csv')
OUT_PATH_PSSM_AUX= os.path.join(EVORATOR_EXEC, 'trained_pssm_model')
PSSM_NAME_AUX = 280722
NAPS_API = os.path.join(EVORATOR_EXEC, 'api_bulk.py')

RESULT_MSG = 'Unresolved error'


CONTAINER_WIDTH = 'width: 80%'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_FONT = 'font-size: 20px'

CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}; {CONTAINER_FONT}'

PROCESSING_MSG = f'<i>{WEBSERVER_NAME.upper()}</i> is now processing your request. This page will be automatically ' \
    f'updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, ' \
    f'several links to the output files will appear below. '

PROGRESS_BAR_ANCHOR = '''<!--progress_bar_anchor-->'''
PROGRESS_BAR_TAG = '''<div class="progress">
        <div class="progress-bar progress-bar-striped active" role="progressbar" style="width:100%">
        </div>
    </div>'''
