import click
# define Python user-defined exceptions
def NNModel_class(input_dim=None, num_classes=2,n_hidden=147,vote=False):
    X_input = Input(input_dim)
    X = Dense(input_dim[0],kernel_initializer='normal', activation='relu', name='fc1')(X_input)
    X = Dense(n_hidden, activation='relu', name='fc2')(X)
    # X = Dense(60, activation='relu', name='fc3')(X)
    # X = Dense(30, activation='relu', name='fc4')(X)
    # X = Dense(num_classes, activation='sigmoid', name='fc_out')(X)
    X = Dense(num_classes, activation='softmax', name='fc_out')(X) # as in preprint
    # X = Dense(num_classes,kernel_initializer='normal', name='fc_out')(X)
    if vote:
        X =Flatten()(X)

    # X = Dense(num_classes, activation='relu', name='fc_out')(X_input)
    model = Model(inputs=X_input, outputs=X)
    # model.compile('adam', 'binary_crossentropy', metrics=['binary_crossentropy'])
    model.compile('adam', 'kullback_leibler_divergence', metrics=['kullback_leibler_divergence']) # as in preprint
    # model.compile('adam', 'mean_squared_error', metrics=['mean_squared_error'])
    return model


def get_consurf_dfs(consurfdb_query, out_path):
        main_link = f"https://consurfdb.tau.ac.il/DB/{consurfdb_query}/consurf_summary.txt"
        # print(main_link)
        url_content = urllib.request.urlopen(main_link).read()
        soup = BeautifulSoup(url_content)
        f=open(out_path,'w')
        f.write(soup.get_text())
        f.close()

class Error(Exception):
    """Base class for other exceptions"""
    pass


class ChainNotFound(Error):
    """Raised when the chain was not found for the input structure"""
    pass

class QueryNotFound(Error):
    """Raised when the structure was not found in ConSurf-DB"""
    pass

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def finalize_html(html_path, error_path, run_number,job_title,results_dir,consurf_output=False,r2_test=None, task='',done_figure_path=''):
    succeeded = not os.path.exists(error_path)
    logging.debug(f'SUCCEEDED = {succeeded}')
    consurf_output_flg = True if consurf_output else False
    r2_test = str(round(r2_test,2)) if consurf_output_flg and task=='PERfIFR' else ''
    if succeeded and consurf_output_flg:
        if task=='PERfIFR':
            edit_success_html(CONSTS, html_path,job_title,results_dir,consurf_output_flg,r2_test,task,done_figure_path)
        elif task=='PERfGR':
            edit_success_html(CONSTS, html_path,job_title,results_dir,consurf_output_flg,task=task)
    elif succeeded and not consurf_output_flg:
        edit_success_html(CONSTS, html_path,job_title,results_dir)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def edit_success_html(CONSTS, html_path, job_title,results_dir,consurf_output=False,r2_test='',task='',done_figure_path=''):

    update_html(html_path, 'RUNNING', 'FINISHED')
    # href = "https://consurf.tau.ac.il/ngl/viewer.php?job=1627817377"
    if consurf_output:
        if task=='PERfIFR':

            append_to_html(html_path, f'''
                <br>
                <h4 class="output_title">Final Results</h4>
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'>
                <li>
                R squared: {r2_test} 
                </li>
                <img id="regression" src="{job_title}_regression.png">
                <br>
                <br>
                <li>
                <a href='{job_title}_features_and_predictions.csv' target='_blank'>Download features and predictions table</a>
                </li>
                <li>
                    View EvoRator results<a href='/ngl/viewer.php?job={os.path.split(results_dir)[-1]}' target='_blank'>
                    <b> with NGL viewer  </b>
                    </a>
                     
                </li>
                <li>
                    View EvoRator results <a href="network.html">
                    <b> projected onto a network representation of your protein </b>
                    </a> 
                     
                </li>
                </div>
    ''')
        elif task=='PERfGR':
            append_to_html(html_path, f'''
                            <br>
                            <h4 class="output_title">Final Results</h4>
                            <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'>
                            <li>
                            <a href='{job_title}_features_and_predictions.csv' target='_blank'>Download features and predictions table</a>
                            </li>
                            <li>
                    View EvoRator results<a href='/ngl/viewer.php?job={os.path.split(results_dir)[-1]}' target='_blank'>
                    <b> with NGL viewer  </b>
                    </a>
                     
                </li>
<li>
                    View EvoRator results <a href="network.html">
                    <b> projected onto a network representation of your protein </b>
                    </a>

                </li>
                            </div>
                ''')

    else:
        append_to_html(html_path, f'''
            <br>
            <h4 class="output_title">Final Results</h4>
            <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'>
            <li>
            <a href='{job_title}_features_and_predictions.csv' target='_blank'>Download features and predictions table</a>
            </li>
            <li>
                    View EvoRator results<a href='/ngl/viewer.php?job={os.path.split(results_dir)[-1]}' target='_blank'>
                    <b> with NGL viewer  </b>
                    </a>
                     
                </li>
<li>
                    View EvoRator results <a href="network.html">
                    <b> projected onto a network representation of your protein </b>
                    </a>

                </li>
            </div>
''')


def edit_failure_html(CONSTS, error_msg, html_path, run_number):
    update_html(html_path, 'RUNNING', 'FAILED')
    append_to_html(html_path,
                   f'<div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3>\n'
                   f'<font color="red">{error_msg}</font></h3><br><br>'
                   f'Please make sure your input is OK and then try to re-run your job or '
                   f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME}%20Run%20Number:%20{run_number}">'
                   f'contact us'
                   f'</a> '
                   f'for further information.<br>'
                   f'</div>\n')


def add_closing_html_tags(html_path, CONSTS, run_number):
    FORMER_MSG = 'EvoRator is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, the output will appear below.'
    update_html(html_path, FORMER_MSG, '')  # remove "web server is now processing your request" message
    update_html(html_path, 'progress-bar-striped active', 'progress-bar-striped')  # stop_progress_bar

    append_to_html(html_path, f'''
            <hr>
                <h4 class=footer>
                    <p align='center'>Questions and comments are welcome! Please
                        <span class="admin_link"> 
                        <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> 
                        </span>
                    </p>
                </h4>
                <div id="bottom_links" align="center"><span class="bottom_link">
                <a href="{CONSTS.WEBSERVER_URL}" target="_blank">Home</a>  
            </span>
        </div>
        <br><br><br>
    </body>
</html>''')

    # have to be last thing that is done
    sleep(2 * CONSTS.RELOAD_INTERVAL)
    update_html(html_path, CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')  # stop refresh


def initialize_html(CONSTS, output_dir_path, html_path):
    path_tokens = output_dir_path.split('/')
    # e.g., "/bioseq/data/results/sincopa/12345678/outputs"
    run_number = path_tokens[path_tokens.index(CONSTS.WEBSERVER_NAME) + 1]

    update_html(html_path, 'QUEUED', 'RUNNING')
    # update_html(html_path, CONSTS.PROGRESS_BAR_ANCHOR, CONSTS.PROGRESS_BAR_TAG)  # add progress bar

    return run_number

@click.command()
@click.argument('pdb_name', type=str)
@click.argument('pdb_file', type=click.Path())
@click.argument('pdb_chain', type=str)
@click.argument('catalytic_sites',type=click.Path(exists=True))
@click.argument('results_dir', type=click.Path(exists=True))
@click.option('--orphan-prediction',type=bool,default=True,show_default=True,help='If False then delta(Consurf_scores,predictions)'
                                                                                  'will be projected onto the structure.'
                                                                                  'For positions for which delta does not exist (i.e. indel)'
                                                                                  'only predictions will be projected'
                                                                                  'If True, predictions will be projected')
@click.option('--trained-regressor',type=click.Path(),help='If orphan_prediction is True then the path trained regressor is required')
@click.option('--consurfdb-query',type=str,default='',show_default=True,help='Provide consurf output file (scores table)')
@click.option('--consurf-output',type=str,default='',show_default=True,help='Provide consurf output file (scores table)')
@click.option('--job-title',type=str,default='',show_default=True,help='Insert job title')
@click.option('--html-path',type=str,default='',show_default=True)
@click.option('--predict',type=bool,default=True,hidden=True)
@click.option('--prediction-task',type=str)
def main(pdb_name,pdb_file,pdb_chain,catalytic_sites,results_dir, orphan_prediction, trained_regressor,
         consurf_output,consurfdb_query,job_title,html_path,predict,prediction_task):
    print('semek')
    file_path = os.path.realpath(__file__)
    print(f'Hi, {file_path}')  # Press Ctrl+F8 to toggle the breakpoint.
    # try:
    if not html_path:
        print('semek')
        if pdb_name:
            subdir_name = pdb_name + '_' + pdb_chain
        else:
            subdir_name = os.path.split(pdb_file)[-1].split(".")[0] + '_' + pdb_chain

        results_dir = os.path.join(results_dir, subdir_name)
        os.makedirs(results_dir, exist_ok=True)

        print(results_dir)
        print('semek')
    else:
        run_number = initialize_html(CONSTS, results_dir, html_path)
        # final_zip_path = f'{os.path.split(results_dir)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'

        os.makedirs(results_dir, exist_ok=True)

    print('semek')
    logging.basicConfig(filename=results_dir + '/' + 'feature_extraction_and_prediction.log',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logging.debug(f'FINAL EVORATOR')
    scripts_dir = CONSTS.EVORATOR_EXEC
    error_path = f'{results_dir}/error.txt'

    logging.debug(f'sklearn version = {sklearn.__version__}')

    if consurf_output:
        with open(consurf_output,'r') as f:
            l = f.readline().strip()
            if not l.startswith('Amino Acid Conservation Scores'):
                edit_failure_html(CONSTS,
                                  "The provided ConSurf table of rates is in the wrong format or the wrong input was provided. An example input that has the correct format can be found <a href=https://evorator.tau.ac.il/3LZMA_consurf_summary.txt target=_blank>here</a>.",
                                  html_path, run_number)
                add_closing_html_tags(html_path, CONSTS, run_number)

    if consurfdb_query:
        consurf_db_file = os.path.join(CONSTS.CONSURFDB_DB, consurfdb_query, f'{consurfdb_query}_consurf_summary.txt')
        logging.debug(f'consurfdb file={consurf_db_file}')
        consurf_output = os.path.join(results_dir, f"{consurfdb_query}_consurf_summary.txt")
        Identical2UniqueDict = open(CONSTS.IDENTICAL_2_UNIQUE_DICT, 'r')
        try:
            cmd = f'ssh bioseq@powerweb1 cp {consurf_db_file} {results_dir}'

            logging.debug(f'copying {consurf_db_file}: {cmd}')
            subprocess.check_output(cmd, shell=True)
        except Exception as e1:
            for line in Identical2UniqueDict:
                if consurfdb_query in line:
                    unique_consurfdb_query = line.strip().split(":")[1]
                    break
            try:
                consurf_db_file = os.path.join(CONSTS.CONSURFDB_DB, unique_consurfdb_query, f'{unique_consurfdb_query}_consurf_summary.txt')
                cmd = f'ssh bioseq@powerweb1 cp {consurf_db_file} {results_dir}'
                logging.debug(f'copying {consurf_db_file}: {cmd}')
                subprocess.check_output(cmd, shell=True)
                consurf_output = os.path.join(results_dir, f"{unique_consurfdb_query}_consurf_summary.txt")

            except Exception as e2:
                logging.debug(f'SUCCEEDED = False')
                logging.debug(e2)
                if html_path:
                    edit_failure_html(CONSTS,
                                      "Your query was not found in ConSurf-DB. The required table of rates can be obtained by running <a href=https://consurf.tau.ac.il target=_blank>ConSurf</a> for your structure of interest",
                                      html_path, run_number)
                    add_closing_html_tags(html_path, CONSTS, run_number)
                    return



    logging.debug(f'174')
    identifier_for_scannet_obtain_pdb_routine = pdb_name + "_0-" + pdb_chain

    if pdb_name:
        pdb_file_name = "pdb" + pdb_name.lower() + ".ent.gz"
        local_gz_file_path = os.path.join(CONSTS.PDB_DIVIDED, pdb_name.lower()[1:3], pdb_file_name)  # pdb file absolute
        final_pdb_file_path = os.path.join(results_dir, pdb_name.upper() + ".pdb")
        if not os.path.exists(local_gz_file_path):
            try:
                print(301)

                pdb_file_single_chain_path, chain_ids, sequence_from_pdb, residue_pdb_index, backbone_coordinates = obtain_pdb.obtain_pdb(identifier_for_scannet_obtain_pdb_routine,results_dir)

                print(pdb_file_single_chain_path)
        # edgelist_path = get_edge_list(identifier_for_scannet_obtain_pdb_routine,chain_ids,backbone_coordinates,residue_pdb_index,results_dir)
        # print(edgelist_path)
        # print(edgelist_path)
        # pdb_file_name = "pdb"+pdb_name.lower()+".ent.gz"
        # local_gz_file_path =  os.path.join(PDB_DIVIDED,pdb_name.lower()[1:3],pdb_file_name)  # pdb file absolute
        # final_pdb_file_path = os.path.join(results_dir,pdb_name.upper()+".pdb")

        # if not os.path.exists(local_gz_file_path):
        #     try:
        #         pdb_query = f'wget -P {results_dir} ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/{pdb_file_name}'
        #         subprocess.check_output(pdb_query,shell=True)
        #         downloaded_gz_file_path = os.path.join(results_dir, pdb_file_name)
        #         local_gz_file_path = downloaded_gz_file_path
        #         subprocess.check_output(f'gunzip {local_gz_file_path}', shell=True)
        #         local_pdb_file_path = os.path.join(*os.path.split(local_gz_file_path)[:-1],
        #                                            ".".join(os.path.split(local_gz_file_path)[-1].split(".")[:2]))
        #     except Exception as e:
        #         logging.debug(f'SUCCEEDED = False')
        #         logging.debug(e)
        #         if html_path:
        #             edit_failure_html(CONSTS,
        #                               "The PDB ID was not found",
        #                               html_path, run_number)
        #             add_closing_html_tags(html_path, CONSTS, run_number)

        # else:
        #     subprocess.check_output(f'cp {local_gz_file_path} {results_dir}', shell=True)
        #     pdb_file_unzipped_path = os.path.join(results_dir, pdb_file_name)
        #     local_pdb_file_path = os.path.join(results_dir,".".join(pdb_file_name.split(".")[:2]))
        #     if not os.path.exists(local_pdb_file_path): # file exists
        #         subprocess.check_output(f'gunzip {pdb_file_unzipped_path}', shell=True)
        #
        # subprocess.check_output(f'cp {local_pdb_file_path} {final_pdb_file_path}',shell=True)
                local_pdb_file_path = pdb_file_single_chain_path

                shutil.copyfile(local_pdb_file_path, final_pdb_file_path)

                print(340,local_pdb_file_path)

            except Exception as e:
                logging.debug(f'SUCCEEDED = False')
                logging.debug(e)
                if html_path:
                    edit_failure_html(CONSTS,
                                      "The PDB ID was not found",
                                      html_path, run_number)
                    add_closing_html_tags(html_path, CONSTS, run_number)
        else:
            logging.debug(f'found in local db')
            logging.debug(f'copy cmd = cp {local_gz_file_path} {results_dir}')
            subprocess.check_output(f'cp {local_gz_file_path} {results_dir}', shell=True)
            pdb_file_unzipped_path = os.path.join(results_dir, pdb_file_name)
            local_pdb_file_path = os.path.join(results_dir, ".".join(pdb_file_name.split(".")[:2]))
            if not os.path.exists(final_pdb_file_path):  # file exists
                logging.debug(f'gunzip cmd = cp {local_gz_file_path} {results_dir}')
                subprocess.check_output(f'gunzip {pdb_file_unzipped_path}', shell=True)
                shutil.copyfile(local_pdb_file_path, final_pdb_file_path)
            chain_ids = [(0, pdb_chain)]
            PDBio.extract_chains(final_pdb_file_path,chain_ids , final_pdb_file_path )
            chains = PDBio.load_chains(file=final_pdb_file_path, chain_ids=chain_ids)[1]
            sequence_from_pdb = PDB_processing.process_chain(chains)[0]

            backbone_coordinates = PDB_processing.process_chain(chains)[2]
            residue_pdb_index = PDB_processing.get_PDB_indices(chains, return_model=True, return_chain=True)
        pdb_input = final_pdb_file_path


        # subprocess.check_output(f'cp {local_pdb_file_path} {final_pdb_file_path}', shell=True)

        # PDBio.extract_chains('dowbloaded', [(0,pdb_chain)], final_pdb_file_path)
        # pdb_input = final_pdb_file_path
        # pdb_input = final_pdb_file_path
        # print(356)
        # _, chain_obj = PDBio.load_chains(file=pdb_input, chain_ids= [ (0, pdb_chain)])
        # print(356)
        # pdb_file_single_chain_path, chain_ids = PDBio.extract_chains(‘my_file.cif’,  [ (0, ‘A’)] , ‘final_file.pdb’ )
        # pdb_file_single_chain_path, chain_ids, sequence, residue_pdb_index, backbone_coordinates = PDB_processing.process_chain(
        #     chain_obj)

    else:
        print(results_dir)
        pdb_input = pdb_file
        final_pdb_file_path = os.path.join(results_dir,subdir_name+'.pdb')
        print(final_pdb_file_path)
        chains = PDBio.load_chains(file=pdb_input, chain_ids=[(0, pdb_chain)])[1]
        sequence_from_pdb = PDB_processing.process_chain(chains)[0]
        backbone_coordinates = PDB_processing.process_chain(chains)[2]
        residue_pdb_index = PDB_processing.get_PDB_indices(chains, return_model=True, return_chain=True)
        PDBio.extract_chains(pdb_input, [(0, pdb_chain)], final_pdb_file_path)
        print(final_pdb_file_path)
        print(382)
        try:
            os.path.exists(pdb_input)
        except:
            raise FileNotFoundError
        try:
        #    cmd = f'module load python/python-anaconda3.7-itaym; python /bioseq/evorator/auxiliaries/pdb_validate.py {pdb_input} > {results_dir}/pdb_val.output'
        #     cmd = f'module load python/python-anaconda3.7-itaym; source activate /groups/pupko/natannag/conda/envs/NatanEnv; python /bioseq/evorator/auxiliaries/pdb_validate.py {pdb_input} > {results_dir}/pdb_val.output'
            cmd = f'/groups/pupko/natannag/conda/envs/NatanEnv/bin/python /bioseq/evorator/auxiliaries/pdb_validate.py {pdb_input} > {results_dir}/pdb_val.output'
            subprocess.check_output(cmd, shell=True)
            pdb_val_out =  os.path.join(results_dir,'pdb_val.output')
            pdb_val_out_f = open(pdb_val_out,'r')
            content= pdb_val_out_f.read()
            pdb_val_out_f.close()
            if "OK" not in content:
                logging.debug('unusual pdb file')

        except:
            logging.debug(f'pdb_val.py failed')
    # except Exception as e:
    #     logging.debug(f'SUCCEEDED = False')
    #     logging.debug(e)
    #     pdb_val_out =  os.path.join(results_dir,'pdb_val.output')
    #     pdb_val_out_f = open(pdb_val_out,'r')
    #     content= pdb_val_out_f.read()
    #     pdb_val_out_f.close()
    #     if "OK" not in content:
    #         edit_failure_html(CONSTS,
    #                           "PDB File was not found or is not in the correct <a href=http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM target=_blank>format</a>.",
    #                           html_path, run_number)
    #         add_closing_html_tags(html_path, CONSTS, run_number)


        # get NAPS  input

    # io = PDBIO()
    # pdb_parser = PDBParser().get_structure(os.path.split(pdb_input)[-1].split(".")[0].upper(), pdb_input)
    # chain_list = [c.get_id() for c in pdb_parser.get_chains()]
    # logging.debug(f'chain_list:{chain_list}')
    #
    # pdb_chain = pdb_chain.upper()
    # for c in pdb_parser.get_chains():
    #     logging.debug(f'{c.get_id()}\n\n')
    #     if c.get_id() == pdb_chain:
    #
    #         io.set_structure(c)
    #         logging.debug(f'{pdb_parser.get_full_id()}\n\n')
    #         pdb_file_single_chain_path =os.path.join(results_dir,str(pdb_parser.get_id()) + "_" + c.get_id() + ".pdb")
            #pdb_file_single_chain_path =os.path.join(results_dir,str(pdb_parser.get_id()).upper() + ".pdb")
            # try:
            #     io.save(pdb_file_single_chain_path)
            # except:
            #     subprocess.check_output(f'cp {pdb_input} {pdb_file_single_chain_path}',shell=True)

            # pdb_parser = PDBParser().get_structure(os.path.split(pdb_input)[-1].split(".")[0].upper(), pdb_input)
            # f = open(pdb_input,'r')
            # line_to_take= [line for line in f.readlines() if line.startswith("CRYST1") or line.startswith("SITE     ") ]
            # f.close()
            #
            # line_to_take= "\n".join(line_to_take)
            #
            # f = open(pdb_file_single_chain_path,'r')
            # f_content = f.read()
            # f.close()
            #
            # f= open(pdb_file_single_chain_path,'w')
            # f.write(line_to_take+f_content)
            # f.close()
            #
            # break
        # else:
        #     logging.debug(f'No chains found for {os.path.split(pdb_input)[-1].split(".")[0].upper()}\n\n')
        #     continue

    logging.debug(f'backbone_coordinates={backbone_coordinates}')
    logging.debug(f'getting edge_list_file')
    # identifier_edgelist_calculator = pdb_name+'_0-'+pdb_chain
    # edgelist_file = extract_pairwise_distances.get_edge_list(identifier_edgelist_calculator,results_dir)
    # logging.debug(f'edge_list_file:{edgelist_file}')
    # Compute_bulk_input =os.path.join(results_dir,'Compute_bulk_input.txt')
    #Compute_bulk_output =os.path.join(results_dir,'Compute_bulk_output.zip')
    # f = open(Compute_bulk_input,'w')
    # f.write(" ".join(os.path.split(pdb_file_single_chain_path)[-1].split(".")[0].split("_")))
    # f.close()
    # subprocess.check_output(f'module list 2> {os.path.join(results_dir, "test.txt")}', shell=True)
    # logging.debug(f'{results_dir}')
#    cmd = f'module load python/python-anaconda3.7-itaym; python {os.path.join(scripts_dir,"Compute_bulk.py")} {Compute_bulk_input} edgelist ca unweighted 7 0 1 {results_dir}'
#     cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv; python {os.path.join(scripts_dir,"Compute_bulk.py")} {Compute_bulk_input} edgelist ca unweighted 7 0 1 {results_dir}'
#     cmd = f'/groups/pupko/natannag/conda/envs/NatanEnv/bin/python {os.path.join(scripts_dir,"Compute_bulk.py")} {Compute_bulk_input} edgelist ca unweighted 7 0 1 {results_dir}'
    # cmd = f'python {os.path.join(scripts_dir,"Compute_bulk.py")} {Compute_bulk_input} edgelist ca unweighted 7 0 1 {results_dir}'
    # logging.debug(f'{cmd}')
    # subprocess.check_output(cmd,shell=True)
    #subprocess.check_output(f'unzip {Compute_bulk_output}',shell=True)
    edgelist_file = extract_pairwise_distances.get_edge_list(identifier_for_scannet_obtain_pdb_routine,chain_ids,backbone_coordinates,residue_pdb_index,results_dir)
    logging.debug(f'{edgelist_file}')


    # Feature extraction
    done_path = os.path.join(results_dir, job_title + "_features_and_predictions.csv")
    done_path_neigh = os.path.join(results_dir, job_title + "_neighbor_map.csv")
    done_path_features_only = os.path.join(results_dir, job_title + "_features.csv")
    done_figure_path = os.path.join(results_dir, job_title + "_regression.png")
    df_merged = extract_features_4_a_single_query.extract_features(edgelist_file, pdb_input,pdb_input, pdb_chain, catalytic_sites,results_dir, consurf_output=consurf_output,job_title=job_title)
    logging.debug(f'raw feature table {df_merged.columns.tolist()}')
    if not predict:
        n2v_flag = False
        if os.path.exists(done_path_features_only) and n2v_flag:
            print(done_path_features_only)
            print("Already exists")
            return
        df_merged['glycosylation'] = df_merged['glycosylation'].str.replace('None glycosylated', 'Non glycosylated')
        df_merged['protein_interaction'] = df_merged['protein_interaction'].replace('missing', np.nan)
        df_merged['structure'] = df_merged['structure'].str.replace('missing', 'loop')
        df_merged['catalysis'] = df_merged['catalysis'].str.replace('None catalytic', 'Non catalytic')
        df_merged.to_csv(done_path_features_only,index=False)
        return


    complete_features = ['2_clique', '3_clique', '4_clique', '5_clique', '6_clique', '7_clique', 'average_neighbor_degree',
                    'betweenness_centrality', 'clustering_coefficient', 'degree_centrality', 'eigenvector_centrality',
                    'graphlet1',
                    'graphlet10', 'graphlet11', 'graphlet12', 'graphlet13', 'graphlet14', 'graphlet15', 'graphlet16',
                    'graphlet17',
                    'graphlet18', 'graphlet19', 'graphlet2', 'graphlet20', 'graphlet21', 'graphlet22', 'graphlet23',
                    'graphlet24',
                    'graphlet25', 'graphlet26', 'graphlet27', 'graphlet28', 'graphlet29', 'graphlet3', 'graphlet30',
                    'graphlet31',
                    'graphlet32', 'graphlet33', 'graphlet34', 'graphlet35', 'graphlet36', 'graphlet37', 'graphlet38',
                    'graphlet39',
                    'graphlet4', 'graphlet40', 'graphlet41', 'graphlet42', 'graphlet43', 'graphlet44', 'graphlet45',
                    'graphlet46',
                    'graphlet47', 'graphlet48', 'graphlet49', 'graphlet5', 'graphlet50', 'graphlet51', 'graphlet52',
                    'graphlet53',
                    'graphlet54', 'graphlet55', 'graphlet56', 'graphlet57', 'graphlet58', 'graphlet59', 'graphlet6',
                    'graphlet60',
                    'graphlet61', 'graphlet62', 'graphlet63', 'graphlet64', 'graphlet65', 'graphlet66', 'graphlet67',
                    'graphlet68',
                    'graphlet69', 'graphlet7', 'graphlet70', 'graphlet71', 'graphlet72', 'graphlet73', 'graphlet8',
                    'graphlet9',
                    'median_rsa_neigh', 'median_wcn_ca_neigh', 'median_wcn_sc_neigh', 'node_degree', 'rsa',
                    'total_A_neigh',
                    'total_Aliphatic_neigh', 'total_Aromatic_neigh', 'total_C_neigh', 'total_Charged_neigh',
                    'total_D_neigh',
                    'total_Diverse_neigh', 'total_E_neigh', 'total_F_neigh', 'total_G_neigh', 'total_H_neigh',
                    'total_Hydrophobic_neigh', 'total_I_neigh', 'total_K_neigh', 'total_L_neigh', 'total_M_neigh',
                    'total_N_neigh',
                    'total_P_neigh', 'total_Polar_neigh', 'total_Q_neigh', 'total_R_neigh', 'total_S_neigh',
                    'total_T_neigh',
                    'total_Tiny_neigh', 'total_V_neigh', 'total_W_neigh', 'total_Y_neigh', 'total_catalytic_neigh',
                    'total_contact_neigh', 'total_disordered_neigh', 'total_glycosylated_neigh',
                    'total_interface_neigh',
                    'total_nan_neigh', 'total_site_neigh', 'total_ss_B_neigh', 'total_ss_E_neigh', 'total_ss_G_neigh',
                    'total_ss_H_neigh', 'total_ss_I_neigh', 'total_ss_P_neigh', 'total_ss_S_neigh', 'total_ss_T_neigh',
                    'total_ss_nan_neigh', 'wcn_ca', 'wcn_sc', 'glycosylation', 'protein_interaction', 'disorder',
                    'binding',
                    'catalysis', 'pdb_aa', 'aa_group_5', 'aa_group_HP', 'structure']
    # if 'pdb_aa' not in df_merged.columns.tolist():
    # logging.debug(df_merged['pdb_aa_x'])
    # logging.debug(df_merged['pdb_aa_y'])
    try:
        df_merged['pdb_aa'] = df_merged['pdb_aa_x'].fillna(df_merged['pdb_aa_y'])
    except:
        logging.debug(df_merged['pdb_aa'])
        logging.debug(df_merged['pdb_aa'].dropna())


    for f in complete_features:
        if f not in df_merged.columns.tolist():
            logging.debug(f'missing feature\n\n {f}')
            df_merged[f] = np.nan

    not_features = ['pdb_position', 'chain_x', 'chain', 'chain_y', 'pdb_aa_x', 'pdb_aa_y'] + [f for f in
                                                                                              df_merged.columns.tolist()
                                                                                              if f.startswith(
            'neighbor_pos_')]
    # categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa_x', 'aa_group_5', 'aa_group_HP','structure']
    categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                            'aa_group_5', 'aa_group_HP', 'structure']
    # if len(chain_list) <= 1 or orphan_prediction:
    #     categorical_features.remove('protein_interaction')
    all_features = [f for f in df_merged.columns.tolist() if f not in not_features]
    numeric_features = [f for f in all_features if f not in categorical_features]
    if not orphan_prediction:
        numeric_features.append('mean_neigh_r4s_score')

    logging.debug(f'numeric\n\n{numeric_features}')
    logging.debug(f'cat\n\n{categorical_features}')
    all_features = [*numeric_features, *categorical_features]

    logging.debug(f'evorator: {254}')
    names_d = {'normalized_score': 'consurf_score',
               'total_ss_nan_neigh': 'total_ss_loop_neigh',
               'total_nan_neigh': 'total_non_interfacing_neigh'}
    names_d_inv = {v: k for k, v in names_d.items()}
    if not orphan_prediction:
        logging.debug(f'evorator: {256}')

        consurf_df = pd.read_csv(consurf_output,sep='\t',skiprows=15,skipfooter=4,header=None)
        # consurf_df = consurf_df[[0,1,2,3,5]]
        logging.debug(f'{consurf_df}')
        logging.debug(f'{consurf_df[[0, 1, 2, 3,4, 5, 7,8,9,10]]}')

        # consurf_df.columns = ['position', 'seq', 'pdb_seq','normalized_score','color']
        #consurf_df = consurf_df[[0, 1, 2, 3,4,5,6,7,8,9,10,11,12,13]]
        if not consurfdb_query:
            consurf_df = consurf_df[[0, 1, 2, 3, 5, 6]]
            # consurf_df = consurf_df[[0, 1, 2, 3, 5, 7]]
        else:
            consurf_df = consurf_df[[0, 1, 2, 3, 5, 6]]
            # consurf_df = consurf_df[[0, 1, 2, 3, 5, 7]]
        consurf_df.columns = ['position', 'seq', 'pdb_seq', 'normalized_score', 'color', 'ci']
        logging.debug(f'{consurf_df.values.tolist()}')

        # ALA117: A

        consurf_df = consurf_df[consurf_df['pdb_seq']!="         -"]
        consurf_df['chain'] =  consurf_df['pdb_seq'].str.split(":").str[-1]
        consurf_df['pdb_position_consurf'] = consurf_df['pdb_seq'].str.split(":").str[0].str.extract('(\d+)')
        consurf_df['pdb_position'] = consurf_df['chain'].astype(str) +  consurf_df['pdb_position_consurf'].astype(str) + "_" + job_title
        if prediction_task == 'PERfIFR':
            logging.debug(f'{consurf_df.columns}')
            logging.debug(f'{consurf_df["ci"]}')
            logging.debug(f'{consurf_df["ci"].str.split(",")}')
            consurf_df['ci_lower'] = consurf_df['ci'].str.split(",").str[0].astype(float)
            consurf_df['ci_upper'] = consurf_df['ci'].str.split(",").str[1].astype(float)
            ci_lower_d = dict(zip(consurf_df['pdb_position'], consurf_df['ci_lower']))
            ci_upper_d = dict(zip(consurf_df['pdb_position'], consurf_df['ci_upper']))

        df_merged = df_merged.merge(consurf_df[['pdb_position','normalized_score','color']],on='pdb_position',how='left')
        neigh_df = df_merged[['pdb_position',*[f for f in df_merged.columns.tolist() if 'neighbor_pos_' in f],'normalized_score']]
        consurf_score_d = dict(zip(neigh_df['pdb_position'],neigh_df['normalized_score']))
        neigh_d = dict(zip(neigh_df['pdb_position'],neigh_df[[f for f in df_merged.columns.tolist() if 'neighbor_pos_' in f]].values.tolist()))
        mean_neigh_r4s_score_d = {}
        for k in neigh_d:
            tmp=[]
            for n in neigh_d[k]:
                tmp.append(consurf_score_d.get(n,np.nan))
            mean_neigh_r4s_score_d[k] = np.nanmean(tmp)
        df_merged['mean_neigh_r4s_score'] = df_merged['pdb_position'].map(mean_neigh_r4s_score_d)

        features4consurf =all_features
        # features4consurf =[f for f in features4consurf if f not in [*[str(i)+"_clique" for i in range(100)], *["graphlet"+str(i) for i in range(1,100)],*[f for f in features4consurf if (df_merged[f] == 0).all()]]]

        logging.debug(f'final features: {features4consurf}')
        logging.debug(f'cat features: {categorical_features}')
        logging.debug(f'numeric_features: {numeric_features}')
        for c in complete_features:
            if c not in df_merged.columns.tolist():
                df_merged[c] = np.nan

        df_train = df_merged[~df_merged['color'].astype(str).str.endswith('*')].dropna(subset=['normalized_score'],axis=0)
        df_test = df_merged[df_merged['color'].astype(str).str.endswith('*')]



        for cat in categorical_features:
            df_train[cat].fillna('missing', inplace=True)
            df_test[cat].fillna('missing', inplace=True)

        X_test  = df_test[features4consurf]
        X_train  = df_train[features4consurf]

        y_test  = df_test['normalized_score']
        y_train  = df_train['normalized_score']

        if X_test.shape[0]==1:
            X_test = X_test.values.reshape(1, -1)
            y_test = y_test.values.reshape(1, -1)

        logging.debug(f'train: {X_train.shape}')
        logging.debug(f'train: {y_train.shape}')
        logging.debug(f'train: {y_train.dropna().shape}')
        logging.debug(f'test: {X_test.shape}')
        logging.debug(f'test: {y_test.shape}')
        regressor = joblib.load(trained_regressor)




        logging.debug(f'test line {320}')

#        regressor.fit(X_train,y_train)

        logging.debug(f'test line {359}')

        y_cv_pred = regressor.predict(X_train)
        r2_test = r2_score(y_train,y_cv_pred)
        logging.debug(f'test R2 cross-validation: {r2_test}')

        if X_test.shape[0] != 0:
            pred = regressor.predict(pd.DataFrame(X_test,columns=features4consurf))
        else:
            pred = []

        final_df_train = pd.concat([df_train[['pdb_position', *features4consurf, 'normalized_score']].reset_index(),pd.Series(y_cv_pred,name='cv_predicted_score')],axis=1)
        final_df_test = pd.concat([df_test[['pdb_position', *features4consurf]].reset_index(), pd.Series(pred, name='predicted_score_4_missing_in_consurf')],axis=1)
        logging.debug(f'final train: {final_df_train.shape}')
        logging.debug(f'final test: {final_df_test.shape}')
        final_df = pd.concat([final_df_train,final_df_test],axis=0)
        logging.debug(f'final df: {final_df.tail(10)}')
        logging.debug(f'train:{df_train.shape}')
        logging.debug(f'test:{df_test.shape}')

        neigh_df.to_csv(done_path_neigh,index=False)
        features4consurf_2_write = [names_d.get(f, f) for f in features4consurf]
        if prediction_task == 'PERfGR':
            final_df['predicted_score'] = final_df['cv_predicted_score'].fillna(final_df['predicted_score_4_missing_in_consurf'])
            final_df = final_df.rename(columns=names_d)
            df_2_write=final_df[['pdb_position', *features4consurf_2_write, 'consurf_score',
                      'predicted_score_4_missing_in_consurf']].drop(columns=['total_non_interfacing_neigh'])
            df_2_write.to_csv(done_path, index=False)
            final_df = final_df.rename(columns=names_d_inv)
        elif prediction_task == 'PERfIFR':
            final_df['predicted_score'] = final_df['cv_predicted_score']
            final_df['ci_lower'] = final_df['pdb_position'].map(ci_lower_d)
            final_df['ci_upper'] = final_df['pdb_position'].map(ci_upper_d)
            final_df['upper_e'] = final_df['ci_upper'].astype(float) - final_df[
                'normalized_score'].astype(float)
            final_df['lower_e'] = final_df['normalized_score'].astype(float) - final_df[
                'ci_lower'].astype(float)
            final_df = final_df.rename(columns=names_d)
            final_df['glycosylation'] = final_df['glycosylation'].str.replace('None glycosylated', 'Non glycosylated')
#            final_df['protein_interaction'] = final_df['protein_interaction'].str.replace('missing', np.nan)
            final_df['protein_interaction'] = final_df['protein_interaction'].replace('missing', np.nan)
            final_df['structure'] = final_df['structure'].str.replace('missing', 'loop')
            final_df['catalysis'] = final_df['catalysis'].str.replace('None glycosylated', 'Non glycosylated')

            final_df = final_df.rename(columns=names_d_inv)

        # final_df=final_df.rename(columns=names_d)
        # final_df[['pdb_position', *features4consurf_2_write, 'consurf_score',
        #       'predicted_score_4_missing_in_consurf', 'predicted_score']].to_csv(done_path, index=False)
        # final_df=final_df.rename(columns=names_d_inv)
    else:
        try:
            df_merged['chain_x'] = df_merged['chain_x'].fillna(df_merged['chain_y'])
            df_merged['pdb_aa_x'] = df_merged['pdb_aa_x'].fillna(df_merged['pdb_aa_y'])
            df_merged = df_merged.rename(columns={'chain_x':'chain', 'pdb_aa_x':'pdb_aa'})
        except:
            pass
        df_merged = df_merged.loc[:, ~df_merged.columns.duplicated()]
        categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                                'aa_group_5', 'aa_group_HP', 'structure']





        for cat in categorical_features:
            logging.debug(f'{cat}')
            df_merged[cat].fillna('missing', inplace=True)

        logging.debug(f'{551}')
        X = df_merged[[*complete_features]]
        X.to_csv(os.path.join(results_dir, 'test_del.csv'))
        logging.debug(f'{X.shape}')
        logging.debug(f'all\n{X.columns[X.isna().any()].tolist()}')
        logging.debug(f'ABC\n{df_merged[[*categorical_features]].columns[df_merged[[*categorical_features]].isna().any()].tolist()}')
        logging.debug(f'#\n{df_merged[[*numeric_features]].columns[df_merged[[*numeric_features]].isna().any()].tolist()}')

        # Load ExtraTreesRegressor (like RandomForest)
        regressor = joblib.load(trained_regressor)

        # Predict
        y_hat = regressor.predict(X)
        df_merged['predicted_score'] = y_hat
        # Write feature set & predictions to file & plot feature importance
        df_merged = df_merged.rename(columns=names_d)
        complete_features_2_write = [names_d.get(f,f) for f in complete_features]
        df_merged['glycosylation'] = df_merged['glycosylation'].str.replace('None glycosylated', 'Non glycosylated')
        df_merged['protein_interaction'] = df_merged['protein_interaction'].replace('missing', np.nan)
        df_merged['structure'] = df_merged['structure'].str.replace('missing', 'loop')
        df_merged['catalysis'] = df_merged['catalysis'].str.replace('None catalytic', 'Non catalytic')

        df_merged = df_merged.rename(columns=names_d_inv)
        final_df = df_merged
    final_df['ResInd'] = final_df['pdb_position'].str.split("_").str[0]
    final_df['ResInd'] = final_df['ResInd'].str.replace('\D+', '')
    final_df['ResInd'] = final_df['ResInd'].astype(int)

    min_bin_neg = final_df['predicted_score'].min()
    max_bin_neg = final_df[final_df['predicted_score'] < 0]['predicted_score'].max()
    diff = max_bin_neg - min_bin_neg
    interval = diff / 4.5
    bins = []
    bins.append(min_bin_neg)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)
    bins.append(bins[-1] + interval)

    logging.debug(f'scores:{bins[0]}')
    logging.debug(f'scores:{bins[1]}')
    logging.debug(f'scores:{interval}')
    logging.debug(f'scores:{bins}')
    data = final_df['predicted_score']
    final_df['Score'] = np.digitize(data, bins[::-1])
    final_df.loc[final_df['Score'] == 0, 'Score'] = 1


    pdbFileOriginal = [f for f in os.listdir(results_dir) if f.endswith(".ent") or f.endswith(".pdb")]
    pdbFileOriginal = [f for f in pdbFileOriginal if "_" not in f][0]

    if not orphan_prediction:
        if prediction_task=='PERfGR':
            final_df['diff'] = final_df['normalized_score'].fillna(final_df['predicted_score'])
            if final_df['diff'].min() < final_df['normalized_score'].min():
                min_bin_neg = final_df['diff'].min()
            else:
                min_bin_neg = final_df['normalized_score'].min()
            if final_df[final_df['diff'] < 0]['diff'].max()> final_df[final_df['normalized_score'] < 0]['normalized_score'].max():
                max_bin_neg = final_df[final_df['diff'] < 0]['diff'].max()
            else:
                max_bin_neg = final_df[final_df['normalized_score'] < 0]['normalized_score'].max()
        elif prediction_task=='PERfIFR':
            final_df['diff'] = final_df['normalized_score'] - final_df['predicted_score']
            # final_df['diff'] = final_df['normalized_score'] # test consistency with consurf
        elif prediction_task=='PERfIFR':
            min_bin_neg = final_df['diff'].min()
            max_bin_neg = final_df[final_df['diff'] < 0]['diff'].max()

        diff = max_bin_neg - min_bin_neg
        interval = diff / 4.5
        bins = []
        bins.append(min_bin_neg)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        bins.append(bins[-1] + interval)
        data = final_df['diff']
        final_df['dScore'] = np.digitize(data, bins[::-1])
        logging.debug(f'{final_df}')
        final_df.loc[final_df['dScore'] == 0, 'dScore'] = 1
        final_df['dScore'] = final_df['dScore'].fillna(0)
        if prediction_task == 'PERfIFR':
            final_df.loc[final_df['normalized_score'].isnull()==True, 'dScore'] = 0

        logging.debug(f'{final_df}')
        logging.debug(f'scores:{final_df[["ResInd", "Score", "predicted_score", "normalized_score"]]}')
        final_df[["ResInd", "dScore"]].sort_values(by="ResInd").to_csv(
            os.path.join(results_dir, 'evorator_diff.scores'),
            index=False)
        final_df[["pdb_position", 'pdb_aa', 'dScore', 'ResInd']].sort_values(by="ResInd")[
            ["pdb_position", 'pdb_aa', 'dScore']].to_csv(os.path.join(results_dir, 'evorator.scores.for.2d'),
                                                        index=False)

        job_info = {"chainId": pdb_chain, "scoresFile": "evorator_diff.scores", "jobTitle": job_title,
                    "pdbFile": pdbFileOriginal}

    else:

        logging.debug(f'scores:{final_df[["ResInd","predicted_score"]]}')
        final_df[["ResInd", "Score"]].sort_values(by="ResInd").to_csv(os.path.join(results_dir, 'evorator.scores'),
                                                                      index=False)
        final_df[["pdb_position",'pdb_aa', 'Score','ResInd']].sort_values(by="ResInd")[["pdb_position",'pdb_aa', 'Score']].to_csv(os.path.join(results_dir, 'evorator.scores.for.2d'),
                                                                      index=False)

        final_df["Score"]=  final_df["predicted_score"].round(2)
        final_df[["ResInd", "Score"]].sort_values(by="ResInd").to_csv(
            os.path.join(results_dir, 'evorator_raw.scores'),
            index=False)

        job_info = {"chainId": pdb_chain, "scoresFile": "evorator.scores", "jobTitle": job_title,
                    "pdbFile": pdbFileOriginal}

        # job_info = {"chainId": pdb_chain, "scoresFile": "evorator_raw.scores", "jobTitle": job_title,
        #             "pdbFile": pdbFileOriginal}


    # final_df[
    #     ['pdb_position', *features, *[f for f in missing_features if f not in features], 'normalized_score',
    #      'predicted_score_train', 'predicted_score_test','predicted_score','Score']].to_csv(done_path, index=False)


    # job_info = {"chainId": pdb_chain, "scoresFile": "evorator.scores", "jobTitle": job_title, "pdbFile": pdbFileOriginal}
    with open(os.path.join(results_dir,'job_info.json'), 'w') as f:
        json.dump(job_info, f)

    # logging.debug(f'error shape={final_df[["lower_e", "upper_e"]].T}')
    # if os.path.exists(done_path) or os.path.exists(done_figure_path):
    final_df.to_csv(done_path)
    if os.path.exists(done_path):
        logging.debug(f'768')

        if pdb_name=='':
            pdb_name = os.path.split(pdb_input)[-1].split(".")[0]
        COORD_FILE_FOR_DRAWING_NETWORK = os.path.join(results_dir,pdb_name.upper() + "_" + pdb_chain + "_xyz.txt")
        EDGELIST_FILE_FOR_DRAWING_NETWORK = os.path.join(results_dir,pdb_name.upper() + "_" + pdb_chain + "_edgelist_draw.txt")
        le = LabelEncoder()
        content_df = pd.read_csv(edgelist_file,header=None,sep='\t')
        n_lines= len(content_df)
        content = content_df.iloc[:, 0].tolist() + content_df.iloc[:, 1].tolist()
        print()

        unique_content = np.unique(content)
        unique_content_ints = [int(re.search(r"\d+",i).group()) for i in unique_content]
        sorter = np.argsort(unique_content_ints)
        unique_content_sorted = unique_content[sorter]
        print(unique_content_sorted)
        print(unique_content.shape)
        le.fit(unique_content_sorted)
        # content_coded = le.transform(content)

        content_coded_1 = le.transform(content_df.iloc[:, 0][sorter]).astype(str)
        content_coded_2 = le.transform(content_df.iloc[:, 1][sorter]).astype(str)

        print(content_coded_1)
        print(content_coded_2)

        with open(COORD_FILE_FOR_DRAWING_NETWORK,'w') as f:

            Calpha_coordinates = backbone_coordinates[:, 2, :]
            for x,y,z in zip(residue_pdb_index,Calpha_coordinates, sequence_from_pdb):
                f.write(''.join(x[1:])+'\t'+str(y[0])+'\t'+str(y[1])+'\t'+str(y[2])+'\t'+z+'\n')


        cmd = f'/groups/pupko/natannag/conda/envs/NatanEnv/bin/python {os.path.join(scripts_dir,"draw_3d_network.py")} {COORD_FILE_FOR_DRAWING_NETWORK} {EDGELIST_FILE_FOR_DRAWING_NETWORK} {os.path.join(results_dir, "evorator.scores.for.2d")} {results_dir}'
#        cmd = f'module load python/python-anaconda3.7-itaym; python {os.path.join(scripts_dir,"draw_network.py")} {results_dir}/Compute_bulk_input/xyz/{pdb_name.upper()}_{pdb_chain}_xyz.txt {results_dir}/Compute_bulk_input/edgelist/{pdb_name.upper()}_{pdb_chain}_edgelist.txt {os.path.join(results_dir, "evorator.scores.for.2d")} {results_dir}'
#        cmd = f'module load python/python-anaconda3.7-itaym; python {os.path.join(scripts_dir,"draw_3d_network.py")} {results_dir}/Compute_bulk_input/xyz/{pdb_name.upper()}_{pdb_chain}_xyz.txt {results_dir}/Compute_bulk_input/edgelist/{pdb_name.upper()}_{pdb_chain}_edgelist.txt {os.path.join(results_dir, "evorator.scores.for.2d")} {results_dir}'
        logging.debug(f'creating network image: {cmd}')
        subprocess.check_output(cmd, shell=True)
        logging.debug(f'773')
        # try:
        if prediction_task=='PERfGR':
            finalize_html(html_path, error_path, run_number, job_title,results_dir,consurf_output,task='PERfGR')
        elif prediction_task=='PERfO':


            # pssm prediction #
            # if prediction_task != 'PERfGR' and prediction_task != 'PERfIFR':
            matrix = matlist.blosum62
            aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                  'Y']
            blosum_d = {}
            for AA1 in aa:
                tmp_vec = []
                for AA2 in aa:
                    tmp_vec.append(matrix.get((AA1, AA2), matrix.get((AA2, AA1))))
                blosum_d[AA1] = tmp_vec
            blosum_cols = ['blosum' + str(i) for i in range(len(blosum_d['A']))]
            # logging.debug(f"{final_df['pdb_aa'].isnull().all()}")
            # logging.debug(f"{final_df['pdb_aa'].tolist()}")
            if final_df['pdb_aa'].unique().tolist() == ['missing']:
                for c in blosum_cols:
                    final_df[c] = np.nan
                final_df['aa_freq_ExPASy'] = np.nan
            else:
                final_df[blosum_cols] = final_df['pdb_aa'].map(blosum_d).apply(pd.Series)

                aa_freqs = {'A': 8.25 / 100, 'Q': 3.93 / 100, 'L': 9.65 / 100, 'S': 6.63 / 100, 'R': 5.53 / 100,
                            'E': 6.72 / 100, 'K': 5.80 / 100, 'T': 5.35 / 100, 'N': 4.05 / 100,
                            'G': 7.08 / 100, 'M': 2.41 / 100, 'W': 1.09 / 100, 'D': 5.46 / 100, 'H': 2.27 / 100,
                            'F': 3.86 / 100, 'Y': 2.92 / 100, 'C': 1.38 / 100, 'I': 5.91 / 100, 'P': 4.73 / 100,
                            'V': 6.86 / 100}
                final_df['aa_freq_ExPASy'] = final_df['pdb_aa'].map(aa_freqs)
            numeric_features.append('aa_freq_ExPASy')
            for c in blosum_cols:
                numeric_features.append(c)
            all_features_pssm = [*numeric_features, *categorical_features]
            logging.debug(all_features_pssm)
            # all_features_pssm = ['2_clique', '3_clique', '4_clique', '5_clique', '6_clique', '7_clique', 'average_neighbor_degree', 'betweenness_centrality', 'clustering_coefficient', 'degree_centrality', 'eigenvector_centrality', 'graphlet1', 'graphlet10', 'graphlet11', 'graphlet12', 'graphlet13', 'graphlet14', 'graphlet15', 'graphlet16', 'graphlet17', 'graphlet18', 'graphlet19', 'graphlet2', 'graphlet20', 'graphlet21', 'graphlet22', 'graphlet23', 'graphlet24', 'graphlet25', 'graphlet26', 'graphlet27', 'graphlet28', 'graphlet29', 'graphlet3', 'graphlet30', 'graphlet31', 'graphlet32', 'graphlet33', 'graphlet34', 'graphlet35', 'graphlet36', 'graphlet37', 'graphlet38', 'graphlet39', 'graphlet4', 'graphlet40', 'graphlet41', 'graphlet42', 'graphlet43', 'graphlet44', 'graphlet45', 'graphlet46', 'graphlet47', 'graphlet48', 'graphlet49', 'graphlet5', 'graphlet50', 'graphlet51', 'graphlet52', 'graphlet53', 'graphlet54', 'graphlet55', 'graphlet56', 'graphlet57', 'graphlet58', 'graphlet59', 'graphlet6', 'graphlet60', 'graphlet61', 'graphlet62', 'graphlet63', 'graphlet64', 'graphlet65', 'graphlet66', 'graphlet67', 'graphlet68', 'graphlet69', 'graphlet7', 'graphlet70', 'graphlet71', 'graphlet72', 'graphlet73', 'graphlet8', 'graphlet9', 'median_rsa_neigh', 'median_wcn_ca_neigh', 'median_wcn_sc_neigh', 'node_degree', 'rsa', 'total_A_neigh', 'total_Aliphatic_neigh', 'total_Aromatic_neigh', 'total_C_neigh', 'total_Charged_neigh', 'total_D_neigh', 'total_Diverse_neigh', 'total_E_neigh', 'total_F_neigh', 'total_G_neigh', 'total_H_neigh', 'total_Hydrophobic_neigh', 'total_I_neigh', 'total_K_neigh', 'total_L_neigh', 'total_M_neigh', 'total_N_neigh', 'total_P_neigh', 'total_Polar_neigh', 'total_Q_neigh', 'total_R_neigh', 'total_S_neigh', 'total_T_neigh', 'total_Tiny_neigh', 'total_V_neigh', 'total_W_neigh', 'total_Y_neigh', 'total_catalytic_neigh', 'total_contact_neigh', 'total_disordered_neigh', 'total_glycosylated_neigh', 'total_interface_neigh', 'total_site_neigh', 'total_ss_B_neigh', 'total_ss_E_neigh', 'total_ss_G_neigh', 'total_ss_H_neigh', 'total_ss_I_neigh', 'total_ss_P_neigh', 'total_ss_S_neigh', 'total_ss_T_neigh', 'total_ss_nan_neigh', 'wcn_ca', 'wcn_sc', 'predicted_score', 'blosum0', 'blosum1', 'blosum2', 'blosum3', 'blosum4', 'blosum5', 'blosum6', 'blosum7', 'blosum8', 'blosum9', 'blosum10', 'blosum11', 'blosum12', 'blosum13', 'blosum14', 'blosum15', 'blosum16', 'blosum17', 'blosum18', 'blosum19', 'aa_freq_ExPASy','glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa', 'aa_group_5', 'aa_group_HP', 'structure']
            print(all_features_pssm)
            # X_test_pssm = final_df[all_features_pssm]
            final_df.to_csv(done_path)
            # cmd_predict_pssm = 'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv; python /groups/pupko/natannag/consurf_n2v/huang/preprocessor_pssm_backup_290722.py /groups/pupko/natannag/consurf_n2v/huang/predict_huang/huang_with_evorator_features_and_predictions.csv /groups/pupko/natannag/consurf_n2v/huang/predict_huang/huang_aa_variance_pssm.csv /groups/pupko/natannag/consurf_n2v/huang/trained_pssm_model 280722 /bioseq/data/results/evorator/1658918346/2LZM_features_and_predictions.csv /bioseq/data/results/evorator/1658918346     test_del'
            cmd_predict_pssm = f"module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv;!@# python {CONSTS.PSSM_PRED_EXE} {CONSTS.HUANG_DATA_PATH} {CONSTS.HUANG_PSSM_PATH} {CONSTS.OUT_PATH_PSSM_AUX} {CONSTS.PSSM_NAME_AUX} {done_path} {results_dir} {job_title}\t{job_title}_pssm_evorator"
            with open(f'{results_dir+"/"+"predict_pssm_cla_test.cmds"}','w') as f:
                f.write(cmd_predict_pssm)
            # subprocess.check_output(f'echo {cmd_predict_pssm} > {results_dir+"/"+"predict_pssm_cla_test.cmds"}',shell=True)
            qsub_predict_pssm = f'ssh bioseq@power9login /bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {results_dir+"/"+"predict_pssm_cla_test.cmds"} {results_dir} -q pupkoweb --verbose > {results_dir}/qsub.log'
            # final_df= preprocessor_pssm.predict_new_case(final_df,CONSTS.HUANG_DATA_PATH,CONSTS.HUANG_PSSM_PATH,
            #                                              CONSTS.OUT_PATH_PSSM_AUX,str(CONSTS.PSSM_NAME_AUX))
            logging.debug(f'pred pssm = {qsub_predict_pssm}')
            subprocess.check_output(qsub_predict_pssm,shell=True)
            # final_df = pd.read_csv(done_path)
            i=0
            while not os.path.exists(results_dir+'/'+job_title+"_done_pssm_pred.flag"):
                if i==0:
                    logging.debug(f'waiting for pssm prediction to finish')
                time.sleep(1)
                i+=1

                # read file
                # logging.debug(f'X_test_pssm shape = {X_test_pssm.shape}\n')
                # logging.debug(f'loading trained pssm model\n')
                # TRAINED_MODEL_PSSM = load_model("/groups/pupko/natannag/consurf_n2v/huang/trained_pssm_model/final_output_ANN_n_hidden_147_trained_model_del.pkl")
                # TRAINED_MODEL_PSSM = load_model(CONSTS.TRAINED_PSSM_ANN)
                # logging.debug(f'trained pssm model loaded\n')
                # logging.debug(f'loading fitted pssm preprocessor\n')
                # PREPROCESSOR_PSSM = pickle.load(open("/groups/pupko/natannag/consurf_n2v/huang/trained_pssm_model/preprocessor_pssm_270722.pkl",'rb'))
                # PREPROCESSOR_PSSM = pickle.load(open(CONSTS.FITTED_PSSM_PREPROCESSOR,'rb'))
                # logging.debug(f'pssm preprocessor loaded\n')
                # logging.debug(f'preprocessing test data\n')
                # X_test_pssm_processed = PREPROCESSOR_PSSM.transform(X_test_pssm)
                # logging.debug(X_test_pssm_processed.shape)
                # logging.debug(f'test data processed\n')
                # logging.debug(f'predicting PSSM\n')
                # f'new features={preprocessor.named_steps["preprocessor"].transformers_[1][1].named_steps["onehot"].get_feature_names(categorical_features)}')

                # PREDICTED_PSSM = TRAINED_MODEL_PSSM.predict(X_test_pssm_processed)
                # for i,aa1 in enumerate(pssm_pred_cols):
                #     final_df[aa1]  = PREDICTED_PSSM[:,i]
                # final_df = pd.read_csv(done_path)

                if os.path.isfile(results_dir + '/' + job_title + "_done_pssm_pred.flag"):

                    pssm_pred_cols = [AA1 + '_pred' for AA1 in aa]
                    final_df = pd.read_csv(done_path)

                    df_2_write = final_df[['pdb_position', *all_features_pssm, 'predicted_score',*pssm_pred_cols]]
                    df_2_write.to_csv(done_path, index=False)
                    finalize_html(html_path, error_path, run_number, job_title,results_dir,consurf_output)
        elif prediction_task == 'PERfIFR':
            # logging.debug(f'{final_df.columns.tolist()}')
            final_df = final_df.rename(columns={'total_ss_nan_neigh':'total_ss_loop_neigh','total_nan_neigh':'total_non_interfacing_neigh','normalized_score':'consurf_score'})
            df_2_write = final_df[
                ['pdb_position', *features4consurf_2_write, 'consurf_score', 'lower_e', 'upper_e',
                 'predicted_score_4_missing_in_consurf', 'predicted_score']].drop(
                columns=['total_non_interfacing_neigh'])
            df_2_write.to_csv(done_path, index=False)

            plt.errorbar(y=final_df["predicted_score"],
                         x=final_df["consurf_score"],
                         xerr=final_df[['lower_e', 'upper_e']].T.values.tolist(),
                         marker='.', ls='', lw=0.5)
            plt.plot([final_df["consurf_score"].min() - 0.5, final_df["consurf_score"].max() + 0.5],
                     [final_df["consurf_score"].min() - 0.5, final_df["consurf_score"].max() + 0.5],
                     color='black', linestyle='--')
            plt.xlabel('ConSurf rate', size=12, fontweight="bold")
            plt.ylabel('EvoRator rate', size=12, fontweight="bold")
            plt.savefig(done_figure_path)
            finalize_html(html_path, error_path, run_number, job_title, results_dir, consurf_output, r2_test,
                          task='PERfIFR', done_figure_path=done_figure_path)
        else:
            logging.debug(f'779')
            finalize_html(html_path, error_path, run_number, job_title,results_dir,consurf_output)

# except Exception as e:
    if html_path:
        edit_failure_html(CONSTS,
                          "",
                          html_path, run_number)
        add_closing_html_tags(html_path, CONSTS, run_number)
        return

        # except Exception as e:
        #     logging.debug(f'SUCCEEDED = False')
        #     logging.debug(e)
        #     logging.debug(f"""results_dir: {results_dir}\n
        #                 pdb_file:{pdb_input}\npdb_chain:{pdb_chain}""")
        #     if html_path:
        #         # error_msg = e.args[-1]
        #         error_msg = e.args
        #         if os.path.exists(error_path):
        #             with open(error_path) as f:
        #                 error_msg = f.read()
        #         edit_failure_html(CONSTS, error_msg, html_path, run_number)
        #         add_closing_html_tags(html_path, CONSTS, run_number)




if __name__ == '__main__':
    import sys

    import os
    # sys.path.append('/bioseq/evorator/auxiliaries')
    import EvoRator.evorator_CONSTANTS as CONSTS  # from /effectidor/auxiliaries
    from time import sleep

    if os.path.exists('/bioseq/evorator'):  # remote run
        sys.path.append('/bioseq/evorator/auxiliaries/')
        sys.path.append('/bioseq/bioSequence_scripts_and_constants/')
        sys.path.append('/groups/pupko/natannag/consurf_n2v/huang')
        sys.path.append('/groups/pupko/natannag/ScanNet_dev/')
        sys.path.append('/groups/pupko/natannag/EvolutionPrediction/')
        #
    # C:\Users\natan\Documents\EvolutionPrediction\extract_pairwise_distances.py
    laptop = False if os.path.exists("/groups/pupko/natannag") else True
    if laptop:
        path2github = "../"
        path2scannet = "C:/Users/natan/Documents/ScanNet_dev/"
        path2evolutionprediction = "C:/Users/natan/Documents/EvolutionPrediction/"
    else:
        path2github = "/groups/pupko/natannag"
        path2scannet = '/groups/pupko/natannag/ScanNet_dev/'
        path2evolutionprediction = '/groups/pupko/natannag/EvolutionPrediction/'
    sys.path.append(path2github)
    sys.path.append(path2scannet)
    sys.path.append(path2evolutionprediction)
    # import evorator_CONSTANTS as CONSTS  # from /bioseq/natan_conservation_webserver/auxiliaries/
    # from GENERAL_CONSTANTS import PDB_DIVIDED  # from /bioseq/bioSequence_scripts_and_constants/
    import extract_pairwise_distances
    from preprocessing import PDBio, PDB_processing

    import obtain_pdb
    import json
    import time
    import pandas as pd
    import numpy as np
    import subprocess
    from Bio.PDB import PDBParser, PDBIO

    import extract_features_4_a_single_query
    import logging
    import matplotlib.pyplot as plt
    from sklearn.metrics import r2_score
    import scipy
    import matplotlib.pyplot as plt
    import json
    from auxiliaries import fail, update_html, append_to_html  # from /effectidor/auxiliaries
    import urllib
    from bs4 import BeautifulSoup
    if not laptop:
        from sklearn.externals import joblib #
    else:
        import joblib
    from Bio.SubsMat import MatrixInfo as matlist
    from tensorflow.python.keras.models import Model, load_model
    import sklearn
    from sklearn.preprocessing import LabelEncoder

    import shutil
    import re
    main()
