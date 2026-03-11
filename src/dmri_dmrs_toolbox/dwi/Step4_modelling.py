"""
Script to analyse GM microstructure by fitting models such as Nexi, Sandi, ...
It uses the SwissKnife module installed in the SwissKnife python environment.

Last changed Jan 2025
@author: Rita O
"""
import os
import pandas as pd
import numpy as np
import subprocess
import shutil
import glob
import fnmatch
from pathlib import Path
from dmri_dmrs_toolbox.misc.bids_structure import create_bids_structure
from dmri_dmrs_toolbox.misc.custom_functions import (
    create_directory,
    remove_folder,
    copy_files_BIDS,
    get_file_in_folder,
    read_numeric_txt,
    modify_units_bvals,
    multiply_by_mask,
    plot_summary_params_model,
    register_outputfits_to_anat,
    calculate_pwd_avg,
    compute_micro_FA,
    get_param_names_model,
    copy_file
)
from importlib.resources import files
from dmri_dmrs_toolbox.misc.atlas_functions import prepare_atlas_labels, create_ROI_mask

def Step4_modelling(cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    cfg['model_list'] = cfg['model_list_GM'] + cfg['model_list_WM']
    # Define path to docker
    docker_path = '/data' 
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in cfg['subj_list']:
        
        print('Modelling ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        subj_data      = subj_data[subj_data['analyse'] == 'y']

        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['sessNo'].unique()) :
          
            print('Working on session ' + str(sess) + '...')
            
            ########################## DELTA-WISE OPERATIONS ##########################      
            filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['sessNo'] == sess) & (subj_data['noBval'] > 1) & (subj_data['acqType'] == 'PGSE')]
            Delta_list = filtered_data['diffTime'].unique().astype(int).tolist()
            
            for Delta in Delta_list:

                # Define bids structure for the processed data
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                #bids_strc_prep.set_param(description='Delta_' + str(Delta) + '_fwd') # repricated: done on each Delta processed separately
                bids_strc_prep.set_param(description='allDelta-allb/Delta_' + str(Delta)) # new: done on each Delta processed together ('allDelta')

                ######## Run DTI and DKI ########  
                 
                # Define BIDS structure for the analysis data
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=f'DTI_DKI_Delta_{Delta}', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
              
                # Make output folder 
                output_path = bids_strc_analysis.get_path()
                
                # Decide which model to run
                run_DTIDKI = 1
                bvals = np.unique(read_numeric_txt(get_file_in_folder(bids_strc_prep, '*bvalsNom.txt'))[0])
                has_low  = np.any(bvals <= 1000)
                has_high = np.any(bvals > 1000)
                
                do_model = []
                if has_high:  # at least one shell above 1000
                    do_model.append('-DKI')
                if has_low:  # all bvals <= 1000
                    do_model.append('-DTI')
                if has_high==0 and has_low==0:  # no valid shells found
                    run_DTIDKI = 0
                    print('No valid shells found for neither DTI nor DKI!')
                do_model = ' '.join(do_model)

                # Just run model if it doesn't exist on the folder yet
                if not os.path.exists(output_path) or cfg['redo_modelling']:
                                        
                    if cfg['use_server_mount']==1:
                        data_path_temp  = os.path.join(os.path.expanduser('~'), 'temp', 'current_sub')
                        input_path      = os.path.join(data_path_temp, 'inputs')
                        mount_path      = data_path_temp
                        remove_folder(mount_path)
                        needs_copy_back = True
                    else:
                        input_path      = os.path.join(output_path, 'inputs')
                        mount_path      = data_path
                        needs_copy_back = False
                    
                    create_directory(input_path)
                    
                    # --- Copy inputs
                    dwi  = copy_files_BIDS(bids_strc_prep, input_path, 'dwi_dn_gc_ec.mif')
                    mask = copy_files_BIDS(bids_strc_prep, input_path, 'mask.nii.gz')
                    
                    # --- Convert to docker paths
                    dwi        = dwi.replace(mount_path, docker_path)
                    mask       = mask.replace(mount_path, docker_path)
                    out_folder = docker_path if cfg['use_server_mount'] else output_path.replace(data_path, docker_path)
                    
                    # --- Run model
                    if run_DTIDKI == 1:
                    
                        extra = '-maxb 7' if cfg['is_alive'] == 'ex_vivo' else ''
                    
                        call = [
                            f'docker run -v {mount_path}:{docker_path} '
                            f'nyudiffusionmri/designer2:v2.0.13 '
                            f'tmi {do_model} {extra} '
                            f'{dwi} {out_folder} -fit_constraints 0,1,0'
                        ]
                    
                        print(' '.join(call))
                        os.system(' '.join(call))
                    
                    # --- Restore local paths
                    bids_strc_analysis.set_param(root=data_path)
                    bids_strc_prep.set_param(root=data_path)
                    output_path = bids_strc_analysis.get_path()
                    
                    # --- Copy back only if needed
                    if needs_copy_back:
                        shutil.copytree(data_path_temp, output_path, dirs_exist_ok=True) 

                if run_DTIDKI==1:
                    # Mask output with brain mask for better visualization
                    for filename in os.listdir(output_path):
                        if filename.endswith(".nii"):
                            multiply_by_mask(os.path.join(output_path, filename), # filename input
                                             os.path.join(output_path,'Output_masked'), # output folder
                                                     os.path.join(output_path,'inputs', 'mask.nii.gz'),cfg) # mask
                            
                    # Plot summary in dwi space
                    bids_strc_prep.set_param(description='allDelta-allb') # new: done on each Delta processed together ('allDelta')
                    plot_summary_params_model(os.path.join(output_path,'Output_masked'), 'DTI_DKI', cfg,bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))
                    
                    # Register to anat space
                    bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                              folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                    register_outputfits_to_anat(os.path.join(output_path,'Output_masked'),
                                                os.path.join(output_path,'Output_in_anat'),
                                            'DTI_DKI',cfg, bids_strc_anat, bids_strc_prep)
                
                    # Plot summary plot in anat space
                    bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']  +'-To-'+cfg['anat_format'], root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                    bids_strc_reg.set_param(base_name='')
                    if cfg['subject_type']=='organoid' and os.path.exists(bids_strc_reg.get_path('mask_organoids.nii.gz')):
                             atlas=bids_strc_reg.get_path(f"atlas_in_{cfg['anat_format']}.nii.gz")
                             atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])
                             mask = create_ROI_mask(atlas, atlas_labels, [], 'organoids', cfg['tpm_thr'], bids_strc_reg)
                             plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), 'DTI_DKI', cfg, bids_strc_anat.get_path(f"{cfg['anat_format']}_bc_brain.nii.gz"), bids_strc_reg.get_path('mask_organoids.nii.gz'))
    
                    else:
                        plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), 'DTI_DKI', cfg, bids_strc_anat.get_path(f"{cfg['anat_format']}_bc_brain.nii.gz"))

                
                ######## Compute PWD for LTE data ######## 
                
                # Create BIDS structures for LTE
                bids_LTE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f'allDelta-allb/Delta_{Delta}')
                bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f'pwd_avg_Delta_{Delta}')
              
                # Create pwd average of LTE, just if it doesn't exist yet
                if not os.path.exists(bids_LTE.get_path()) or cfg['redo_modelling']:
                    create_directory(bids_LTE.get_path())
                    calculate_pwd_avg(get_file_in_folder(bids_LTE_temp,'*dwi_dn_gc_ec.nii.gz'),
                                      get_file_in_folder(bids_LTE_temp,'*bvalsNom.txt'),
                                      get_file_in_folder(bids_LTE_temp,'*bvalsEff.txt'),
                                      bids_LTE.get_path(),
                                      np.nan)
             

            ######## Compute PWD for STE data (if exists) ######## 
            
            # Create BIDS structures for STE
            bids_STE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg')
            
            # Create pwd average of STE, just if it doesn't exist yet
            if os.path.exists(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz')):                
                 
                  if not os.path.exists(bids_STE.get_path()) or cfg['redo_modelling']:
                      create_directory(bids_STE.get_path())
                      calculate_pwd_avg(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz'),
                                        bids_STE_temp.get_path('bvalsNom.txt'),
                                        bids_STE_temp.get_path('bvalsEff.txt'),
                                        bids_STE.get_path(),
                                        np.nan)
                      
            ######## Compute PWD for STE data registered to LTE (if exists) ######## 
             
            # Create BIDS structures for STE
            bids_STE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                           folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                           folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg_in_LTE')
            bids_STE_reg      = create_bids_structure(subj=subj, sess=sess, datatype='registration', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='STE-To-LTE_'+ "allDelta-allb")
            bids_STE_reg.set_param(base_name='')
            
            # Create pwd average of STE, just if it doesn't exist yet
            if get_file_in_folder(bids_STE_reg,'*dn_gc_topup.nii.gz'):                
                  
                   if not os.path.exists(bids_STE.get_path()) or cfg['redo_modelling']:
                       create_directory(bids_STE.get_path())
                       calculate_pwd_avg(get_file_in_folder(bids_STE_reg,'*dn_gc_topup.nii.gz'),
                                         bids_STE_temp.get_path('bvalsNom.txt'),
                                         bids_STE_temp.get_path('bvalsEff.txt'),
                                         bids_STE.get_path(),
                                         np.nan)
                          
            
            ######## Compute MicroFA if data exists ########  
            if get_file_in_folder(bids_STE_reg,'*dn_gc_topup.nii.gz'):   
               
                # 1. Computing microFA in matlab 
                      # Define BIDS structure
                # bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                #               folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA_matlab')
                # output_path = bids_STE.get_path()
                # bids_STE_reg      = create_bids_structure(subj=subj, sess=sess, datatype='registration', root=cfg['data_path'] , 
                #               folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='STE-To-LTE_'+ "allDelta-allb")
                # bids_STE_reg.set_param(base_name='')
                # bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                #               folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
                # bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                #              folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f"allDelta-allb/Delta_{cfg['LTEDelta_for_microFA']}")
                # header        =  get_file_in_folder(bids_LTE,'*mask.nii.gz')
                
                     # Just run model if it doesn't exist on the folder yet
                # if not os.path.exists(output_path) or cfg['redo_modelling']:
                #     mdm_matlab(bids_LTE, bids_STE, bids_STE_reg, header, output_path, cfg['code_path2'], cfg['toolboxes'],low_b=True)

                # 2. Computing microFA in python
                    # Define BIDS structure
                bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f"pwd_avg_Delta_{cfg['LTEDelta_for_microFA']}")
                bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                               folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg_in_LTE')
                bids_LTE_aux  = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f"allDelta-allb/Delta_{cfg['LTEDelta_for_microFA']}")
                mask          =  get_file_in_folder(bids_LTE_aux,'*mask.nii.gz')
                bids_STE_out      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                              folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA')
                output_path = bids_STE_out.get_path()
                
                    # Just run model if it doesn't exist on the folder yet
                if not os.path.exists(output_path) or cfg['redo_modelling']:
                    compute_micro_FA(bids_LTE, bids_STE, mask, output_path)
                
                  
            ########################## MODEL-WISE OPERATIONS ##########################       
            for model in cfg['model_list']:
                
                print('Working with ' + model + '...')

                if 'Nexi' in model or 'Smex' in model or 'Sandix' in model:  
                    data_used = 'allDelta-allb'
                elif 'Sandi' in model: # lowest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['sessNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmin')()
                    #data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'   # depricated
                    data_used = 'allDelta-allb/Delta_'+str(int(filtered_data['diffTime'][ind_folder]))   # new
                elif model=='SMI' or model=='SMI_wSTE': # largest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['sessNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmax')()
                    #data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'   # depricated
                    data_used = 'allDelta-allb/Delta_'+str(int(filtered_data['diffTime'][ind_folder]))   # new

                print(f'Using data: {data_used}')
                # Define bids structure 
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=model)
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=data_used, root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                
                # Define paths   
                output_path = bids_strc_analysis.get_path()
                input_path = os.path.join(output_path,'inputs')

                # Just run model if it doesn't exist on the folder yet
                if not os.path.exists(output_path) or cfg['redo_modelling']:
                   
                    # Make output folder 
                    create_directory(input_path)

                    # Copy necessary files for analysis
                    dwi         = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.nii.gz')
                    big_delta   = copy_files_BIDS(bids_strc_prep,input_path,'DiffTime.txt')
                    small_delta = copy_files_BIDS(bids_strc_prep,input_path,'DiffDuration.txt')
                    bvals       = copy_files_BIDS(bids_strc_prep,input_path,'bvalsNom.txt')
                    sigma       = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_sigma.nii.gz')
                    mask        = copy_files_BIDS(bids_strc_prep,input_path,'mask_dil.nii.gz')
        
                    # Get diffusion duration (assumes the same value for all acquisitions)
                    #small_delta = np.loadtxt(small_delta)[0]
             
                    # Modify units of bvals for NEXI          
                    new_bvals = bvals.replace('.txt','_units.txt')
                    modify_units_bvals(bvals, new_bvals )
            
                    # Copy necessary files for analysis 
                    if 'Nexi'in model or 'Sandi' in model or 'Smex' in model or 'Sandix' in model:  
                        bids_strc_lowb = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description="allDelta-lowb", root=data_path, 
                                                    folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                        sigma     = copy_files_BIDS(bids_strc_lowb,input_path,'dwi_dn_sigma.nii.gz')
                    elif model=='SMI':
                        input_file =  copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                        others = '' 
                    elif model=='SMI_wSTE':
                        bids_strc_STE = create_bids_structure(subj=subj, sess=sess, datatype="dwi_STE", root=data_path, 
                                                    folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
                        STE        = copy_files_BIDS(bids_strc_STE,input_path,'dwi_dn_gc_topup.mif').replace(data_path,docker_path)
                        LTE        = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                        input_file = LTE + ',' + STE
                        others     = '-echo_time 51,51 -bshape 1,0 -compartments EAS,IAS -debug'
               
                    # Run SwissKnife models
                    if 'Nexi'in model or model=='Sandi' or 'Smex' in model or 'Sandix' in model:  
                        
                        # Define arguments 
                        args = [model, 
                                output_path, 
                                dwi,  
                                new_bvals, 
                                big_delta,  
                                small_delta, 
                                sigma,
                                mask,
                                cfg['is_alive'],
                                '--debug']
                                                  
                        # Run script
                        if "wmicroFA" in model:
                            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA')
                            uFA = os.path.join(bids_STE.get_path(),'Uaniso.nii')
                            args.insert(-1, uFA)  
                            script_path = files("dmri_dmrs_toolbox.misc").joinpath("auxiliar_modelling.py")
                            command = [
                                cfg["conda_exe"], "run", "-n", "SwissKnife_exp",
                                "python", str(script_path)
                            ] + args
                            subprocess.run(command, check=True)
                            
                        elif "mrsinformed" in model:
                            bids_mrs      = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=f'dmrs-to-allDelta-allb', root=cfg['data_path'], 
                                           folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                            args[7] = get_file_in_folder(bids_mrs,'*voxel_mrs.nii.gz')

                            # mrs_radius_s = os.path.join(bids_mrs.get_path(),'mrs_radius_s.nii')
                            args.insert(-1, '10.5')  # sub-01
                            args.insert(-1, '10')  # sub-01
                            script_path = files("dmri_dmrs_toolbox.misc").joinpath("auxiliar_modelling.py")
                            command = [
                                cfg["conda_exe"], "run", "-n", "SwissKnife_exp",
                                "python", str(script_path)
                            ] + args
                            subprocess.run(command, check=True)
                        else:
                            script_path = files("dmri_dmrs_toolbox.misc").joinpath("auxiliar_modelling.py")
                            command = [
                                cfg["conda_exe"], "run", "-n", "SwissKnife",
                                "python", str(script_path)
                            ] + args
                            subprocess.run(command, check=True)
            
                        
                    # Run Designer models
                    elif model=='SMI' or model=='SMI_wSTE':  
                        # Define arguments 
                        args = [model, 
                                output_path.replace(data_path,docker_path), 
                                input_file,  
                                mask.replace(data_path,docker_path),  
                                sigma.replace(data_path,docker_path), 
                                data_path,
                                others]
                    
                        # Run script
                        script_path = files("dmri_dmrs_toolbox.misc").joinpath("auxiliar_modelling.py")
                        command = [
                            cfg["conda_exe"], "run", "-n", "pipeline",
                            "python", str(script_path)
                        ] + args
                        subprocess.run(command, check=True)
                        
                    # Run matlab models
                    elif model=='Sandi_MP':  
   
                        # Put in the format Sandi toolbox wants
                        bvecs       = copy_files_BIDS(bids_strc_prep,input_path,'bvecsRotated.txt')
                        for file in [sigma, dwi, bvals, bvecs, mask]:
                            src = Path(file)
                            sub = [p for p in src.parts if p.startswith("sub-")][0]
                            ses = [p for p in src.parts if p.startswith("ses-")][0]
                            acq = "acq-01"
                            run = "run-01"
                            if file==sigma:
                                ending='noisemap.nii.gz'
                            elif file==dwi:
                                ending='dwi.nii.gz'
                            elif file==mask:
                                ending='mask.nii.gz'
                            elif file==bvals:
                                ending='dwi.bval'
                            elif file==bvecs:
                                ending='dwi.bvec'
                            new_name = f"{sub}_{ses}_{acq}_{run}_desc-preproc_{ending}"
                            
                            dst_dir = src.parent.parent / "derivatives" / "preprocessed" / sub / ses
                            create_directory(dst_dir)
                          
                            dst = dst_dir / new_name
                            shutil.copy2(src, dst)   # copy, do not move

                        big_delta_val = float(np.loadtxt(big_delta)[0])
                        small_delta_val = float(np.loadtxt(small_delta)[0])
                        
                        # # rename the folder inputs → preprocessed
                        # inputs_dir = src.parent
                        # preproc_dir = inputs_dir.parent / "derivatives"
                        # inputs_dir.rename(preproc_dir)
                        
                        # Matlab command
                        sandi_folder = src.parent.parent
                        matlab_cmd = (
                            "try, "
                            f"addpath(genpath('{os.path.join(cfg['toolboxes'], 'SANDI')}')); "
                            f"SANDI_batch_analysis('{sandi_folder}', {big_delta_val}, {small_delta_val}, []); "
                            "catch, exit(1), end, exit(0);"
                        )
                        cmd = [
                             "matlab", "-nodisplay", "-nosplash", "-nodesktop",
                             "-r", matlab_cmd
                        ]
                        # Run matlab command
                        subprocess.run(cmd)
                        
                        # Put the files in the normal format
                        src_dir = sandi_folder / "derivatives" / 'SANDI_analysis' / sub / ses / 'SANDI_Output'  
                        dst_dir = sandi_folder 
                        
                        pattern_map = {
                            "SANDI-fit_Din.nii.gz": "sandi_di.nii.gz",
                            "SANDI-fit_De.nii.gz": "sandi_de.nii.gz",
                            "SANDI-fit_fneurite.nii.gz": "sandi_fneurite.nii.gz",
                            "SANDI-fit_fsoma.nii.gz": "sandi_fsoma.nii.gz",
                            "SANDI-fit_Rsoma.nii.gz": "sandi_rs.nii.gz",
                        }

                        for pattern, new_name in pattern_map.items():
                            matches = glob.glob(str(src_dir / pattern))
                        
                            if not matches:
                                print(f"No file found for pattern: {pattern}")
                                continue
                        
                            src = Path(matches[0])
                        
                            dst = dst_dir / new_name
                        
                            shutil.move(src, dst)
                        
                            print(f"Moved: {src} -> {dst}")

                # Mask output for better visualization
                patterns, lims, maximums = get_param_names_model(model,cfg['is_alive'])
                for filename in os.listdir(output_path):
                    if any(fnmatch.fnmatch(filename, pattern) for pattern in patterns):
                        multiply_by_mask(os.path.join(output_path, filename), # filename input
                                         os.path.join(output_path,'Output_masked'), # output folder
                                         os.path.join(output_path,'inputs', 'mask_dil.nii.gz'),cfg) # mask
                # Plot summary plot in dwi space
                bids_strc_prep.set_param(description='allDelta-allb') # new: done on each Delta processed together ('allDelta')
                plot_summary_params_model(os.path.join(output_path,'Output_masked'), model, cfg, bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))
                
                # Register to anat space
                bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['prep_foldername'])                
                register_outputfits_to_anat(os.path.join(output_path,'Output_masked'),
                                            os.path.join(output_path,'Output_in_anat'),
                                            model,cfg, bids_strc_anat, bids_strc_prep)
                
                # Plot summary plot in anat space
                if cfg['subject_type']=='organoid':
                         bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']  +'-To-'+cfg['anat_format'], root=data_path, 
                                                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                         bids_strc_reg.set_param(base_name='')
                         atlas=bids_strc_reg.get_path(f"atlas_in_{cfg['anat_format']}.nii.gz")
                         atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])
                         mask = create_ROI_mask(atlas, atlas_labels, [], 'organoids', cfg['tpm_thr'], bids_strc_reg)
                         plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), model, cfg, bids_strc_anat.get_path(f"{cfg['anat_format']}_bc_brain.nii.gz"), bids_strc_reg.get_path('mask_organoids.nii.gz'))

                else:
                    plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), model, cfg, bids_strc_anat.get_path(f"{cfg['anat_format']}_bc_brain.nii.gz"))

           
            
            


