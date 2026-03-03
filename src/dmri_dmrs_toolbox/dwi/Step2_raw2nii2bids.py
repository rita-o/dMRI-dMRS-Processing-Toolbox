"""
Script to convert raw data to nifti format and organize it in BIDS format
It uses the Dicomifier module, installed in the Dicomifier python environment

Last changed Jan 2025
@author: Rita O
"""

import os
import re
import json
import sys
import pandas as pd

from dmri_dmrs_toolbox.misc.custom_functions import (
    create_directory, raw_to_nifti, copy_file, extract_methods, plot_bvals
)
from dmri_dmrs_toolbox.misc.bids_structure import create_bids_structure


def Step2_raw2nii2bids(cfg):
    
    data_path       = cfg['data_path']     
    scan_list       = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in cfg['subj_list']:
    
        print('Converting to nifti ' + subj + '...')
    
        # Extract data for this subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        
        # Generate paths and convert data 
        raw_path        = os.path.join(data_path, 'raw_data', list(subj_data['raw_data_folder'].unique())[0]) 
        nifti_path      = os.path.join(data_path, 'nifti_data', 'unsorted', subj)
        create_directory(nifti_path)
        raw_to_nifti(raw_path, nifti_path, cfg)
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['sessNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='nifti_data', workingdir='sorted')
            
            # Index of scans for this session 
            study_indx  = subj_data.index[subj_data['sessNo'] == sess].tolist()
            
            ###### SCAN-WISE OPERATIONS ######
            for scn_ctr in study_indx:
                
                # scan folder number
                scan_no = subj_data['scanNo'][scn_ctr]
                
                if subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'PGSE':
        
                    method_path   = os.path.join(raw_path,str(scan_no), 'method')
                    # Get paths and directories
                    if subj_data['phaseDir'][scn_ctr] == 'fwd':
                        nii_path    = os.path.join(nifti_path,str(scan_no) + '_1_' + subj_data['acqSeq'][scn_ctr])
                    elif subj_data['phaseDir'][scn_ctr] == 'rev':
                        with open(os.path.join(raw_path,str(scan_no), 'acqp'), 'r') as f:
                            for line in f:
                                if '##$ACQ_scan_name=' in line: 
                                    match=re.search(r'\((.*?)\)',next(f))
                                    ref_name=match.group(1)[1:] 
                    
                        if 'V3.5' in subj_data['PV'][scn_ctr]:
                            nii_path    = os.path.join(nifti_path,str(scan_no) + '_1_' + 'ADJ_REVPE_E' + ref_name)
                        elif  'V1.1' in subj_data['PV'][scn_ctr]:
                            nii_path    = os.path.join(nifti_path,str(scan_no) + '_1_' + subj_data['acqSeq'][scn_ctr])
                    bids_strc.set_param(datatype='dwi',description='Delta_'+str(int(subj_data['diffTime'][scn_ctr]))+'_'+subj_data['phaseDir'][scn_ctr])

                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('dwi.nii.gz')])
                    extract_methods(method_path, bids_strc, 'PGSE')
                    plot_bvals(bids_strc)
        
                elif subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'STE':
        
                    method_path   = os.path.join(raw_path,str(scan_no), 'method')
                    # Get paths and directories
                    if subj_data['phaseDir'][scn_ctr] == 'fwd':
                        nii_path    = os.path.join(nifti_path,str(scan_no) + '_1_' + subj_data['acqSeq'][scn_ctr])
                    elif subj_data['phaseDir'][scn_ctr] == 'rev':
                        with open(os.path.join(raw_path,str(scan_no), 'acqp'), 'r') as f:
                            for line in f:
                                if '##$ACQ_scan_name=' in line: 
                                    match=re.search(r'\((.*?)\)',next(f))
                                    ref_name=match.group(1)[1:] 
                        nii_path    = os.path.join(nifti_path,str(scan_no) + '_1_' + 'ADJ_REVPE_E' + ref_name)
                    
                    bids_strc.set_param(datatype='dwi_STE',description='STE_'+ subj_data['phaseDir'][scn_ctr])
        
                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('dwi.nii.gz')])
                    extract_methods(method_path, bids_strc, 'STE', cfg)
                    plot_bvals(bids_strc)
    
                elif subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'T2W':
        
                    # Get paths and directories
                    nii_path    = os.path.join(nifti_path,str(scan_no) + '_1_' + subj_data['acqSeq'][scn_ctr])
                    bids_strc.set_param(datatype='anat',description=None)

                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('T2w.nii.gz')])
        
                elif subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'B0': # probably not working
        
                    # Get paths and directories
                    nii_path    = os.path.join(nifti_path,str(scan_no)  + '_1_' + subj_data['acqSeq'][scn_ctr])
                    bids_strc.set_param(datatype='B0',description=None)

                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('B0.nii.gz')])
        
            #shutil.rmtree(os.path.join(data_path, 'nifti_data', 'unsorted'))


