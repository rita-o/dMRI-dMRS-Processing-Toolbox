#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare groups microstrcture
@author: localadmin
"""

import os
import sys
import matplotlib.pyplot as plt
from dmri_dmrs_toolbox.dmrs.dmrsmodel import DMRSModel

subj_list = [f'sub-{i:02d}' for i in [3,4,5,6,7,8,10,11,12,14,15]]    # list of subjects to analyse [8,10,11,12,14,15]]#

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = "/media/localadmin/DATA/data/CTD/"          # path to where the data from the cohort is
cfg['prep_foldername']      = 'preprocessed'    # name of the preprocessed folder (keep 'preprocessed' as default)
cfg['analysis_foldername']  = 'analysis'        # name of the analysis folder (keep 'analysis' as default)
cfg['common_folder']        = "/home/localadmin/Software/dMRI_dMRS_toolbox/common/"  # path to the common folder with files needed throught the pipeline
cfg['scan_list_name']       = 'ScanList_CTD.xlsx'   # name of the excel file containing the metadata of the cohort
cfg['model_list']           =  ['cylinder']
cfg['metabolites']          = ['NAA+NAAG','Glu','Ins','GPC+PCho','Cr+PCr','Tau','Gln']              # metabolites for analysis


from src.dmri_dmrs_toolbox.misc.bids_structure import *
from src.dmri_dmrs_toolbox.misc.custom_functions import *

from scipy.optimize import curve_fit
import pandas as pd
import glob
import copy
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu
from pathlib import Path

output_folder = Path(cfg['data_path'])/'results'
output_folder.mkdir(parents=True, exist_ok=True)

scan_list   = pd.read_excel(os.path.join(cfg['data_path'], cfg['scan_list_name']))

def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx - 1)*np.var(x, ddof=1) + (ny - 1)*np.var(y, ddof=1)) / dof)



######## MODEL-WISE OPERATIONS ########
dmrsmodel = DMRSModel()
for model in cfg['model_list']:
    print(f'Working with {model}...')
    params = dmrsmodel.model_specs[model].param_names
    nb_params = len(params)

    df_all_data = pd.DataFrame()

    for subj in cfg['subj_list']:
        print('Working on subject ' + subj + '...')
        # Extract data for subject
        subj_data = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        data_path = cfg['data_path']
        group = subj_data['group'].iloc[0]
        if '?' in group:
            group = 'no group yet'

        # List of acquisition sessions
        sess_list = [x for x in list(subj_data['sessNo'].unique()) if not math.isnan(x)]  # clean NaNs

        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            print('Working on session ' + str(sess) + '...')

            bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=data_path,
                                                       folderlevel='derivatives', workingdir=cfg['analysis_foldername'],
                                                       description=model)

            for metab in cfg['metabolites']:
                fit_params_path = Path(bids_strc_analysis.get_path())/"csvs"/f'fit_parameters_{metab}_{model}.csv'
                if fit_params_path.is_file():
                    fit_parameters = pd.read_csv(fit_params_path)
                else:
                    fit_parameters = pd.DataFrame(
                        [np.nan]*nb_params,
                        columns = params
                    )
                    group = "no data"
                fit_parameters['Metabolite'] = metab
                fit_parameters["Subject"] = subj
                fit_parameters["Group"] = group
                fit_parameters["Session"] = int(sess)

                df_all_data = pd.concat([df_all_data, fit_parameters], ignore_index=True)

    ######## FOR ALL THE DATA OF THIS MODEL ########
    plot_rows = []

    for _, row in df_all_data.iloc[::2].iterrows(): # only take mean values for now
        metab = row['Metabolite']
        group = row['Group']
        session = row['Session']
        subject = row['Subject']

        for param in params:
            values = np.ravel(row[param])
            for val in values:
                plot_rows.append({
                    "Parameter": param,
                    "Metabolite": metab,
                    "Group": group,
                    "Session": session,
                    "Subject": subject,
                    "Model": model,
                    "Value": float(val)
                })
    df_long = pd.DataFrame(plot_rows)

    ## Create labels
    group_counts = df_long[['Group', 'Subject']].drop_duplicates().groupby('Group').count()['Subject']

    desired_order = ['WT', 'KI', 'HTZ', 'no group yet', 'no data']
    HUE_ORDER = [g for g in desired_order if g in group_counts.index]

    df_long['Group'] = (df_long['Group']
                        .astype(str)
                        .str.strip())
    df_long['Group'] = pd.Categorical(df_long['Group'], categories=desired_order, ordered=True)

    ## Plot
    print('Plotting')
    params = np.array(params)
    plotting_params = params[np.where(params!="amp")]

    fig, axes = plt.subplots(len(plotting_params), 1, figsize=(15, 15))
    fig.subplots_adjust(wspace=0.05, hspace=0.11, top=0.95, bottom=0.1, left=0.05, right=0.95)

    lims = dmrsmodel.model_specs[model].bounds
    for i, (ax, param) in enumerate(zip(axes, params)):
        lim

        sns.boxplot(
            data=df_long[df_long['Parameter'] == param],
            x='Metabolite',
            y='Value',
            hue='Group',
            hue_order=HUE_ORDER,
            ax=ax,
            #split=False,
            #inner='quart',
            #cut=0
        )

        metabs = df_long['Metabolite'].unique()
        for metab in metabs:
            data_roi = df_long[df_long['Metabolite'] == metab]
            group_WT = data_roi[data_roi['Group'] == 'WT']['Value']
            group_KI = data_roi[data_roi['Group'] == 'KI']['Value']

            # Choose test (use Mann-Whitney U if non-normal, or t-test if assumptions met)
            # if not group_WT.empty and not group_KI.empty:
            #     result = mannwhitneyu(group_WT, group_KI, alternative='two-sided')  # non-parametric

            #     pval = result.pvalue
            #     effect_size = cohen_d(group_WT.values, group_KI.values)
            #     if pval < 0.05 and abs(effect_size) >= 0.3:
            #         # Get median y for positioning *
            #         y_max = max(data_roi['Value'].max(), lim[1]) * 1.02
            #         xpos = list(rois).index(roi)
            #         ax.text(x=xpos, y=lim[1]*0.9, s='*', ha='center', va='bottom', fontsize=16, color='red')

        if i == len(axes) - 1:
            ax.set_xlabel("Metabolite", fontsize=11)
            ax.tick_params(axis='x', rotation=30)
            handles, lgn = ax.get_legend_handles_labels()
            new_lgn = []
            for l in lgn:
                count = group_counts[l]
                new_lgn.append(f"{l} (n={count})")
            ax.legend(handles, new_lgn, loc='upper right')
        else:
            ax.set_xlabel("")
            ax.set_xticklabels([])
            ax.get_legend().remove()

        ax.set_ylim(lim[0], lim[1] * 1.1)
        ax.set_ylabel(dmrsmodel.model_specs[model].param_names_latex, fontsize=11)
        ax.tick_params(axis='y', labelsize=11)
        ax.tick_params(axis='x', labelsize=11)

    fig.suptitle(f"{dmrsmodel.model_specs[model].label}", fontsize=16)
    plt.savefig(os.path.join(output_folder, f'{model}_group_comparison.png'))
    plt.show()



