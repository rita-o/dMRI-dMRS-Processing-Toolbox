# dMRI-dMRS Processing Toolbox

This toolbox is intended for researchers working with advanced diffusion MRI (dMRI), in particular *multi-shell* and *multi–diffusion-time* acquisitions; and diffusion MRS (dMRS) data with multi–diffusion-time and multi–b-value acquired with the *SPECIAL* sequence, on Bruker or Siemens scanners.

It provides **preprocessing and analysis pipelines** for both modalities:

  🔵 **Diffusion MRI (dMRI)**: detailed documentation for the dMRI processing pipeline is available [here](README_dMRI.md).  
  🟡 **Diffusion MRS (dMRS)**: detailed documentation for the dMRS processing pipeline is available [here](README_dMRS.md).  

> 🛠️ **This toolbox is a continuous work in progress.**  
> Please pull the latest changes frequently.  
> If you encounter any issues or have questions, **let us know so we can improve it.**  
> **Contact**: ana.veiga-de-oliveira@chuv.ch

 <br> 

## 🚀 Quick Start  

1. Clone the repository and install prerequisites as described below.
2. Prepare the cohort Excel file (see `common/example_study.xlsx`). (⚠️ Note the change of columns names in the cohort file in Dec/2025).
3. Put raw data from the scanner under:
     <pre>
   folder_study_name (name of your project's folder)
       └── raw_data  
         └── raw_data_folder (name of the folder created in the MRI scanner)
      </pre>
4. (Optional!)  
   For ROI-based analysis of dMRI data, prepare an atlas and/or tissue probability map (TPM) and save it in `common/atlas` (see notes [here](README_dMRI.md)).
5. Choose the appropriate processing script:  

   🔵 For **dMRI processing** (inside folder `processing_dwi`):
      - **Main_rat.py** — For rodent dMRI data  
      - **Main_human.py** — For human dMRI data (not extensively tested)  
      - **Main_organoid.py** — For organoid dMRI data   
   
   🟡 For **dMRS processing** (inside folder `processing_dmrs`):
      - **Main.py** — For rodent dMRI data  
8. Edit the configuration structure `cfg` at the start of each script to customize the parameters for the analysis.
9. Run individual steps (e.g., `StepX`) (recommended) or the full pipeline.
   
 <br>
 
 
## 📌 Notes

### 1. Dataset Metadata

Each dataset must include an Excel file with **cohort metadata** (e.g., subject ID, group, scan date).  
Some columns will be automatically filled during Step1 of the pipeline.  
An example file is provided in the `common` folder.The following columns must be pre-filled manually before running the script:
  
> - **study_name**: Assigned study name (e.g., `sub-01`, `sub-02`, …)  
> - **raw_data_folder**: Name of the raw data folder saved in the MRI system  
>   &nbsp;  • Must match the folder containing this subject’s raw data  
>   &nbsp;  • Should match the name used in the methods file
> - **other_name_tag**: (optional) Other name tag, like the rat number or batch...
> - **sex**: (optional) Sex of the subject (`F` or `M`) 
> - **group**: (optional) Group tag (e.g.,`1` or `2`, or `WT`, `disease`) 
> - **scanNo**: Folder number of raw imaging data (integer)  
> - **acqType**: Acquisition type (`T2W`, `PGSE`, `STE`, `SPECIAL`); pay attention to capital letters. 
> - **sessNo**: Session number (usually `1`, unless it’s a rescan)  
> - **Reorient**:  Data collected on a Bruker scanner is typically in the orientation:  `x: L→R`, `y: P→A`, `z: I→S`. To match standard atlas orientations, it is recommended to reoriented dMRI data to: `x: L→R`, `y: S→I`, `z: A→P` (This corresponds to axis flipping as: `x −z y`)  
>   &nbsp;      This standard orientation allows easier integration with online atlases and tools. *(Required for dMRI data; not required for dMRS.)*
> - **VoxMidHem**: voxel of the mid coronal plane in dwi space to then define left and right hemispheres. If you don't know or don't care set it to zero and ignore the results of the dMRI metrics plots left and right. *(Required for dMRI data; not required for dMRS.)*
> - **anat_thr**: intensity threshold to be used as initial guess for the creation of an anatomical brain mask for dMRI processing. *(Required for dMRI data; not required for dMRS.)*
> - **Notes**: Notes regarding that specidic subejct/acquisition.
> - **analyse**: 'y' (yes) or 'n' (no) if that row of data is to be analyzed or not (for example if there are repeated/bad scans put that column to 'y' only on the one you want to keep).
> - **TM**: For dMRS analysis only! Mixing time in ms.
> - **dMRS_acq_type**: For dMRS analysis only! Options are: 'water' or 'metab', if the data was acquired for water or metabolites, respectively. 
> - **coil_type**: For dMRS analysis only! Options are: 'rat' for animals scanned with rat cryo prob; or 'mouse' for animals scanned with mouse cryo probe

### 2. `common/` Folder

The `common` folder contains **shared resources** required across multiple processing pipelines:

- Pre-configured **Anaconda environments** (you can dowload them and add to you Anaconda)
- Example **ScanList_example.xlsx**
- dMRI: **STE sequence b-values and b-vecs**, which cannot be retrieved from the methods file  
  (plus placeholder/fake b-vectors needed for parts of the analysis)
- dMRI: **Atlas folders** for ROI-based analysis of dMRI data ([see read me](README_dMRI.md)). Atlases used: Rodents: *WHS_v4*; Humans: *DKT*, *Juelich*. ⚠️ Atlas folders are too large for GitHub. Please **contact us** if you need access to what we usually use.
- dMRI: **.cnf topup files** with parameter details for topup analysis of dMRI data.
- dMRS: **mrs_basis_sets**: basis datasets for dMRS quantification

### 3. Common Python Scripts

Several Python modules support multiple components of the pipeline:

- **`auxiliar_modelling.py`**  
  Auxiliary functions for modeling routines (e.g., fitting options such as `nexi`).

- **`bids_structure.py`**  
  Utilities for organizing data in **BIDS** format.

- **`custom_functions.py`**  
  General-purpose helper functions used throughout the pipeline.

- **`atlas_functions.py`**  
  Atlas-specific utilities: preparing and resampling atlases; parsing atlas label files; mapping ROI names to atlas region IDs. When adding a **new atlas**, this module must be updated accordingly.

 <br>
 
## SOFTWARE REQUIREMENTS

This package is written in Python and uses Conda to manage multiple environments.

The pipeline relies on:

1. Command-line tools available in the system PATH, so they can be called from any directory without their fullpath (e.g fsleyes). 
2. Other tools installed in dedicated Conda environments.
3. Matlab runtime if dealing with dMRS data.

### 1. Command-line tools
The following tools must be installed and accessible from the command line:

- [**RATS_MM**](https://iibi.uiowa.edu/rats-rodent-brain-mri) Add to your system's `PATH` after installation. Need only for brain extraction of rodent data if this option is chosen (available options: RATS, UNET). Need for dMRI data analysis.

- [**ANTs (Advanced Normalization Tools)**](https://github.com/ANTsX/ANTs) Add to your system's `PATH` after installation. Needed throughout preprocessing. Need for dMRI data analysis.

- [**FSL (FMRIB Software Library)**](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/) Add to your system's `PATH` after installation. Needed throughout preprocessing. Need for dMRI data analysis.

- [**MRtrix3**](https://www.mrtrix.org/) Needed throughout preprocessing. Need for dMRI data analysis.

- [**DESIGNER**](https://nyu-diffusionmri.github.io/DESIGNER-v2/) Needed for denoising and DTI/DKI fitting. Need for dMRI data analysis.
  
- [**LCModel**](https://s-provencher.com/lcmodel.shtml) Needed for dMRS metabolite quantification. It does **not** need to be added to the system `PATH`, as the scripts directly call the executable from its installation directory.

 <br> 

### 2. Tools installed in Conda environments
Each of the following tools should be installed in its own Conda environment:

- [**Dicomifier**](https://github.com/lamyj/dicomifier) Environment name: `Dicomifier`; Purpose: Conversion of Bruker data to NIfTI. Only needed for dMRI data acquired with Bruker scanner - on rodents or organoids for example.
 
- [**dcm2niix**](https://github.com/rordenlab/dcm2niix) Environment name: `niix2bids`; Purpose: Conversion of Siemens data to NIfTI. Only needed for dMRI data acquired with human Siemens scanner.

- [**SwissKnife**](https://github.com/QuentinUhl/graymatter_swissknife) Environment name: `SwissKnife`; Purpose: Apply microstructural models to the dMRI data and dMRS data. Needed to apply NEXI, SANDI or SMEX on dMRI data and to fit sphere and stick models to dMRS data.
  
- [**ANTS**](https://github.com/ANTsX/ANTsPy) Environment name: `ants`; Purpose: python interface to ANTs. Note: although ANTs is installed and accessible from the command line, this Conda environment provides the Python API and additional utilities required for generating a NIfTI representation of the MRS voxel when dMRS data are present.
 
- [**FSL MRS**](https://open.win.ox.ac.uk/pages/fsl/fsl_mrs/) Environment name: `fsl_mrs`; Purpose: analyze dMRS data

- [**RodentSkullStrip UNET**](https://github.com/CAMRIatUNC/RodentMRISkullStripping) Environment name: `RodentSkullStrip`; Purpose: skull strip of rodent data with U-NET. Need only for brain extraction of rodent data if this option is chosen (available options: RATS, UNET)
  
- **MATLAB**: 1) Required if using the MATLAB-based denoising options, with the [**MPPCA**](https://github.com/Neurophysics-CFIN/MP-PCA-Denoising) and [**tMPPCA**](https://github.com/Neurophysics-CFIN/Tensor-MP-PCA) toolboxes. If MATLAB is not available, denoising can instead be performed using **MRtrix** or **DESIGNER** (see `Step3.py`).  2) The pipeline previously also relied on [**md-dmri-master**](https://github.com/markus-nilsson/md-dmri/tree/master) and [**SPM12**](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) to compute *microscopic FA (µFA)* when *STE* data were acquired. These dependencies are now commented out in the scripts, as a *Python implementation of the µFA computation* has been integrated.


### 3. MATLAB runtime

The dMRS processing codes are provided as *compiled MATLAB executables*.  
To run them, you **do not need a MATLAB license**, but you must install the MATLAB Runtime (R2025a): https://ch.mathworks.com/products/compiler/matlab-runtime.html


 <br> 
 
## AUTHORS

Authors: Rita Oliveira & Malte Brammerloh

Supervisor: Ileana Jelescu

Microstructure Mapping Lab (mic-map),
Department of Radiology,
Lausanne University Hospital and University of Lausanne (CHUV),
Rue Pépinet 3, 1003 Lausanne, Switzerland

Email: ana.veiga-de-oliveira@chuv.ch, malte.brammerloh@chuv.ch 

Last updated: January 2026

