"""
Atlas handling utilities
========================

This module centralizes everything that is atlas-specific:
- how an atlas file is prepared to match our anatomical space
- how atlas label files are parsed into a usable table/dict
- how ROI names used in the pipeline map to atlas region IDs (including ROI grouping)
    
# Definitions:
    - An ATLAS is defined as a file containing labeled regions and must be located in
    folders matching the pattern `*Atlas_*`.

    - A TPM (tissue probability map) atlas is another type of atlas, consisting of a 5D volume
    containing tissue probability maps. By convention, the dimensions correspond to:
        1) GM (gray matter)
        2) WM (white matter)
        3) CSF
    The 4th and 5th volumes typically correspond to skull and non-brain tissues, but
    this may vary depending on the TPM atlas used.
    TPMs must be located in folders matching the pattern `*TPM_*`.


# What you need to do:
When a NEW ATLAS is added, please update the three functions below 
(add a new `if atlas_name == "..."` condition where needed):

1) prepare_atlas()
       An ATLAS must include:
          * an anatomical template image ('*template_brain*')
          * an atlas image with integer region labels ('*atlas*')
       A TPM atlas must include: 
          * an anatomical template image ('*template_brain*')
          * an map of tissue probabilities ('*TPM*')
   Then you can adapt the atlas to better match our anatomical data type / reference space. Typical steps:
   - reorient (e.g., to RAS)
   - resample to the anatomical voxel size/grid (nearest-neighbour for label images)
   - crop/pad (if required)
   - optional cleanup (remove unwanted label IDs, merge IDs, etc.)

2) prepare_atlas_labels()
       An ATLAS must include:
           * a label file mapping region IDs to region names ('*label*')
   Adapt the way to load the labels. There are several atlas label-file format 
   (txt/csv/xml/json, separators, headers). If none of the existing conditions fits 
   your atlas label-file format, add a new one.

3) create_ROI_mask()
   Define how ROI strings used in the main pipeline map to atlas region IDs/names for ROI-bases analysis.
   Grouped ROIs should be defined here (e.g., "CTX" = ["M1","S1","V1"] or
   "Thalamus" = [ID1, ID2, ID3]). Keep atlas IDs and atlas region names contained
   here (do not scatter atlas-specific mappings across the pipeline).

    
Last updated: Jan 2026
@author: Rita O
"""

import os
import re
import shutil
import numpy as np
import nibabel as nib
import nibabel.processing as nip
import pandas as pd
import glob as glob
       
##### ATLAS HANDLE #####


def prepare_atlas(atlas_name, atlas_folder, atlas_type):
    import glob

    if atlas_name== 'Atlas_WHS_v4' and atlas_type=='atlas':
        
        # Define atlas 
        atlas      = glob.glob(os.path.join(atlas_folder, atlas_name, '*atlas.nii.gz'))[0]
        template   = glob.glob(os.path.join(atlas_folder, atlas_name, '*template_brain.nii.gz'))[0]
       
        # If atlas wasn't refined before, do it
        if (not os.path.exists(atlas.replace('.nii.gz', '_crop_lowres.nii.gz'))) and (not os.path.exists(template.replace('.nii.gz', '_crop_lowres.nii.gz'))):

            # Remove extra regions to make it more like brain only
            img  = nib.load(atlas)
            data_a = img.get_fdata()
            data_a[data_a == 42] = 0
            data_a[data_a == 41] = 0
            data_a[data_a == 45] = 0
            data_a[data_a == 76] = 0
            img  = nib.load(template)
            data_t = img.get_fdata()           
            mask = (data_a != 0).astype(np.uint8)
            data_t = data_t*mask
            nib.save(nib.Nifti1Image(data_a, img.affine), atlas.replace('.nii.gz', '_crop.nii.gz'))
            nib.save(nib.Nifti1Image(data_t, img.affine), template.replace('.nii.gz', '_crop.nii.gz'))
            
            for image in (atlas,template):
             
             # Crop template/atlas - otherwise too much data to register
             img  = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
             data = img.get_fdata()
             masked_data = np.zeros_like(data)
             masked_data[:, 230:840, :] = data[:, 230:840, :]
             nib.save(nib.Nifti1Image(masked_data, img.affine), image.replace('.nii.gz', '_crop.nii.gz'))
             
             # Downsample template/atlas to avoid segmentation faults
             input_img = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
             if image==atlas:
                 resampled_img = nip.resample_to_output(input_img, [0.05, 0.5, 0.05],order=0)
             elif image==template:
                 resampled_img = nip.resample_to_output(input_img, [0.05, 0.5, 0.05])
             nib.save(resampled_img,  image.replace('.nii.gz', '_crop_lowres.nii.gz')) 
           
        # Define atlas 
        atlas      = atlas.replace('.nii.gz', '_crop_lowres.nii.gz')
        template   = template.replace('.nii.gz', '_crop_lowres.nii.gz')
        
    elif atlas_name == 'TPM_C57Bl6'  and atlas_type=='TPM':
    
        # Define TPM 
        atlas    = glob.glob(os.path.join(atlas_folder, atlas_name, '*TPM_C57Bl6_n30.nii'))[0]
        template = glob.glob(os.path.join(atlas_folder, atlas_name, '*C57Bl6_T2_n10_template_brain.nii'))[0]
        
        # If atlas wasn't refined before, do it
        if (not os.path.exists(atlas.replace('.nii.gz', '_rescaled.nii.gz'))) and (not os.path.exists(template.replace('.nii.gz', '_rescaled.nii.gz'))):

            for image in (atlas,template):
                
                # Correct scale
                img = nib.load(image)
                data = img.get_fdata()
                affine = img.affine.copy()
                affine[:3, :3] /= 10  # Correct for the scale factor
                corrected_img = nib.Nifti1Image(data, affine, img.header)
                
                # Save the rescaled image
                nib.save(corrected_img, image.replace('.nii', '_rescaled.nii'))
                
                # Separate TPMs into different files
                if image==atlas:
                    
                    # get path of rescaled image
                    input_path = image.replace('.nii', '_rescaled.nii')
                    out_path = image.replace('.nii', '_rescaled_vol_')
                    
                    call = [f'fslsplit',
                            f'{input_path}',
                            f'{out_path} -t']
                    print(' '.join(call))
        
                    os.system(' '.join(call))
    
    
                #     # Resample each 3D volume of the 4D TPM atlas
                #     data = input_img.get_fdata()
                #     affine = input_img.affine
                #     header = input_img.header
            
                #     resampled_volumes = []
                #     for dim in range(data.shape[-1]):
                #        vol_3d = nib.Nifti1Image(data[..., dim], input_img.affine)
                #        resampled_vol = nib.processing.resample_to_output(vol_3d, voxel_sizes=[0.08, 0.5, 0.08], order=3)
                #        nib.save(resampled_img, image.replace('.nii', '_rescaled_lowres.nii'))
        
                #        resampled_volumes.append(resampled_vol.get_fdata())
            
                #     resampled_data = np.stack(resampled_volumes, axis=-1)
                #     resampled_img = nib.Nifti1Image(resampled_data, resampled_vol.affine, resampled_vol.header)
            
                # elif image==template:
                #     # Resample the 3D template directly
                #     resampled_img = nibabel.processing.resample_to_output(input_img, [0.08, 0.5, 0.08])
             
                # # Save final lowres version
                # nib.save(vol_3d, image.replace('.nii', '_rescaled_lowres.nii'))
               
        # Define TPM 
        atlas      = atlas.replace('.nii', '_rescaled.nii')
        template   = template.replace('.nii', '_rescaled.nii')
  
    elif (atlas_name== 'Atlas_postnatal_P24' or atlas_name== 'Atlas_postnatal_P40' or atlas_name== 'Atlas_postnatal_P80')  and atlas_type=='atlas':

        # Define atlas 
        atlas      = glob.glob(os.path.join(atlas_folder, atlas_name, '*atlas.nii.gz'))[0]
        template   = glob.glob(os.path.join(atlas_folder, atlas_name, '*template_brain.nii.gz'))[0]
        
        # If atlas wasn't refined before, do it
        if (not os.path.exists(atlas.replace('.nii.gz', '_crop_lowres.nii.gz'))) and (not os.path.exists(template.replace('.nii.gz', '_crop_lowres.nii.gz'))):
            
            for image in (atlas,template):
                
                # Crop template/atlas - otherwise too much data to register
                img  = nib.load(image)
                data = img.get_fdata()
                masked_data = np.zeros_like(data)
                masked_data[:, 190:1075, :] = data[:, 190:1075, :]
                nib.save(nib.Nifti1Image(masked_data, img.affine), image.replace('.nii.gz', '_crop.nii.gz'))
                
                # Downsample template/atlas to avoid segmentation faults
                input_img = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
                if image==atlas:
                   resampled_img = nip.resample_to_output(input_img, [0.1, 0.1, 0.1],order=0)
                elif image==template:
                   resampled_img = nip.resample_to_output(input_img, [0.1, 0.1, 0.1])
                nib.save(resampled_img,  image.replace('.nii.gz', '_crop_lowres.nii.gz')) 
      
        # Define atlas 
        atlas      = atlas.replace('.nii.gz', '_crop_lowres.nii.gz')
        template   = template.replace('.nii.gz', '_crop_lowres.nii.gz')

    elif atlas_name== 'Atlas_WHS_mouse' and atlas_type=='atlas':

         # Define atlas 
         atlas      = glob.glob(os.path.join(atlas_folder, atlas_name, '*atlas.nii.gz'))[0]
         template   = glob.glob(os.path.join(atlas_folder, atlas_name, '*template_brain.nii.gz'))[0]
         
         # If atlas wasn't refined before, do it
         if (not os.path.exists(atlas.replace('.nii.gz', '_crop_lowres.nii.gz'))) and (not os.path.exists(template.replace('.nii.gz', '_crop_lowres.nii.gz'))):
             
             for image in (atlas,template):
                 
                 # Crop template/atlas - otherwise too much data to register
                 img  = nib.load(image)
                 data = img.get_fdata()
                 masked_data = np.zeros_like(data)
                 masked_data[:, 59:739, :] = data[:, 59:739, :]
                 nib.save(nib.Nifti1Image(masked_data, img.affine), image.replace('.nii.gz', '_crop.nii.gz'))
                 
                 # Downsample template/atlas to avoid segmentation faults
                 input_img = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
                 if image==atlas:
                    resampled_img = nip.resample_to_output(input_img, [0.1, 0.1, 0.1],order=0)
                 elif image==template:
                    resampled_img = nip.resample_to_output(input_img, [0.1, 0.1, 0.1])
                 nib.save(resampled_img,  image.replace('.nii.gz', '_crop_lowres.nii.gz')) 
       
         # Define atlas 
         atlas      = atlas.replace('.nii.gz', '_crop_lowres.nii.gz')
         template   = template.replace('.nii.gz', '_crop_lowres.nii.gz')
         
    # If no adjustments need to be made in atlas it just reads the files
    elif atlas_type=='atlas':
   
       # Define atlas 
       atlas      = glob.glob(os.path.join(atlas_folder, atlas_name, '*atlas*'))[0]
       template   = glob.glob(os.path.join(atlas_folder, atlas_name, '*template_brain*'))[0]
   
    # If no adjustments need to be made in atlas it just reads the files
    elif atlas_type=='TPM':
     
        # Define TPM 
        atlas_list    = glob.glob(os.path.join(atlas_folder, atlas_name, '*TPM*'))
        atlas = [a for a in atlas_list if 'vol' not in os.path.basename(a)][0]
        template = glob.glob(os.path.join(atlas_folder,atlas_name, '*template_brain*'))[0]
        
        # Separate TPMs into different files
        out_path = atlas.replace('.nii', '_vol_')     
        call = [f'fslsplit',
                f'{atlas}',
                f'{out_path} -t']
        print(' '.join(call))
        os.system(' '.join(call))
      
    return atlas, template

def prepare_atlas_labels(atlas_name, atlas_label_path):
    import glob
    import random
    
    #atlas_label_path = glob.glob(os.path.join(atlas_folder, atlas_name, '*label*'))[0]

    # Handle Atlas_WHS_v4 which have the labels in a specific formar 
    if atlas_name== 'Atlas_WHS_v4' or atlas_label_path.endswith('.label'):
        atlas_labels = pd.read_csv(atlas_label_path, sep=r'\s+', skiprows=14, header=None,
                                   names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], quotechar='"')
     
    # Handle DKT-style .txt file
    elif atlas_label_path.endswith('.txt') and atlas_name== 'Atlas_DKT':
        with open(atlas_label_path, 'r') as f:
            content = f.read()
    
        # Extract label entries like [1002, "left caudal anterior cingulate"]
        matches = re.findall(r'\[\s*(\d+)\s*,\s*"([^"]+)"\s*\]', content)
    
        labels = []
        for match in matches:
            index = int(match[0])
            name = match[1].strip()
    
            # Assign random RGB colors
            R, G, B = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
            A, VIS, MSH = 1, 1, 0
    
            labels.append({
                'IDX': index,
                'R': R,
                'G': G,
                'B': B,
                'A': A,
                'VIS': VIS,
                'MSH': MSH,
                'LABEL': name
            })
    
        atlas_labels = pd.DataFrame(labels)
        atlas_labels.sort_values(by='IDX', inplace=True)
    
    # Handle DKT-style .txt file
    elif atlas_label_path.endswith('.txt') and 'Atlas_postnatal' in atlas_name:
        labels = []

        with open(atlas_label_path, 'r') as f:
            for line in f:
                line = line.strip()
    
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
    
                parts = line.split()
    
                # Expect at least: IDX LABEL R G B A
                if len(parts) < 6:
                    continue
    
                index = int(parts[0])
                name  = parts[1]
    
                R = int(parts[2])
                G = int(parts[3])
                B = int(parts[4])
                A = int(parts[5])
    
                VIS = 1
                MSH = 0
    
                labels.append({
                    'IDX': index,
                    'R': R,
                    'G': G,
                    'B': B,
                    'A': A,
                    'VIS': VIS,
                    'MSH': MSH,
                    'LABEL': name
                })
    
        atlas_labels = pd.DataFrame(labels)
        atlas_labels.sort_values(by='IDX', inplace=True)
    
    elif atlas_label_path.endswith('.ilf'):
       import xml.etree.ElementTree as ET

       def hex_to_rgb(hex_color: str):
           """'#RRGGBB' -> (R,G,B). Returns random if malformed."""
           if not hex_color:
               return (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
           hex_color = hex_color.strip()
           if hex_color.startswith('#'):
               hex_color = hex_color[1:]
           if len(hex_color) != 6:
               return (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
           return (int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16))

       tree = ET.parse(atlas_label_path)
       root = tree.getroot()

       labels = []

       # In .ilf, labels are nested <label .../> inside <structure>
       for lbl in root.findall(".//structure//label"):
           id_str = lbl.attrib.get("id", None)
           if id_str is None:
               continue

           try:
               index = int(id_str)
           except ValueError:
               continue

           # Prefer 'name', fallback to 'abbreviation'
           name = (lbl.attrib.get("name") or lbl.attrib.get("abbreviation") or "").strip()

           # color="#ccbcb0" etc.
           color_hex = lbl.attrib.get("color", None)
           R, G, B = hex_to_rgb(color_hex)

           A = 1
           VIS = 1
           MSH = 0

           labels.append({
               'IDX': index,
               'R': R,
               'G': G,
               'B': B,
               'A': A,
               'VIS': VIS,
               'MSH': MSH,
               'LABEL': name
           })

       atlas_labels = pd.DataFrame(labels)

       # De-duplicate in case repeated IDs appear (keep first)
       if not atlas_labels.empty:
           atlas_labels.sort_values(by='IDX', inplace=True)
           atlas_labels.drop_duplicates(subset=['IDX'], keep='first', inplace=True)
    
    # If none of the previous, assumes there is an xml file    
    else:
        import xml.etree.ElementTree as ET

        # Load the XML file
        tree = ET.parse(atlas_label_path)
        root = tree.getroot()
        
        # Extract label data
        labels = []
        for label in root.findall(".//label"):
            
            if 'index' in label.attrib:
                # Format 1: attribute-based
                index = int(label.attrib['index'])
                name = label.text.strip()
            else:
                # Format 2: element-based
                index_elem = label.find('index')
                name_elem = label.find('name')
                if index_elem is not None:
                    index = int(index_elem.text.strip())
                else:
                    continue  # Skip if no index
                name = name_elem.text.strip() if name_elem is not None else ''
        
            # Assign random RGB colors
            R = random.randint(0, 255)
            G = random.randint(0, 255)
            B = random.randint(0, 255)
        
            # Constants for alpha, visibility, mesh
            A = 1
            VIS = 1
            MSH = 0
        
            # for some reason this atlas have the labels shifted by 1 value, so we need to correct
            if 'Atlas_Juelich' in atlas_name: 
               index = index +1
               
            labels.append({
                'IDX': index,
                'R': R,
                'G': G,
                'B': B,
                'A': A,
                'VIS': VIS,
                'MSH': MSH,
                'LABEL': name
            })
        
        # Convert to DataFrame
        df = pd.DataFrame(labels)
        
        # Optional: Sort by index
        df.sort_values(by='IDX', inplace=True)

        atlas_labels = df


    return atlas_labels

def make_atlas_manual_organoid(organoid_mask,folder,label_path, cfg):
    
    atlas = organoid_mask.replace('_mask.nii.gz', '_atlas.nii.gz')
    
    candidate_masks = glob.glob(os.path.join(folder, '*organoid?_mask.nii.gz'))
    pattern = re.compile(r'organoid([A-Z])_mask\.nii\.gz$')
    organoid_masks = [f for f in candidate_masks if pattern.search(os.path.basename(f))]

    temp_files = []
    label_lines = []
    color_base = (255, 105, 180)  # Hot pink base

    for idx, mask_path in enumerate(organoid_masks, start=1):
        temp_file = mask_path.replace('.nii.gz', f'_temp.nii.gz')
        temp_files.append(temp_file)

        # Multiply mask by its unique label index
        call = ['fslmaths', mask_path, '-mul', str(idx), temp_file]
        os.system(' '.join(call))

        # Add a label line
        r, g, b = color_base  # You could vary color if needed
        match = re.search(r'(organoid[A-Za-z])_mask\.nii\.gz$', os.path.basename(mask_path))
        label = match.group(1) if match else 'unknown'
        label_lines.append(f'   {idx}   {r}  {g}  {b}        1  1  0    "{label}"\n')

    # Merge all the labeled masks
    exe = os.path.join(cfg["fsl_path"], "fslmaths")
    merged_command = [exe, temp_files[0]]
    for temp_file in temp_files[1:]:
        merged_command += ['-add', temp_file]
    merged_command.append(atlas)
    os.system(' '.join(merged_command))

    # Remove temporary files
    for temp_file in temp_files:
        os.remove(temp_file)

    # Write the label file
    header = """\
    ################################################
    # ITK-SnAP Label Description File
    # File format: 
    # IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
    # Fields: 
    #    IDX:   Zero-based index 
    #    -R-:   Red color component (0..255)
    #    -G-:   Green color component (0..255)
    #    -B-:   Blue color component (0..255)
    #    -A-:   Label transparency (0.00 .. 1.00)
    #    VIS:   Label visibility (0 or 1)
    #    MSH:   Label mesh visibility (0 or 1)
    #  LABEL:   Label description 
    ################################################
    """

    with open(label_path, 'w', encoding='utf-8') as f:
        f.write(header + '\n')
        f.writelines(label_lines)
        
   
        
#def make_atlas_label_organoid(organoid_mask, non_organoid_mask,label_path):
    
    # temp_non_organoid = non_organoid_mask.replace('.nii.gz', '_temp.nii.gz')
    # call = ['fslmaths',
    #          f'{non_organoid_mask}',
    #          '-mul', '2',
    #          f'{temp_non_organoid}'
    #          ]
  
    # os.system(' '.join(call))
    
    
    # atlas = organoid_mask.replace('_mask.nii.gz', '_atlas.nii.gz')
    # call = [f'fslmaths',
    #         f'{organoid_mask}',
    #         f'-add {temp_non_organoid}',
    #         f'{atlas}']

    # os.system(' '.join(call))
    # remove_file(temp_non_organoid)
    
    # header = """\
    #             ################################################
    #             # ITK-SnAP Label Description File
    #             # File format: 
    #             # IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
    #             # Fields: 
    #             #    IDX:   Zero-based index 
    #             #    -R-:   Red color component (0..255)
    #             #    -G-:   Green color component (0..255)
    #             #    -B-:   Blue color component (0..255)
    #             #    -A-:   Label transparency (0.00 .. 1.00)
    #             #    VIS:   Label visibility (0 or 1)
    #             #    MSH:   Label mesh visibility (0 or 1)
    #             #  LABEL:   Label description 
    #             ################################################
    #             """

    # label_line = ['   1   255  105  180        1  1  0    "organoids"\n'
    #               '   2   255  105  180        1  1  0    "medium"\n']

    # with open(label_path, "w", encoding="utf-8") as f:
    #     f.write(header)
    #     f.write("\n")
    #     f.writelines(label_line)        
    
def create_ROI_mask(atlas, atlas_labels, TPMs, ROI, tpm_thr, bids_strc_reg):
 
     # Define tpms and threshold at 0.9
     atlas_shape = nib.load(atlas).shape
     
     # GM
     tpm_GM = [f for f in TPMs if 'GM' in f and f and os.path.exists(f)]
     tmp_GM = nib.load(tpm_GM[0]).get_fdata() > tpm_thr if tpm_GM else np.ones(atlas_shape, dtype=bool)
    
     # WM
     tpm_WM = [f for f in TPMs if 'WM' in f and f and os.path.exists(f)]
     tmp_WM = nib.load(tpm_WM[0]).get_fdata() > tpm_thr if tpm_WM else np.ones(atlas_shape, dtype=bool)
    
     # CSF
     tpm_CSF = [f for f in TPMs if 'CSF' in f and f and os.path.exists(f)]
     tmp_CSF = nib.load(tpm_CSF[0]).get_fdata() > tpm_thr if tpm_CSF else np.ones(atlas_shape, dtype=bool)
 

     if 'Atlas_WHS_v4' in atlas:
         # Define ROI labels for each ROI asked
         roi_definitions = {
            'hippocampus': ['hippocampus','Cornu ammonis 3','Cornu ammonis 2','Cornu ammonis 1','Dentate gyrus'],
            'M1': ['Primary motor area'],
            'M2': ['Secondary motor area'],
            'S1': ['Primary somatosensory'],
            'S2': ['Secondary somatosensory'],
            'V1': ['Primary visual'],
            'V2': ['Secondary visual'],
            'Parietal': ['Parietal'],
            'CC': ['corpus callosum'],
            'Amygdala': ['Amygdaloid'],
            'AC': ['commissure'],
            'Auditory': ['auditory'],
            'Cingulate': ['Cingulate','Prelimbic'],
            'Cochlea': ['Cochlea'],
            'Entorhinal': ['entorhinal'],
            'Caudade_Putamen':['Ventral striatal region','Caudate putamen'],
            'GP': ['Globus_Pallidus'],
            'Pallidum': ['Ventral pallidum','Bed nucleus of the stria terminalis'],
            'SN': ['Substantia nigra'],
            'Hypothalamus': ['Hypothalamic','Zona incerta, dorsal part','Zona incerta, ventral part','Zona incerta, caudal part'],
            'Inferior colliculus': ['Inferior colliculus, external cortex','Inferior colliculus, dorsal cortex','Inferior colliculus, central nucleus'],
            'Superior_colliculus': ['superior colliculus'],
            'MGN': ['medial geniculate body, ventral','medial geniculate body, dorsal'],
            'Nucleus accumbens':['Nucleus accumbens'],
            'Thal': ['thalamic nucleus', 'Pineal gland','pretectothalamic lamina','Pregeniculate nucleus','thalamic nuclear'],
            'CSF': ['Ventricular'],
            'Cereb GM': ['erebellum','cerebellar'],
            'Cereb WM': ['erebellum','cerebellar'],
            'insula': ['insula'],
            'VTA': ['Ventral tegmental area'],
            'WB': ['whole brain']
        } 
     elif 'Atlas_Juelich' in atlas:
        
         roi_definitions = {
            'hippocampus': ['ippocampus'],
            'V1': ['V1'],
            'V2': ['V2'],
            'premotor': ['Premotor'],
            'parietal': ['parietal'],
            'S1': ['Primary somatosensory'],
            'S2': ['Secondary somatosensory'],
            'M1': ['Primary motor'],
            'Broca': ['Broca'],
            'CC': ['callosal'],
            'WB': ['whole brain']
        } 
     elif 'Atlas_Neuromorphometrics' in atlas:
        
         roi_definitions = {
            'frontal': ['frontal cortex', 'frontal gyrus'],
            'precentral': ['precentral gyrus'],
            'postcentral': ['postcentral gyrus'],
            'occipital': ['occipital gyrus'],
            'parietal': ['parietal lobule'],
            'temporal': ['temporal lobe','temporal gyrus'],
            'WB': ['whole brain']
        } 
         
     elif 'Atlas_DKT' in atlas:
        
         roi_definitions = {
            'frontal': ['frontal'],
            'precentral': ['precentral'],
            'postcentral': ['postcentral'],
            'cuneus': [' cuneus'],
            'occipital': ['occipital'],
            'parietal': ['parietal'],
            'temporal': ['temporal'],
            'WB': ['whole brain']
        } 
     
     elif 'Atlas_postnatal_P24' in atlas  or 'Atlas_postnatal_P40' in atlas:
       
        roi_definitions = {
           'Isocortex': ['Isocortex'],
           'SN': ['Substantia_Nigra'],
           'Cereb GM': ['Cerebellum'],
           'Cereb WM': ['Cerebellum'],
           'Cerebellum': ['Cerebellum'],
           'Pallidum': ['Pallidum'],
           'Hypothalamus': ['Hypothalamus'],
           'Hippocampal Form': ['Hippocampal_Formation'],
           'Striatum':['Striatum'],
           'CC': ['Corpus_Callosum'],

       } 
        
     elif 'Atlas_WHS_mouse' in atlas:
      
        roi_definitions = {
        # Cortex / telencephalon
        'Isocortex': [
            'Cerebral cortex', 'Cerebral cortex'.replace(' ', '_'),
            'Cerebral cortex'
        ],
        'Olfactory_Bulb': [
            'olfactory bulb', 'olfactory_bulb'
        ],

        # Hippocampal system
        'Hippocampal_Formation': [
            'hippocampal formation', 'hippocampal_formation',
            'fimbria', 'fornix'
        ],

        # Basal ganglia / ventral striatum / basal forebrain
        'Basal_Ganglia': [
            'striatum',
            'caudate putamen', 'caudate_putamen',
            'caudate nucleus', 'caudate_nucleus',
            'globus pallidus', 'globus_pallidus',
            'nucleus accumbens', 'nucleus_accumbens',
            'septal nucleus', 'septal_nucleus',
            'basal forbrain', 'basal_forbrain'
        ],

        # Diencephalon
        'Thalamus': ['thalamus'],
        'Hypothalamus': ['hypothalamus'],
        'Epithalamus_Pineal': [
            'Epithalumus', 'epithalamus',   # (note: typo in the file: "Epithalumus")
            'pineal gland', 'pineal_gland'
        ],

        # Midbrain
        'Midbrain': [
            'midbrain', 'midbrain unsegmented', 'midbrain_unsegmented',
            'superior colliculus', 'superior_colliculus',
            'inferior colliculus', 'inferior_colliculus',
            'periaqueductal gray', 'periaqueductal_gray',
            'substantia nigra', 'substantia_nigra',   # appears twice in WHS labels (IDs 24 and 28)
            'interpeduncular nucleus', 'interpeduncular_nucleus',
            'fasciculus retroflexus', 'fasciculus_retroflexus',
            'un-named 01', 'un_named_01'
        ],

        # Hindbrain
        'Hindbrain': [
            'cerebellum',
            'pons',
            'medulla'
        ],

        # Major white matter bundles
        'WM_Tracts': [
            'corpus callosum', 'corpus_callosum',
            'internal capsule', 'internal_capsule',
            'anterior commissure', 'anterior_commissure'
        ],

        # Ventricles / CSF
        'Ventricles': [
            'lateral ventricles', 'lateral_ventricles',
            '3rd ventricle', '3rd_ventricle',
            '4th ventricle', '4th_ventricle'
        ],

        # Cranial nerves / sensory structures
        'Cranial_Nerves_Sensory': [
            'sensory system', 'sensory_system',
            'cranial nerves', 'cranial_nerves',
            'optic nerve', 'optic_nerve',
            'optic chiasm', 'optic_chiasm',
            'optic tract', 'optic_tract',
            'trigeminal', 'trigeminal tract', 'trigeminal_tract',
            'acoustic nerve', 'acoustic_nerve',
            'inner ear', 'inner_ear'
        ],

        # Spinal
        'Spinal_Cord': ['spinal cord', 'spinal_cord'],
    }
             
     elif 'anat_space_organoids' in atlas:
         
          roi_definitions = {
             'organoids': ['organoid'],
             'medium': ['medium'],

         } 
          
     else:
         print('>> Error: you did not defined which regions corresponds to each ROI. Please add a condition for that atlas in create_ROI_mask')
        
     # Load the atlas data
     template = nib.load(atlas)
     atlas_data = template.get_fdata()
      
     # Find matching indices for ROI
     if ROI=='WB':
         if  'Atlas_WHS_v4' in atlas:
            label_mask = ~atlas_labels["LABEL"].str.contains("olfactory", regex=True)
            match_idx = atlas_labels.loc[label_mask & (atlas_labels["IDX"] != 0), "IDX"].to_numpy()
         else:
            match_idx = atlas_labels["IDX"].to_numpy()[atlas_labels["IDX"].to_numpy() != 0]
     elif ROI == 'Thal':
         #  exclude "hypothalamic"
         ind_list = atlas_labels["LABEL"].str.contains('thalamic nucleus', regex=True) & ~atlas_labels["LABEL"].str.contains('Hypothalamic', regex=True) & ~atlas_labels["LABEL"].str.contains('unspecified', regex=True)
         match_idx = atlas_labels["IDX"][ind_list].to_numpy()
     else:
         ind_list = atlas_labels["LABEL"].str.contains("|".join(roi_definitions[ROI]), regex=True)
         match_idx = atlas_labels["IDX"][ind_list].to_numpy()
     
     # Create the mask for the ROI
     mask_indexes = np.isin(atlas_data, match_idx)
     masked_data = (mask_indexes > 0).astype(np.uint8)
    
     # Save the mask
     masked_img = nib.Nifti1Image(masked_data, affine=template.affine, header=template.header)
     nib.save(masked_img, bids_strc_reg.get_path(f'mask_{ROI}.nii.gz'))
    
     # Multiply by tissue probability map
     if ROI=='CC': # is WM
         masked_data = masked_data*tmp_WM
     elif ROI=='CSF':
         masked_data = masked_data*tmp_CSF
     elif ROI=='Cereb GM':
         masked_data = masked_data*tmp_GM 
     elif ROI=='Cereb WM':
         masked_data = masked_data*tmp_WM 
     elif ROI=='WB':
         masked_data = masked_data* (tmp_GM | tmp_WM | tmp_CSF)
     else:
         masked_data = masked_data*tmp_GM

     # Save the mask
     masked_img = nib.Nifti1Image(masked_data, affine=template.affine, header=template.header)
     nib.save(masked_img, bids_strc_reg.get_path(f'mask_{ROI}.nii.gz'))


     return masked_data
 