"""Template script for the generation of an atlas. Note that the script
has to be renamed to match the name of the atlas (e.g. allen_mouse.py)
"""

__version__ = "0"  # will be used to set minor version of the atlas

from pathlib import Path
import zipfile
import imio
import pandas as pd
import scipy
import numpy as np
from bg_atlasgen.wrapup import wrapup_atlas_from_data

from bg_atlasapi import utils
#from bg_atlasgen.mesh_utils import create_region_mesh, Region
from bg_atlasgen.wrapup import wrapup_atlas_from_data
from bg_atlasapi.structure_tree_util import get_structures_tree

def create_atlas(working_dir, resolution):
    """Function to generate source data for an atlas.

    Parameters
    ----------
    working_dir : Path object
        Path where atlas will be created.
    resolution :
        Resolution of the atlas, in um.

    Returns
    -------
    Path object
        Path to the final compressed atlas file.

    """

    ATLAS_NAME = "oldenburg_blackcap"
    SPECIES = "Sylvia atricapilla"
    ATLAS_LINK = ""
    CITATION = "unpublished"
    ORIENTATION = "rai"
    ATLAS_FILE_URL = "https://www.dropbox.com/sh/qxvf1fwxdq96vot/AADRRBaWCef9xct4McwHpwkWa?dl=1"

    # Temporary folder for  download:
    download_dir_path = working_dir / "downloads"
    download_dir_path.mkdir(exist_ok=True)
    atlas_files_dir = download_dir_path / "atlas_files"

    # Download atlas_file
    utils.check_internet_connection()
    destination_path = download_dir_path / "atlas_download.zip"
    print("Downloading files")
    utils.retrieve_over_http(ATLAS_FILE_URL, destination_path)
    zf = zipfile.ZipFile(destination_path, "r")
    zf.extractall(atlas_files_dir)

    template_file = atlas_files_dir / "101010_ds_SW_BC74white_220217_120749_10_10_ch04_chan_4_blue_raw_oriented.nii"
    structures_file = atlas_files_dir / "structures.csv"

    # Katrin & Isabelle are currently annotating separately so there 2 sets of files.
    # Katrin also made one for the whole brain. ITSnap can't have overlapping labels so it has be by itself.
    # But actually can't have overlapping structures.
    annotations_file_kh = atlas_files_dir / "BC74white_100um_annotations_KH_210722.nii"
    annotations_file_kh_brain = atlas_files_dir / "BC74white_brain_outline_correction_KH_230522.nii"
    annotations_file_im = atlas_files_dir / "BC74white_100um_annoations_IM_200522.nii.gz"

    itksnap_label_file_kh = atlas_files_dir / "Label_descriptions_BC_brainatlas_KH_21072022.txt"
    itksnap_label_file_kh_brain = atlas_files_dir / "Label_description_brain_outline_KH_230522.txt"
    itksnap_label_file_im = atlas_files_dir / "Label_description_BC74white_IM_200522.txt"

    # do stuff to create the atlas

    # ---------------------------------------------------------------------------- #
    #                             REFERENCE VOLUME                                 #
    # ---------------------------------------------------------------------------- #
    native_resolution = 10  #um
    scale = native_resolution / resolution
    scaling = (scale, scale, scale)

    print("Loading template")
    template_volume = imio.load_any(template_file)
    # Normalise with clipping to remove noise and fix contrast.
    print("Cleaning template")
    template_volume = clean_norm(template_volume, dtype=np.dtype(np.uint16), clean=True)
    # imio.load_any does not support scaling of nii, so need to do it here.
    print("Scaling template")
    template_volume = scipy.ndimage.zoom(template_volume,
                                         scaling,
                                         order=1,
                                         mode='nearest')

    # ---------------------------------------------------------------------------- #
    #                             ANNOTATED VOLUME                                 #
    # ---------------------------------------------------------------------------- #

    # Annotation is 100um in AP.
    print("Loading annotations")
    annotated_volume_100_kh = imio.load_any(annotations_file_kh, (1, 1, 1))
    annotated_volume_100_kh_brain = imio.load_any(annotations_file_kh_brain, (1, 1, 1))
    annotated_volume_100_im = imio.load_any(annotations_file_im, (1, 1, 1))

    # Combine KH & IM's annotations (volumes).
    n_structs_kh = np.sum(np.unique(annotated_volume_100_kh) > 0)
    n_structs_kh_brain = np.sum(np.unique(annotated_volume_100_kh_brain) > 0)
    if n_structs_kh_brain > 1:
        raise Exception("There was more than 1 structure in the outline volume")
    n_structs_im = np.sum(np.unique(annotated_volume_100_im) > 0)
    # Offset the outline structure by the number made by Katrin.
    annotated_volume_100_kh_brain[annotated_volume_100_kh_brain > 0] = annotated_volume_100_kh_brain[annotated_volume_100_kh_brain > 0] + n_structs_kh
    # Offset all Isabelle's structures by the number made by Katrin, + 1 for outline
    annotated_volume_100_im[annotated_volume_100_im > 0] = annotated_volume_100_im[annotated_volume_100_im > 0] + n_structs_kh + 1

    # Combine.
    print("Combining annotations")
    # There could be some overlap between KH & IM structures. Take KH when this happens.
    annotated_volume_100_im[annotated_volume_100_kh > 0] = 0
    annotated_volume_100_comb = annotated_volume_100_kh + annotated_volume_100_im
    # Combine the whole brain. Just use it to fill in the gaps, because you can't have overlapping structures.
    annotated_volume_100_kh_brain[annotated_volume_100_comb > 0] = 0
    annotated_volume_100_comb = annotated_volume_100_comb + annotated_volume_100_kh_brain

    # Resample to required resolution in AP with nearest neighbour so as to not mess up the labels.
    native_annotated_ap_res = 100  # um
    annotated_ap_scale = (native_annotated_ap_res / resolution)
    print("Scaling annotations")
    annotated_volume_comb = scipy.ndimage.zoom(annotated_volume_100_comb,
                                              (scale, annotated_ap_scale, scale),
                                              order=0,
                                              mode='nearest')


    # Create a new empty annoated volume to be filled out later.
    annotated_volume = np.full_like(annotated_volume_comb, 0)

    # ---------------------------------------------------------------------------- #
    #                             STRUCTURES HIERARCHY                             #
    # ---------------------------------------------------------------------------- #

    print("Loading ITKSnap data")
    # Load the ITKSnap description file to get the colours.
    n_rows_header = 14
    df_itksnap_kh = pd.read_csv(itksnap_label_file_kh,
                                delim_whitespace=True,
                                skiprows=n_rows_header,
                                names=["IDX", "-R-", "-G-", "-B-", "-A-", "VIS", "MESH-VIS", "LABEL"])
    # Drop first row, it's just an ITKSnap place holder for clear labels.
    df_itksnap_kh = df_itksnap_kh.iloc[1:, :]

    if n_structs_kh != df_itksnap_kh.shape[0]:
        raise Exception("The number of KH labels in the volume did not match the number of labels in the descriptions.")

    df_itksnap_kh_brain = pd.read_csv(itksnap_label_file_kh_brain,
                                        delim_whitespace=True,
                                        skiprows=n_rows_header,
                                        names=["IDX", "-R-", "-G-", "-B-", "-A-", "VIS", "MESH-VIS", "LABEL"])
    # Drop first row, it's just an ITKSnap place holder for clear labels.
    df_itksnap_kh_brain = df_itksnap_kh_brain.iloc[1:, :]

    df_itksnap_im = pd.read_csv(itksnap_label_file_im,
                                delim_whitespace=True,
                                skiprows=n_rows_header,
                                names=["IDX", "-R-", "-G-", "-B-", "-A-", "VIS", "MESH-VIS", "LABEL"])
    # Drop first row, it's just an ITKSnap place holder for clear labels.
    df_itksnap_im = df_itksnap_im.iloc[1:, :]

    if n_structs_im != df_itksnap_im.shape[0]:
        raise Exception("The number of IM labels in the volume did not match the number of labels in the descriptions.")

    # Offset the outline structure by the number made by Katrin
    df_itksnap_kh_brain["IDX"] = df_itksnap_kh_brain["IDX"] + n_structs_kh
    # Offset all Isabelle's structures by the number made by Katrin, + 1 for outline
    df_itksnap_im["IDX"] = df_itksnap_im["IDX"] + n_structs_kh + 1

    # Combine KH & IM's annotations (ITKSnap files)
    df_itksnap = pd.concat([df_itksnap_kh, df_itksnap_kh_brain, df_itksnap_im])

    # Make sure there are no duplicates.
    if np.sum(df_itksnap.duplicated(subset="LABEL")) > 0:
        raise Exception("There were duplicate labels across KH & IM ITKSnap descriptions")

    root_id = 1  # id of the root structure

    # Parse region names & hierarchy
    df = pd.read_csv(structures_file)

    # split by "/" and convert list of strings to list of ints
    df["structure_id_path"] = (
        df["structure_id_path"]
            .str.split(pat="/")
            .map(lambda x: [int(i) for i in x[1:-1]])
    )
    structures = df.to_dict("records")
    for structure in structures:
        if structure["id"] == root_id:
            # root doesn't have a parent or color.
            structure["structure_id_path"] = [root_id]
            structure.update({"rgb_triplet": [255, 255, 255]})
        else:
            # Structures have a parent and a color.
            structure["structure_id_path"].append(structure["id"])
            struc_index = df_itksnap["LABEL"] == structure['acronym']
            struc_red = int(df_itksnap.loc[struc_index, "-R-"].values[0])
            struc_blue = int(df_itksnap.loc[struc_index, "-B-"].values[0])
            struc_green = int(df_itksnap.loc[struc_index, "-G-"].values[0])
            structure.update({"rgb_triplet": [struc_red, struc_blue, struc_green]})

            # The label index in ITKSnap does not match the index from structures.
            # So look up the ITKSnap index and set it the structures index.
            itksnap_index = int(df_itksnap.loc[struc_index, "IDX"].values[0])
            annotated_volume[annotated_volume_comb == itksnap_index] = int(structure["id"])





    # ---------------------------------------------------------------------------- #
    #                             MESHES                                           #
    # ---------------------------------------------------------------------------- #

    meshes_dict = {}  # dictionary of files with region meshes

    # # Mesh creation
    # closing_n_iters = 2
    # decimate_fraction = 0.2
    # smooth = False  # smooth meshes after creation
    # start = time.time()
    # if PARALLEL:
    #
    #     pool = mp.Pool(mp.cpu_count() - 2)
    #
    #     try:
    #         pool.map(
    #             create_region_mesh,
    #             [
    #                 (
    #                     meshes_dir_path,
    #                     node,
    #                     tree,
    #                     labels,
    #                     rotated_annotations,
    #                     ROOT_ID,
    #                     closing_n_iters,
    #                     decimate_fraction,
    #                     smooth,
    #                 )
    #                 for node in tree.nodes.values()
    #             ],
    #         )
    #     except mp.pool.MaybeEncodingError:
    #         pass  # error with returning results from pool.map but we don't care
    # else:
    #     for node in track(
    #             tree.nodes.values(),
    #             total=tree.size(),
    #             description="Creating meshes",
    #     ):
    #         create_region_mesh(
    #             (
    #                 meshes_dir_path,
    #                 node,
    #                 tree,
    #                 labels,
    #                 rotated_annotations,
    #                 ROOT_ID,
    #                 closing_n_iters,
    #                 decimate_fraction,
    #                 smooth,
    #             )
    #         )

    # ---------------------------------------------------------------------------- #
    #                          ADDITIONAL REFERENCES                               #
    # ---------------------------------------------------------------------------- #

    # Put here additional reference stacks
    # (different genotypes, filtered volumes, etc)
    additional_references = dict()

    # ---------------------------------------------------------------------------- #
    #                              WRAP UP ATLAS                                   #
    # ---------------------------------------------------------------------------- #

    output_filename = wrapup_atlas_from_data(
        atlas_name=ATLAS_NAME,
        atlas_minor_version=__version__,
        citation=CITATION,
        atlas_link=ATLAS_LINK,
        species=SPECIES,
        resolution=(resolution,) * 3,  # if isotropic - highly recommended
        orientation=ORIENTATION,
        root_id=root_id,
        reference_stack=template_volume,
        annotation_stack=annotated_volume,
        structures_list=structures,
        meshes_dict=meshes_dict,
        working_dir=working_dir,
        additional_references=additional_references,
        hemispheres_stack=None,
        cleanup_files=False,
        compress=False,
    )

    return output_filename

def clean_norm(image, dtype, clean=True, clip_prc_lo=0.3, clip_prc_hi=99.7):

    if clean:
        # Remove low and high clipping pixels.
        clip_lo = np.percentile(image, clip_prc_lo)
        clip_hi = np.percentile(image, clip_prc_hi)
        image[image < clip_lo] = clip_lo
        image[image > clip_hi] = clip_hi

    image = norm_type(image, dtype)

    return image

def norm_type(vol, dtype):
    vol_min = np.min(vol)
    vol_max_ = np.max(vol)

    vol = ((vol - vol_min) / (vol_max_ - vol_min)) \
          * np.iinfo(dtype).max

    return vol.astype(dtype, copy=False)

# To test stuff locally:
if __name__ == "__main__":
    resolution = 50  # some resolution, in microns

    # Generated atlas path:
    bg_root_dir = Path.home() / "brainglobe_workingdir" / "blackcap"
    bg_root_dir.mkdir(exist_ok=True, parents=True)

    create_atlas(bg_root_dir, resolution)
