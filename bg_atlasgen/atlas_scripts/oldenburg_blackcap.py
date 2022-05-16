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
    # utils.retrieve_over_http(ATLAS_FILE_URL, destination_path)
    #
    # zf = zipfile.ZipFile(destination_path, "r")
    # zf.extractall(atlas_files_dir)

    template_file = atlas_files_dir / "101010_ds_SW_BC74white_220217_120749_10_10_ch04_chan_4_blue_raw_oriented.nii"
    structures_file = atlas_files_dir / "structures.csv"
    itksnap_label_file = atlas_files_dir / "Label_descriptions_BC74white_KH_12042022.txt"
    annotations_file = atlas_files_dir / "BC74white_100um_annotations_120422.nii"

    # do stuff to create the atlas

    # ---------------------------------------------------------------------------- #
    #                             REFERENCE VOLUME                                 #
    # ---------------------------------------------------------------------------- #
    native_resolution = 10 #um
    scale = native_resolution / resolution
    scaling = (scale, scale, scale)

    template_volume = imio.load_any(template_file)
    # Normalise with clipping to remove noise and fix contrast.
    template_volume = clean_norm(template_volume, dtype=np.dtype(np.uint16), clean=False)
    # imio.load_any does not support scaling of nii, so need to do it here.
    template_volume = scipy.ndimage.zoom(template_volume,
                                         scaling,
                                         order=1,
                                         mode='nearest')

    # ---------------------------------------------------------------------------- #
    #                             ANOTATED VOLUME                                  #
    # ---------------------------------------------------------------------------- #

    # Annotation is 100um in AP.
    annotated_volume_100 = imio.load_any(annotations_file, (1, 1, 1))

    # Resample to required resolution in AP with nearest neighbour so as to not mess up the labels.
    native_annotated_ap_res = 100  # um
    annotated_ap_scale = (native_annotated_ap_res / resolution)
    annotated_volume = scipy.ndimage.zoom(annotated_volume_100,
                                          (scale, annotated_ap_scale, scale),
                                          order=0,
                                          mode='nearest')

    # ---------------------------------------------------------------------------- #
    #                             STRUCTURES HIERARCHY                             #
    # ---------------------------------------------------------------------------- #


    # Load the ITKSnap description file to get the colours.
    n_rows_header = 14
    df_itksnap = pd.read_csv(itksnap_label_file,
                             delim_whitespace=True,
                             skiprows=n_rows_header,
                             names=["IDX", "-R-", "-G-", "-B-", "-A-", "VIS", "MESH-VIS", "LABEL"])

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
    #structures[0000]["structure_id_path"] = [root_id]
    for structure in structures:
        if structure["id"] == root_id:
            # root doesn't have a parent or color.
            structure["structure_id_path"] = [root_id]
            structure.update({"rgb_triplet": [255, 255, 255]})
        else:
            # Structures have a parent and a color.
            structure["structure_id_path"].append(structure["id"])
            struc_index = df_itksnap["LABEL"] == structure['name']
            struc_red = int(df_itksnap.loc[struc_index, "-R-"].values[0])
            struc_blue = int(df_itksnap.loc[struc_index, "-B-"].values[0])
            struc_green = int(df_itksnap.loc[struc_index, "-G-"].values[0])
            structure.update({"rgb_triplet": [struc_red, struc_blue, struc_green]})




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

def clean_norm(image, dtype, clean=True, clip_prc_lo=0.1, clip_prc_hi=99.9):

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
