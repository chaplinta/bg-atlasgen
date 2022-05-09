"""Template script for the generation of an atlas. Note that the script
has to be renamed to match the name of the atlas (e.g. allen_mouse.py)
"""

__version__ = "0"  # will be used to set minor version of the atlas

from pathlib import Path
import zipfile
import imio
import pandas as pd

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

    ATLAS_NAME = "blackcap"
    SPECIES = "Sylvia atricapilla"
    ATLAS_LINK = ""
    CITATION = ""
    ORIENTATION = "asr"
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

    template_file = atlas_files_dir / "1001010_ds_SW_BC74white_220217_120749_10_10_ch04_chan_4_blue_raw_oriented.nii"
    structures_file = atlas_files_dir / "structures.csv"
    itksnap_label_file = atlas_files_dir / "Label_descriptions_BC74white_KH_12042022.txt"
    annotations_file = atlas_files_dir / "BC74white_100um_annotations_120422.nii"

    # do stuff to create the atlas

    # ---------------------------------------------------------------------------- #
    #                             REFERENCE VOLUME                                 #
    # ---------------------------------------------------------------------------- #
    scaling = (1, 1, 1)
    template_volume = imio.load_any(
        template_file,
        scaling[1],
        scaling[2],
        scaling[0],
    )

    # ---------------------------------------------------------------------------- #
    #                             ANOTATED VOLUME                                  #
    # ---------------------------------------------------------------------------- #

    annotated_volume = None  # volume with structures annotations

    # ---------------------------------------------------------------------------- #
    #                             STRUCTURES HIERARCHY                             #
    # ---------------------------------------------------------------------------- #

    root_id = 1  # id of the root structure

    # Parse region names & hierarchy
    df = pd.read_csv(structures_file)

    # split by "/" and convert list of strings to list of ints
    df["structure_id_path"] = (
        df["structure_id_path"]
            .str.split(pat="/")
            .map(lambda x: [int(i) for i in x[1:-1]])
    )
    #df["structure_id_path"] = df["structure_id_path"].map(lambda x: x[:-1])
    structures = df.to_dict("records")
    structures[0000]["structure_id_path"] = root_id
    for structure in structures:
        structure.update({"rgb_triplet": [255, 255, 255]})
        # root doesn't have a parent
        if structure["id"] != root_id:
            structure["structure_id_path"].append(structure["id"])


    # ---------------------------------------------------------------------------- #
    #                             MESHES                                           #
    # ---------------------------------------------------------------------------- #

    meshes_dict = None  # dictionary of files with region meshes

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
        compress=True,
    )

    return output_filename


# To test stuff locally:
if __name__ == "__main__":
    resolution = 10  # some resolution, in microns

    # Generated atlas path:
    bg_root_dir = Path.home() / "brainglobe_workingdir" / "blackcap"
    bg_root_dir.mkdir(exist_ok=True, parents=True)

    create_atlas(bg_root_dir, resolution)
