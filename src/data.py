"""Data handling

A module which locates files and filenames, 
and returns filename descriptions of all project input data
along with pathnames where results should be stored

"""
# imports
import os
import sys
import pandas as pd
import glob
from pathlib import Path
from typing import Dict, TypedDict, List, Tuple
from datetime import datetime


# Helper functions
def get_date_time(
) -> str:
    """Gets current date.
    
    Module obtains current date and returns it as a string variable. 
    This variable can then be used later when outputting results into 
    new folders labelled by date."""
    
    now = datetime.now() # gets current date

    # converts to YYYY-MM-DD format, e.g., 2022-10-31
    dt_string = now.strftime(f"%Y-%m-%d")
    
    return dt_string


def make_dir(path: Path
) -> None:
    """Creates directory if it does not exist."""
    if not os.path.exists(path):
        os.mkdir(path)


# Folder paths
def get_path(study: str,
            folder: str,
) -> tuple[str, str]:
    """Gets data/results folderpath of chosen study.
    
     Args:
        study: Study name of interest (e.g., 'SixTestCompounds').
        folder: Folder of interest (e.g., 'results').

    Returns:
        Folderpath for data or results folder (as selected by user)
        related to study of choice in string format.
    """
    cwd = os.getcwd()
    path = f"{cwd}\\..\\{folder}\\{study}\\"
    
    return path


# Data filenames & filepaths
def get_files(study: str,
                dataset: str
) -> tuple[list, list]:
    """Gets list of filepaths and filenames.
    
    Searches dataset folder within chosen study and returns
    list of filepaths and filenames contained within.

    Args:
        study: Study name of interest (e.g., 'SixTestCompounds').
        dataset: Dataset of interest (e.g., '01_signals').
    """
    study_data_path = get_path(study, 'data')
    dataset_path = f"{study_data_path}{dataset}"
    files = glob.glob(dataset_path + '*\\*')
    filenames = [Path(path).stem for path in files]

    return files, filenames


# File metadata
def get_metadata(filename: str,
                files: Path
) -> TypedDict('metadata', {'drug': str, 
                            'site': str,
                            'subject': int,
                            'day': int,
                            'signals': pd.DataFrame
                            }
                ):
    """Extracts study descriptions from filenames.
    
    This function allows the user to extract metadata from
    a filename along with the corresponding dataframe contained
    within the file. Works only when filename is formatted as 
    a string containing study descriptors (metadata) separated 
    by underscores, i.e.,
    filename = "compound_site_RatNumber_dayNumber_dataType"
    e.g., "Asunaprevir_E_Rat2_2_Signals"

    Args:
        filename: File name of interest.
        files: File of interest.

    Returns:
        A dict mapping keys to the corresponding metadata and
        signal data contained in csv file.
    """
    filename_elements = filename.split('_')
    metadata = {}
    metadata['drug'] = filename_elements[0]
    metadata['site'] = filename_elements[1]
    metadata['subject'] = int(filename_elements[2][3:])
    metadata['day'] = int(filename_elements[3])
    metadata['signals'] = pd.read_csv(files)

    return metadata


# Subdirectory paths
def get_subdir(parentdir: str,
                subdir_name: str
) -> str:
    """Gets subdirectory path.

    Args:
        parentdir: Parent directory where subdirectory
            should be stored.
        subdir_name: Name of subdirectory.

    Returns:
        Subdirectory path name in string format.    
    """
    subdir_path = f"{parentdir}\\{subdir_name}\\"
    make_dir(subdir_path)
    
    return subdir_path


# Results folders
def get_results_folder(study: str,
                       folder_name: str,
                       subfolder_name: str,
                       extra_subfolder_name: str,
                       filename: str,
                       file_type: str
) -> str:
    
    dt_string = get_date_time()    
    parent_dir = get_path(study, 'results')
    date_folder = get_subdir(parent_dir, dt_string)
    output_folder = get_subdir(date_folder, folder_name)
    
    if (subfolder_name==None) & (extra_subfolder_name==None):
        filepath = f"{output_folder}\\{filename}"
        
    elif (subfolder_name!=None) & (extra_subfolder_name==None):
        subfolder = get_subdir(output_folder, subfolder_name)
        filepath = f"{subfolder}\\{filename}"
        
    elif (subfolder_name!=None) & (extra_subfolder_name!=None):
        subfolder = get_subdir(output_folder, subfolder_name)
        extra_subfolder = get_subdir(subfolder, extra_subfolder_name)
        filepath = f"{extra_subfolder}\\{filename}"
        
    save_name = f"{filepath}.{file_type}"
    
    return save_name
