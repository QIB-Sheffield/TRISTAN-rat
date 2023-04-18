"""Data handling.

A module which organises the datasets by locating files and filenames,
returning filename descriptions of all project input data, along with
pathnames where results should be stored. Additionally cleans data, e.g.,
removing errors due to model non-convergence and datasets with insufficient
sample numbers.
"""
# imports
import os
import sys
import pandas as pd
import glob
from pathlib import Path
from datetime import datetime
if sys.version_info >= (3, 8):
    from typing import TypedDict, Tuple
else:
    from typing import Tuple
    from typing_extensions import TypedDict


# Helper functions
def get_date_time(
) -> str:
    """Gets current date.

    Module obtains current date and returns it as a string variable.
    This variable can then be used for outputting results into new
    folders labelled by date.
    """
    now = datetime.now()  # gets current date

    # converts to YYYY-MM-DD format, e.g., 2022-10-31
    dt_string = now.strftime("%Y-%m-%d")

    return dt_string


def make_dir(path: Path
             ) -> None:
    """Creates directory if it does not exist."""
    if not os.path.exists(path):
        os.mkdir(path)


# Folder paths
def get_path(study: str,
             folder: str
             ) -> Tuple[str, str]:
    """Gets data/results folder path of chosen study.

     Args:
        study: Study name of interest (e.g., 'SixTestCompounds').
        folder: Folder of interest (e.g., 'results').

    Returns:
        Folder path (in string format) for data or results folder
        (as selected by user) related to study of choice.
    """
    cwd = os.getcwd()
    path = f"{cwd}\\{folder}\\{study}\\"

    return path


# Data filenames & filepaths
def get_files(study: str,
              dataset: str
              ) -> Tuple[list, list]:
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
                 ) -> TypedDict('metadata', {'substudy': str,
                                             'drug': str,
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
    metadata['substudy'] = filename_elements[0]
    metadata['drug'] = filename_elements[1]
    metadata['site'] = filename_elements[2]
    metadata['subject'] = int(filename_elements[3][3:])
    metadata['day'] = int(filename_elements[4])
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
    subdir_path = f"{parentdir}{subdir_name}\\"
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
    """Obtains results folder and subfolders.

    Args:
        study: Study name of interest (e.g., 'SixTestCompounds').
        folder_name: Results folder of interest (e.g., '02_effect_sizes').
        subfolder_name: Results subfolder of interest (e.g., 'figures'). If
            no subfolders are required use None.
        extra_subfolder_name: Extra results subfolder of interest (e.g.,
            'per_Rat'). If no extra subfolders are required use None.
        filename: Chosen file name.
        file_type: File format, e.g., 'csv'.

    Returns:
        Directory path name (in string format) for storing results.
    """
    make_dir('results')
    parent_dir = get_path(study, 'results')
    make_dir(parent_dir)
    dt_string = get_date_time()
    date_folder = get_subdir(parent_dir, dt_string)
    output_folder = get_subdir(date_folder, folder_name)

    if (subfolder_name is None) & (extra_subfolder_name is None):
        filepath = f"{output_folder}\\{filename}"

    elif (subfolder_name is not None) & (extra_subfolder_name is None):
        subfolder = get_subdir(output_folder, subfolder_name)
        filepath = f"{subfolder}\\{filename}"

    elif (subfolder_name is not None) & (extra_subfolder_name is not None):
        subfolder = get_subdir(output_folder, subfolder_name)
        extra_subfolder = get_subdir(subfolder, extra_subfolder_name)
        filepath = f"{extra_subfolder}\\{filename}"

    save_name = f"{filepath}.{file_type}"

    return save_name


# Data cleaning
def remove_data_errors(parameter_data: pd.DataFrame,
                       study: str
                       ) -> pd.DataFrame:
    """Cleans data.

    Removes computational fitting errors in biomarker prediction caused by
    non-convergence of the model fitting algorithm, i.e., the output of
    boundary values instead of true values, e.g., gadoxetate extraction
    fraction, E=100%.

    Args:
        parameter_data: DataFrame containing estimated parameter variables.
        study: Study name of interest (e.g., 'SixTestCompounds').

    Returns:
        Cleaned DataFrame.
    """
    data_pivoted = pd.pivot_table(parameter_data,
                                  values='Value',
                                  columns=['Symbol'],
                                  index=['Substudy',
                                         'Site',
                                         'Drug',
                                         'Rat',
                                         'Day'])

    # Remove computational fitting errors based on subjects where gadoxetate
    # extraction fraction, E is close or equal to 100% (i.e., >= 99%)
    fit_errors = data_pivoted[data_pivoted['E'] >= 99.95]
    fit_errors_removed = (data_pivoted[~data_pivoted
                                       .index.isin(fit_errors.index)])

    # Save index metadata for computational fitting errors
    save_name = get_results_folder(study,
                                   '01_model_outputs',
                                   None,
                                   None,
                                   'fit_errors',
                                   'txt')
    with open(save_name, "w") as output:
        output.write(str(list([fit_errors.index])))

    cleaned_parameter_data = fit_errors_removed.stack().reset_index()
    cleaned_parameter_data.rename(columns={0: 'Value'}, inplace=True)

    return cleaned_parameter_data


def remove_insufficient_data(parameter_data: pd.DataFrame,
                             study: str
                             ) -> pd.DataFrame:
    """Removes data with insufficient number of observations.

    Removes cases where only one acquisition per subject is present
    (i.e., day 1 or day 2 data are missing), or whole substudies
    where data for only 1 subject is present.

    Args:
        parameter_data: DataFrame containing estimated parameter variables.
        study: Study name of interest (e.g., 'Reproducibility').

    Returns:
        Cleaned DataFrame.
    """
    data_pivoted = pd.pivot_table(parameter_data,
                                  values='Value',
                                  columns=['Symbol'],
                                  index=['Substudy',
                                         'Drug',
                                         'Site',
                                         'Fstrength',
                                         'Site_Fstrength',
                                         'Time_period',
                                         'Rat',
                                         'Day'])
    # Remove subjects with missing acquisition on day 1 or day 2
    missing_days_removed = (data_pivoted[data_pivoted
                                         .groupby(['Substudy',
                                                   'Site',
                                                   'Rat'])
                                         .transform('count') > 1]
                            .dropna())
    missing_days_removed_clean = missing_days_removed.stack().reset_index()
    missing_days_removed_clean.rename(columns={0: 'Value'}, inplace=True)

    # Count number of substudies with only 1 subject
    counts = missing_days_removed_clean.groupby(['Symbol',
                                                 'Substudy',
                                                 'Site',
                                                 'Drug',
                                                 'Day'])['Rat'].count()

    # Remove substudies with only 1 subject
    substudy_to_drop = []
    for i in counts.index:
        if counts[i] == 1:
            substudy_to_drop.append(i[1])
    insufficient_data_removed = (missing_days_removed_clean[~missing_days_removed_clean['Substudy']
                                                            .isin(substudy_to_drop)])
    insufficient_data_removed.reset_index(drop=True, inplace=True)

    return insufficient_data_removed
