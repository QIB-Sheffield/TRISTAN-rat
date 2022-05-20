import os
import pandas as pd

def results_path():

    return "C:\\Users\\steve\\Dropbox\\Projects\\IMI-TRISTAN\\WP2\\pytristan\\Results\\"

def description(filename):

    x = filename.split('_')
    y = {}
    y['drug'] = x[0]
    y['site'] = x[1]
    y['subject'] = int(x[2][3:])
    y['day'] = int(x[3])
    return y

def example():

    filename = "Asunaprevir_MSD_Rat2_2_Signals"
    filename = "Asunaprevir_MSD_Rat5_1_Signals"
    file = "C:\\Users\\steve\\Dropbox\\Projects\\IMI-TRISTAN\\WP2\\pytristan\\SixTestCompounds\\Asunaprevir\\" + filename +".csv"
    df = pd.read_csv(file)
    return df, filename

def folder(study):

    path = "C:\\Users\\steve\\Dropbox\\Projects\\IMI-TRISTAN\\WP2\\pytristan\\Data" + study

    files = [item.path for item in _scan_tree(path) if item.is_file()]
    filenames = [os.path.basename(path) for path in files]
    filenames = [os.path.splitext(path)[0] for path in filenames]
    return files, filenames

def _scan_tree(directory):
    """Helper function: yield DirEntry objects for the directory."""

    for entry in os.scandir(directory):
        if entry.is_dir(follow_symlinks=False):
            yield from _scan_tree(entry.path)
        else:
            yield entry


