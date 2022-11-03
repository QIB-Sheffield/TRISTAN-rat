"""Main module

A module which 

"""
import os
import sys
import configparser
import argparse
import ast
from pathlib import Path
import itertools
import data
import models
import effect_sizes
import plots
from rat import TristanRat


def main(study: str
) -> None:
    """Import & prepare data from TK modelling and output statistical summaries

    Args:
        study: study of interest (e.g., 'SixTestCompounds')
    """

    files, filenames = data.get_files(study, '01_signals')
    signals = models.split_groups(files, filenames)
    all_parameters = models.fit_data(study, filenames, files, signals, TristanRat)

    subject_range = models.get_num_subjects(signals)
    for curve in ['Delta R1 Liver (s-1)', 'Delta R1 Liver fit (s-1)', 'Delta R1 Spleen (s-1)']:
        models.get_signal_averages(signals, subject_range, curve)

    for drug, day in list(itertools.product(signals.keys(), [1, 2])):
        plots.get_deltaR1_plots(signals, drug, 'Liver', study, is_fitted=True, YLIM=(0, 5))
        plots.get_deltaR1_plots(signals, drug, 'Liver', study, is_fitted=False, YLIM=(0, 5))
        plots.get_deltaR1_plots(signals, drug, 'Spleen', study, is_fitted=False, YLIM=(0, 1))
        
    # Pivot dataframe
    site_names = {'Bosentan':'Bosentan_2mg', 'BosentanHigh':'Bosentan_high', 
    'Cyclosporine':'Ciclosporin'} # create dictionary to differentiate to update site names

    all_parameters.replace({'Drug':site_names}, inplace=True)

    all_parameters['Site_drug'] = all_parameters['Site'] + ' ' + all_parameters['Drug']

    all_parameters_cleaned = effect_sizes.remove_data_errors(all_parameters, study)

    variables = ['Drug','Symbol','Site']
    params = ['Ktrans', 'kbh', 'khe']

    # save effect_sizes
    effect_sizes.save_effect_sizes(all_parameters_cleaned, params, variables, study)
        
    # Pairplots    
    for i in params:
        plots.pairplots(all_parameters, str(i), study)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Import & prepare data from TK modelling and output statistical summaries',
        usage='python main.py --study study_name')
    parser.add_argument('--study', required=True, help='Study name')
    args = parser.parse_args()
    main(args.study)