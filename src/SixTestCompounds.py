"""SixTestCompounds analysis main script.

A script which can be called from the command line to perform
the tracer kinetic modelling and generate all resulting reports
and figures for a specific study of interest.
"""
from pathlib import Path
import os
import shutil
import argparse
import itertools
import data
import signals
import analyses
import plots
from rat import TristanRat


def main(study: str
         ) -> None:
    """Import & prepare data from TK modelling & output statistical summaries.

    Args:
        study: Study name of interest (e.g., 'SixTestCompounds')
    """
    # Get original files and filenames
    orig_files, orig_filenames = data.get_files(study, '01_signals')
    # Modify filenames to include study number
    # to comply with updated functionality used in Reproducibility study
    renamed_folder = os.path.join(Path(orig_files[0]).parent.parent,
                                  '01_signals_renamed')
    data.make_dir(renamed_folder)
    for n in orig_files:
        if Path(n).stem.split('_')[0]=='Asunaprevir':
            shutil.copyfile(n, os.path.join(renamed_folder, '5_' + Path(n).stem + '.csv'))
        elif Path(n).stem.split('_')[0]=='Pioglitazone':
            shutil.copyfile(n, os.path.join(renamed_folder, '6_' + Path(n).stem + '.csv'))
        elif Path(n).stem.split('_')[0]=='Ketoconazole':
            shutil.copyfile(n, os.path.join(renamed_folder, '7_' + Path(n).stem + '.csv'))
        elif Path(n).stem.split('_')[0]=='Cyclosporine':
            shutil.copyfile(n, os.path.join(renamed_folder, '8_' + Path(n).stem + '.csv'))
        elif Path(n).stem.split('_')[0]=='Bosentan':
            shutil.copyfile(n, os.path.join(renamed_folder, '10_' + Path(n).stem + '.csv'))
        elif Path(n).stem.split('_')[0]=='BosentanHigh':
            shutil.copyfile(n, os.path.join(renamed_folder, '9_' + Path(n).stem + '.csv'))
        elif Path(n).stem.split('_')[0]=='Rifampicin':
            shutil.copyfile(n, os.path.join(renamed_folder, '12_' + Path(n).stem + '.csv'))
    # Get updated files and filenames
    files, filenames = data.get_files(study, '01_signals_renamed')
    # Split control and treatment groups
    signal_dict = signals.split_groups(files, filenames)
    # Fit data and get all estimated parameter variables
    all_parameters = signals.fit_data(study,
                                      filenames,
                                      files,
                                      signal_dict,
                                      TristanRat)

    # Get time curve averages per drug and per day
    subject_list = signals.get_subject_list(signal_dict)
    for curve in ['Delta R1 Liver (s-1)', 'Delta R1 Liver fit (s-1)',
                  'Delta R1 Spleen (s-1)']:
        signals.get_average_curves(signal_dict, subject_list, curve)
        
    fits = signal_dict
    # Plot average delta R1 time curves per drug and per day
    for substudy, day in list(itertools.product(fits.keys(), [1, 2])):
        try:
            print(f"{substudy}, Liver fit: Saving average deltaR1 plot")
            # For fitted liver data
            plots.get_deltaR1_plots(fits,
                                    substudy,
                                    'Liver',
                                    study,
                                    is_fitted=True,
                                    YLIM=(-1.5, 4.5))
            # For observed liver data only
            print(f"{substudy}, Liver: Saving average deltaR1 plot")
            plots.get_deltaR1_plots(fits,
                                    substudy,
                                    'Liver',
                                    study,
                                    is_fitted=False,
                                    YLIM=(-1.5, 4.5))
            # For observed spleen data only
            print(f"{substudy}, Spleen: Saving average deltaR1 plot")
            plots.get_deltaR1_plots(fits,
                                    substudy,
                                    'Spleen',
                                    study,
                                    is_fitted=False,
                                    YLIM=(-0.25, 1))
        except KeyError:
            continue

    # Convert substudy string number labels to integers
    all_parameters['Substudy'] = all_parameters['Substudy'].astype(int)

    # Create dictionary to rename sites into more comprehensive format
    site_names = {'Bosentan': 'Bosentan_2mg',
                  'BosentanHigh': 'Bosentan_high',
                  'Cyclosporine': 'Ciclosporin'}
    all_parameters.replace({'Drug': site_names}, inplace=True)

    # Create extra column storing combined site and drug labels
    all_parameters['Site_drug'] = (all_parameters['Site'] +
                                   ' ' + all_parameters['Drug'])

    # Quality control
    # Remove missing data
    # and computational fitting errors from all estimated parameter data
    print("Removing computational fitting errors and missing data")
    all_parameters_cleaned = data.remove_data_errors(all_parameters,
                                                     study)
    # Remove subjects with insufficient number of observations
    all_parameters_cleaned = (data
                              .remove_insufficient_data(all_parameters_cleaned,
                                                        study))

    # Create list of condition variables to group by
    variables = ['Drug', 'Symbol', 'Site']
    # Create list of biomarkers (parameters) of interest
    params = ['Ktrans', 'kbh', 'khe']

    # Obtain effect size summaries and save as csv
    print("Calculating average effect sizes")
    # Get statistical summary for saline-saline data
    single_subject, overall = (analyses.get_retest_results(study,
                                                           'effects',
                                                           all_parameters_cleaned.query("Symbol in @params")))

    # Plot biomarker distributions between Day 1 and Day 2 per rat
    print("Plotting individual biomarker \
              distributions between Day 1 and Day 2")
    plots.pairplots(study,
                    'effects',
                    all_parameters_cleaned.query("Symbol in @params"),
                    'Day',
                    'Rat',
                    'Drug',
                    'Symbol',
                    'rocket',
                    95,
                    ylabels=['$K_{trans}$', '$k_{he}$', '$k_{bh}$'],
                    sharey='row')
    print("Done!")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Import & prepare data \
                                     from TK modelling and output \
                                     statistical summaries',
                                     usage='python main.py --study study_name')
    parser.add_argument('--study', required=True, help='Study name of \
                                                        interest, e.g., \
                                                        SixTestCompounds')
    args = parser.parse_args()
    main(args.study)
