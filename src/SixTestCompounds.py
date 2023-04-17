"""SixTestCompounds analysis main script.

A script which can be called from the command line to perform
the tracer kinetic modelling and generate all resulting reports
and figures for a specific study of interest.
"""
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
    # Get files and filenames
    files, filenames = data.get_files(study, '01_signals')
    # Split control and treatment groups
    signal_dict = signals.split_groups(files, filenames)
    # Fit data and get all estimated parameter variables
    all_parameters = signals.fit_data(study, filenames, files,
                                     signal_dict, TristanRat)

    # Get time curve averages per drug and per day
    subject_list = signals.get_subject_list(signal_dict)
    for curve in ['Delta R1 Liver (s-1)', 'Delta R1 Liver fit (s-1)',
                  'Delta R1 Spleen (s-1)']:
        signals.get_average_curves(signal_dict, subject_list, curve)

    # Update dictionary keys for average delta R1 plots
    fits = signal_dict
    fits['G2 Ciclosporin'] = fits.pop('Cyclosporine')
    fits['G2 Rifampicin'] = fits.pop('Rifampicin')
    fits['D Ketoconazole'] = fits.pop('Ketoconazole')
    fits['E Asunaprevir'] = fits.pop('Asunaprevir')
    fits['E Pioglitazone'] = fits.pop('Pioglitazone')
    fits['G1 Bosentan_2mg'] = fits.pop('Bosentan')
    fits['G1 Bosentan_high'] = fits.pop('BosentanHigh')

    # Plot average delta R1 time curves per drug and per day
    for drug, day in list(itertools.product(fits.keys(), [1, 2])):
        print(f"{drug}, Liver fit: Saving average deltaR1 plot")
        # For fitted liver data
        plots.get_deltaR1_plots(fits, drug, 'Liver', study,
                                is_fitted=True, YLIM=(-1.5, 4.5))
        # For observed liver data only
        print(f"{drug}, Liver: Saving average deltaR1 plot")
        plots.get_deltaR1_plots(fits, drug, 'Liver', study,
                                is_fitted=False, YLIM=(-1.5, 4.5))
        # For observed spleen data only
        print(f"{drug}, Spleen: Saving average deltaR1 plot")
        plots.get_deltaR1_plots(fits, drug, 'Spleen', study,
                                is_fitted=False, YLIM=(-0.25, 1))

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
    all_parameters_cleaned = analyses.remove_data_errors(all_parameters,
                                                             study)

    # Create list of condition variables to group by
    variables = ['Drug', 'Symbol', 'Site']
    # Create list of biomarkers (parameters) of interest
    params = ['Ktrans', 'kbh', 'khe']

    # Obtain effect size summaries and save as csv
    print("Calculating average effect sizes")
    analyses.save_effect_sizes(all_parameters_cleaned, params,
                                   variables, study)

    # Plot biomarker distributions between Day 1 and Day 2 per rat
    for biomarker in params:
        print(f"{biomarker}: Plotting individual biomarker distributions between Day 1 and Day 2")
        plots.pairplots(all_parameters, str(biomarker), study)

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
