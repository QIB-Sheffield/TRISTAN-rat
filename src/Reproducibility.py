"""Reproducibility analysis main script.

A script which can be called from the command line to perform
the tracer kinetic modelling and generate all resulting reports
and figures for a specific study of interest.
"""
import argparse
import itertools
import numpy as np
import pandas as pd
import data
import signals
import analyses
import plots
from rat import TristanRat


def main(study: str
         ) -> None:
    """Import & prepare data from TK modelling & output statistical summaries.

    Args:
        study: Study name of interest (e.g., 'Reproducibility')
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

    fits = signal_dict

    for substudy, day in list(itertools.product(fits.keys(), [1, 2])):
        try:
            # For fitted liver data
            print(f"{substudy}, Liver fit: Saving average deltaR1 plot")
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

    # CLEAN DATA ######
    # Convert substudy string number labels to integers
    all_parameters['Substudy'] = all_parameters['Substudy'].astype(int)

    # Quality control
    # Remove computational fitting errors from all estimated parameter data
    print("Removing computational fitting errors")
    all_parameters_cleaned = data.remove_data_errors(all_parameters,
                                                     study)

    # DATA PREP/CURATION ####
    # create dictionary storing field strengths per site
    field_strengths = {'E': '7', 'D': '4.7', 'G1': '7', 'G2': '4.7'}
    # create column for field strengths by duplicating site labels
    all_parameters_cleaned['Fstrength'] = all_parameters_cleaned.loc[:, 'Site']
    # replace site labels with field strengths
    all_parameters_cleaned.replace({'Fstrength': field_strengths},
                                   inplace=True)

    # create column for field strength & site labels combined
    all_parameters_cleaned['Site_Fstrength'] = (all_parameters_cleaned['Site']
                                                .astype(str).str[0] + ' (' + all_parameters_cleaned['Fstrength']
                                                .astype(str) + 'T)')

    # Combine G1 & G2 into same site G
    # create dictionary assigning G to G1 & G2 labels
    G = {'G1': 'G', 'G2': 'G'}
    all_parameters_cleaned.replace({'Site': G}, inplace=True)

    # create dictionary storing time periods per substudy
    time_periods = {1: 0, 2: 0, 3: 0, 4: 0,
                    5: 12, 6: 13, 7: 19, 8: 2,
                    9: 15, 10: 10, 11: 16, 12: 0,
                    13: 0}
    # create column for time periods by duplicating substudy labels
    all_parameters_cleaned['Time_period'] = (all_parameters_cleaned
                                             .loc[:, 'Substudy'])
    # replace substudy labels with time periods
    all_parameters_cleaned.replace({'Time_period': time_periods}, inplace=True)
    all_parameters_cleaned['Time_period'] = np.where((all_parameters_cleaned['Substudy'] == 13) & (all_parameters_cleaned['Fstrength'] == '7'),
                                                     14, all_parameters_cleaned['Time_period'])

    # Sort dataframe
    (all_parameters_cleaned
     .sort_values(['Substudy', 'Site', 'Fstrength'],
                  ascending=[True, True, True],
                  inplace=True))

    # GET STATISTICS
    # Create list of biomarkers (parameters) of interest
    params = ['Ktrans', 'kbh']
    # Split data into variance of interest groups
    # across substudies
    substudy = all_parameters_cleaned.query('Symbol in @params and Day==1')
    # across sites (at specific field strengths)
    site = all_parameters_cleaned.query('Substudy<=4 and Symbol in @params \
                                        and Day==1')
    # across field strengths (at specific site G)
    fieldstrength = all_parameters_cleaned.query('Substudy<=4 and Site=="G" \
                                                 and Symbol in @params and \
                                                 Day==1')
    # across field strengths (at specific centre and subject)
    fieldstrength_pairwise = all_parameters_cleaned.query('Substudy==13 and \
                                                          Symbol in @params \
                                                          and Rat!=[4, 5, 6]')
    # across time periods (at specific sites and field strengths)
    timeperiod = all_parameters_cleaned.query('Symbol in @params and Day==1')

    # Get statistical summaries per group
    group_summaries = pd.concat([analyses
                                 .summarise_group(substudy,
                                                  'substudy',
                                                  ['Substudy'],
                                                  None),
                                 analyses
                                 .summarise_group(site,
                                                  'site',
                                                  ['Site'],
                                                  ['Fstrength']),
                                 analyses
                                 .summarise_group(fieldstrength,
                                                  'fieldstrength',
                                                  ['Fstrength'],
                                                  ['Site']),
                                 analyses
                                 .summarise_group(fieldstrength_pairwise,
                                                  'fieldstrength_pairwise',
                                                  ['Fstrength'],
                                                  ['Site']),
                                 analyses
                                 .summarise_group(timeperiod,
                                                  'timeperiod',
                                                  ['Time_period'],
                                                  ['Site_Fstrength'])])
    print("Performing reproducibility analysis")
    # Get deconstructed reproducibility uncertainties
    reproducibility_deconstructed = pd.concat([analyses
                                               .deconstruct_reproducibility(group_summaries,
                                                                            'substudy'),
                                               analyses
                                               .deconstruct_reproducibility(group_summaries,
                                                                            'site'),
                                               analyses
                                               .deconstruct_reproducibility(group_summaries,
                                                                            'fieldstrength'),
                                               analyses
                                               .deconstruct_reproducibility(group_summaries,
                                                                            'fieldstrength_pairwise'),
                                               analyses
                                               .deconstruct_reproducibility(group_summaries,
                                                                            'timeperiod')])

    # Benchmarks (mean of substudy means +/- 95%CI)
    benchmarks = (reproducibility_deconstructed
                  .query("Variables=='substudy'")[['reproducibility',
                                                   'mean',
                                                   'CI95']])
    # Calculate absolute detection limits
    absolute = []
    for biomarker in params:
        (absolute
         .append(analyses
                 .detect_absolute(benchmarks['reproducibility'][biomarker],
                                  benchmarks['mean'][biomarker],
                                  benchmarks['CI95'][biomarker])[0]))
    benchmarks['absolute'] = absolute

    print("Plotting biomarker distributions per variance group")
    # Plot group distributions
    plots.plot_distributions(study,
                             'substudy',
                             substudy,
                             'Substudy',
                             None,
                             benchmarks,
                             None)
    plots.plot_distributions(study,
                             'site',
                             site,
                             'Site',
                             'Fstrength',
                             benchmarks,
                             None)
    plots.plot_distributions(study,
                             'fieldstrength',
                             fieldstrength,
                             'Fstrength',
                             'Site',
                             benchmarks,
                             None)
    plots.plot_distributions(study,
                             'fieldstrength_pairwise',
                             fieldstrength_pairwise,
                             'Fstrength',
                             'Site',
                             benchmarks,
                             None)
    # Plot site-fieldstrength combinations according to ascending time period
    sites_list = ['G (4.7T)', 'E (7T)', 'G (7T)', 'D (4.7T)']
    plots.plot_distributions(study,
                             'timeperiod',
                             timeperiod,
                             'Time_period',
                             'Site_Fstrength',
                             benchmarks,
                             sites_list)

    # Get list of unique constants for each variable and
    # store in dictionary
    constants = {}
    for variable in group_summaries['Variables'].unique():
        constants[variable] = analyses.get_constants(group_summaries, variable)

    # Calculate one-way ANOVA p-values for each variance
    # group in deconstructed reproducibility table
    one_way_anovas = []
    # substudy
    for biomarker in params:
        one_way_anovas.append(analyses
                              .get_anova_oneway(substudy,
                                                biomarker,
                                                None,
                                                None,
                                                'Substudy'))
    # site
    for biomarker, constant in list(itertools
                                    .product(params,
                                             constants['site'])):
        one_way_anovas.append(analyses
                              .get_anova_oneway(site,
                                                biomarker,
                                                'Fstrength',
                                                constant,
                                                'Site'))
    # fieldstrength
    for biomarker in params:
        one_way_anovas.append(analyses
                              .get_anova_oneway(fieldstrength,
                                                biomarker,
                                                'Site',
                                                'G',
                                                'Fstrength'))
    # fieldstrength pairwise
    for biomarker in params:
        one_way_anovas.append(analyses
                              .get_anova_oneway(fieldstrength_pairwise,
                                                biomarker,
                                                'Site',
                                                'G',
                                                'Fstrength'))
    # timeperiod
    for biomarker, constant in list(itertools
                                    .product(params,
                                             constants['timeperiod'])):
        one_way_anovas.append(analyses.get_anova_oneway(timeperiod,
                                                        biomarker,
                                                        'Site_Fstrength',
                                                        constant,
                                                        'Time_period'))
    # Append p-values to deconstructed reproducibility table
    reproducibility_deconstructed['pvalue'] = one_way_anovas
    # Save table of deconstructed reproducbility uncertainties
    save_name = data.get_results_folder(study,
                                        '02_analyses',
                                        'reproducibility',
                                        None,
                                        'reproducibility_deconstructed',
                                        'csv')
    reproducibility_deconstructed.to_csv(save_name)

    print("Performing repeatability and rifampicin-effect analysis")
    # Extract retest data (i.e., substudies 1-4)
    retest_data = (data
                   .remove_insufficient_data(all_parameters_cleaned
                                             .query("Symbol in @params \
                                                    and Substudy<=4"),
                                             study))

    # Separate saline-rifampicin subjects
    rifampicin_data = (retest_data
                       .query("Drug=='Rifampicin'")[['Substudy',
                                                     'Site',
                                                     'Symbol',
                                                     'Rat']]
                                                     .merge(retest_data))

    # and saline-saline subjects
    combined = retest_data.merge(rifampicin_data,
                                 on=['Substudy',
                                     'Site',
                                     'Symbol',
                                     'Rat'],
                                 how='outer',
                                 indicator=True)
    saline_data = combined[combined._merge != 'both'].dropna(axis=1)
    saline_data.columns = saline_data.columns.str.split('_x').str[0]
    saline_data.drop(saline_data.columns[-1], axis=1, inplace=True)
    saline_data.sort_values('Substudy', inplace=True)
    saline_data.reset_index(drop=True, inplace=True)

    # Get statistical summary for saline-saline data
    single_subject_saline, overall_saline = (analyses
                                             .get_retest_results(study,
                                                                 'saline-saline',
                                                                 saline_data))
    # Get statistical summary for saline-rifampicin data
    single_subject_rifampicin, overall_rifampicin = (analyses
                                                     .get_retest_results(study,
                                                                         'saline-rifampicin',
                                                                         rifampicin_data))

    # Calculate relative detection limits
    benchmarks['repeatability'] = (overall_saline['repeatability']['mean']
                                   .values)
    benchmarks.eval("relative = repeatability * sqrt(2)", inplace=True)

    # Save benchmarks, uncertainties and absolute/relative detection limits
    save_name = data.get_results_folder(study,
                                        '02_analyses',
                                        None,
                                        None,
                                        'benchmarks',
                                        'csv')
    benchmarks.to_csv(save_name)

    # plot saline-saline data
    print("saline-saline: Plotting individual biomarker distributions between Day 1 and Day 2")
    plots.pairplots(study,
                    'saline-saline',
                    saline_data,
                    'Day',
                    'Substudy',
                    'Symbol',
                    None,
                    'rocket',
                    95,
                    ylabels=['$K_{trans}$', '$k_{bh}$'])
    # plot saline-rifampicin data
    print("saline-rifampicin: Plotting individual biomarker distributions between Day 1 and Day 2")
    plots.pairplots(study,
                    'saline-rifampicin',
                    rifampicin_data,
                    'Day',
                    'Substudy',
                    'Symbol',
                    None,
                    'rocket',
                    95,
                    ylabels=['$K_{trans}$', '$k_{bh}$'])

    # MIXED ANOVA (saline retest data)
    print("Performing mixed ANOVA for saline-saline retest data")
    saline_variances = {}
    for biomarker in params:
        saline_variances[biomarker] = (analyses
                                       .get_mixed_anova(saline_data,
                                                        biomarker))

    variance_components = pd.DataFrame(data=saline_variances['Ktrans'][1],
                                       index=saline_variances['Ktrans'][0],
                                       columns=['Ktrans'])
    variance_components['kbh'] = saline_variances['kbh'][1]

    # Save mixed ANOVA results
    save_name = data.get_results_folder(study,
                                        '02_analyses',
                                        None,
                                        None,
                                        'mixed_anova',
                                        'csv')
    variance_components.to_csv(save_name)

    print("Done!")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Import & prepare data \
                                     from TK modelling and output \
                                     statistical summaries',
                                     usage='python main.py --study study_name')
    parser.add_argument('--study', required=True, help='Study name of \
                                                        interest, e.g., \
                                                        Reproducibility')
    args = parser.parse_args()
    main(args.study)
