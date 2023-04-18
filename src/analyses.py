"""Data analysis.

A module which provides statistical summaries of the data
calculating effect sizes, repeatability, and reproducibility
errors, and performing ANOVA and Student's t-test.
"""
# imports
import pandas as pd
import pingouin as pg
from scipy import stats
import data
from typing import Tuple


def get_relative_error(variability: int,
                       mean: int
                       ) -> int:
    """Calculates relative error.

    A function to calculate the relative error across a set of values,
    at the level of 95% confidence (i.e., 1.96 x standard deviation).

    Args:
        variability: Integer value of variability of interest of
            a dataset (e.g., standard deviation or standard error).
        mean: Integer mean value a dataset.

    Returns:
        Relative error as integer value.
    """

    relative_error = 1.96 * 100 * (variability / mean)

    return relative_error


def summarise_group(group_data: pd.DataFrame,
                    variable_name: str,
                    variables: list,
                    constants: list
                    ) -> pd.DataFrame:
    """Summarises data group statistics.

    Provides an overview of data for a variable group of
    interest. for a variable of interest, e.g., across time
    periods (variables) at a specific centre and field strength
    (constants).

    Args:
        group_data: DataFrame containing parameter values for a
            variable group of interest.
        variable_name: Name of variable of interest.
        variables: List of variables of interest
            (e.g., ['Time_period']).
        constants: List of factors to keep constant
            (e.g., ['Site', 'Fstrength']).

    Returns:
        Dataframe containing group means, confidence intervals, and
        between-subject reproducibility and relative errors.
    """
    if constants is None:
        group_means = (group_data
                       .groupby(['Symbol'] + variables)['Value']
                       .agg(['mean', 'std', 'sem']))
        group_means['Constants'] = (group_means
                                    .index.
                                    get_level_values(variables[0]))
    else:
        group_means = (group_data
                       .groupby(['Symbol'] + variables + constants)['Value']
                       .agg(['mean', 'std', 'sem']))
        group_means['Constants'] = (group_means
                                    .index.get_level_values(constants[0]))

    group_means['between_subject_reproducibility'] = get_relative_error(group_means['std'],
                                                                        group_means['mean'])
    group_means['between_subject_relative_error'] = get_relative_error(group_means['sem'],
                                                                       group_means['mean'])

    group_means['CI95'] = group_means['sem'].mul(1.96)
    group_means['Variables'] = variable_name

    return group_means


def deconstruct_reproducibility(group_summaries: pd.DataFrame,
                                variable: str
                                ) -> pd.DataFrame:
    """Calculates deconstructed reproducibility per variable group.

    Provides a summary of reproducibility errors for a variable group of
    interest, e.g., across time periods (variables) at a specific centre
    and field strength (constants).

    Args:
        group_summaries: DataFrame containing statistical summaries for a
            variable group of interest.
        variable: Name of variable of interest (e.g., ['Time_period']).

    Returns:
        Dataframe containing group means, confidence intervals, average
        between-subject reproducibility and relative errors, and
        between-substudy reproducibility errors.
    """

    if variable == 'substudy':
        grouper = ['Symbol', 'Variables']
    else:
        grouper = ['Symbol', 'Constants', 'Variables']

    reproducibility_deconstructed = (group_summaries
                                     .query("Variables==@variable")
                                     .groupby(grouper)['mean']
                                     .agg(['mean', 'std', 'sem']))
    reproducibility_deconstructed['reproducibility'] = get_relative_error(reproducibility_deconstructed['std'],
                                                                          reproducibility_deconstructed['mean'])
    reproducibility_deconstructed[['between_subject_reproducibility',
                                   'between_subject_sem']] = (group_summaries
                                                              .query("Variables==@variable")
                                                              .groupby(grouper)['between_subject_reproducibility']
                                                              .agg(['mean', 'sem']))

    reproducibility_deconstructed['CI95'] = (reproducibility_deconstructed['sem']
                                             .mul(1.96))

    return reproducibility_deconstructed


def detect_absolute(uncertainty: int,
                    benchmark_mean: int,
                    benchmark_CI95: int
                    ) -> int:
    """Calculates absolute detection limit.

    A function to calculate the absolute detection limit of any given
    biomarker, defined as the smallest change producing a 95% tolerance
    interval (TI) on the absolute value of the biomarker which does not
    overlap with the 95% CI of the saline benchmark value. The 95% TI in
    this case was estimated by applying the reproducibility error in saline.

    Args:
        uncertainty: Integer value of reproducibility error in saline.
        mean: Integer mean absolute saline benchmark value of biomarker.

    Returns:
        Absolute detection limit as integer value.
    """
    # uncertainty = uncertainty on x due to reproducibility error
    # lowerCI95 = mean_benchmark - 1.96*(SEM_benchmark)
    # x + uncertainty < lowerCI95
    # x = possible drug dose detected [mL/min/mL]

    lowerCI95 = benchmark_mean - benchmark_CI95

    # solve for x: x + uncertainty*x < lowerCI95
    x = lowerCI95/(1 + uncertainty/100)

    # absolute change detectable
    absolute = (benchmark_mean - x) / benchmark_mean

    return absolute


def get_constants(group_summaries: pd.DataFrame,
                  variable: str
                  ) -> list:
    """Gets list of constants per variance group.

    Creates a list of unique factors that are constant within a variable group
    of interest.

    Args:
        group_summaries: DataFrame containing statistical summaries for a
            variable group of interest.
        variable: Name of variable of interest (e.g., ['Time_period']).

    Returns:
        List of unqiue contant factors for variance group.
    """
    constants_list = (group_summaries
                      .query("Variables==@variable")['Constants']
                      .unique())

    return constants_list


def get_anova_oneway(group_data: pd.DataFrame,
                     biomarker: str,
                     constant_name: str,
                     constant: str,
                     variables: list):
    """Calculates one-way ANOVA p-values per variance group.

    Performs one-way ANOVA to help understand the relative impact
    of different variables on reproducibility, (i.e., according to
    stratification of Day 1 data by either centre, field strength,
    or time period).

    Args:
        group_data: DataFrame containing parameter values for a
            variable group of interest.
        biomarker: Name of biomarker of interest (e.g., 'Ktrans')
        constant_name: Column name containing values for constants
        constant: Name of factor to keep constant (e.g., 'Site').
        variables: List of variables across which to calculate ANOVA
            (e.g., 'Time_period').

    Returns:
        Statistical p-value result as integer value calculated from
        one-way ANOVA.
    """
    if constant is not None:
        aov_data = group_data[(group_data[constant_name] == constant) & (group_data['Symbol'] == biomarker)]
    else:
        aov_data = group_data[group_data['Symbol'] == biomarker]

    aov = pg.anova(data=aov_data, dv='Value', between=variables)
    pvalue = aov['p-unc'].values[0]

    return pvalue


def get_subject_repeatability(cleaned_parameter_data: pd.DataFrame
                              ) -> pd.DataFrame:
    """Calculates single-subject between-day repeatability.

    Calculates repeatability error of administered test compounds on the
    biomarkers of interest between Day 1 and Day 2 data for a single subject.

    Args:
        cleaned_paramter_data: Cleaned DataFrame containing estimated
            parameter variables.

    Returns:
        A dataframe containing single-subject between-day mean,
        standard deviation, percentage change, and repeatability
        error per biomarker per substudy.
    """
    subject_data = cleaned_parameter_data.copy()
    subject_data['pct_change'] = (subject_data
                                  .pct_change(axis=1).mul(100)[2])
    subject_data.reset_index(inplace=True)
    subject_data[['between_day_mean',
                  'between_day_std']] = (subject_data[[1, 2]]
                                         .agg(['mean', 'std'], axis=1))
    subject_data['repeatability'] = get_relative_error(subject_data['between_day_std'],
                                                       subject_data['between_day_mean'])

    return subject_data


def get_summary_stats(subject_data: pd.DataFrame,
                      cleaned_parameter_data: pd.DataFrame,
                      variables: list
                      ) -> pd.DataFrame:
    """Gets substudy-average statistical data summary.

    Calculates substudy-average summary statistics conveying the observed
    effects of administered test compounds on the biomarkers of interest,
    i.e., between day 1 and day 2, average subject values per substudy are
    calculated for gadoxetate percent change, effect size + 95% confidence
    intervals (CI95), and paired T-test p-values.

    Args:
        subject_data: DataFrame containing single-subject repeatability
            summary.
        cleaned_parameter_data: DataFrame containing estimated parameter
            variables.
        variables: List of condition variables to group by for statistical
            summary (e.g., ['drug', 'site']).

    Returns:
        A dataframe containing substudy simple, standardised, and
        percentage-change effect sizes, mean values, repeatability errors,
        p-value results from Student's t-test, and between-subject variation.
    """
    rats = subject_data.copy()
    global_std = cleaned_parameter_data.stack().groupby(variables).std()

    rats_avg = (rats.groupby(variables)[[1, 2, 'repeatability']]
                .agg(['mean', 'std', 'sem']))
    rats_avg['repeatability', 'CI95'] = (rats_avg['repeatability']['sem']
                                         .mul(1.96))

    rats_avg['simple_effect_size'] = (rats.groupby(variables)[1]
                                      .mean() -
                                      rats.groupby(variables)[2]
                                      .mean())
    rats_avg['percentage_effect_size'] = (rats_avg['simple_effect_size'] / (rats
                                                                            .groupby(variables)[1]
                                                                            .mean())
                                                                            .mul(100))

    rats_avg['standardised_effect_size'] = rats_avg['simple_effect_size'] / global_std

    rats.dropna(inplace=True)
    rats_avg['p-value'] = (rats.groupby(variables)
                           .apply(lambda df: stats
                                  .ttest_rel(df[1], df[2])[1]))
    for i in [1, 2]:
        rats_avg[i,
                 'between_subject_variation'] = get_relative_error(rats_avg[i]['std'],
                                                                   rats_avg[i]['mean'])
        rats_avg[i, 'CI95'] = rats_avg[i]['sem'].mul(1.96)

    rats_avg['between_subject_variation_average'] = (rats_avg
                                                     .groupby(level=1, axis=1)
                                                     .mean()['between_subject_variation'])
    rats_avg['between_subject_variation_CI95'] = 1.96*(rats_avg
                                                       .groupby(level=1, axis=1)
                                                       .sem()['between_subject_variation'])

    rats_avg.sort_index(axis=1, inplace=True)

    return rats_avg


def get_retest_results(study: str,
                       dataset_name: str,
                       cleaned_parameter_data: pd.DataFrame
                       ) -> Tuple[pd.DataFrame,
                                  pd.DataFrame]:
    """Gets retest results statistical data summary.

    Calculates both single-subject and substudy average repeatability
    and effect-size statistics and saves both in tabular csv format.

    Args:
        study: Study name of interest (e.g., 'Reproducibility')
        dataset_name: Name of dataset for saving results
            (e.g., 'saline_data')
        cleaned_paramter_data: Cleaned DataFrame containing estimated
            parameter variables.

    Returns:
        A tuple containing two dataframes with statistical summaries for
        repeatability and effect size statistics on single-subject and
        substudy average levels, respectively.
    """
    # Pivot single-subject data
    data_per_subject = pd.pivot_table(cleaned_parameter_data,
                                      values='Value',
                                      index=['Substudy',
                                             'Symbol',
                                             'Site',
                                             'Fstrength',
                                             'Time_period',
                                             'Rat'],
                                      columns='Day')
    # Get mean substudy values
    data_per_substudy = (cleaned_parameter_data
                         .groupby(['Substudy',
                                   'Symbol',
                                   'Day'])['Value'].mean().unstack())

    # Get single-subject repeatability/effect-size statistics
    subject_stats = get_subject_repeatability(data_per_subject)
    subject_stats_overall = get_summary_stats(subject_stats, data_per_subject,
                                              ['Substudy', 'Symbol'])
    # Get substudy average repeatability/effect-size statistics
    substudy_stats = get_subject_repeatability(data_per_substudy)
    substudy_stats_overall = get_summary_stats(substudy_stats,
                                               data_per_substudy, ['Symbol'])

    # Save single-subject retest statistics
    save_name_singlesubject = data.get_results_folder(study,
                                                      '02_analyses',
                                                      'repeatability',
                                                      None,
                                                      dataset_name + '_singlesubject',
                                                      'csv')
    subject_stats_overall.to_csv(save_name_singlesubject)

    # Save substudy-average retest statistics
    save_name_overall = data.get_results_folder(study,
                                                '02_analyses',
                                                'repeatability',
                                                None,
                                                dataset_name + '_overall',
                                                'csv')
    substudy_stats_overall.to_csv(save_name_overall)

    return subject_stats_overall, substudy_stats_overall


def get_mixed_anova(saline_data: pd.DataFrame,
                    biomarker: str
                    ) -> Tuple[list, list]:
    """Calculates variance contributions from mixed ANOVA.

    Performs a mixed ANOVA on saline retest dataset and determines
    relative contributions of variability (i.e., from between-subjects,
    between-centres, between-days, and residual mean squared errors).

    Args:
        saline_data: DataFrame containing retest data parameter values
            for a biomarker of interest after administration of saline
            on Day 1 and Day 2 scans.
        biomarker: Name of biomarker of interest (e.g., 'Ktrans')

    Returns:
        A tuple of two lists containing column names and percentage
        values for contribution of each source of variance to total
        variance, respectively.
    """
    retest_saline_data = saline_data.query('Symbol==@biomarker')
    # Create column denoting Site-Rat combinations
    retest_saline_data['Site_Rat'] = (retest_saline_data['Site']
                                      .astype(str)) + (retest_saline_data['Rat']
                                                       .astype(str))

    # Perform two-way ANOVA to get mean squared error (MS)
    # attributed to Day, Substudy, and interaction term (Day *Substudy)
    # and to get sum of squared error (SS) for Between-subjects (Residual)
    twoway_aov = pg.anova(data=retest_saline_data,
                          dv='Value',
                          between=['Day', 'Substudy'])

    MSsubstudy = twoway_aov.query("Source=='Substudy'")['MS'].values[0]
    MSday = twoway_aov.query("Source=='Day'")['MS'].values[0]
    MSinteraction = (twoway_aov
                     .query("Source=='Day * Substudy'")['MS'].values[0])

    # Perform mixed ANOVA to get degrees of freedom 2 (DF2)
    # for calculating MS Between-subjects
    # MSbetween_subjects = SSbetween_subjects / DF2
    # DF2 = (total #of observations) - (#substudies - 1) - (#days -1)
    # => 9 -(3-1)-(2-1) = 6
    mixed_aov = pg.mixed_anova(data=retest_saline_data,
                               dv='Value',
                               within='Day',
                               subject='Site_Rat',
                               between='Substudy')

    DF2 = mixed_aov['DF2'][0]
    MSbetweensubjects = (twoway_aov
                         .query("Source=='Residual'")['SS'] / DF2).values[0]

    # Calculate mean squared error (MS) (i.e., total variation)
    MStot = MSsubstudy + MSday + MSinteraction + MSbetweensubjects

    variance_names = ['reproducibility_between_substudies',
                      'repeatability_between_day',
                      'interaction_effect',
                      'between_subject_variation']
    variance_values = [(MSsubstudy / MStot)*100,
                       (MSday / MStot)*100,
                       (MSinteraction / MStot)*100,
                       (MSbetweensubjects / MStot)*100]

    return variance_names, variance_values
