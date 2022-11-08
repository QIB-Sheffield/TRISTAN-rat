"""Data analysis.

A module which cleans missing and erroneous data and
provides statistical summaries of the data.
"""
# imports
import pandas as pd
from scipy import stats
import data
from typing import Tuple


def remove_data_errors(parameter_data: pd.DataFrame,
                       study: str
                       ) -> pd.DataFrame:
    """Cleans data.

    Removes computational fitting errors in biomarker prediction caused by
    non-convergence of the model fitting algorithm, i.e., the output of
    boundary values instead of true values, e.g., gadoxetate extraction
    fraction, E=100%. Additionally removes cases where only one acquisition
    per subject is present (i.e., day 1 or day 2 data are missing).

    Args:
        parameter_data: DataFrame containing estimated parameter variables.
        study: Study name of interest (e.g., 'SixTestCompounds').

    Returns:
        Cleaned DataFrame.
    """
    data_pivoted = pd.pivot_table(parameter_data,
                                  values='Value',
                                  columns=['Symbol'],
                                  index=['Drug',
                                         'Site',
                                         'Site_drug',
                                         'Rat',
                                         'Day'])

    # Remove computational fitting errors based on subjects where gadoxetate
    # extraction fraction, E is close or equal to 100% (i.e., >= 99%)
    fit_errors = data_pivoted[data_pivoted['E'] >= 99.95]
    fit_errors_removed = (data_pivoted[~data_pivoted
                                       .index.isin(fit_errors.index)])

    # Save index metadata for computational fitting errors
    save_name = data.get_results_folder(study,
                                        '02_effect_sizes',
                                        None,
                                        None,
                                        'fit_errors',
                                        'txt')
    with open(save_name, "w") as output:
        output.write(str(list([fit_errors.index])))

    # Remove subjects with missing acquisition on day 1 or day 2
    missing_days_removed = (fit_errors_removed[fit_errors_removed
                                               .groupby(['Site',
                                                         'Drug',
                                                         'Rat'])
                                               .transform('count') > 1]
                            .dropna())
    cleaned_parameter_data = missing_days_removed.stack().reset_index()
    cleaned_parameter_data.rename(columns={0: 'Value'}, inplace=True)

    return cleaned_parameter_data


def get_stats(cleaned_parameter_data: pd.DataFrame,
              variables: list
              ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Gets statistical data summary.

    Calculates summary statistics conveying the observed effects of
    administered test compounds on the biomarkers of interest, i.e.,
    between day 1 and day 2, individual subject and average subject
    values per site are calculated for gadoxetate percent change, effect
    size + 95% confidence intervals (CI95), and paired T-test p-values.

    Args:
        rats: Cleaned DataFrame containing estimated parameter variables.
        variables: List of condition variables to group by for statistical
            summary (e.g., ['drug', 'site']).

    Returns:
        A tuple (rats, rats_avg), where rats is a DataFrame containing
        statistical summaries per rat, whereas rats_avg is a DataFrame
        containing statistical summaries for average rat values per drug.
    """
    rats = cleaned_parameter_data.copy()
    rats['pct_change'] = rats.pct_change(axis=1).mul(100)[2]
    rats['diff'] = rats[1] - rats[2]

    rats_avg = (rats.groupby(level=variables)[[1, 2, 'pct_change']]
                .agg(['mean', 'sem']))

    rats_avg['simple_effect_size'] = (rats.groupby(level=variables)[1]
                                      .mean() -
                                      rats.groupby(level=variables)[2]
                                      .mean())

    rats_avg['lowerCI95'] = (rats_avg['simple_effect_size'] -
                             rats.groupby(level=variables)['diff']
                             .sem().mul(1.96))

    rats_avg['upperCI95'] = (rats_avg['simple_effect_size'] +
                             rats.groupby(level=variables)['diff']
                             .sem().mul(1.96))

    rats.dropna(inplace=True)
    rats_avg['p-value'] = (rats.groupby(level=variables)
                           .apply(lambda df: stats
                                  .ttest_rel(df[1], df[2])[1]))

    rats_avg.sort_index(axis=1, inplace=True)

    return rats, rats_avg


def save_effect_sizes(cleaned_parameter_data: pd.DataFrame,
                      params: list,
                      variables: list,
                      study: str
                      ) -> None:
    """Saves average effect size summaries per dug.

    Obtains statistical summaries for specific MRI biomarkers of interest
    and saves output as csv in results folder.

    Args:
        cleaned_parameter_data: Cleaned DataFrame containing estimated
            parameter variables.
        params: MRI biomarkers of interest.
        variables: List of condition variables to group by for statistical
            summary (e.g., ['drug', 'site']).
        study: Study name of interest (e.g., 'SixTestCompounds').
        """
    data_pivoted = (pd
                    .pivot_table(cleaned_parameter_data[cleaned_parameter_data['Symbol']
                                                        .isin(params)],
                                 values='Value',
                                 index=variables + ['Rat'],
                                 columns='Day',
                                 margins=True,
                                 margins_name='mean')[:-1])

    # Get tables
    _, per_drug = get_stats(data_pivoted.query('Symbol in @params'), variables)
    effect_sizes = per_drug[['simple_effect_size',
                             'lowerCI95',
                             'upperCI95',
                             'p-value']].swaplevel(axis=1)
    effect_sizes.reset_index(inplace=True)
    effect_sizes.sort_index(axis=1, inplace=True)
    effect_sizes_pivoted = (effect_sizes
                            .pivot_table(columns=['Symbol'],
                                         index='Drug').swaplevel(axis=1))
    effect_sizes_pivoted.sort_index(axis=1, inplace=True)

    save_name = data.get_results_folder(study,
                                        '02_effect_sizes',
                                        None,
                                        None,
                                        'effect_sizes',
                                        'csv')
    effect_sizes_pivoted.to_csv(save_name)
