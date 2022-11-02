"""Statistical analysis

A module which 

"""


# effect_sizes.py
import pandas as pd
import numpy as np
from scipy import stats
import data

def remove_data_errors(effect_size_data: pd.DataFrame,
                        study: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    data_pivoted = pd.pivot_table(effect_size_data, values = 'Value', columns=['Symbol'],
                   index=['Drug', 'Site', 'Site_drug', 'Rat', 'Day'])
    fit_errors = data_pivoted[data_pivoted['E']>=95]
    clean_data = data_pivoted[~data_pivoted.index.isin(fit_errors.index)]
    no_missing_days = clean_data[clean_data.groupby(['Site', 'Drug', 'Rat']).transform('count')>1].dropna()
    no_missing_days = no_missing_days.stack().reset_index()
    no_missing_days.rename(columns={0:'Value'}, inplace=True)
    
    save_name = data.get_results_folder(study, '02_effect_sizes', None, None, 'fit_errors', 'txt')
    with open(save_name, "w") as output:
        output.write(str(list([fit_errors.index])))
    
    return no_missing_days


def get_stats(rats: pd.DataFrame,
            variables: list
) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    rats['pct_change'] = rats.pct_change(axis=1).mul(100)[2]
    rats['diff'] = rats[1] - rats[2]
    
    rats_avg = rats.groupby(level=variables)[[1, 2, 'pct_change']].agg(['mean', 'sem'])
    rats_avg['simpleES'] = rats.groupby(level=variables)['diff'].mean()
    rats_avg['simpleES CI95'] = rats.groupby(level=variables)['diff'].sem().mul(1.96)
    rats.dropna(inplace=True)
    rats_avg['p-value'] = rats.groupby(level=variables).apply(lambda df: stats.ttest_rel(df[1], df[2])[1])
    
    rats_avg.sort_index(axis=1, inplace=True)
            
    return rats, rats_avg


def save_effect_sizes(effect_size_data: pd.DataFrame,
                        params: list,
                        variables: list,
                        study: str):

        data_pivoted = pd.pivot_table(effect_size_data[effect_size_data['Symbol'].isin(params)], 
                                            values = 'Value',
                                            index=variables + ['Rat'],
                                            columns='Day',
                                            margins=True,
                                            margins_name='mean')[:-1]

        ## Get tables
        perRat, perDrug = get_stats(data_pivoted.query('Symbol in @params'), variables)
        effect_sizes = perDrug[['simpleES', 'simpleES CI95', 'p-value']].swaplevel(axis=1)
        effect_sizes.reset_index(inplace=True)
        effect_sizes.sort_index(axis=1, inplace=True)
        effect_sizes_pivoted = effect_sizes.pivot_table(columns=['Symbol'], index='Drug').swaplevel(axis=1)
        effect_sizes_pivoted.sort_index(axis=1, inplace=True)

        save_name = data.get_results_folder(study, '02_effect_sizes', None, None, 'effect_sizes', 'csv')
        effect_sizes_pivoted.to_csv(save_name)