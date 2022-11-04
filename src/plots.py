"""Data visualisation.

A module which creates plots for visualising data distributions.
"""
# imports
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import data
import sys
if sys.version_info >= (3, 8):
    from typing import Dict, TypedDict, List, Tuple
else:
    from typing import Dict, List, Tuple
    from typing_extensions import TypedDict


# Helper function
def extract_curves(signals: dict,
                    drug: str,
                    time_curve: str
) -> pd.DataFrame:
    """Extract observed or fitted deltaR1 time curve data.

    Args:
        signals: Dictionary containing all observed and fitted 
            liver and spleen data for all rats and all acquistions.
        drug: Test compound of interest.
        time_curve: DeltaR1 time curve of interest (e.g., 'Liver' or 
            'Liver fit').
    Returns:
        DataFrame containing extracted deltaR1 time curve data.
    """
    extracted_curves = pd.concat([signals[drug][1][f"Average Delta R1 {time_curve} (s-1)"],
                                signals[drug][2][f"Average Delta R1 {time_curve} (s-1)"]['Average deltaR1 (s-1)']], 
                                axis=1)
    extracted_curves.columns = ['Time (s)', 'Control', 'Treatment']
    
    extracted_curves.dropna(inplace = True)
    
    return extracted_curves

# Define default settings for plots
sns.set_theme(style="white", color_codes=True)
sns.set(font_scale=3)
sns.set_style('whitegrid')


def get_signal_plots(study: str,
                     filename: str,
                     time_series: np.ndarray,
                     fitted_signal: np.ndarray,
                     sample_times: np.ndarray,
                     observed_signal: np.ndarray,
                     metadata: dict
) -> None:
    """Plots individual MRI signal intensity curves.
    
    Plots fitted vs. observed gadoxetate signal-time activity curves for the
    liver for each individual subject and save to results folder.

    Args:
        study: Study name of interest (e.g., 'SixTestCompounds').
        filename: File name of interest containing MRI signal-time data in
            csv format.
        time_series: MRI time series data (seconds).
        fitted_signal: MRI fitted liver signal data (a.u.).
        sample_times: Sampled MRI time series data (seconds).
        observed_signal: MRI observed liver signal data (a.u.).
        metadata: A dict mapping keys to the corresponding metadata and
            signal data contained in csv file.
    """
    plt.figure(figsize=(12,8))
    # Create the plot
    plt.plot(time_series, fitted_signal, color='green', label='Simulated Liver Signal')
    plt.plot(sample_times, observed_signal, marker="x", color='green', linewidth=0, label='Observed Liver Signal')
    plt.xlabel("Time (sec)")
    plt.ylabel('Signal (a.u.)')
    plt.ylim(bottom=0, top=16)
    plt.title(metadata['drug'] + ' (Rat ' + str(metadata['subject']) + ', Day ' + str(metadata['day']) + ')')
    plt.legend(loc='best')
    save_name = data.get_results_folder(study, '01_model_outputs', 'figures', 'per_rat', f"fig_{filename}", 'png')
    plt.savefig(save_name)
    plt.close()


def get_deltaR1_plots(signals: dict,
                      drug: str,
                      ROI: str,
                      study: str,
                      *, is_fitted: bool = False,
                      YLIM: tuple
) -> None:
    """Plots average MRI delta R1 curves per drug.
    
    Plots average gadoxetate deltaR1-time activity curves within a region
    of interest for each test compound administered.

    Args:
        signals: Dictionary containing all observed and fitted 
            liver and spleen data for all rats and all acquistions.
        drug: Test compound of interest.
        ROI: Region of interest (e.g., 'Liver' or 'Spleen').
        study: Study name of interest (e.g., 'SixTestCompounds').
        is_fitted: Bool which allows user to define whether signal time.
            curve of interest has been fitted. Default value = False and
            should be used only for unfitted signal time curve data.
        YLIM: Tuple containing upper and lower y-axis limits for the plot.
    """
    observed = extract_curves(signals, drug, ROI)
    
    plt.figure(figsize=(12,8))
    plt.errorbar(data=observed,
                    x="Time (s)", y='Control',
                    yerr=observed['Control'].std(),
                    color='royalblue',
                    fmt="o", elinewidth=1, capthick=0.01, capsize=0.01, alpha=0.95,
                    label = "Control - observed")

    plt.errorbar(data=observed,
                    x="Time (s)", y='Treatment',
                    yerr=observed['Treatment'].std(), color='darkorange',
                    fmt="x", elinewidth=1, capthick=0.01, capsize=0.01, alpha=0.95,
                    label = "Treatment - observed")
        
    if is_fitted==True:
        fitted = extract_curves(signals, drug, f"{ROI} fit")
        g = sns.lineplot(data=fitted, x="Time (s)", y='Control',
                            ls='-', label = "Control - fitted")

        g = sns.lineplot(data=fitted, x="Time (s)", y='Treatment',
                            ls='-', label = "Treatment - fitted")
                            
        fig_name = f"{drug}_{ROI}_deltaR1_fitted"
    else:
        g = sns.lineplot(data=observed, x="Time (s)", y='Control', ls='',)
        fig_name = f"{drug}_{ROI}_deltaR1"

    g.set_title(f"{drug} group mean {ROI} gadoxetate profiles in control and inhibitory phases \n (error bars represent standard deviation)")
    g.set_xlabel("Time (s)", weight='bold')
    g.set_ylabel("\u0394 $R_{1}$ ($s^{-1}$)", weight='bold')

    g.set(ylim=YLIM)
    g.legend(loc='best')
    
    save_name = data.get_results_folder(study, '01_model_outputs', 'figures', 'per_drug', fig_name, 'png')
    plt.savefig(save_name)
    plt.close()


def pairplots(effect_size_data: pd.DataFrame, 
                biomarker: str,
                study: str
) -> None:
    """Plots paired data distributions per biomarker.
    
    Plots estimated values for MRI biomarker of interest with connecting line
    between day 1 (baseline saline) and day 2 (follow-up test compound)
    datapoints, highlighting observed effect sizes caused by compounds.

    Args:
        effect_size_data: DataFrame containing statistical effect size 
            summaries.
        biomarker: MRI biomarker of interest (e.g., 'Ktrans').
        study: Study name of interest (e.g., 'SixTestCompounds').
    """
    g = sns.catplot(data=effect_size_data[effect_size_data['Symbol']==biomarker],
                        x='Day', y='Value', sharey=True, 
                        hue='Rat', col='Drug', col_wrap=4,
                        kind="point",
                        height=8, aspect=1.2,
                        legend=False,
                        errorbar=None)

    (g.set_xticklabels(["Control", "Treatment"], weight='bold')
     .set_titles("{col_name}", weight='bold')
     .despine(left=False))
    
    if biomarker=='Ktrans':
        (g.set_axis_labels("", "$K_{trans}$ (mL/min/mL)", weight='bold')
         .set(ylim=(0, 1.4)))
    elif biomarker=='kbh':
        (g.set_axis_labels("", "$k_{bh}$ (mL/min/mL)", weight='bold')
        .set(ylim=(0, 0.35)))
    elif biomarker=='khe':
        (g.set_axis_labels("", "$k_{he}$ (mL/min/mL)", weight='bold')
         .set(ylim=(0, 15)))
    
    for ax in g.axes.flatten():
        ax.tick_params(labelleft=True, labelbottom=True)
        
    plt.legend(title='Rat', bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    g.fig.tight_layout()
    save_name = data.get_results_folder(study, '02_effect_sizes', 'figures', None, f"{biomarker}_plot", 'png')
    plt.savefig(save_name)
    plt.close()
