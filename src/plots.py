"""Data visualisation.

A module which creates plots for visualising data distributions.
"""
# imports
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import data


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
    extracted_curves = (pd.
                        concat([signals[drug][1][f"Average Delta R1 {time_curve} (s-1)"],
                                signals[drug][2][f"Average Delta R1 {time_curve} (s-1)"]['Average deltaR1 (s-1)']],
                               axis=1))
    extracted_curves.columns = ['Time (s)', 'Control', 'Treatment']

    extracted_curves.dropna(inplace=True)

    return extracted_curves


# Define default settings for plots
sns.set_theme(style="white", color_codes=True)
sns.set(style="whitegrid", font_scale=2, rc={"savefig.dpi": 300})


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
    plt.figure(figsize=(12, 8))
    # Create the plot
    plt.plot(time_series, fitted_signal, color='green',
             label='Simulated Liver Signal')
    plt.plot(sample_times, observed_signal, marker="x", color='green',
             linewidth=0, label='Observed Liver Signal')
    plt.xlabel("Time [sec]")
    plt.ylabel('Signal [a.u.]')
    plt.ylim(bottom=0, top=16)
    plt.title(metadata['drug'] + ' (Rat ' + str(metadata['subject']) +
              ', Day ' + str(metadata['day']) + ')')
    plt.legend(loc='best')
    plt.tight_layout()

    save_name = data.get_results_folder(study,
                                        '01_model_outputs',
                                        'figures',
                                        'per_rat',
                                        f"fig_{filename}",
                                        'png')
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
    observed['Time (min)'] = observed['Time (s)']/60

    plt.rcParams["axes.labelsize"] = 50
    plt.rcParams["axes.titlesize"] = 50
    plt.rcParams["axes.labelweight"] = 'bold'
    plt.rcParams["axes.titleweight"] = 'bold'
    plt.rcParams["font.weight"] = 'bold'
    plt.rcParams['savefig.dpi'] = 300
    plt.rc('axes', linewidth=2)
    plt.rc('xtick', labelsize=40) 
    plt.rc('ytick', labelsize=40) 

    plt.figure(figsize=(12, 10))
    # use set_position
    ax = plt.gca()
    #ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('k')
    #ax.spines['right'].set_color('none')
    #ax.spines['bottom'].set_position('zero')

    ax.axhline(y=0, lw=3, color='k')
    
    plt.errorbar(data=observed,
                 x="Time (min)", y='Control',
                 yerr=observed['Control'].std(),
                 color='#00008B',
                 fmt="o", elinewidth=2, capthick=3, capsize=6,
                 markersize=12, markeredgewidth=2, alpha=0.95,
                 label="Control - observed")

    plt.errorbar(data=observed,
                 x="Time (min)", y='Treatment',
                 yerr=observed['Treatment'].std(), color='#8B2323',
                 fmt="v", elinewidth=2, capthick=3, capsize=6,
                 markersize=12, markeredgewidth=2, alpha=0.95,
                 label="Treatment - observed")

    if is_fitted is True:
        fitted = extract_curves(signals, drug, f"{ROI} fit")
        fitted['Time (min)'] = fitted['Time (s)']/60
        g = sns.lineplot(data=fitted, x="Time (min)", y='Control',
                         ls='-', linewidth=4, label="Control - fitted")

        g = sns.lineplot(data=fitted, x="Time (min)", y='Treatment',
                         ls='-', linewidth=4, label="Treatment - fitted")

        fig_name = f"{drug}_{ROI}_deltaR1_fitted"
    else:
        g = sns.lineplot(data=observed, x="Time (min)", y='Control', ls='',)
        fig_name = f"{drug}_{ROI}_deltaR1"

    plt.suptitle(f"Group mean {ROI} gadoxetate profiles in control and inhibitory phases \n (error bars represent standard deviation)")
    g.set_title(f"{drug}", weight='bold')
    g.set_xlabel("Time [min]", weight='bold')
    g.set_ylabel("\u0394 $R_{1}$ [$s^{-1}$]", weight='bold')

    g.set(ylim=YLIM)
    g.set(xlim=(0,30))
    g.get_legend().remove()
    #g.legend(loc='best', framealpha=1)
    plt.tight_layout()

    save_name = data.get_results_folder(study,
                                        '01_model_outputs',
                                        'figures',
                                        'per_drug',
                                        fig_name,
                                        'png')
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
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams["axes.labelsize"] = 50
    plt.rcParams["axes.titlesize"] = 50
    plt.rcParams["axes.labelweight"] = 'bold'
    plt.rcParams["axes.titleweight"] = 'bold'
    plt.rcParams["font.weight"] = 'bold'
    plt.rc('axes', linewidth=2)
    plt.rc('xtick', labelsize=40) 
    plt.rc('ytick', labelsize=40) 
    plt.rcParams["lines.linewidth"] = 4
    plt.rcParams['lines.markersize'] = 12

    # use set_position
    ax = plt.gca()
    ax.spines['left'].set_color('k')
    ax.spines['bottom'].set_color('k')
    
    g = sns.catplot(data=effect_size_data[effect_size_data['Symbol'] == biomarker],
                    x='Day', y='Value', sharey=True,
                    hue='Rat', col='Drug', col_wrap=4,
                    kind="point",
                    height=8, aspect=1.2,
                    legend=False,
                    errorbar=None)

    (g.set_xticklabels(["Control", "Treatment"], weight='bold')
     .set_titles("{col_name}", weight='bold')
     .despine(left=False))

    if biomarker == 'Ktrans':
        (g.set_axis_labels("", "$K^{trans}$ [mL/min/mL]", weight='bold')
         .set(ylim=(0, 1.4)))
    elif biomarker == 'kbh':
        (g.set_axis_labels("", "$k_{bh}$ [mL/min/mL]", weight='bold')
         .set(ylim=(0, 0.35)))
    elif biomarker == 'khe':
        (g.set_axis_labels("", "$k_{he}$ [mL/min/mL]", weight='bold')
         .set(ylim=(0, 15)))

    for ax in g.axes.flatten():
        ax.tick_params(labelleft=True, labelbottom=True)

    plt.legend(title='Rat', bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    g.fig.tight_layout()
    save_name = data.get_results_folder(study,
                                        '02_effect_sizes',
                                        'figures',
                                        None,
                                        f"{biomarker}_plot",
                                        'png')
    plt.savefig(save_name)
    plt.close()
