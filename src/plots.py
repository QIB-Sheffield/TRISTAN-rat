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
                   substudy: str,
                   time_curve: str
                   ) -> pd.DataFrame:
    """Extract observed or fitted deltaR1 time curve data.

    Args:
        signals: Dictionary containing all observed and fitted
            liver and spleen data for all rats and all acquistions.
        substudy: Substudy of interest.
        time_curve: DeltaR1 time curve of interest (e.g., 'Liver' or
            'Liver fit').
    Returns:
        DataFrame containing extracted deltaR1 time curve data.
    """
    extracted_curves = (pd.
                        concat([signals[substudy][1][f"Average Delta R1 {time_curve} (s-1)"],
                                signals[substudy][2][f"Average Delta R1 {time_curve} (s-1)"]['Average deltaR1 (s-1)']],
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
    plt.title(metadata['substudy'] + ' (Rat ' + str(metadata['subject']) +
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
                      substudy: str,
                      ROI: str,
                      study: str,
                      *, is_fitted: bool = False,
                      YLIM: tuple
                      ) -> None:
    """Plots average MRI delta R1 curves per substudy.

    Plots average gadoxetate deltaR1-time activity curves within a region
    of interest for each test compound administered.

    Args:
        signals: Dictionary containing all observed and fitted
            liver and spleen data for all rats and all acquistions.
        substudy: Substudy of interest.
        ROI: Region of interest (e.g., 'Liver' or 'Spleen').
        study: Study name of interest (e.g., 'SixTestCompounds').
        is_fitted: Bool which allows user to define whether signal time.
            curve of interest has been fitted. Default value = False and
            should be used only for unfitted signal time curve data.
        YLIM: Tuple containing upper and lower y-axis limits for the plot.
    """
    observed = extract_curves(signals, substudy, ROI)
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
    ax.spines['left'].set_color('k')

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
        fitted = extract_curves(signals, substudy, f"{ROI} fit")
        fitted['Time (min)'] = fitted['Time (s)']/60
        g = sns.lineplot(data=fitted, x="Time (min)", y='Control',
                         ls='-', linewidth=4, label="Control - fitted")

        g = sns.lineplot(data=fitted, x="Time (min)", y='Treatment',
                         ls='-', linewidth=4, label="Treatment - fitted")

        fig_name = f"{substudy}_{ROI}_deltaR1_fitted"
    else:
        g = sns.lineplot(data=observed, x="Time (min)", y='Control', ls='',)
        fig_name = f"{substudy}_{ROI}_deltaR1"

    plt.suptitle(f"Group mean {ROI} gadoxetate profiles in control and \
                 \n inhibitory phases \
                 \n (error bars represent standard deviation)")
    g.set_title(f"{substudy}", weight='bold')
    g.set_xlabel("Time [min]", weight='bold')
    g.set_ylabel("\u0394 $R_{1}$ [$s^{-1}$]", weight='bold')

    g.set(ylim=YLIM)
    g.set(xlim=(0, 30))
    g.get_legend().remove()
    plt.tight_layout()

    save_name = data.get_results_folder(study,
                                        '01_model_outputs',
                                        'figures',
                                        'per_substudy',
                                        fig_name,
                                        'png')
    plt.savefig(save_name)

    plt.close()


def plot_distributions(study: str,
                       variable_name: dict,
                       group_data: pd.DataFrame,
                       variable: str,
                       constant: str,
                       benchmarks: pd.DataFrame,
                       order: list
                       ) -> None:
    """Biomarker distribution plots per variance group of interest.

    Plots MRI biomarker distribution plots for a variable of interest,
    e.g., across time periods (variable) at a specific centre (constant).
    Overlaid with lines reprenting biomarker benchmark values derived as
    the mean +/-95% CI across all substudies.

    Args:
        study: Study name of interest (e.g., 'Reproducibility').
        variable_name: Name of variable of interest.
        group_data: DataFrame containing parameter values for a variable
        of interest.
        variable: Variable of interest (e.g., 'Time_period').
        constant: Factor to keep constant (e.g., 'Site').
        benchmarks: DataFrame containing biomarker benchmark values.
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

    ylabels = ['$K_{trans}$', '$k_{bh}$']
    means = benchmarks['mean']
    lower_cis = means - benchmarks['CI95']
    upper_cis = means + benchmarks['CI95']

    if constant is None:
        g = sns.catplot(data=group_data,
                        x=variable,
                        y='Value',
                        col='Symbol',
                        kind='point',
                        capsize=0.2,
                        sharey=False,
                        join=False,
                        height=14,
                        aspect=1.2,
                        color='k',
                        ci=95)

        (g.set_titles("")
         .axes[0, 1].set(ylim=([0, 0.4])))

        for i in range(len(ylabels)):
            g.axes[0, i].set_ylabel(f"{ylabels[i]} [mL/min/mL]")
            g.axes[0, i].axhline(means[i], color='blue', ls=':')
            g.axes[0, i].axhline(lower_cis[i], color='red', ls='--')
            g.axes[0, i].axhline(upper_cis[i], color='red', ls='--')

    else:
        g = sns.catplot(data=group_data,
                        x=variable,
                        y='Value',
                        col=constant,
                        row='Symbol',
                        kind='point',
                        capsize=0.2,
                        sharey='row',
                        sharex=False,
                        join=False,
                        height=12,
                        aspect=0.72,
                        color='k',
                        ci=95,
                        col_order=order)

        g.axes[1, 0].set(ylim=([0, 0.4]))

        for i in range(int(len(group_data.groupby(['Symbol'] + [constant])
                               .count().index)/2)):
            # Ktrans reference values
            g.axes[0, i].axhline(means[0], color='blue', ls=':')
            g.axes[0, i].axhline(lower_cis[0], color='red', ls='--')
            g.axes[0, i].axhline(upper_cis[0], color='red', ls='--')
            # kbh reference values
            g.axes[1, i].axhline(means[1], color='blue', ls=':')
            g.axes[1, i].axhline(lower_cis[1], color='red', ls='--')
            g.axes[1, i].axhline(upper_cis[1], color='red', ls='--')

    g.axes[0, 0].set(ylim=([0, 1.5]))
    g.fig.tight_layout()

    save_name = data.get_results_folder(study,
                                        '02_analyses',
                                        'figures',
                                        'reproducibility',
                                        f"fig_{variable_name}",
                                        'png')
    plt.savefig(save_name)
    plt.close()


def pairplots(study: str,
              variable_name: dict,
              pair_data: pd.DataFrame,
              x: str,
              hue: str,
              col: str,
              row: str,
              palette: str,
              error: int,
              ylabels: list,
              sharey: 'str'
              ) -> None:
    """Plots paired data distributions per biomarker.

    Plots estimated values for MRI biomarker of interest with connecting line
    between day 1 (baseline saline) and day 2 (follow-up test compound)
    datapoints, highlighting observed effect sizes caused by compounds.

    Args:
        study: Study name of interest (e.g., 'SixTestCompounds').
        variable_name: Name of variable of interest.
        pair_data: DataFrame containing parameter data.
        x: Variable to plot along x-axis.
        hue: Variable for hue settings.
        col: Variable to plot along column subplots.
        row: Variable to plot along row subplots.
        palette: List or string defining colour palette for plot
        (e.g., 'rocket')
        error: Error interval for error bars (e.g., 95 for 95% CI)
        ylabels: List of labels for y-axes.
    """
    g = sns.catplot(data=pair_data,
                    x=x,
                    y='Value',
                    hue=hue,
                    col=col,
                    row=row,
                    kind="point",
                    sharey=sharey,
                    palette=palette,
                    height=8,
                    aspect=1,
                    legend=False,
                    ci=error)
    
    if study=='Reproducibility':
        (g.set_titles("")
         .axes[0, 1].set(ylim=([0, 0.4])))
        g.axes[0, 0].set(ylim=([0, 1.5]))
        for i in range(len(ylabels)):
            g.axes[0, i].set_ylabel(f"{ylabels[i]} [mL/min/mL]")
    else:
        g.set_titles(template='{col_name}')
        for i in range(len(ylabels)):
            g.axes[i, 0].set_ylabel(f"{ylabels[i]} [mL/min/mL]")

    for ax in g.axes.flatten():
        ax.tick_params(labelleft=True, labelbottom=True)

    plt.legend(title=hue, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    g.fig.tight_layout()
    save_name = data.get_results_folder(study,
                                        '02_analyses',
                                        'figures',
                                        'repeatability',
                                        variable_name,
                                        'png')
    plt.savefig(save_name)
    plt.close()
