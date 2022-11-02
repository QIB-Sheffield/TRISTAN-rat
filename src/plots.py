"""Data visualisation

A module which 

"""
# visualisations.py
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import data


# Define default settings for plots
sns.set_theme(style="white", color_codes=True)
sns.set(font_scale=3)
sns.set_style('whitegrid')

def get_signal_plots(study: str,
                     filename: str,
                     time: np.ndarray,
                     signal: np.ndarray,
                     sample_times: np.ndarray,
                     liver_data: np.ndarray,
                     metadata: dict
) -> None:
    
    plt.figure(figsize=(12,8))
    # Create the plot
    plt.plot(time, signal, color='green', label='Simulated Liver Signal')
    plt.plot(sample_times, liver_data, marker="x", color='green', linewidth=0, label='Observed Liver Signal')
    plt.xlabel("Time (sec)")
    plt.ylabel('Signal (a.u.)')
    plt.ylim(bottom=0, top=16)
    plt.title(metadata['drug'] + ' (Rat ' + str(metadata['subject']) + ', Day ' + str(metadata['day']) + ')')
    plt.legend(loc='best')
    save_name = data.get_results_folder(study, '01_model_outputs', 'figures', 'per_rat', f"fig_{filename}", 'png')
    plt.savefig(save_name)
    plt.close()


def pairplots(effect_size_data: pd.DataFrame, 
            biomarker: str,
            study: str
) -> None:
    
    g = sns.catplot(data=effect_size_data[effect_size_data['Symbol']==biomarker], x='Day', y='Value', hue='Rat', col='Drug',
                kind="point", sharey=True, col_wrap=4,
                height=8, aspect=1.2, legend=False, ci=None)

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


def split_fits(signals: dict,
               drug: str,
               time_curve: str
) -> pd.DataFrame:
    
    df = pd.concat([signals[drug][1][f"Average Delta R1 {time_curve} (s-1)"],
                   signals[drug][2][f"Average Delta R1 {time_curve} (s-1)"]['Average deltaR1 (s-1)']], 
                           axis=1)
    df.columns = ['Time (s)', 'Control', 'Treatment']
    
    df.dropna(inplace = True)
    
    return df
    

def get_deltaR1_plots(signals: dict,
                      drug: str,
                      ROI: str,
                      study: str,
                      *, is_fitted: bool = False,
                      YLIM: tuple
) -> None:
    
    observed = split_fits(signals, drug, ROI)
    
    plt.figure(figsize=(12,8))
    plt.errorbar(data=observed, x="Time (s)", y='Control', yerr=observed['Control'].std(), color='royalblue',
                             fmt="o", elinewidth=1, capthick=0.01, capsize=0.01, alpha=0.95, label = "Control - observed")
    plt.errorbar(data=observed, x="Time (s)", y='Treatment', yerr=observed['Treatment'].std(), color='darkorange',
                             fmt="x", elinewidth=1, capthick=0.01, capsize=0.01, alpha=0.95, label = "Treatment - observed")
        
    if is_fitted==True:
        fitted = split_fits(signals, drug, f"{ROI} fit")
        g = sns.lineplot(data=fitted, x="Time (s)", y='Control', ls='-', label = "Control - fitted")
        g = sns.lineplot(data=fitted, x="Time (s)", y='Treatment', ls='-', label = "Treatment - fitted")
        fig_name = f"{drug}_{ROI}_deltaR1_fitted"
    else:
        g = sns.lineplot(data=observed, x="Time (s)", y='Control', ls='',)
        fig_name = f"{drug}_{ROI}_deltaR1"

    # g.set_title(f"{drug} group mean gadoxetate profiles in control and inhibitory phases \n (error bars represent standard deviation)")
    g.set_title(f"{drug}", weight='bold')
    g.set_xlabel("Time (s)", weight='bold')
    g.set_ylabel("\u0394 $R_{1}$ ($s^{-1}$)", weight='bold')

    g.set(ylim=YLIM)
    g.legend(loc='best')
    
    save_name = data.get_results_folder(study, '01_model_outputs', 'figures', 'per_drug', fig_name, 'png')
    plt.savefig(save_name)
    plt.close()
