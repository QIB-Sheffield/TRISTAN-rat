import statistics as stats
import math 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import data
from rat import TristanRat

def plot_time_curves():

    rat = TristanRat()
    rat.plot_whole_body_concentrations()
    rat.plot_signals()
    rat.plot_data()
    rat.plot_fit()

def fit_simulated():

    rat = TristanRat()
    rat.plot_fit()

def set_sequence():

    rat = TristanRat()
    rat.SNR0 = 100
    rat.FA = 50
    rat.TR = 2.0
    rat.plot_fit()

def noise_sensitivity():

    SNR0_range = range(5, 100, 1)
    E = []
    Th = []
    rat = TristanRat()
    for SNR0 in SNR0_range:
        print(SNR0)
        rat.SNR0 = SNR0
        rat.simulate_data()
        rat.fit()
        E.append(rat.liver_variables[0])
        Th.append(rat.liver_variables[1])
        rat.initialize_variables()
    plt.plot(SNR0_range, [rat.E]*len(E), color='red')
    plt.plot(SNR0_range, E, marker="x", color='red', linewidth=0)
    plt.show()
    plt.plot(SNR0_range, [rat.Th]*len(Th), color='red')
    plt.plot(SNR0_range, Th, marker="x", color='red', linewidth=0)
    plt.show()

def simulate_measurement():

    studyrat = TristanRat()
    studyrat.E = 0.1
    studyrat.Th = 5*60
    studyrat.SNR0 = 10
    studyrat.initialize_variables()
    studyrat.calculate_signals()
    studyrat.simulate_measurement()

    rat = TristanRat()
    rat.set_spleen_data(studyrat.spleen_sampling_times, studyrat.spleen_data)
    rat.set_liver_data(studyrat.liver_sampling_times, studyrat.liver_data)
    rat.fit()
    
    plt.plot(studyrat.t, studyrat.spleen_signal, color='blue')
    plt.plot(studyrat.spleen_sampling_times, studyrat.spleen_data, marker="x", color='blue', linewidth=0)
    plt.plot(rat.t, rat.spleen_signal, color='red')

    plt.plot(studyrat.t, studyrat.liver_signal, color='blue')
    plt.plot(studyrat.liver_sampling_times, studyrat.liver_data, marker="x", color='blue', linewidth=0)
    plt.plot(rat.t, rat.liver_signal, color='green')

    plt.show()

    print("Exact E: ", studyrat.E)
    print("Fitted E: ", rat.liver_variables[0])
    print("Exact Th: ", studyrat.Th)
    print("Fitted Th: ", rat.liver_variables[1])

def fit_data_example():

    df_data, filename = data.example()
    ts = df_data["Time (s)"].values

    rat = TristanRat()
    rat.dt = ts[1] - ts[0]
    rat.set_spleen_data(ts, df_data["Spleen (a.u.)"].values)
    rat.set_liver_data(ts, df_data["Liver (a.u.)"].values)
    rat.fit()

    plt.plot(rat.t, rat.spleen_signal, color='red')
    plt.plot(rat.spleen_sampling_times, rat.spleen_data, marker="x", color='red', linewidth=0)
    plt.plot(rat.t, rat.liver_signal, color='green')
    plt.plot(rat.liver_sampling_times, rat.liver_data, marker="x", color='green', linewidth=0)
    save_file = data.results_path() + 'fig_' + filename + ".png"
    plt.savefig(save_file)

    df_results = pd.DataFrame({"Time fit (s)": rat.t})
    df_results["Spleen fit (a.u.)"] = rat.spleen_signal
    df_results["Liver fit (a.u.)"] = rat.liver_signal
    df_output = pd.concat([df_data, df_results], axis=1)
    save_file = data.results_path() + 'fit_' + filename + ".csv"
    try:
        df_output.to_csv(save_file)
    except:
        print("Can't write to file ", save_file)
        print("Please close the file before saving data")
    save_file = data.results_path() + 'par_' + filename + ".csv"
    try:
        rat.export_variables().to_csv(save_file)
    except:
        print("Can't write to file ", save_file)
        print("Please close the file before saving data")


def report_six_test_compounds(all, study):

    save_path = data.results_path(study) + 'all_' + "SixTestCompounds" + '_'

    figdrug, axdrug = plt.subplots(2, 3)
    row = 0
    col = 0

    for drug in all.Drug.unique():
        label = []
        score = []
        for symbol in all.Symbol.unique():

            # Determine before and after values
            all_drug = all[(all.Drug == drug) & (all.Symbol == symbol)]
            data_pre = all_drug[all_drug.Day == 1]
            data_post = all_drug[all_drug.Day == 2]
            x = []
            y = []
            z = []
            n = []
            for rat in data_pre.Rat.values:
                if rat in data_post.Rat.values:
                    rat_data_pre = data_pre[data_pre.Rat == rat]
                    rat_data_post = data_post[data_post.Rat == rat]
                    xv = rat_data_pre.Value.values[0]
                    yv = rat_data_post.Value.values[0]
                    n.append(rat)
                    x.append(xv)
                    y.append(yv)
                    z.append(100*(yv-xv)/xv)
            zmean = stats.mean(z)
            zsdev = stats.pstdev(z)
            zerr = 1.96 * zsdev/math.sqrt(len(z))
            label.append(symbol)
            if zerr != 0:
                score.append(abs(zmean)/zerr)
            else:
                score.append(0)

            # plot before vs after per drug and per parameter
            fig, ax = plt.subplots()
            plt.plot(x, y, 'bx', linewidth=0)
            for i, rat in enumerate(n):
                plt.text(x[i], y[i], str(rat), 
                    ha='left', va='top',
                    transform = ax.transData)
            x.append(0)
            plt.plot(x, x, color='black', linewidth=1)
            text = 'Drug effect (%): ' + str(round(zmean,0)) + '+/- ' + str(round(zerr,0))
            plt.text(
                0.0, 1.0, text, 
                ha='left', va='top',
                transform = fig.transFigure)
            plt.savefig(save_path + drug + '_' + symbol + ".png")
            plt.clf()
            plt.close()
        
        # plot significance score per drug
        ind = np.arange(len(label))

    #    plt.bar(ind, score)
    #    plt.xlabel("Parameter")
    #    plt.ylabel('Effect size')
    #    plt.title(drug)
    #    plt.xticks(ind, label)
    #    plt.plot(ind, 1+ind*0, color='red', linewidth=1)
    #    plt.savefig(save_path + '1.' + drug + '_' + ".png")
    #    plt.clf()
    #    plt.close()

        axdrug[row,col].bar(label, score)
        if row==1:
            axdrug[row,col].set(xlabel="Parameter")
        if col==0:
            axdrug[row,col].set(ylabel = 'Effect size')
        axdrug[row,col].set_title(drug)
        axdrug[row,col].plot(ind, np.ones(len(label)), color='red', linewidth=1)
        axdrug[row,col].margins(0)
        axdrug[row,col].set_xticklabels(label, rotation=90)

        if col < 2:
            col += 1
        else:
            col = 0
            row += 1

    figdrug.tight_layout()
    figdrug.savefig(save_path + 'all_effect_sizes' + ".png")
#    figdrug.clf()
#    figdrug.close()


def fit_study(study):
    
    all_vars = None
    files, filenames = data.folder(study)
    for f, file in enumerate(files):

        print("Fitting ", 100*f/len(files), "%")
        print(filenames[f])
        descr = data.description(filenames[f])

        # Perform the fit
        df_data = pd.read_csv(file)
        ts = df_data["Time (s)"].values
        rat = TristanRat()
        rat.dt = ts[1] - ts[0]
        rat.dose = 0.0075
        if descr['site'] == 'MSD':
            rat.tstart = 4.0*60 + 45  # no info - check       
            rat.tduration = 30  # no info - check
            rat.field_strength = 7.0                              
            rat.FA = 30   # From DICOM header - check
            rat.TR = 5.8/1000
        elif descr['site'] == 'Sanofi': 
            rat.tstart = 4.0*60 + 45         
            rat.tduration = 30  
            rat.field_strength = 7.0                                 
            rat.FA = 20   
            rat.TR = 5.8/1000
        elif descr['site'] == 'Bayer':
            rat.tstart = 4.0*60 + 45 + 7        
            rat.tduration = 22  
            rat.field_strength = 4.7                                
            rat.FA = 20   
            rat.TR = 5.8/1000
        elif descr['site'] == 'Antaros':
            rat.tstart = 3.0*60 + 45 + 7               
            rat.tduration = 15  
            rat.field_strength = 7.0                  
            rat.FA = 30   
            rat.TR = 36.0/1000 # check
        rat.set_spleen_data(ts, df_data["Spleen (a.u.)"].values)
        rat.set_liver_data(ts, df_data["Liver (a.u.)"].values)
        # rat.fit_direct()
        # rat.fit()
        # rat.fit_joint()
        rat.fit_standard()

        # Create the plot
        save_file = data.results_path(study) + 'fig_' + filenames[f] + ".png"
        plt.plot(rat.t, rat.spleen_signal, color='red')
        plt.plot(rat.spleen_sampling_times, rat.spleen_data, marker="x", color='red', linewidth=0)
        plt.plot(rat.t, rat.liver_signal, color='green')
        plt.plot(rat.liver_sampling_times, rat.liver_data, marker="x", color='green', linewidth=0)
        plt.xlabel("Time (sec)")
        plt.ylabel('Signal (a.u.)')
        # plt.ylim(bottom=0)
        plt.ylim(bottom=0, top=16)
        plt.title(descr['drug'] + ' (Rat ' + str(descr['subject']) + ', Day ' + str(descr['day']) + ')')
        plt.savefig(save_file)
        plt.clf()

        # Save time curves
        save_file = data.results_path(study) + 'fit_' + filenames[f] + ".csv"
        df_results = pd.DataFrame({"Time fit (s)": rat.t})
        df_results["Spleen fit (a.u.)"] = rat.spleen_signal
        df_results["Liver fit (a.u.)"] = rat.liver_signal
        df_output = pd.concat([df_data, df_results], axis=1)
        try:
            df_output.to_csv(save_file)
        except:
            print("Can't write to file ", save_file)
            print("Please close the file before saving data")

        # Save parameter file
        vars = rat.export_variables()
        name = pd.DataFrame({"Data file": [file]*vars.shape[0]})
        drug = pd.DataFrame({"Drug": [descr['drug']]*vars.shape[0]})
        site = pd.DataFrame({"Site": [descr['site']]*vars.shape[0]})
        subj = pd.DataFrame({"Rat": [descr['subject']]*vars.shape[0]})
        day = pd.DataFrame({"Day": [descr['day']]*vars.shape[0]})
        vars = pd.concat([name, drug, site, subj, day, vars], axis=1)

        # Add to main output
        if all_vars is None:
            all_vars = vars
        else:
            all_vars = pd.concat([all_vars, vars], axis=0)

    save_file = data.results_path(study) + 'all_' + study + ".csv"
    try:
        all_vars.to_csv(save_file)
    except:
        print("Can't write to file ", save_file)
        print("Please close the file before saving data")

    if study == "SixTestCompounds":
        report_six_test_compounds(all_vars, study)


def run_all():

    plot_time_curves()
    fit_simulated()
    set_sequence()
    noise_sensitivity()
    simulate_measurement()

# plot_time_curves()
# fit_simulated()
# set_sequence()
# noise_sensitivity()
# simulate_measurement()
# fit_data_example()
fit_study("MultipleDosing")

# run_all()