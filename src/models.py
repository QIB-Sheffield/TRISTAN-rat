"""Tracer kinetic modelling.

A module which excutes the tracer kinetic modelling upon a study
of interest using the TristanRat class, outputting fitted data and 
estimated parameter variables in tabular and graphical formats.
"""
# imports
import os
import sys
import pandas as pd
import numpy as np
import itertools
import math
import data
import plots
from rat import TristanRat
if sys.version_info >= (3, 8):
    from typing import Dict, TypedDict, List, Tuple
else:
    from typing import Dict, List, Tuple
    from typing_extensions import TypedDict


def split_groups(files: list,
               filenames: list
) -> Dict[str, Dict[str, Dict[int, pd.DataFrame]]]:
    """Splits signal data into control and treatment groups.

    A function to separate control (day 1) and treatment (day 2)
    time curve data for each respective compound administered to 
    rats.

    Args:
        files: List of files where signal data are stored.
        filenames: List of filenames for each time curve.

    Returns:
        A nested dictionary storing time curves for respective
        compounds ('drug') administered to rats. Within this, there
        are two nested dictionaries, dividing time curves into control
        (baseline, day '1') and treatment (follow-up, day '2') groups.
    """
    time_curves = {}
    for n, file in enumerate(files):
        metadata = data.get_metadata(filenames[n], file)

        time_curves[metadata['drug']] = {}
        time_curves[metadata['drug']][1] = {}
        time_curves[metadata['drug']][2] = {}
        
    for n, file in enumerate(files):
        metadata = data.get_metadata(filenames[n], file)
        time_curves[metadata['drug']][metadata['day']].update({metadata['subject']:metadata['signals']})
        
    return time_curves


def convert_to_deltaR1(combined_signals: pd.DataFrame,
                       time_curve: str,
                       R10: float,
                       FA: int,
                       TR: float,
                       signal_model: np.float64,
                       signals: Dict[str, Dict[str, Dict[int, pd.DataFrame]]],
                       metadata: dict
) -> None:
    """Converts signal to delta R1.

    A function to convert MRI signal (S(t)) time curve data to ∆R1
    data. The functions ∆R1,L(t) and ∆R1,S(t) measure the change 
    in relaxation rate induced by the gadoxetate in liver (L) and 
    spleen (S), respectively, as a function of time. They are derived 
    from the measured signals S(t) and precontrast relaxation rate R10 in 
    liver and spleen, respectively. The TRISTAN rat model uses fixed 
    literature-based values for the precontrast relaxation rates R10.
    
    Args:
        combined_signals: Dataframe containing observed and fitted 
            liver and spleen data for single rat from one 
            acquisition.
        time_curve: Time curve of interest, e.g., 'Liver fit'.
        R10: Precontrast relaxation rate corresponding to time
            curve of interest.
        FA: MRI sequence flip angle (degrees).
        TR: MRI sequence Repetitition Time (sec).
        signal_model: MRI signal model for time curve of interest.
        signals: Dictionary containing all observed and fitted 
            liver and spleen data for all rats and all acquistions.
        metadata: Metadata corresponding to specific rat and
            acquisition of interest.
    """
    S = combined_signals[time_curve]
    S0 = S[:4].mean() # S0 = precontrast signal (measured by averaging precontrast signals S(t))
    
    # Deriving relaxation rates from signals ->
    # Functions deltaR1(t) are derived from measured signals S(t) and precontrast R10
    # by inverting the signal model for a spoiled gradient echo sequence in the steady state:
    cFA = math.cos(FA*math.pi/180) # cos of FA in radians
    X = signal_model*(S/S0)
    R1 = (-1/TR)*np.log((1 - X)/(1 - cFA*X))

    dR1 = R1 - R10
    dR1 = dR1[np.logical_not(np.isnan(dR1))]
    
    signals[metadata['drug']][metadata['day']][metadata['subject']][f"Delta R1 {time_curve.replace(' (a.u.)', '')} (s-1)"] = dR1


def get_subject_list(signals: dict
) -> list:
    """Get subject index list per drug and day.
    
    Args:
        signals: Dictionary containing all observed and fitted 
            liver and spleen data for all rats and all acquistions.
    
    Returns:
        List of index combinations for each drug, day, and subject in
        study of interest.
    """
    subject_list = []
    for drug, day in list(itertools.product(signals.keys(), [1, 2])):
        subject_list.append([drug, day, list(signals[drug][day].keys())])
        
    return subject_list


def get_average_curves(signals: dict,
                        subject_list: list,
                        time_curve: str
) -> None:
    """Extracts time curve averages per drug and day.

    Calculates time curve averages over all subjects and stores
    result in new key within the 'signals' dictionary.

    Args:
        signals: Dictionary containing all observed and fitted 
            liver and spleen data for all rats and all acquistions.
        subject_list: List of index combinations for each drug, day, and 
            subject in study of interest.
        time_curve: Time curve of interest, e.g., 'Liver fit'.    
    """
    for i in subject_list:
        drug = i[0]
        subjectRange = i[2]
        day = i[1]
        signal = 0
        num_subjects = len(subjectRange)
        for subject in range(subjectRange[0], num_subjects):
            signal = signal + signals[drug][day][subject][time_curve]
            time_observed = signals[drug][day][subject]['Time (s)']
            time_fit = signals[drug][day][subject]['Time fit (s)']
            
        signal_average = signal/num_subjects
          
        if 'fit' in time_curve:
            average_signal = pd.DataFrame({'Time (s)': time_fit,
                              'Average deltaR1 (s-1)': signal_average})
        else:
            average_signal = pd.DataFrame({'Time (s)': time_observed,
                              'Average deltaR1 (s-1)': signal_average})

        signals[drug][day]['Average ' + time_curve] = average_signal


def fit_data(study: str,
            filenames: list,
            files: list,
            signals: Dict[str, Dict[str, Dict[int, pd.DataFrame]]],
            model: str
) -> pd.DataFrame:
    """Fits liver time_curve data.
    
    Args:
        study. Study name of interest (e.g., 'SixTestCompounds).
        filenames: List of filenames where MRI signal data are contained.
        files: List of files where MRI signal data are contained.
        signals: Dictionary containing all observed and fitted liver and 
            spleen data for all rats and all acquistions.
        model: Tracer kinetic model used for fitting MRI signal data.
        
    Returns:
        DataFrame containing estimated parameter variables.    
    """
    all_vars = None
    for n, file in enumerate(files):
        print("Fitting ", np.round(100*n/len(files),2), "%")
        print(filenames[n])
        rat = model()
        metadata = data.get_metadata(filenames[n], file)
        # Perform the fit
        signal_df = signals[metadata['drug']][metadata['day']][metadata['subject']]
        ts = signal_df["Time (s)"].values 
        rat.dt = ts[1] - ts[0]
        rat.dose = 0.0075

        # Assign site-specific criteria (from DICOM headers or relayed by sites)
        if metadata['site'] == 'E':
            rat.tstart = 4.0*60 + 45     
            rat.tduration = 30
            rat.field_strength = 7.0                              
            rat.FA = 20
            rat.TR = 5.8/1000
        elif metadata['site'] == 'G2': 
            rat.tstart = 4.0*60 + 45         
            rat.tduration = 30  
            rat.field_strength = 4.7                                 
            rat.FA = 20   
            rat.TR = 5.8/1000
        elif metadata['site'] == 'G1': 
            rat.tstart = 4.0*60 + 45         
            rat.tduration = 30  
            rat.field_strength = 7.0                                 
            rat.FA = 20   
            rat.TR = 5.8/1000
        elif metadata['site'] == 'D':
            rat.tstart = 4.0*60 + 45 + 7        
            rat.tduration = 22  
            rat.field_strength = 4.7                                
            rat.FA = 20   
            rat.TR = 5.8/1000

        R10L = rat.R10L
        R10S = rat.R10S
        liver_signal_model = rat._signal(R10L, 1)
        spleen_signal_model = rat._signal(R10S, 1)
        FA = rat.FA # actual Flip Angle (known sequence parameter - for rat data, assumption is made that actual FA = nominal FA)
        TR = rat.TR # Repetition Time (known sequence parameter)

        rat.set_liver_data(ts, signal_df["Liver (a.u.)"].values)
        rat.set_spleen_data(ts, signal_df["Spleen (a.u.)"].values)
        rat.fit_standard()
        
        plots.get_signal_plots(study, filenames[n], 
                                rat.t, rat.liver_signal, rat.liver_sampling_times,
                                rat.liver_data, metadata)

        fitted_signals_df = pd.DataFrame({"Time fit (s)": rat.t})
        fitted_signals_df["Spleen fit (a.u.)"] = rat.spleen_signal
        fitted_signals_df["Liver fit (a.u.)"] = rat.liver_signal
        combined_signals = pd.concat([signal_df, fitted_signals_df], axis=1)
        signals[metadata['drug']][metadata['day']][metadata['subject']] = combined_signals
        col_names = combined_signals.columns
        liver_curves = [x for x in col_names if 'Liver' in x]
        spleen_curves = [x for x in col_names if 'Spleen' in x]

        for curve in liver_curves:
            convert_to_deltaR1(combined_signals, curve, R10L, FA, TR, liver_signal_model, signals, metadata)
        for curve in spleen_curves:
            convert_to_deltaR1(combined_signals, curve, R10S, FA, TR, spleen_signal_model, signals, metadata)

        save_name = data.get_results_folder(study, '01_model_outputs',
                                       'relaxation_rates_and_signals', None, f"fit_{filenames[n][:-8]}", 'csv')
        combined_signals.to_csv(save_name)

        # Create DataFrame for storing estimated parameters
        vars = rat.export_variables()
        name = pd.DataFrame({"Data file": [file]*vars.shape[0]})
        drug = pd.DataFrame({"Drug": [metadata['drug']]*vars.shape[0]})
        site = pd.DataFrame({"Site": [metadata['site']]*vars.shape[0]})
        subj = pd.DataFrame({"Rat": [metadata['subject']]*vars.shape[0]})
        day = pd.DataFrame({"Day": [metadata['day']]*vars.shape[0]})
        vars = pd.concat([name, drug, site, subj, day, vars], axis=1)

        # Add to main output
        if all_vars is None:
            all_vars = vars
        else:
            all_vars = pd.concat([all_vars, vars], axis=0)

    save_name = data.get_results_folder(study, '01_model_outputs', None, None, 'all_parameters', 'csv')
    try:
        all_vars.to_csv(save_name)
    except:
        print("Can't write to file ", save_name)
        print("Please close the file before saving data")
        
    return all_vars
