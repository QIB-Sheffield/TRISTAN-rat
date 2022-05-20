import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit


def expconv(T, time, a):
    """Convolve a 1D-array with a normalised exponential.

    expconv() uses an efficient and accurate numerical formula to calculate the convolution,
    as detailed in the appendix of Flouri et al., Magn Reson Med, 76 (2016), pp. 998-1006.

    Arguments
    ---------
    a : numpy array
        the 1D array to be convolved.
    time : numpy array
        the time points where the values of ca are defined
        these do not have to to be equally spaced.
    T : float
        the characteristic time of the the exponential function.
        time and T must be in the same units.

    Returns
    -------
    a numpy array of the same shape as ca.

    Example
    -------
    coming soon..

    """
    if T==0: return a

    n = len(time)
    f = np.zeros(n)
    x = (time[1:n] - time[0:n-1])/T
    da = (a[1:n] - a[0:n-1])/x
    E = np.exp(-x)
    E0 = 1-E
    E1 = x-E0
    add = a[0:n-1]*E0 + da*E1
    for i in range(0,n-1):
        f[i+1] = E[i]*f[i] + add[i]      
    return f


def propagate_2cxm(t, ca, KP, KE, KB):
    """Calculate the propagators for the individual compartments in the 2CXM 
    
    For details and notations see appendix of 
    Sourbron et al. Magn Reson Med 62:672–681 (2009)

    Arguments
    ---------

    t : numpy array
        time points (sec) where the input function is defined
    ca : numpy array
        input function (mmol/mL)
    KP : float
        inverse plasma MTT (sec) = VP/(FP+PS)
    KE : float
        inverse extracellular MTT (sec) = VE/PS
    KB : float
        inverse blood MTT (sec) = VP/FP

    Returns
    -------
    cp : numpy array
        concentration in the plasma compartment (mmol/mL)
    ce : numpy array
        concentration in the extracellular compartment (mmol/mL)

    Examples
    --------
    coming soon..

    """

    KT = KP + KE
    sqrt = math.sqrt(KT**2-4*KE*KB)

    Kpos = 0.5*(KT + sqrt)
    Kneg = 0.5*(KT - sqrt)

    cpos = expconv(1/Kpos, t, ca)
    cneg = expconv(1/Kneg, t, ca)

    Eneg = (Kpos - KB)/(Kpos - Kneg)

    cp = (1-Eneg)*cpos + Eneg*cneg
    ce = (cneg*Kpos - cpos*Kneg) / (Kpos -  Kneg) 

    return cp, ce



class TristanRat():

    # Internal time resolution & acquisition time
    _tstep = 0.5             # sec
    _tmax = 40*60.0            # Total acquisition time (sec)

    # Experimental variables
    tstart = 4.6*60         # Injection start (sec)
    tduration = 30          # Injection duration (sec)
    dose = 0.0075           # Injection dose (mmol)
    field_strength = 4.7    # Field strength (T)
    dt = 60                 # Sampling duration (sec)
    FA = 20                 # Flip angle (degrees)
    TR = 5.0/1000           # Repetition time (sec)
    SNR0 = 20               # SNR at baseline (simulation only)
    
    # Spleen parameters
    veS = 0.314         # Spleen extracellular volume (mL/mL)
    S0spleen = 250      # Baseline signal (a.u.) 
                        # Randomly chosen

    # Liver parameters
    Fb = 3.61/60        # Blood flow (mL/sec/mL)
                        # https://doi.org/10.1021/acs.molpharmaceut.1c00206
    E = 0.4             # Gadoxetate extraction fraction
    veL = 0.230         # Liver extracellular volumen (mL/mL)
    Th = 20*60          # Hepatocellular mean transit time (sec)
    vh = 0.722          # Hepatocellular volume fraction (mL/mL)
    S0liver = 200       # Baseline signal (a.u.) 
                        # Randomly chosen

    # Whole body parameters
    Hct = 0.418         # Hematocrit
                        # Cremer et al, J Cereb Blood Flow Metab 3, 254-256 (1983)
    VL = 8.47           # Liver volume (mL) 
                        # Scotcher et al 2021, DOI: 10.1021/acs.molpharmaceut.1c00206
                        # Supplementary material, Table S2 
    GFR = 0.023         # Glomerular Filtration Rate (mL/sec) 
                        # https://doi.org/10.1152/ajprenal.1985.248.5.F734
    P = 0.172           # Permeability-surface for the whole body, with permeability for gadoxetate (mL/sec).
                        # Estimated from rat repro study data using PBPK model
                        # Table 3 in Scotcher et al 2021
                        # DOI: 10.1021/acs.molpharmaceut.1c00206
    VB = 15.8           # Whole body blood volume (mL) 
                        # 0.06 X BW + 0.77, Assuming body weight (BW) = 250 g
                        # Lee and Blaufox. Blood volume in the rat. 
                        # J Nucl Med. 1985 Jan;26(1):72-6.
    VE = 30             # Whole body extracellular volume (mL) 
                        # All tissues, including liver
                        # Derived from Supplementary material, Table S2
                        # Scotcher et al 2021
                        # DOI: 10.1021/acs.molpharmaceut.1c00206

    @property
    def K(self):
        return self.GFR + self.Ktrans * self.VL

    @property
    def Fp(self):
        return (1-self.Hct) * self.Fb

    @property
    def VP(self):
        return (1-self.Hct) * self.VB 

    @property
    def Ktrans(self):
        return self.E * self.Fp 

    @property
    def rp(self):
        field = math.floor(self.field_strength)
        if field == 4.0: return 6.4     # relaxivity of blood in Hz/mM ;assume spleen relaxivity is the same
        if field == 7.0: return 6.2     # relaxivity of blood in Hz/mM ;assume spleen relaxivity is the same
        if field == 9.0: return 6.1     # relaxivity of blood in Hz/mM ;assume spleen relaxivity is the same

    @property
    def rh(self):
        field = math.floor(self.field_strength)
        if field == 4.0: return 7.6     # relaxivity of hepatocytes in Hz/mM
        if field == 7.0: return 6.0     # relaxivity of hepatocytes in Hz/mM
        if field == 9.0: return 6.1     # relaxivity of hepatocytes in Hz/mM

    @property
    def R10L(self):
        field = math.floor(self.field_strength)
        if field == 4.0: return 1.281     # liver R1 in 1/sec (Changed from 1.285 on 06/08/2020)
        if field == 7.0: return 1.109     # liver R1 in 1/sec (Changed from 0.8350 on 06/08/2020)
        if field == 9.0: return 0.920     # per sec - liver R1 (https://doi.org/10.1007/s10334-021-00928-x)

    @property
    def R10S(self):
        field = math.floor(self.field_strength)
        if field == 4.0: return 0.631     # spleen R1 in 1/sec (Changed from 0.7458 on 23/07/2020)
        if field == 7.0: return 0.611     # spleen R1 in 1/sec (Changed from 0.6313 on 23/07/2020)
        if field == 9.0: return 0.600     # spleen R1 in 1/sec

    @property
    def t(self):
        return np.arange(0, self._tmax+self._tstep, self._tstep)

    @property
    def J(self):
        tend = self.tstart + self.tduration
        Jmax = self.dose/self.tduration
        t_inject = (self.t > self.tstart) & (self.t < tend)
        J = np.zeros(self.t.shape)
        J[np.argwhere(t_inject)] = Jmax
        return J

    @property
    # t = 0 is the start of data acquisition
    # sample points are defined at the center of the sampling interval
    def ts(self): 
        n = math.floor((self.t[-1] - self.t[0])/self.dt)
        ts = self.dt/2 + np.arange(n)*self.dt
        return ts

    def __init__(self):
        
        self.initialize_variables()

    def initialize_variables(self):

        self.initialize_spleen_variables()
        self.initialize_liver_variables()

    def initialize_spleen_variables(self):

        self.spleen_variables = np.array([self.K, self.P, self.VP, self.VE, self.S0spleen])

    def initialize_liver_variables(self):

        self.liver_variables = np.array([self.E, self.Th, self.S0liver, self.R10L])

    def export_variables(self):

        names = [
            "Gadoxetate extraction fraction", 
            "Hepatic plasma clearance rate",
            "Hepatocellular uptake rate", 
            "Biliary excretion rate", 
            "Glomerular Filtration Rate", 
            "Whole-body permeability-surface area",
            "Whole-body plasma volume", 
            "Whole-body extracelullar volume", 
            "Area under liver signal enhancement", 
        #    "Liver relaxation time", 
        ]
        symbols = [
            "E", 
            "Ktrans",
            "khe", 
            "kbh", 
            "GFR", 
            "PS",
            "Vp", 
            "Ve", 
            "AUC",
        #    "T10L", 
        ]
        units = [
            "%", 
            "mL/min/mL", 
            "mL/min/mL", 
            "mL/min/mL",
            "mL/min", 
            "mL/min",
            "mL", 
            "mL", 
            "(a.u.)*sec", 
        #    "sec", 
        ]
        values = [
            self.liver_variables[0] * 100,      # E
            self.Fp*self.liver_variables[0] * 60, # EF
            self.Fp*self.liver_variables[0]/(1-self.liver_variables[0]) * 60,   # FE/(1-E)
            self.vh / self.liver_variables[1] * 60,  # vh/Th
        #    self.liver_variables[4] * 60, # EF
        #    self.liver_variables[4] * self.VL * 60, # EF*VL
            (self.spleen_variables[0] - self.Fp*self.liver_variables[0]*self.VL) * 60, # K - EF*VL
            self.spleen_variables[1] * 60, 
            self.spleen_variables[2], 
            self.spleen_variables[3], 
            np.sum(self.liver_data-self.liver_data[0]) * self.dt, 
        #    1/self.liver_variables[3], 
        ]
        df1 = pd.DataFrame({"Variable": names})
        df2 = pd.DataFrame({"Symbol": symbols})
        df3 = pd.DataFrame({"Units": units})
        df4 = pd.DataFrame({"Value": values})
        return pd.concat([df1, df2, df3, df4], axis=1)

    def _signal(self, R1, S0):

        E = np.exp(-self.TR*R1)
        cFA = math.cos(self.FA*math.pi/180)
        return S0 * (1-E) / (1-cFA*E)

    def calculate_whole_body_concentrations(self):

        K = self.spleen_variables[0]
        P = self.spleen_variables[1]
        VP = self.spleen_variables[2]
        VE = self.spleen_variables[3]

        KP = (K + P)/VP
        KE = P/VE
        KB = K/VP

        self.cp, self.ce = propagate_2cxm(self.t, self.J/K, KP, KE, KB)

        self.cp *= 1000     # (mM)
        self.ce *= 1000     # (mM)

    def calculate_spleen_signal(self):

        S0 = self.spleen_variables[4]

        self.calculate_whole_body_concentrations()
        R1 = self.R10S + self.rp*self.veS*self.cp
        self.spleen_signal = self._signal(R1, S0)

    def calculate_liver_signal(self):

        E = self.liver_variables[0]
        Th = self.liver_variables[1]
        S0 = self.liver_variables[2]
        R10 = self.liver_variables[3]
    #    R10 = self.R10L
        
        #Te = self.veL*(1-E)/self.Fp

        X = self.rp * self.veL * (1-E)
        Y = self.rh * self.Fp * E
        R1 = R10 + X*self.cp + Y*Th*expconv(Th, self.t, self.cp)
        # R1 = R10 + X*expconv(Te, self.t, self.cp) + Y*Th*expconv(Th, self.t, self.cp)
        self.liver_signal = self._signal(R1, S0)

    def calculate_signals(self):

        self.calculate_spleen_signal()
        self.calculate_liver_signal()

    def _sample(self, ts, S):

        Ss = np.zeros(len(ts))
        for k, tk in enumerate(ts):
            tacq = (self.t > tk-self.dt/2) & (self.t < tk+self.dt/2)
            Ss[k] = np.average(S[np.argwhere(tacq)])
        return Ss 

    def _measure(self, ts, S):

        Ss = self._sample(ts, S)
        noise = np.random.normal(loc=0, scale=Ss[0]/self.SNR0, size=len(Ss))
        return Ss + noise

    def simulate_measurement(self):

        data = self._measure(self.ts, self.spleen_signal)
        self.set_spleen_data(self.ts, data)

        data = self._measure(self.ts, self.liver_signal)
        self.set_liver_data(self.ts, data)
    
    def set_spleen_data(self, ts, data):

        self._tmax = ts[-1] + self.dt/2
        self.spleen_sampling_times = ts
        self.spleen_data = data

    def set_liver_data(self, ts, data):

        self._tmax = ts[-1] + self.dt/2
        self.liver_sampling_times = ts
        self.liver_data = data
        
    def simulate_data(self):

        self.calculate_signals()
        self.simulate_measurement()

    def _fit_spleen_func(self, ts, *params):

        self.spleen_variables = np.array(params)
        self.calculate_spleen_signal()
        return self._sample(self.spleen_sampling_times, self.spleen_signal)

    def fit_spleen(self):

        self.initialize_spleen_variables()
        S0 = np.mean(self.spleen_data[0:4]) / self._signal(self.R10S, 1)
        self.spleen_variables[4] = S0
        variables = self.spleen_variables[:4]
        try:
            variables, _ = curve_fit(
                lambda t, K, P, VP, VE: self._fit_spleen_func(t, K, P, VP, VE, S0), 
                self.spleen_sampling_times, 
                self.spleen_data, 
                p0 = variables, 
                bounds = (0, np.inf),
            )
            self.spleen_variables[:4] = variables
        except:
            pass
        self.calculate_spleen_signal()

    def _spleen_signal(self, ts, cp):

        S0 = self.spleen_variables[4]
        R1 = self.R10S + self.rp*self.veS*cp
        return self._signal(R1, S0)

    def fit_spleen_direct(self):

        self.initialize_spleen_variables()
        S0 = np.mean(self.spleen_data[0:4]) / self._signal(self.R10S, 1)
        self.spleen_variables[4] = S0
        cp_sampled = []
        for i, ts in enumerate(self.spleen_sampling_times):
            cp, _ = curve_fit(
                self._spleen_signal, 
                [ts], 
                [self.spleen_data[i]], 
                p0 = [0.1], 
                bounds = (0, np.inf),
            )
            cp_sampled.append(cp[0])
        self.cp = np.interp(self.t, self.spleen_sampling_times, cp_sampled)
        R1 = self.R10S + self.rp*self.veS*self.cp
        self.spleen_signal = self._signal(R1, S0)

    def fit_direct(self):

        self.fit_spleen_direct()
        self.fit_liver()

    def _fit_liver_func(self, ts, *params):

        self.liver_variables = np.array(params)
        self.calculate_liver_signal()
        return self._sample(self.liver_sampling_times, self.liver_signal)

    def fit_liver(self):

        self.initialize_liver_variables()
        R10 = self.liver_variables[3]
        S0 = np.mean(self.liver_data[:4]) / self._signal(R10, 1)
        self.liver_variables[2] = S0
        try:
            variables = self.liver_variables[:2]
            variables, _ = curve_fit(
                lambda t, E, Th: self._fit_liver_func(t, E, Th, S0, R10), 
                self.liver_sampling_times, 
                self.liver_data, 
                p0 = variables, 
                bounds = (0, [1, np.inf]), 
            )
            self.liver_variables[:2] = variables
        except:
            pass
        self.calculate_liver_signal()

    def fit(self):

        self.fit_spleen()
        self.fit_liver()

    def fit_standard(self):

        self.initialize_spleen_variables()
        S0 = np.mean(self.spleen_data[0:4]) / self._signal(self.R10S, 1)
        self.spleen_variables[4] = S0
        self.calculate_spleen_signal()
        self.fit_liver()

    def _fit_joint_func(self, ts, *params):

        spleen_data = self._fit_spleen_func(self.spleen_sampling_times, *params[:5])
        liver_data = self._fit_liver_func(self.liver_sampling_times, *params[5:])
        return np.concatenate((spleen_data, liver_data))

    def fit_joint(self):

        self.fit() # initialize
    #    self.initialize_spleen_variables()
    #    S0_estimate = self.spleen_data[0] / self._signal(self.R10S, 1)
    #    self.spleen_variables[4] = S0_estimate
    #    self.initialize_liver_variables()
    #    S0_estimate = self.liver_data[0] / self._signal(self.R10L, 1)
    #    self.liver_variables[2] = S0_estimate
        S0s = self.spleen_variables[4]
        S0l = self.liver_variables[2]
        R10L = self.liver_variables[3]
        variables = np.concatenate((self.spleen_variables[:4], self.liver_variables[:2]))
        variables, _ = curve_fit(
            lambda t, K, P, VP, VE, E, Th: self._fit_joint_func(t, K, P, VP, VE, S0s, E, Th, S0l, R10L), 
            np.concatenate((self.spleen_sampling_times, self.liver_sampling_times)), 
            np.concatenate((self.spleen_data, self.liver_data)), 
            p0 = variables, 
            bounds = (0, [np.inf, np.inf, np.inf, np.inf, 1, np.inf]), 
        )
        self.spleen_variables[:4] = variables[:4]
        self.liver_variables[:2] = variables[4:]

    def plot_whole_body_concentrations(self):

        self.calculate_whole_body_concentrations()

        plt.plot(self.t, self.cp, color='red')
        plt.plot(self.t, self.ce, color='green')
        plt.show()

    def plot_signals(self):

        self.calculate_signals()

        plt.plot(self.t, self.spleen_signal, color='red')
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.show()

    def plot_data(self):

        self.simulate_data()

        plt.plot(self.t, self.spleen_signal, color='red')
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.plot(self.spleen_sampling_times, self.spleen_data, marker="x", color='red', linewidth=0)
        plt.plot(self.liver_sampling_times, self.liver_data, marker="x", color='blue', linewidth=0)
        plt.show()

    def plot_spleen_fit(self):

        self.simulate_data()
        plt.plot(self.t, self.spleen_signal, color='blue')
        plt.plot(self.spleen_sampling_times, self.spleen_data, marker="x", color='red', linewidth=0)
        self.fit_spleen()
        plt.plot(self.t, self.spleen_signal, color='red')
        plt.show()

    def plot_liver_fit(self):

        self.simulate_data()
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.plot(self.liver_sampling_times, self.liver_data, marker="x", color='red', linewidth=0)
        self.fit_spleen()
        self.fit_liver()
        plt.plot(self.t, self.liver_signal, color='red')
        plt.show()

    def plot_fit(self):

        self.simulate_data()
        plt.plot(self.t, self.spleen_signal, color='blue')
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.plot(self.spleen_sampling_times, self.spleen_data, marker="x", color='red', linewidth=0)
        plt.plot(self.liver_sampling_times, self.liver_data, marker="x", color='green', linewidth=0)
        self.fit()
        plt.plot(self.t, self.spleen_signal, color='red')
        plt.plot(self.t, self.liver_signal, color='green')
        plt.show()
