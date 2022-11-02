"""TRISTAN model for dynamic gadoxetate-enhanced MRI in rats.

A class TristanRat which defines the tracer kinetic modelling equations 
and associated default variables used in the preclinical dynamic gadoxetate-
enhanced MR imaging work of the IMI-TRISTAN WP2 project.

"""
# imports
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from typing import Dict, Tuple


# Functions
def expconv(T: float, 
            time: np.ndarray, 
            a: np.ndarray
) -> np.ndarray:
    """Convolves a 1D-array with a normalised exponential.

    Uses an efficient and accurate numerical formula to calculate the convolution,
    as detailed in the appendix of Flouri et al., Magn Reson Med, 76 (2016), pp. 998-1006.

    Args:
        T: The characteristic time of the the exponential function.
            time and T must be in the same units.
        time: The time points where the values of ca are defined.
            These do not have to to be equally spaced.
        a: The 1D array to be convolved.

    Returns: 
        The convolved array.
        this is the same shape as ca.
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


def propagate_2cxm(t: np.ndarray,
                    ca: np.ndarray,
                    KP: float,
                    KE: float,
                    KB: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculates the propagators for the individual compartments in the 2CXM.
    
    For details and notations see appendix of 
    Sourbron et al. Magn Reson Med 62:672–681 (2009).

    Args:
        t: time points (sec) where the input function is defined
        ca: input function (mmol/mL)
        KP: inverse plasma MTT (sec) = VP/(FP+PS)
        KE: inverse extracellular MTT (sec) = VE/PS
        KB: inverse blood MTT (sec) = VP/FP

    Returns:
        A tuple (cp, ce), where cp is the concentration in the plasma compartment, 
        and ce is the concentration in the extracellular compartment.
        Both are in mmol/mL.
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
    """The TRISTAN-rat model.
    
    A class describing the tracer kinetic model used in the IMI-TRISTAN WP2 
    preclinical studies. Default values for tracer kinetic dosing, MRI scanning
    parameters and rat physiology are provided, with additional references cited
    in block or inline comments beside variables.  

    Attributes:
        _tstep: an integer count of the MRI sequence internal time resolution (sec).
        _tmax: an integer count of the MRI sequence total acquisition time (sec).
        tstart: integer value for the start time (sec) of contrast agent injection.
        tduration: integer value for the duration (sec) of contrast agent injection.
        dose: integer value for the dose (mmol) of contrast agent administered.
        field_strength: integer value for the MRI field strength (T) used.
        dt: integer count of the MRI sequence sampling duration (sec).
        FA: integer value for the MRI sequence flip angle (degrees).
        TR: integer count of the MRI sequence repetition time (sec).
        SNR0: integer value for the signal-to-noise ratio (SNR) at baseline.
            Simulation only.
        veS: integer value for the rat spleen extracellular volume (mL/mL).
        S0spleen: integer value for the rat spleen baseline signal (a.u.).
            Randomly chosen.
        Fb: integer describing the rate of blood flow (mL/sec/mL) in rat liver.
        E: integer value for the rat gadoxetate extraction fraction (%).
        veL: integer value for the rat liver extracellular volume (mL/mL).
        Th: integer value for the rat hepatocellular mean transit time (sec).
        vh: integer value for the rat hepatocellular volume fraction (mL/mL).
        S0liver: integer value for the rat liver baseline signal (a.u.).
            Randomly chosen.
        Hct: integer value for Hematocrit in rat.
        VL: integer value for the rat liver volume (mL).
        GFR: integer describing the rate of glomerular filtration in rat (mL/sec).
        P: integer value of permeability-surface for the whole rat body, with 
            permeability for gadoxetate (mL/sec).
        VB: integer value for the rat whole body blood volume (mL).
        VE: integer value for the rat whole body extracellular volume (mL).
    """
    
    # Internal time resolution & acquisition time
    _tstep = 0.5
    _tmax = 40*60.0

    # Experimental variables
    tstart = 4.6*60
    tduration = 30
    dose = 0.0075
    field_strength = 4.7
    dt = 60
    FA = 20
    TR = 5.0/1000
    SNR0 = 20
    
    # Spleen parameters
    veS = 0.314
    S0spleen = 250

    # Liver parameters
    Fb = 2.27/60        # https://doi.org/10.1021/acs.molpharmaceut.1c00206
                        # (Changed from 3.61/60 on 07/03/2022)
	                    # From Brown the cardiac output of rats is 110.4 mL/min (table 3-1) ~ 6.62L/h
	                    # From table 3-4, sum of hepatic artery and portal vein blood flow is 17.4% of total cardiac output ~ 1.152 L/h
	                    # Mass of liver is 9.15g, with density of 1.08 kg/L, therefore ~8.47mL
                        #  9.18g refers to the whole liver, i.e. intracellular tissue + extracellular space + blood
 	                    # Dividing 1.152L/h for 8.47mL we obtain ~2.27 mL/h/mL liver
	                    # Calculation done with values in Table S2 of our article lead to the same results

    E = 0.4
    veL = 0.230
    Th = 20*60
    vh = 0.722
    S0liver = 200

    # Whole body parameters
    Hct = 0.418         # Cremer et al, J Cereb Blood Flow Metab 3, 254-256 (1983)
    VL = 8.47           # Scotcher et al 2021, DOI: 10.1021/acs.molpharmaceut.1c00206
                        # Supplementary material, Table S2 
    GFR = 0.023         # https://doi.org/10.1152/ajprenal.1985.248.5.F734
    P = 0.172           # Estimated from rat repro study data using PBPK model
                        # Table 3 in Scotcher et al 2021
                        # DOI: 10.1021/acs.molpharmaceut.1c00206
    VB = 15.8           # 0.06 X BW + 0.77, Assuming body weight (BW) = 250 g
                        # Lee and Blaufox. Blood volume in the rat. 
                        # J Nucl Med. 1985 Jan;26(1):72-6.
    VE = 30             # All tissues, including liver.
                        # Derived from Supplementary material, Table S2
                        # Scotcher et al 2021
                        # DOI: 10.1021/acs.molpharmaceut.1c00206

    @property
    def K(self):
        """"Total excretion rate (mL/min/mL)."""
        return self.GFR + self.Ktrans * self.VL

    @property
    def Fp(self):
        """Rat liver plasma flow (mL/min/mL)."""
        return (1-self.Hct) * self.Fb

    @property
    def VP(self):
        """Rat whole body plasma volume (mL)."""
        return (1-self.Hct) * self.VB 

    @property
    def Ktrans(self):
        """Rat hepatic plasma clearance rate (mL/min/mL)."""
        return self.E * self.Fp 

    @property
    def rp(self):
        """Relaxivity of rat blood (Hz/mM) depending on MRI field strength used."""
        field = math.floor(self.field_strength)
        if field == 4.0: return 6.4     # relaxivity of blood in Hz/mM ;assume spleen relaxivity is the same
        if field == 7.0: return 6.2     # relaxivity of blood in Hz/mM ;assume spleen relaxivity is the same
        if field == 9.0: return 6.1     # relaxivity of blood in Hz/mM ;assume spleen relaxivity is the same

    @property
    def rh(self):
        """Relaxivity of rat hepatocytes (Hz/mM) depending on MRI field strength used."""
        field = math.floor(self.field_strength)
        if field == 4.0: return 7.6     # relaxivity of hepatocytes in Hz/mM
        if field == 7.0: return 6.0     # relaxivity of hepatocytes in Hz/mM
        if field == 9.0: return 6.1     # relaxivity of hepatocytes in Hz/mM

    @property
    def R10L(self):
        """Precontrast rat liver relaxation rate (1/sec) depending on MRI field strength used."""
        field = math.floor(self.field_strength)
        if field == 4.0: return 1.281     # liver R1 in 1/sec (Changed from 1.285 on 06/08/2020)
        if field == 7.0: return 1.109     # liver R1 in 1/sec (Changed from 0.8350 on 06/08/2020)
        if field == 9.0: return 0.920     # per sec - liver R1 (https://doi.org/10.1007/s10334-021-00928-x)

    @property
    def R10S(self):
        """Precontrast rat spleen relaxation rate (1/sec) depending on MRI field strength used."""
        field = math.floor(self.field_strength)
        if field == 4.0: return 0.631     # spleen R1 in 1/sec (Changed from 0.7458 on 23/07/2020)
        if field == 7.0: return 0.611     # spleen R1 in 1/sec (Changed from 0.6313 on 23/07/2020)
        if field == 9.0: return 0.600     # spleen R1 in 1/sec

    @property
    def t(self):
        """1D time series array spanning the total length of the MRI acquisition."""
        return np.arange(0, self._tmax+self._tstep, self._tstep)

    @property
    def J(self):
        """Gadoxetate influx (mmol/sec)."""
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
        """Sampling time (sec).""" 
        n = math.floor((self.t[-1] - self.t[0])/self.dt)
        ts = self.dt/2 + np.arange(n)*self.dt
        return ts

    def __init__(self):
        """Initializes class variables."""
        
        self.initialize_variables()

    def initialize_variables(self):
        """Initializes rat spleen and liver variables."""

        self.initialize_spleen_variables()
        self.initialize_liver_variables()

    def initialize_spleen_variables(self):
        """Initializes rat spleen variables into array format."""

        self.spleen_variables = np.array([self.K, self.P, self.VP, self.VE, self.S0spleen])

    def initialize_liver_variables(self):
        """Initializes rat liver variables into array format."""

        self.liver_variables = np.array([self.E, self.Th, self.S0liver, self.R10L])

    def export_variables(self
    ) -> pd.DataFrame:
        """Exports estimated parameter variables.
        
        A function to export all estimated parameter variables and return them
        in DataFrame format after fitting the data with the TRISTAN-rat tracer 
        kinetic model.
        """
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
            (self.spleen_variables[0] - self.Fp*self.liver_variables[0]*self.VL) * 60, # K - EF*VL
            self.spleen_variables[1] * 60, 
            self.spleen_variables[2], 
            self.spleen_variables[3], 
            np.sum(self.liver_data-self.liver_data[0]) * self.dt,
        ]
        df1 = pd.DataFrame({"Variable": names})
        df2 = pd.DataFrame({"Symbol": symbols})
        df3 = pd.DataFrame({"Units": units})
        df4 = pd.DataFrame({"Value": values})
        return pd.concat([df1, df2, df3, df4], axis=1)

    def _signal(self, 
                R1: np.ndarray,
                S0: np.float64
    ) -> np.ndarray:
        """Calculates MRI signal intensity."""
        E = np.exp(-self.TR*R1)
        cFA = math.cos(self.FA*math.pi/180)
        return S0 * (1-E) / (1-cFA*E)

    def calculate_whole_body_concentrations(self
    ) -> np.ndarray:
        """Simulates whole body concentrations for the 2cxm.
        
        Returns:
            An array of the plasma compartment concentration, cp (mmol/mL).    
        """
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

        return self.cp

    def calculate_spleen_signal(self
    ) -> None:
        """Calculates MRI signal intensity in rat spleen."""
        S0 = self.spleen_variables[4]

        self.calculate_whole_body_concentrations()
        R1 = self.R10S + self.rp*self.veS*self.cp
        self.spleen_signal = self._signal(R1, S0)

    def calculate_liver_signal(self
    ) -> None:
        """Calculates MRI signal intensity in rat liver."""
        E = self.liver_variables[0]
        Th = self.liver_variables[1]
        S0 = self.liver_variables[2]
        R10 = self.liver_variables[3]

        X = self.rp * self.veL * (1-E)
        Y = self.rh * self.Fp * E
        R1 = R10 + X*self.cp + Y*Th*expconv(Th, self.t, self.cp)
        self.liver_signal = self._signal(R1, S0)

    def calculate_signals(self
    ) -> None:
        self.calculate_spleen_signal()
        self.calculate_liver_signal()

    def _sample(self, 
                ts: np.ndarray, 
                S: np.ndarray
    ) -> np.ndarray:
        """Sample a pseudo-continuous MRI signal at given sampling times."""
        Ss = np.zeros(len(ts))
        for k, tk in enumerate(ts):
            tacq = (self.t > tk-self.dt/2) & (self.t < tk+self.dt/2)
            Ss[k] = np.average(S[np.argwhere(tacq)])
        return Ss 

    def _measure(self,
                ts: np.ndarray,
                S: np.ndarray
    ) -> np.ndarray:
        """Returns noisy sampled MRI signal."""
        Ss = self._sample(ts, S)
        noise = np.random.normal(loc=0, scale=Ss[0]/self.SNR0, size=len(Ss))
        return Ss + noise

    def simulate_measurement(self):
        """Simulates MRI measurement."""
        data = self._measure(self.ts, self.spleen_signal)
        self.set_spleen_data(self.ts, data)

        data = self._measure(self.ts, self.liver_signal)
        self.set_liver_data(self.ts, data)
    
    def set_spleen_data(self,
                        ts: np.ndarray,
                        data: np.ndarray
    ) -> None:
        """Assigns rat spleen data to be fitted."""
        self._tmax = ts[-1] + self.dt/2
        self.spleen_sampling_times = ts
        self.spleen_data = data

    def set_liver_data(self,
                        ts: np.ndarray,
                        data: np.ndarray
    ) -> None:
        """Assigns rat liver data to be fitted."""
        self._tmax = ts[-1] + self.dt/2
        self.liver_sampling_times = ts
        self.liver_data = data
        
    def simulate_data(self):
        """Simulates rat liver and spleen data/signals."""
        self.calculate_signals()
        self.simulate_measurement()

    def _fit_spleen_func(self,
                        ts: np.ndarray,
                        *params: np.ndarray
    ) -> np.ndarray:
        """Fits sample rat spleen data."""
        self.spleen_variables = np.array(params)
        self.calculate_spleen_signal()
        return self._sample(self.spleen_sampling_times, self.spleen_signal)

    def fit_spleen(self):
        """Fits rat spleen data."""
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

    def _spleen_signal(self,
                        ts: np.ndarray,
                        cp: np.ndarray
    ) -> np.ndarray:
        """Returns rat spleen signal."""
        S0 = self.spleen_variables[4]
        R1 = self.R10S + self.rp*self.veS*cp
        return self._signal(R1, S0)

    def fit_spleen_direct(self):
        """Directly calculates and fits rat spleen data."""
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

    def _fit_liver_func(self,
                        ts: np.ndarray,
                        *params: np.ndarray
    ) -> np.ndarray:
        """Fits sample rat liver data."""
        self.liver_variables = np.array(params)
        self.calculate_liver_signal()
        return self._sample(self.liver_sampling_times, self.liver_signal)

    def fit_liver(self):
        """Fits rat liver data."""
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
        """Fits data using rat spleen and liver data."""
        self.fit_spleen()
        self.fit_liver()

    def fit_standard(self):
        """Fits data using a standardised cp(t).
        
        Considering the difficulty in reliably measuring cp(t) in rats, this 
        functions performs the tracer kinetic modelling by implementing a 
        standardised cp(t) derived from a simplified 2cxm of the rat 
        circulation, as defined by the propagate_2cxm() and 
        calculate_whole_body_concentrations() functions.
        """
        self.initialize_spleen_variables()
        S0 = np.mean(self.spleen_data[0:4]) / self._signal(self.R10S, 1)
        self.spleen_variables[4] = S0
        self.calculate_spleen_signal()
        self.fit_liver()

    def _fit_joint_func(self,
                        ts: np.ndarray,
                        *params: np.ndarray
    ) -> np.ndarray:
        """Returns fitted sample spleen and liver data combined."""
        spleen_data = self._fit_spleen_func(self.spleen_sampling_times, *params[:5])
        liver_data = self._fit_liver_func(self.liver_sampling_times, *params[5:])
        return np.concatenate((spleen_data, liver_data))

    def fit_joint(self):
        """Jointly fits spleen and liver data."""
        self.fit() # initialize
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
        """Plots simulated whole body concentrations for the 2cxm."""
        self.calculate_whole_body_concentrations()

        plt.plot(self.t, self.cp, color='red')
        plt.plot(self.t, self.ce, color='green')
        plt.show()

    def plot_signals(self):
        """Plots calculated rat spleen and liver signals."""
        self.calculate_signals()

        plt.plot(self.t, self.spleen_signal, color='red')
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.show()

    def plot_data(self):
        """Plots sample rat spleen and liver data."""
        self.simulate_data()

        plt.plot(self.t, self.spleen_signal, color='red')
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.plot(self.spleen_sampling_times, self.spleen_data, marker="x", color='red', linewidth=0)
        plt.plot(self.liver_sampling_times, self.liver_data, marker="x", color='blue', linewidth=0)
        plt.show()

    def plot_spleen_fit(self):
        """Plots fitted sample spleen data."""
        self.simulate_data()
        plt.plot(self.t, self.spleen_signal, color='blue')
        plt.plot(self.spleen_sampling_times, self.spleen_data, marker="x", color='red', linewidth=0)
        self.fit_spleen()
        plt.plot(self.t, self.spleen_signal, color='red')
        plt.show()

    def plot_liver_fit(self):
        """Plots fitted sample liver data."""
        self.simulate_data()
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.plot(self.liver_sampling_times, self.liver_data, marker="x", color='red', linewidth=0)
        self.fit_spleen()
        self.fit_liver()
        plt.plot(self.t, self.liver_signal, color='red')
        plt.show()

    def plot_fit(self):
        """Plots fitted sample spleen and liver data together."""
        self.simulate_data()
        plt.plot(self.t, self.spleen_signal, color='blue')
        plt.plot(self.t, self.liver_signal, color='blue')
        plt.plot(self.spleen_sampling_times, self.spleen_data, marker="x", color='red', linewidth=0)
        plt.plot(self.liver_sampling_times, self.liver_data, marker="x", color='green', linewidth=0)
        self.fit()
        plt.plot(self.t, self.spleen_signal, color='red')
        plt.plot(self.t, self.liver_signal, color='green')
        plt.show()
