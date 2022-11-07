# TRISTAN-RAT
This archive contains an assortment of scripts used for analysing dynamic 
gadoxetate-enhanced magnetic resonance imaging (MRI) signal data in rat
liver as part of the work conducted by 
[work package 2 (WP2) of the IMI-TRISTAN TRISTAN EU project](https://www.imi-tristan.eu/liver).

These scripts accompany the manuscript "Use of physiologically based 
pharmacokinetic and tracer kinetic models for prospective prediction and 
quantification of hepatic transporters drug-drug interaction via in vivo imaging 
in rats" by Nicola Melillo, Daniel Scotcher, J. Gerry Kenna, Claudia Green, 
Catherine D. G. Hines, Iina Laitinen, Paul D. Hockings, Kayode Ogungbenro, 
Ebony R. Gunwhy, Steven Sourbron, John C. Waterton, Gunnar Schuetz, and 
Aleksandra Galetin (hereafter referred to as the 'Six Test Compounds' study).

This paper has been submitted for publication in the open access journal 
[Pharmaceutics from MDPI](https://www.mdpi.com/journal/pharmaceutics).

# Getting started
## Project environment and dependencies
A working Python environment is required to run the scripts contained within 
this repository. All required Python dependencies are specified within the 
`requirements.txt` file located within the root directory of this project.

In this project, `conda` virtual environments were used to manage the project 
dependencies in isolation. Anaconda may be installed within the user's directory 
without causing conflicts with a system's Python installation, therefore it is 
recommended to set up a working environment by downloading the `conda` package 
manager from [Anaconda's Python distribution](https://www.anaconda.com/download/).

**Note**
<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #8a6d3b;; background-color: #fcf8e3; border-color: #faebcc;">
The following steps assume that Anaconda has already been installed and that 
commands are run from a Windows OS. If replicating from a different OS, please 
adapt commands to the appropriate related invocation [(Some examples here)](https://kinsta.com/blog/python-commands/).
</div>

## Project setup and installation
1. Open Git Bash in an interface of your choice and navigate to an empty project directory by inputting the following command:

        cd <path_to_project_name>

3. Clone the project from GitHub:

        git clone https://github.com/<account_name>/<project_name>.git .

   Alternatively, download the whole archive locally as a zip file. 
   [More details on cloning a repository may be found here.](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)

4. From the project root directory, run the following command to create a separate virtual environment:

        conda create --name <project_name> python=3.7

5. Activate the virtual environment:
    
        conda activate <project_name>

6. Install required packages:

        pip install -r requirements.txt

## Running the scripts
### `data` directory
The MRI signal data acquired in this study are provided in csv format and can 
be found in the `data/SixTestCompounds/01_signals` folder.

Each csv file contains the liver and spleen region of interest (ROI) time 
curve data for a single rat at a specific site and day, after administration of
the test compound of interest.

The corresponding filename for each file is formatted as a string containing 
study descriptors (metadata) separated by underscores, i.e.,
filename = `compound_site_RatNumber_dayNumber_dataType`,
e.g., `Asunaprevir_E_Rat2_2_Signals`.

### `src` directory
All necessary source code for reproducing the results presented in the 
manuscript are located within the `src` directory of this project. This
consists of 4 modules:
```
data.py
models.py
effect_sizes.py
plots.py
```
1 class:
```
rat.py
```
and 1 script used to generate all outputs from the modules:
```
main.py
```

### `results` directory
Results may be found within the `results` directory and are generated from the
code by running the main script using the following command from within the 
project's root directory:

    python src/main.py --study SixTestCompounds

This will create a results folder with the following structure:
```
results/
|
|---- StudyName (e.g., 'SixTestCompounds')/
|    |
|    |---- DateOfAnalysis YYYY-MM-DD (e.g., '2022-11-04')/
|    |     |
|    |     |---- 01_model_outputs/
|    |     |    |---- figures/
|    |     |    |       |---- per_drug/
|    |     |    |       |---- per_rat/
|    |     |    |---- relaxation_rates_and_signals
|    |     |    |---- all_parameters.csv
|    |     |---- 02_effect_sizes/
|    |     |    |---- figures/
|    |     |    |---- effect_sizes.csv
|    |     |    |---- fit_errors.txt
```
As the tracer kinetic model used in this study produces estimated parameter
variables, modelling outputs may vary slightly between different iterations.
Therefore, upon each execution of the code, a top-level directory named after 
the date the analysis was conducted is created for storing the results from
that particular execution. For reference, the results and figures shown in this 
manuscript were created using the outputs contained in 
`results/SixTestCompounds/2022-09-01`.

`01_model_ouputs` contains all outputs generated as result of the tracer kinetic
model fitting. Within this, plotted signal time curves for each acquistion per 
rat may be found in `figures/per_rat`, while average deltaR1 curves per drug may
be found in `figures/per_drug`. All estimated parameter variables from the fitting
are stored in `all_parameters.csv`. The folder `relaxation_rates_and_signals` 
contains the fitted MRI signal data in a similar format to the original csv data 
in `data/01_signals`. Additional columns have been included showing the MRI signals 
converted to R1.

`02_effect_sizes` contains results of the statistical analysis performed on the 
tracer kinetic modelling output. This is summarised in tabular format in 
`effect_sizes.csv`, while graphical distributions are provided within the `figures/`
folder. `fit_errors.txt` contains an ID list for any computational fitting errors
found during quality control of the tracer kinetic modelling.

### `notebooks` directory
An additional Jupyter notebook `TristanRat_examples.ipynb` has been provided for 
interacting with some of the functions and methods of the TristanRat class 
described in the `src/rat.py` script. To try some of these examples, run the command:

    jupyter notebook

This will start a server and open the Jupyter interface within the operating 
system's default browser. In the contents, navigate to the notebook 
`notebooks/TristanRat_examples.ipynb` and select to open. 

Each example in the notebook is provided in a separate cell. Pressing 
`Shift + Enter` on the keyboard executes the code within a cell and prints the 
output below. Feel free to interact with the examples by adjusting the code and
exploring the different outputs obtained as a result.

## Contributors

|Name     |  GitHub Handle   | Email    |
|---------|------------------|----------|
|[Ebony R. Gunwhy](https://github.com/EbonyGunwhy)  | [@EbonyGunwhy](https://github.com/EbonyGunwhy)     | e.gunwhy@sheffield.ac.uk   |
|[Steven Sourbron](https://github.com/plaresmedima) | [@plaresmedima](https://github.com/plaresmedima)   | s.sourbron@sheffield.ac.uk |

## Funding
The research leading to these results received funding from the [Innovative Medicines 
Initiative](https://www.imi.europa.eu/) 2 Joint Undertaking under grant agreement No 
116106. This Joint Undertaking receives support from the [European Union’s Horizon 2020](https://research-and-innovation.ec.europa.eu/funding/funding-opportunities/funding-programmes-and-open-calls/horizon-2020_en) research and innovation programme 
and [EFPIA](https://www.efpia.eu/).
