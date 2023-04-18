# TRISTAN-RAT

This archive contains an assortment of scripts used for analysing dynamic 
gadoxetate-enhanced magnetic resonance imaging (MRI) signal data in rat
liver as part of the work conducted by 
[work package 2 (WP2) of the IMI-TRISTAN TRISTAN EU project](https://www.imi-tristan.eu/liver).

These scripts accompany the following studies:

* `Six Test Compounds` study (1)
* `Reproducibility` study (manuscript in preparation)


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


> **_PLEASE NOTE:_** The following steps assume that Anaconda has already been installed and that 
>       commands are run from a Windows OS. If replicating from a different OS, please 
>       adapt commands to the appropriate related invocation ([Some examples here](https://kinsta.com/blog/python-commands/)).


## Project setup and installation

1. Open Command Prompt in an interface of your choice, create an empty project
directory and navigate into it by inputting the following commands:

        mkdir TRISTAN-rat
        cd TRISTAN-rat

3. Clone the TRISTAN-rat repository from GitHub into the newly created
directory:

        git clone https://github.com/QIB-Sheffield/TRISTAN-rat .

   Alternatively, download the whole archive locally as a zip file. 
   [More details on cloning a repository may be found here.](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)

4. From the project root directory, run the following command to create a separate virtual environment:

        conda create --name <environment_name> python=3.7

5. Activate the virtual environment:
    
        conda activate <environment_name>

6. Install required Python dependencies:

        pip install -r requirements.txt


## Running the scripts

### `src` directory

All necessary source code for reproducing the results from each 
study are located within the `src` directory of this project. In
its current form, the majority of scripts in this directory are
intended to be used specifically for reproducing these results
along with the data accompanying these studies. If the user wishes
to analyse their own DCE-MRI data with the TRISTAN-rat model and
define their own outputs, please see section `Interacting with the
scripts` for guidelines on how to do this.

### `Six Test Compounds` study
Please download the `data.zip` and `results.zip` folders from
https://zenodo.org/record/7506968#.Y7cNCdXP1PY
(DOI: 10.5281/zenodo.7506968) and extract the contents. Please
make sure the extracted contents are stored within the project's
root directory.

The results for this study are generated using the
`SixTestCompounds.py` script located in the `src` directory.
From within the project's root directory, run the script using
the following command:

    python src/SixTestCompounds.py --study SixTestCompounds

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
that particular execution. For reference, the results and figures presented in
the accompanying manuscript to the `Six Test Compounds` study were created using the
outputs contained in `results/SixTestCompounds/2022-09-01`.

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

### `Reproducibility` study

Please download the `data.zip` and `results.zip` folders from
https://zenodo.org/record/7838397#.ZD3jlvzMJPY
(DOI: 10.5281/zenodo.7838397) and extract the contents. Please
make sure the extracted contents are stored within the project's
root directory.

The results for this study are generated using the
`Reproducibility.py` script located in the `src` directory.
From within the project's root directory, run the script using
the following command:

    python src/Reproducibility.py --study Reproducibility

This will create a results folder with the following structure:
```
results/
|
|---- StudyName (e.g., 'Reproducibility')/
|    |
|    |---- DateOfAnalysis YYYY-MM-DD (e.g., '2023-04-14')/
|    |     |
|    |     |---- 01_model_outputs/
|    |     |    |---- figures/
|    |     |    |       |---- per_substudy/
|    |     |    |       |---- per_rat/
|    |     |    |---- relaxation_rates_and_signals
|    |     |    |---- all_parameters.csv
|    |     |    |---- fit_errors.txt
|    |     |---- 02_analyses/
|    |     |    |---- figures/
|    |     |    |---- repeatability/
|    |     |    |---- reproducibility
|    |     |    |---- benchmarks.csv
|    |     |    |---- mixed_anova.csv
```

For reference, the results and figures presented in the accompanying manuscript
were created using the outputs contained in `results/SixTestCompounds/2023-04-17`.

## Interacting with the scripts

The script `rat.py` contains the TRISTAN-rat model used for dynamic
gadoxetate-enhanced MRI in rats. Within this script, all tracer
kinetic modelling equations and associated default variables used
in the preclinical dynamic gadoxetate-enhanced MR imaging work of
the IMI-TRISTAN WP2 project are defined. The interested user may employ
the funtions in this script to develop their own analysis workflow using
their own datasets, by creating similar scripts to those shown in the
TRISTAN-rat `src` directory.

### Input data requirements

The models outlined in `rat.py` may be used to fit liver and spleen region
of interest (ROI) time curve data provided by the user. This data must be
provided in a csv format to match the data contained inside the
`data/SixTestCompounds/01_signals` subfolder. Here, each csv file contains the
liver and spleen ROI time curve data for a single rat at a specific site and
day, after administration of the test compound of interest.

### `notebooks` directory

An additional Jupyter notebook `TristanRat_examples.ipynb` has been provided
for interacting with some of the functions and methods of the TristanRat class 
described in the `src/rat.py` script.

> **_PLEASE NOTE:_** The functions developed in this notebook are *examples* only
> to demonstate how the functions and methods defined in `rat.py` may be used for
> simulations and for fitting DCE-MRI datasets. It is up to the user to build their
> own functions in a similar manner shown in the examples to employ the TRISTAN-rat
> model in their own analysis workflow.

To try some of these examples, run the command:

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


## *References*
*1. Melillo N, Scotcher D, Kenna JG, Green C, Hines CDG, Laitinen I, et al.*
*Use of In Vivo Imaging and Physiologically-Based Kinetic Modelling to Predict*
*Hepatic Transporter Mediated Drug–Drug Interactions in Rats. Pharmaceutics*
*2023;15(3):896. doi: 10.3390/pharmaceutics15030896*
