# TRISTAN-RAT
This archive contains an assortment of scripts used for analysing dynamic 
gadoxetate-enhanced magnetic resonance imaging (MRI) signal data in rat
liver as part of the work conducted by 
[work package 2 (WP2) of the IMI-TRISTAN TRISTAN EU project](https://www.imi-tristan.eu/liver).

These scripts accompany the manuscript "Use of physiologically based 
pharmacokinetic and tracer kinetic models for prospective prediction and 
quantification of hepatic transporters drug-drug interaction via in vivo imaging 
in rats" (hereafter referred to as the 'Six Test Compounds' study).

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

## Project setup and installation
1. Open a Git Bash in an interface of your choice and navigate to an empty project directory by inputting the following command:

        cd <path_to_project_name>

3. Clone the project from GitHub:

        git clone https://github.com/<account_name>/<project_name>.git .

Alternatively, download the whole archive locally as a zip file. 
[More details on cloning a repository may be found here:](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)

4. From the project root directory, run the following command to create a separate virtual environment:

        conda create --name <project_name> python=3.7

5. Activate the virtual environment:
    
        conda activate <project_name>

6. Install required packages:

        pip install -r requirements.txt

## Running the scripts
The MRI signal data acquired in this study are provided in csv format and can 
be found in the `01_signals` folder contained within the project `data` 
directory.

All necessary source code for reproducing the results presented in the 
manuscript are located within the `src` directory of this project.

Results may be found within the `results` directory and are generated from the
code by running the following command from within the project's root directory:

    python src/main.py --study SixTestCompounds

A detailed metadata further describing the contents of each directory may be
found within respective folders.

An additional Jupyter notebook [TristanRat_examples.ipynb] has been provided for 
interacting with some of the functions and methoods of the TristanRat class 
described in the `src/rat.py`. To try some of these examples, run the command:

    jupter notebook

This will start a server and open the Jupyter interface within the operating 
system's default browser. In the contents, navigate to the notebook in
`notebooks/TristanRat_examples.ipynb` and select to open. 

Each example in the notebook is provided in a separate cell. Pressing 
`Shift + Enter` on the keyboard executes the code within a cell and prints the 
output below. Feel free to interact with the examples by adjusting the code and
exploring the different outputs obtained as a result.

## Contributors

|Name     |  GitHub Handle   | Email    |
|---------|------------------|----------|
|[Ebony R. Gunwhy](https://github.com/EbonyGunwhy)  | @EbonyGunwhy    | e.gunwhy@sheffield.ac.uk   |
|[Steven Sourbron](https://github.com/plaresmedima) | @plaresmedima   | s.sourbron@sheffield.ac.uk |
