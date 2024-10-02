# Vasp Balsam workflow for magnetic materials

> **Note:** This documentation is written assuming the you are running on Polaris. This workflow is adaptable to other systems, but particulars such as setting up the balsam database, creating the conda environment, and the particulars of commands may need to be changed. In addition, the python scripts may need to be adapted to the new environment, such as changing Balsam app script preambles and running commands. Check the new compute environment's documentation to see what needs to be changed in the running script

> **Note:** Things marked in <> mean that you should substitute it whatever your own thing is called.

## Balsam Documentaion

User guide and documentation for Balsam can be found here:
https://balsam.readthedocs.io/en/latest/user-guide/installation/

Check the user guide for additional information on Balsam. The setup tutorial here is mostly taken out of the user guide. 

# Getting started with Balsam and VASP

Before getting started, you need to actually be able to access balsam and VASP on Polaris. You can check if you have permission by running the `groups` command on Polaris. You should see something like `MI2Dmaterials balsam-oauth vasp6` in the output. If you don't have the balsam group you will not be able to run those balsam. To get access to those, you need to send a support ticket to the ALCF help desk.

Check here for instructions:
https://docs.alcf.anl.gov/polaris/workflows/balsam/

## Setting up the Conda Environment

Run the following commands to set up a conda environment with the dependencies needed the run the workflow with:

```
module use /soft/modulefiles
module load conda
conda create -n <conda environment name> python=3.9
conda activate <conda environment name>
pip install --pre balsam
conda install pymatgen
conda install jarvis-tools
```

> **Note**: Any python>=3.9 should work fine for this. If you have any additional code that needs to run in this environment, make sure the python type is compatible. If you are adding balsam applications to run in addition, note that you can also define a different python environment for that script to run in.

Now that the conda environment is created, whenever you relog into Polaris, you will need to run the following to use balsma again:

```
module use /soft/modulefiles
module load conda
conda activate <conda environment name>
```

I personally have the those commands `balsam_env.sh` in my home directory that I call with `source ~/balsam_env.sh` whenever I login.

## Logging in to Balsam

Once you have a working conda environment with balsam inside, you will need to log into balsam. Run the following

```
balsam login
```

Follow the printed link to the login page and login using your ALCF account. 

Once you have logged in you will stay logged in unless you have been idle for a while. I believe the logout period is ~2 weeks or so. If Balsam tells you that you are not logged in after a bit, just follow the steps to login again.

## Creating a Balsam Site

The Balsam site is where all your jobs and defined Balsam applications live. Once you have created a conda environment and are logged in you can now create them.

To do so run: 

```
balsam site init <SITE-PATH>
```

Now give the following inputs to Balsam
- Enter the name of your site. **Note:** The name you put here defines the `site_name` which is a unique identifier for the site that will be needed later, so make sure you know the sites name.
- Navigate to the polaris setup option with the up and down keys and press enter
- The name of the primary allocation is `MI2Dmaterials`, so type that and press enter.

Now navigate to the folder that has the site in it (\<SITE-PATH\>) and start the site:

```
cd <SITE-PATH>
balsam site start
```

This will start a process on the login node that will handle Balsam command inputs, run pre/postprocessing on jobs, and handle queue submissions. To check if the balsam site is running, you can run `balsam site ls`. If the site is not running, just run `balsam site start` again in the site folder. This process will expire after ~2 weeks if left alone that long, so just make sure to start it again if its down and you want to use the site.

If you need to restart the site for whatever reason, you can run `balsam site stop` and then start it again with `balsam site start`. If Balsam tells you that you cannot stop the site because it lives on a different login node, you can delete the process file in the balsam site folder (\<something\>.pid) and the site will automatically stop.

# Extra Balsam Commands

## Listing jobs

To check the status of jobs you may run commands such as:
```
balsam job ls
balsam job ls --by-state
balsam job ls --tag <tag key>=<tag value>
```
see balsam docs for more details

# Creating the workflow

## Creating inputs

To create the balsam workflow you first need to setup an input file directory for the jobs to be created from. This directory should have this format:

```
input_files/
├── Mn1Se1S1
│   ├── KPOINTS
│   ├── POSCAR
|   └── POTCAR
├── Mn1Te1S1
│   ├── KPOINTS
│   ├── POSCAR
|   └── POTCAR
...
```

Here each folder within the input directory should be the name of the material whose inputs are contatined in that folder.

### Substitutions

One way to create the inputs of many structures at once is to create substitutions based on the structure of a base input file. To do this use the `gen_inputs_from_subs.py` script.

To do this create a template directory containing a KPOINTS and POSCAR file of your base material. Note that the POSCAR file should contain the maximum substitutions you want to make. This POSCAR does not need to be relaxed or even for a real material, it justs needs to have a different element on each site that may be different in one of the decorated structures you want to make. In my case working with the Janus oxyhalides the POSCAR looks like this:

```
Cr2O2Cl2 monolayer
1
   3.900   0.000   0.000
  -0.000  3.211   0.000
   0.000   0.000   21.342
Cr O Cl Br
2  2 1  1
Cartesian
       1.95011589       0.00000000       9.81283132
      -0.00000000       1.60564807      11.52917711
       1.95011589       1.60564807      11.08065485
       0.00000000      -0.00000000      10.26136249
       1.95011589       1.60564807       8.11185433
       0.00000000       0.00000000      13.23016203
```

Now you need to modify the `spec.yaml` file to define the substitutions you want to make. This file has 3 parts:

1. **base**: This defines the base structure whose sites are substituted. `sites` defines the site labelling for each element in the template POSCAR file. Elements belonging to that site are considered valid substitutions for that element. `elems` controls which element is at that site in the base structure. For instance if you are making TMDs and your base structure is MnSe2, your base would look like:

```
base:
  name: MnSe2
  sites:
    - A
    - X
  elems:
    - Mn
    - Se
```
2. either **mats** or **subs**: This section has two modes, material (*mats*) mode and substitution (*subs*) mode. In *subs* mode, you define which elements are valid substitutions for each site. The script will create an input for every possible material For instance for the TMDs, you might say Ti, V, Cr, Mn, and Fe are valid substitutions for the A sites and so this section will look like:

```
subs:
  A: ["Ti","V","Cr","Mn","Fe"]
  X: ["S","Se","Te"]
```

In *mats* mode, you explicitly define a list of decorations to the base structure you want to make. This is in the form of a list of elements in the same order as `base.sites`. For instance VSe2 would be `["V","Se"]`. The full format looks like:

```
mats:
  - ["V","Se"]
  - ["Cr","Te"]
  - ["Cr","S"]
  - ["Ti","Se"]
```

3. **pseudos**: This section controls the pseudopotentials used for each structure. The format for this is `<element>: <pseudo>`. This looks like:

```
pseudos:
  H: H
  He: He
  Li: Li_sv
  ...
```

The example spec.yaml has the VASP recommended PBE pseudos in it as of 9/16/2024, but I recommend checking here for the most up to date ones:  https://www.vasp.at/wiki/index.php/Available_pseudopotentials#Standard_potentials

To run the `gen_inputs_from_subs.py` run the following command:

```
python gen_inputs_from_subs.py -o <OUT_DIR> -t <TEMPLATE_DIR> -s <SPEC.YAML FILE> -p <PSEUDOPOTENTIAL_DIR> 
```

This will generate new input files at OUT_DIR.

## Creating Jobs

Now that you have inputs, you need to generate the BALSAM jobs. First you may need to adjust the workflow parameters.

### Defining 'magnetic' atoms:

Depending on the structure, which atoms carry the primary magnetic moment may vary. To address this, `vasp_workflow.py` contains a list of atoms that should be considered magnetic by the workflow. By default this is just all the transition metals. For your purposes you may need to change this. To do this modify the `vasp_workflow.py` script and look for the list called `all_mag_atoms`. Adjust the elements in that list as needed.

### Modifying INCARs:

In `vasp_workflow.py` the common INCAR paramaters is controlled by the dictionary `incar_common` to modify shared incar parameters, change the dictionary there following the format `<PARAM> = "<VALUE>", ...`. To change parameters for only the initial job, change the line `parent_incar = Incar.from_dict(dict(...))` using the same format. To change the parameters of the children jobs, change the line `ic = Incar.from_dict(dict(...))`.

### Running the workflow script

Once you have adjusted the workflow to your liking, run the follow command: 

```
python vasp_workflow.py -s <SITE_NAME> -U <DFT+U value> -S <supercell> [-n <# nodes per calc> -r <RPN> -t <TPR> -g <GPR>]
```

> **Note:** vasp_workflow.py must be in the same folder as parse_outcar.py to work correctly

Here supercell is the supercell used in the FM vs AFM relaxation calculations after an initial rough relaxation job. This should be in the format `[a,b,c]` where the supercell is an $a \times b \times c$ supercell. Here `-n`, `-r`, `-t`, and `-g` are optional parameters that define the job running parameters. The defaults for these should generally work fine for most structures with $\leq 20$-ish atoms in the unit cell. You may need to do some testing to see what works for you on these. The defaults on these are:

- nodes per calculation: 1
- ranks per node (RPN): 4
- threads per rank (TPR): 16
- GPUs per rank (GPR): 1

The BALSAM database will now be populated with jobs. Check to see if this worked with `balsam job ls`.

## Running Jobs

To launch the BALSAM workflow run the following command:

```
balsam queue submit -n <# nodes> -t <time in minutes> -q <queue> -A MI2Dmaterials --site <site_name> --job-mode mpi
```

For example to submit to the debug queue running on two nodes, run:

```
balsam queue submit -n 2 -t 60 -q debug -A MI2Dmaterials --site janus_db --job-mode mpi
```

You may want to check to see if this is submitted with `balsam queue ls` and `qstat -u $USER`. Note that it takes a couple seconds for balsam to realize that it should submit to the queue. Once the balsam submitter is running you can run `balsam job ls --state RUNNING` to check which jobs are running. I like to use the script `monitor.sh` to be track of the status of the queue with the command `watch ./monitor.sh`

Congradulations, you now have a working balsam site with magnetic calculations in it. To add jobs for more materials you can run `vasp_workflow.py` on the same site.
