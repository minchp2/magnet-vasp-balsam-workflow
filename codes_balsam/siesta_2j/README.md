# Balsam codes for running Siesta_2J workflow

This is an extension of the magnetic moment workflow. In summary, it will take the relaxed CONTCAR and magnetic configuration from the lowest energy calculation among its parents and will run a series of jobs to calculate the Heisenberg model parameters for the structure using the Green's Function method via Siesta and TB2J.

There are two stages of jobs that are run. First, for each axis, a siesta scf calculation is run with the magnetic moments parallel to that axis. Second, TB2J is used to extract each the heisenberg model parameters, averaging over each axis.

## Preliminaries

To make the workflow able to run, you need to have a working siesta code on your system. Secondly, you need to have set up a balsam site + workflow with the VASP workflow. To set this up see the README file in the head folder.

Now, in `siesta_app.py` replace `SIESTA_BIN` with the path to your executable for siesta and replace `janus_db` in the balsam site name with the name of your balsam site.

> For more information regarding how to use Siesta+TB2J, see the [TB2J documentation](https://tb2j.readthedocs.io/en/latest/src/siesta.html)

## Using the workflow

To load the apps, run `python siesta_app.py`. This will define the apps `Siesta2J` and `TB2J` in the balsam environment. Once these are loaded, you will not need to run this again in the balsam site unless you delete the **apps** with `balsam app rm`. As a note, deleting jobs does not delete the apps.

Once the apps are defined, you can now add jobs using add_siesta_calcs.py. To add new jobs, run:

```
add_siesta_calcs.py -s <site_name> --system System1 System2 ... --parents FM AFM1 --kmesh 7 7 1 --kpoints 11 11 1 
```

 - Replace <site_name> with the name of your balsam site
 - Replace the list System1, etc. with the list of systems you want to add jobs for
 - If you have more than one AFM configuration in the collinear calculations, you can add their names after --parents to add it to the list of configurations that the job considers for picking the ground state. 
 - --kmesh defines the kmesh used by TB2J
 - --kpoints defines the kpoint grid used for siesta.

 Note that siesta will take parameters from the parent incar such as magnetic moments and the value of +U.

## Output

Once these jobs have finished, the output from TB2J can be found in the folder `<data folder>/<system>/TB2J`