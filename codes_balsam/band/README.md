# Balsam codes for running band structure workflow

This is an extension of the magnetic moment workflow. In summary, it will take the relaxed CONTCAR and magnetic configuration from the lowest energy calculation among its parents and will run a VASP scf job with spins-orbit coupling off. Then a VASP nscf band structure calculation will run along with the set kpath.

## Preliminaries

To make the workflow able to run, you need to have set up a balsam site + workflow with the VASP workflow. To set this up, see the README file in the head folder.

Now, in `add_band_calcs.py` replace `VASP_BIN_NCL` and `VASP_BIN_STD` with the path to your executables for vasp and replace `janus_db` in the balsam site name with the name of your balsam site.

## Using the workflow

To add the band structure calculations for a set of systems, run the following

```
add_band_calcs.py -s <site_name> --system System1 System2 ... --parents FM AFM1 --kpoints <KPOINTS file> --kpath <KPATH file> 
```

 - Replace <site_name> with the name of your balsam site
 - Replace the list System1, etc. with the list of systems you want to add jobs for
 - If you have more than one AFM configuration in the collinear calculations, you can add their names after --parents to add it to the list of configurations that the job considers for picking the ground state. 
 - if you want to just use the same kpoints as the parent, replace the --kpoints flag with --keep_parent_kpoints. Otherwise balsam will use the kpoints from the KPOINTS file passed to the kpoints argument

To modify the VASP running parameters, you can set these in add_band_calcs.py. Just change the arguments set in `SCF_incar_overwrite` and `NSCF_incar_overwrite`.