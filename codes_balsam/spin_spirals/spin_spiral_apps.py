from balsam.api import ApplicationDefinition, Job, Site
from pymatgen.io.vasp import Incar, Kpoints, Vasprun
from jarvis.io.vasp.inputs import Poscar
from os.path import abspath, expanduser, exists, isdir, join, dirname
import numpy as np
import shutil
import argparse

def check_contcar(contcar_path):
    if exists(contcar_path):
        if sum(1 for line in open(contcar_path) if line.strip()) > 1:
            return True
    return False

def load_vasp_files(job: Job, site: Site) -> Incar:
    if not job.parent_ids:
        parent_workdir = None
        print("No parents, reading inputs from input directory!")
        source_poscar = job.data["source_poscar"]
        source_potcar = job.data["source_potcar"]
        source_incar = job.data["source_incar"]
        assert exists(source_poscar)
        assert exists(source_potcar)
        assert exists(source_incar)
        shutil.copy(source_poscar,"POSCAR")
        shutil.copy(source_potcar,"POTCAR")
        incar = Incar.from_file(source_incar)
        print("INCAR, POSCAR, and POTCAR read successfully from inputs")

    else:
        print("parent jobs found, reading inputs from lowest energy parent!")
        parents = job.parent_query()
        energies = [p.data["energy"] for p in parents]
        energies = [(en if type(en)!=list else en[-1]) for en in energies]
        ground_ind = np.argmin(energies)
        print(ground_ind)
        parent = parents[int(ground_ind)]
        print(f"parent job {parent.tags['name']} had the lowest energy!")
        parent_workdir = parent.resolve_workdir(site.path / "data")
        contcar = join(parent_workdir, 'CONTCAR')
        assert exists(contcar), f'Could not find parent file: {contcar}'
        assert check_contcar(contcar), f'Parent CONTCAR {contcar} is empty!'

        print(f"Reading Poscar from {contcar}")
        old_poscar = Poscar.from_file(contcar)
        atoms = old_poscar.atoms
        comment = old_poscar.comment
        new_poscar = Poscar(atoms,comment=comment)
        new_poscar.write_file("POSCAR")
        print(f"reading POTCAR from {parent_workdir}")
        shutil.copy(join(parent_workdir,'POTCAR'),'POTCAR')

        incar = Incar.from_file(join(parent_workdir,'INCAR'))
    
    if job.data['kpoints_source']==".PARENT.":
        assert parent_workdir, "kpoint source was .PARENT. but no parent job found"
        shutil.copy(join(parent_workdir,'KPOINTS'),'KPOINTS')
    else:
        shutil.copy(job.data['kpoints_source'],'KPOINTS')

    return incar

def preprocess_scf_vasp(job: Job):
    print("preprocessing!!!!")

    site = Site.objects.get(name=job.app.site)

    # If CONTCAR and OUTCAR already exists, send to postprocessor
    if check_contcar('CONTCAR') and exists('OUTCAR'):
        job.state="RUN_DONE"
        job.state_data={"reason":"job finished before creation"}
        job.save()
        return

    incar = load_vasp_files(job,site)
    
    print(f"Writing scf INCAR file")
    incar.update(job.data['incar_overwrite'])
    incar.write_file("INCAR")

    job.state = "PREPROCESSED"
    job.save()

def handle_scf_vasp_done(job: Job):
    print("Handling RUN_DONE for scf..check for convergence with PyMatGen")
    vasp_run = Vasprun("vasprun.xml")

    if not hasattr(vasp_run, "converged_electronic") or vasp_run.converged_electronic is None:
        print("Pymatgen: Vasprun does not have a converged_electronic!")
        job.state = "FAILED"
        job.state_data = {"reason": "could not find converged_electronic"}
        job.save()
        return

    e_converged = vasp_run.converged_electronic #used to determien electronic convergence when NSW=0
    if e_converged:
        print("calculations is converged with NSW=0")
        energy = vasp_run.final_energy
        dat = job.data
        dat.update(dict(energy=energy))
        job.data = dict(**dat)
        job.state = "POSTPROCESSED"
        job.save()
        return
    else:
        print("Calculation did not converge; marking RUN_ERROR and invoking error handler")
        job.state = "FAILED"
        job.state_data = {"reason":'VASP returned 0; calc did not really finish'}
        job.save()
#        handle_scf_error(job)

def attempt_restart(job, max_retry=4):
    dat = job.data
    retry_count = dat.get("retry_count", 0)
    if retry_count < max_retry:
        retry_count += 1
        job.dat = {**dat, "retry_count": retry_count}
        job.state = "RESTART_READY"
        job.state_data = {"reason": f"restart attempt {retry_count} out of {max_retry}"}
        job.save()
        print("Will retry job: marked RESTART_READY")
        can_retry = True
    else:
        job.state = "FAILED"
        job.state_data = {"reason":f"reached maximum retry attempts of {max_retry}"}
        job.save()
        print("Reached max retry attempts: marked FAILED")
        can_retry = False
    return can_retry

def handle_scf_vasp_timeout(job: Job):
    print("Handling RUN_TIMEOUT!")
    can_retry = attempt_restart(job)
    if not can_retry:
        print("Exceed max retries; marked job as FAILED")
        return

def handle_scf_vasp_error(job: Job):
    job.state = "FAILED"
    job.state_data = {"reason":"job failed"}
    job.save()

VASP_PREAMBLE = '''
module load PrgEnv-nvhpc
module add cray-libsci
module load craype-accel-nvidia80

NVROOT=${NVIDIA_PATH}

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NVROOT/compilers/extras/qd/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/applications/vasp/aol-libs/3.2/amd-blis/lib/ILP64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/applications/vasp/aol-libs/3.2/amd-libflame/lib/ILP64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/applications/vasp/aol-libs/3.2/amd-fftw/lib

# Change to working directory
echo ${PBS_O_WORKDIR}

export MPICH_GPU_SUPPORT_ENABLED=1
'''

VASP_NCL_BIN = '/soft/applications/vasp/vasp.6.4.3/bin/vasp_ncl'

class ScfVasp(ApplicationDefinition):
    site = "altermagnet_db"
    command_template = f"{VASP_NCL_BIN} > vasp.out"

    def shell_preamble(self):
        return VASP_PREAMBLE

    def preprocess(self):
        preprocess_scf_vasp(self.job)

    def postprocess(self):
        handle_scf_vasp_done(self.job)
    
    def handle_timeout(self):
        handle_scf_vasp_timeout(self.job)

    def handle_error(self):
        handle_scf_vasp_error(self.job)

ScfVasp.sync()