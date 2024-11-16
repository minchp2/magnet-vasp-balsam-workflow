from balsam.api import ApplicationDefinition, Job, Site
from pymatgen.io.vasp import Incar, Kpoints, Vasprun
from jarvis.io.vasp.inputs import Poscar
from os.path import abspath, expanduser, exists, isdir, join, dirname
import numpy as np
import shutil
import argparse

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

VASP_STD_BIN = '/soft/applications/vasp/vasp.6.4.3/bin/vasp_std'
VASP_NCL_BIN = '/soft/applications/vasp/vasp.6.4.3/bin/vasp_ncl'

# overwrite these from parent INCAR, everything else is copied from parent
SCF_incar_overwrite = dict(
    ISTART = "0",
    ICHARG = "2",
    LWAVE = True,
    LCHARG = True,
    LORBIT = "11",
    PREC = "Accurate",
    ENCUT = "500",
    EDIFF = "1.0E-8",
    NELM = "100",
    NSW = "0",
    IBRION = "-1"
)

NSCF_incar_overwrite = dict(
    ISTART=1,
    ICHARG = 11,
)

def preprocess_band_scf_vasp(job: Job):
    print("preprocessing!!!!")

    site = Site.objects.get(name=job.app.site)
    from vasp_workflow import check_contcar

    # If CONTCAR and OUTCAR already exists, send to postprocessor
    if check_contcar('CONTCAR') and exists('OUTCAR'):
        job.state="RUN_DONE"
        job.state_data={"reason":"job finished before creation"}
        job.save()
        return

    # Get POSCAR from lowest energy parent
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
    
    print(f"Writing INCAR file")
    incar = Incar.from_file(join(parent_workdir,'INCAR'))
    incar.update(job.data['incar_overwrite'])
    incar.write_file("INCAR")

    shutil.copy(join(parent_workdir,'POTCAR'),'POTCAR')
    if job.data['kpoints_source']==".PARENT.":
        shutil.copy(join(parent_workdir,'KPOINTS'),'KPOINTS')
    else:
        shutil.copy(job.data['kpoints'],'KPOINTS')

    job.state = "PREPROCESSED"
    job.save()

def handle_band_scf_vasp_done(job: Job):
    from parse_outcar import parse_outcar
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
        calc_results = parse_outcar()
        dat = job.data
        dat.update(calc_results)
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

def handle_band_nscf_vasp_done(job: Job):
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
        ### POSTPROCESSING FOR NSCF BAND?????????
        ### stuff
        job.state = "POSTPROCESSED"
        job.save()
        return
    else:
        print("Calculation did not converge; marking RUN_ERROR and invoking error handler")
        job.state = "FAILED"
        job.state_data = {"reason":'VASP returned 0; calc did not really finish'}
        job.save()

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

def handle_band_vasp_timeout(job: Job):
    print("Handling RUN_TIMEOUT!")
    can_retry = attempt_restart(job)
    if not can_retry:
        print("Exceed max retries; marked job as FAILED")
        return

def handle_band_vasp_error(job: Job):
    job.state = "FAILED"
    job.state_data = {"reason":"job failed"}
    job.save()

class BandScfVasp(ApplicationDefinition):
    site = "janus_db"
    command_template = f"{VASP_STD_BIN} > vasp.out"

    def shell_preamble(self):
        return VASP_PREAMBLE

    def preprocess(self):
        preprocess_band_scf_vasp(self.job)

    def postprocess(self):
        handle_band_scf_vasp_done(self.job)
    
    def handle_timeout(self):
        handle_band_vasp_timeout(self.job)

    def handle_error(self):
        handle_band_vasp_error(self.job)

BandScfVasp.sync()

class BandNscfVasp(ApplicationDefinition):
    site = "janus_db"
    command_template = f"{VASP_STD_BIN} > vasp.out"
    def shell_preamble(self):
        return VASP_PREAMBLE

    def preprocess(self):
        incar = Incar.from_file('INCAR')
        incar.update(self.job.data['incar_overwrite'])
        incar.write_file("INCAR")
        shutil.copy(self.job.data['kpath_source'],'KPOINTS')
        self.job.state = "PREPROCESSED"
        self.job.save()

    def postprocess(self):
        handle_band_nscf_vasp_done(self.job)
    
    def handle_timeout(self):
        handle_band_vasp_timeout(self.job)

    def handle_error(self):
        handle_band_vasp_error(self.job)

BandNscfVasp.sync()

NNODES = 1
RPN = 4
TPR = 16 # 64/RPN  updated on 3.12.2019   
GPR = 1

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--site_name', type=str, required=True)
    parser.add_argument('-n', '--nodes-per-calculation', type=int, default=NNODES)
    parser.add_argument('-r', '--ranks-per-node', type=int, default=RPN)
    parser.add_argument('-t', '--threads-per-rank', type=int, default=TPR)
    parser.add_argument('-g', '--gpus-per-rank',type=int,default=GPR)
    parser.add_argument('--system', nargs='*',type=str,required=True, help='list of systems to add band calculations to')
    parser.add_argument('--parents', nargs='*',type=str,required=True, help='name of parents for Band SCF calc')
    parser.add_argument('--kpoints', type=str, default=".PARENT.", help='KPOINTS file for SCF calculations')
    parser.add_argument('--keep_parent_kpoints',action='store_true', help='Will keep parent kpoints if flag is present')
    parser.add_argument('--kpath',type=str,required=True, help='KPATH file for NSCF calculation')

    return parser

def validate_args(args):
    print(args)
    
    if args.keep_parent_kpoints:
        args.kpoints=".PARENT."
        print("using Kpoints from parent")
    elif args.kpoints!=".PARENT.":
        assert exists(args.kpoints), f"{args.kpoints} does not exist"
    else:
        assert False, f"No Kpoints Given!"
    assert exists(args.kpath), f"{args.kpath} does not exist"

    assert args.nodes_per_calculation >= 1
    assert args.ranks_per_node >= 1
    assert args.threads_per_rank >= 1
    assert args.gpus_per_rank >= 1

def main(args,system):

    site_name = args.site_name
    site = Site.objects.get(name=site_name)
    
    print(BandScfVasp._site_id)

    COMMON_VASP_JOB_PARAMS = dict(
        num_nodes = args.nodes_per_calculation,
        ranks_per_node = args.ranks_per_node,
        threads_per_rank = args.threads_per_rank,
        threads_per_core=16,
        gpus_per_rank= args.gpus_per_rank,
        launch_params = {"cpu_bind": "depth"},
    )

    parents = Job.objects.filter(tags={"system":system})
    if any([p.tags["name"]=="band_scf" for p in parents]):
        print(f"band scf job for {system} already created!")
        band_scf = parents.filter(tags={"name":"band_scf"}).first()
    else:
        parent_ids = [p.id for p in parents if p.tags["name"] in args.parents]
        print(f"Creating band scf job for {system} with parent ids {parent_ids}")
    
        band_scf = Job(**COMMON_VASP_JOB_PARAMS,
                    app_id = BandScfVasp.__app_id__, 
                    site_name = BandScfVasp.site,
                    tags={"system":system,"name":"band_scf"},
                    workdir=f"{system}/band",
                    parent_ids=parent_ids)
        dat = band_scf.data
        dat['incar_overwrite']=SCF_incar_overwrite
        dat['kpoints_source']=args.kpoints
        band_scf.data = {**dat}
        band_scf.save()
        print(band_scf.site_id)

    print(f"Creating band nscf job for {system} with parent scf job {band_scf.id}")
    band_nscf = Job(app_id=BandNscfVasp.__app_id__,
                    site_name=BandNscfVasp.site,
                    num_nodes=args.num_nodes,
                    ranks_per_node=args.ranks_per_node,
                    threads_per_rank=args.threads_per_rank,
                    threads_per_core=args.threads_per_core,
                    gpus_per_rank=args.gpus_per_rank,
                    launch_params= {"cpu_bind":"depth"},
                    tags={"system":system,"name":"band_nscf"},
                    workdir=f"{system}/band",
                    parent_ids=[band_scf.id])
    dat = band_nscf.data
    # overwrite incar params from scf calculation
    dat['incar_overwrite']=NSCF_incar_overwrite
    dat['kpath_source']=args.kpath
    band_nscf.data = {**dat}
    band_nscf.save()
    print(band_nscf.site_id)

if __name__ == "__main__":
    from argparse import Namespace
    #print(' call parser \n')
    parser = get_parser()
    args = parser.parse_args()
    validate_args(args)
    #print('\n ******validate_args******* \n')
    for system in args.system:
        main(args,system=system)