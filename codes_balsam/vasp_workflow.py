#!/usr/bin/env python
import argparse
import shutil
from balsam.api import ApplicationDefinition, Job, Site
from os.path import abspath, expanduser, exists, isdir, join, dirname
from jarvis.io.vasp.inputs import Poscar
from pymatgen.io.vasp import Incar, Kpoints
from jarvis.tasks.vasp import vasp as jarvis_vasp

import os
import shutil
from balsam.api import ApplicationDefinition, Job, Site
from os.path import abspath, expanduser, exists, isdir, join, dirname
from jarvis.io.vasp.inputs import Potcar, Poscar
from pymatgen.io.vasp import Incar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun
import sys

path = dirname(__file__)

def check_contcar(contcar_path):
    if exists(contcar_path):
        if sum(1 for line in open(contcar_path) if line.strip()) > 1:
            return True
    return False

def check_outcar(outcar_path):
    if not exists(outcar_path): return False
    print("Examining OUTCAR for convergence...")
    with open(outcar_path) as fp:
        for line in fp:
            if 'reached required accuracy' in line: return True
    return False

def preprocess_vasp(job: Job):
    print("preprocessing!!!!")

    site = Site.objects.get(name=job.app.site)

    # If CONTCAR and OUTCAR already exists, send to postprocessor
    if check_contcar('CONTCAR') and exists('OUTCAR'):
        job.state="RUN_DONE"
        job.state_data={"reason":"job finished before creation"}
        return

    # generate input files
    incar = Incar.from_dict(job.data["incar"])
    kpoints = Kpoints.from_dict(job.data['kpoints'])
    incar.write_file("INCAR")
    kpoints.write_file("KPOINTS")

    # If I am a parent job; I need to copy over poscar + potcar
    if job.tags["name"]  == 'initial':
        assert not job.parent_ids
        source_poscar = job.data["source_poscar"]
        source_potcar = job.data["source_potcar"]
        assert exists(source_poscar)
        assert exists(source_potcar)
        shutil.copy(source_poscar,"POSCAR")
        shutil.copy(source_potcar,"POTCAR")
        print("POSCAR and POTCAR successfully loaded for initial job: OK")
        job.state = "PREPROCESSED"
        return

    # Otherwise, get POSCAR from my parent's CONTCAR
    parents = job.parent_query()
    parent = parents.first()
    
    parent_workdir = parent.resolve_workdir(site.path / "data")
    contcar = join(parent_workdir, 'CONTCAR')
    assert exists(contcar), f'Could not find parent file: {contcar}'
    assert check_contcar(contcar), f'Parent CONTCAR {contcar} is empty!'

    print(f"Reading Poscar from {contcar}")
    old_poscar = Poscar.from_file(contcar)
    atoms = old_poscar.atoms
    comment = old_poscar.comment
    if "supercell" in job.data:
        atoms = atoms.make_supercell(job.data["supercell"])
        print(f"Creating a {'x'.join([str(i) for i in job.data['supercell']])} supercell")
    new_poscar = Poscar(atoms,comment=comment)
    new_poscar.write_file("POSCAR")
    shutil.copy(join(parent_workdir,'POTCAR'),'POTCAR')

    job.state = "PREPROCESSED"


def check_finished(outcar_handle):
    for line in outcar_handle:
        if 'reached required accuracy' in line: return True
    return False

def INCAR_replacer(**kwargs):
    '''Each keyword argument represents an INCAR parameter to replace.

    For example, calling INCAR_replacer(ISYM=0) will cause any line
    starting with ISYM to be deleted, and a single replacement line
    "ISYM = 0" to be added to the file in its place.
    '''
    with open('INCAR') as fp:
        incar_lines = fp.readlines()
    output_lines = []
    added_keys = []
    print("Making INCAR modifications...")

    for i, line in enumerate(incar_lines):
        _line = line.strip().upper()
        for key, value in kwargs.items():
            if _line.startswith(key.upper()):
                print("DELETE:  ", line.strip())
                if key not in added_keys:
                    output_lines.append(f'{key} = {value}\n')
                    added_keys.append(key)
                    print("ADD:  ", output_lines[-1].strip())
                break
        else:
            output_lines.append(line)

    for key, value in kwargs.items():
        if key not in added_keys:
            output_lines.append(f'\n{key} = {value}\n')
            print("ADD:  ", output_lines[-1].strip())
            added_keys.append(key)
    incar_text = ''.join(output_lines)
    with open('INCAR', 'w') as fp:
        fp.write(incar_text+'\n')


def ZBRENT_fixer(job):
    '''Copy CONTCAR to POSCAR if exists'''
    if check_contcar('CONTCAR'):
        print('copy CONTCAR to POSCAR')
        shutil.copy('CONTCAR', 'POSCAR')
        print("Fixed! Copied CONTCAR to POSCAR.")
        return True
    else:
        print('CONTCAR not present')
        return False


def PRICEL_fixer(job):
    '''Modify INCAR: set SYMPREC = 1e-8 ; set ISYM = 0'''
    INCAR_replacer(SYMPREC="1.0E-8", ISYM=0)
    fixed = True
    return fixed


def EDDDAV_fixer(job):
    '''Delete CHGCAR file and modify INCAR: set ALGO = All'''
    if os.path.isfile('CHGCAR'):
        os.remove('CHGCAR')
        print("Deleted CHGCAR file")
    INCAR_replacer(ALGO="All")
    fixed = True
    return fixed


def attempt_restart(job, max_retry=4):
    dat = job.data
    retry_count = dat.get("retry_count", 0)
    if retry_count < max_retry:
        retry_count += 1
        job.dat = {**dat, "retry_count": retry_count}
        job.state = "RESTART_READY"
        job.state_data = {"reason": f"restart attempt {retry_count} out of {max_retry}"}
        print("Will retry job: marked RESTART_READY")
        can_retry = True
    else:
        job.state = "FAILED"
        job.state_data = {"reason":f"reached maximum retry attempts of {max_retry}"}
        print("Reached max retry attempts: marked FAILED")
        can_retry = False
    return can_retry


def handle_vasp_relax_done(job: Job):
    print("Handling RUN_DONE..check for convergence")
    if not os.path.exists('OUTCAR'):
        print("Could not find OUTCAR file!  Marking FAILED.")
        job.state = "FAILED"
        job.state_data = {"reason": "Could not find OUTCAR"}
        return

    with open('OUTCAR') as fp:
        converged = check_finished(fp)
    if converged:
        print("Calculation is actually converged")
        sys.path.append(path)
        from parse_outcar import parse_outcar
        calc_results = parse_outcar()
        dat = job.data
        dat.update(calc_results)
        job.data = dict(**dat)
        job.state = "POSTPROCESSED"
        return

    else:
        print("Calculation did not converge; marking RUN_ERROR and invoking error handler")
        job.state = "RUN_ERROR"
        job.state_data = {"reason": 'VASP returned 0; calc did not really finish ==called handle_done()=='}
        job.save()
        handle_vasp_error(job)

def error_scan(outfile):
    # define a dictionary of (ErrorMessage, error_handler_function) pairs
    HANDLERS = {
        'ZBRENT: fatal error in bracketing' : ZBRENT_fixer,
        'internal error in subroutine PRICEL' : PRICEL_fixer,
        'Error EDDDAV: Call to ZHEGV failed' : EDDDAV_fixer,
        'POSMAP internal error: symmetry equivalent atom not found' : PRICEL_fixer,
    }
    # Scan the outfile, line-by-line for a known Error Message
    for line in outfile:
        for err_msg, handler_func in HANDLERS.items():
            if err_msg in line:
                print("Detected error type:", err_msg)
                return err_msg, handler_func
    # If we didnt catch one of the known error types, FAILED
    raise RuntimeError("The error was not recognized as one of my HANDLERS")


def handle_vasp_error(job):
    outfile_name = "OUTCAR"
    try:
        with open(outfile_name) as fp:
            error_type, error_handler = error_scan(fp)
    except RuntimeError as e:
        print("Could not detect error from OUTCAR:", e)
        outfile_name = job.name + ".out"
        print("Checking the standard out/err file", outfile_name)
        with open(outfile_name) as fp:
            error_type, error_handler = error_scan(fp)

    # error_handler is a function that takes *job* argument and returns True if
    # the error was fixed:
    success = error_handler(job)
    if success:
        attempt_restart(job)
    else:
        print("The error handler for", error_type, "failed.")
        job.state = "FAILED"
        job.state_data = {"reason":"could not succesfully handle error type "+error_type}


def handle_vasp_relax_timeout(job):
    print("Handling RUN_TIMEOUT!")
    if check_contcar('CONTCAR'):
        print("Overwriting POSCAR with CONTCAR")
        shutil.copy('CONTCAR', 'POSCAR')
    else:
        print("There was no CONTCAR from which to continue. Restarting with same POSCAR")
    can_retry = attempt_restart(job)
    if not can_retry:
        print("Exceed max retries; marked job as FAILED")
        return
    

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

class StandardVasp(ApplicationDefinition):
    site = "janus_db"
    command_template = f"{VASP_STD_BIN} > vasp.out"

    def shell_preamble(self):
        return VASP_PREAMBLE

    def preprocess(self):
        preprocess_vasp(self.job)

    def postprocess(self):
        handle_vasp_relax_done(self.job)
    
    def handle_timeout(self):
        handle_vasp_relax_timeout(self.job)

    def handle_error(self):
        handle_vasp_error(self.job)

StandardVasp.sync()

class NoncollinearVasp(ApplicationDefinition):
    site = "janus_db"
    command_template = f"{VASP_NCL_BIN} > vasp.out"

    def shell_preamble(self):
        return VASP_PREAMBLE

    def preprocess(self):
        preprocess_vasp(self.job)

    def postprocess(self):
        handle_vasp_relax_done(self.job)
    
    def handle_timeout(self):
        handle_vasp_relax_timeout(self.job)

    def handle_error(self):
        handle_vasp_error(self.job)

NoncollinearVasp.sync()

NNODES = 1
RPN = 4
TPR = 16 # 64/RPN  updated on 3.12.2019   
GPR = 1

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_directory', nargs='+', help='<Required> Path to directory containing input files')
    parser.add_argument('-s','--site_name', type=str, required=True)
    parser.add_argument('-n', '--nodes-per-calculation', type=int, default=NNODES)
    parser.add_argument('-r', '--ranks-per-node', type=int, default=RPN)
    parser.add_argument('-t', '--threads-per-rank', type=int, default=TPR)
    parser.add_argument('-g', '--gpus-per-rank',type=int,default=GPR)
    parser.add_argument('-U', '--dft-uval',type=float,default=3.)
    parser.add_argument('-S', '--supercell',type=list,default=[1,1,1])
    return parser

def validate_args(args):
    print(args)
    start_dir = abspath(expanduser(args.input_dir))
    assert isdir(start_dir), f'Calculation input dir {start_dir} does not exist'
    args.start_dir = start_dir

    args.poscar = join(start_dir, 'POSCAR')
    args.potcar = join(start_dir, 'POTCAR')
    args.kpoints = join(start_dir, 'KPOINTS')
    args.contcar = join(start_dir, 'CONTCAR')
    args.outcar = join(start_dir, 'OUTCAR')

    assert exists(args.poscar), f"{args.poscar} does not exist"
    assert exists(args.potcar), f"{args.potcar} does not exist"
    args.contcar_exists = check_contcar(args.contcar)
    args.outcar_finished = check_outcar(args.outcar)

    assert args.nodes_per_calculation >= 1
    assert args.ranks_per_node >= 1
    assert args.threads_per_rank >= 1
    assert args.gpus_per_rank >= 1

# def bootstrap_apps(site_name):

#     apps_by_name = ApplicationDefinition.load_by_site(site_name)
#     StdVasp = StandardVasp
#     NCLVasp = NoncollinearVasp

#     StdVasp.sync()
#     NCLVasp.sync()
#     # if "StandardVasp" in apps_by_name:
#     #     StdVasp = apps_by_name["StandardVasp"]
#     # else:
#     #     StandardVasp.sync()
#     # if "NoncollinearVasp" in apps_by_name:
#     #     NCLVasp = apps_by_name["NoncollinearVasp"]
#     # else:
#     #     NoncollinearVasp.sync()
    
#     return StdVasp, NCLVasp

def get_mag_ions(atoms):
    """List all magnetic atoms in the Atoms object."""
    all_mag_elements = [
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Ru",
        "Ir",
        "Rh",
        "Os",
        "Rb",
        "Sc",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "Zn",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Ag"
    ]
    els = atoms.elements
    mag_ions = list(set(all_mag_elements).intersection(set(els)))
    return mag_ions

def main(args):

    site_name = args.site_name
    site = Site.objects.get(name=site_name)

    with open(args.poscar) as fp:
        poscar_comment = fp.readline()
    system_str = ''.join(poscar_comment.split())

    if Job.objects.filter(tags={"system":system_str}).count()>0:
        raise RuntimeError(
            f"workflow for {system_str} already exists! Remove all matching jobs if"
            " you really want to add them again"
        )
    
    poscar = Poscar.from_file(args.poscar)
    atoms = poscar.atoms

    # StandardVasp, NoncollinearVasp = bootstrap_apps(site_name)

    print(StandardVasp._site_id)

    COMMON_VASP_JOB_PARAMS = dict(
        app_id = StandardVasp.__app_id__, 
        site_name = StandardVasp.site, 
        num_nodes = args.nodes_per_calculation,
        ranks_per_node = args.ranks_per_node,
        threads_per_rank = args.threads_per_rank,
        threads_per_core=16,
        gpus_per_rank= args.gpus_per_rank,
        launch_params = {"cpu_bind": "depth"},
    )

    incar_common = dict(
        System = system_str,
        ISTART = "1",
        LREAL = False,
        LWAVE = False,
        ADDGRID = False,
        ICHARG = "2",
        GGA = "PE",
        IVDW = "0",
        LORBIT = "11",
        PREC = "Accurate",
        ALGO = "Normal",
        NPAR = "8",
        KPAR = "4",
        ENCUT = "500.0",
        SYMPREC = "1.0E-6",
        NSW = "500",
        IBRION = "2",
        POTIM = "0.3",
        ISIF = "4",
        ISMEAR = "0",
        SIGMA = "0.03",
        ISPIN = "2",
        SAXIS = "0 0 1",
        GGA_COMPAT = ".FALSE.",
        AMIX = "0.2",
        BMIX = "0.00001",
        AMIX_MAG = "0.8",
        BMIX_MAG = "0.00001",
        MAXMIX = "40"
    )

    exclude_params = ["AMIX","BMIX","AMIX_MAG","BMIX_MAG"]
    LDAU_params = {k:v for k,v in jarvis_vasp.add_ldau_incar(atoms=atoms,Uval=args.dft_uval).items() if k not in exclude_params}
    incar_common.update(LDAU_params)
    
    mag_order = jarvis_vasp.MagneticOrdering(atoms)
    mag_ions = get_mag_ions(atoms)
    magmom = [4.5 if el in mag_ions else 0.5 for el in atoms.elements]
    parent_incar = Incar.from_dict(dict(**incar_common,MAGMOM=magmom,LNONCOLLINEAR=False,LSORBIT=False,EDIFFG="-1e-02",EDIFF="1e-04"))
    parent_kpoints = Kpoints.from_file(args.kpoints)
    parent_job = Job(**COMMON_VASP_JOB_PARAMS,
                     tags={"system":system_str,"name":"initial"},
                     workdir=f"{system_str}/initial")
    dat = parent_job.data
    dat['incar']=parent_incar.as_dict()
    dat['kpoints']=parent_kpoints.as_dict()
    dat['source_poscar']=args.poscar
    dat['source_potcar']=args.potcar
    parent_job.data = {**dat}
    parent_job.save()
    print(parent_job.site_id)

    magmom_list,sup_atoms = mag_order.get_unique_magnetic_structures(atoms,supercell_dim=args.supercell,magnetic_ions=mag_ions,magmom=4.5)
    for i,mm in enumerate(magmom_list):
        print(i)
        config_name = "FM" if i==0 else f"AFM{i}"
        magmom = [m if el in mag_ions else 0.5 for m,el in zip(mm,sup_atoms.elements)]
        ic = Incar.from_dict(dict(**incar_common,MAGMOM=magmom,LNONCOLLINEAR=False,LSORBIT=False,EDIFFG="-1e-03",EDIFF="1e-06"))
        kp = Kpoints.from_file(args.kpoints)
        job = Job(**COMMON_VASP_JOB_PARAMS,
                  tags = {"system":system_str, "name":f"{config_name}_relax"},
                  workdir = f"{system_str}/{config_name}_relax",
                  parent_ids = [parent_job.id]
                )
        dat = job.data
        dat['incar']=ic.as_dict()
        dat['kpoints']=kp.as_dict()
        dat['supercell']=args.supercell
        job.data = {**dat}
        job.save()

if __name__ == "__main__":
    from argparse import Namespace
    #print(' call parser \n')
    parser = get_parser()
    args = parser.parse_args()
    #print('\n ******validate_args******* \n')
    for input_dir in args.input_directory:
        n_args = Namespace(**vars(args),input_dir=input_dir)
        validate_args(n_args)
        main(n_args)
