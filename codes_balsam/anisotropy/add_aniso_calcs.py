from balsam.api import Site, Job
import argparse
from os.path import exists, abspath, join, expanduser, isdir
from ase import Atoms
from ase.data import covalent_radii
import ase.io
import numpy as np
from scipy.spatial.transform import Rotation as R

NNODES = 1
RPN = 4
TPR = 16 # 64/RPN   
GPR = 1

# overwrite these from parent INCAR, everything else is copied from parent
SCF_incar_overwrite = dict(
    ISTART = 0,
    ICHARG = 2,
    ALGO = "All",
    LWAVE = True,
    LCHARG = True,
    LORBIT = 11,
    PREC = "Accurate",
    ENCUT = 500,
    EDIFF = 1.0E-6,
    NELM = 100,
    NELMIN = 6,
    NSW = 0,
    ISYM = -1,
    IBRION = -1,
    LNONCOLLINEAR = False
)

NSCF_incar_overwrite = dict(
    ISTART=0,
    ICHARG = 11,
    LCHARG = False,
    LWAVE = False,
    LSORBIT = True,
    LNONCOLLINEAR = True
)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--site_name', type=str, required=True)
    parser.add_argument('--system', nargs='*',type=str,required=True, help='list of systems to add band calculations to')
    parser.add_argument('--parents', nargs='*',type=str,required=True, help='name of parents for SCF calc')
    parser.add_argument('--kpoints', type=str, default=".PARENT.", help='KPOINTS file for SCF calculations')
    parser.add_argument('--keep_parent_kpoints',action='store_true', help='Will keep parent kpoints if flag is present')

    run_group = parser.add_argument_group('running arguments')
    run_group.add_argument('-N', '--nodes-per-calculation', type=int, default=NNODES)
    run_group.add_argument('-R', '--ranks-per-node', type=int, default=RPN)
    run_group.add_argument('-T', '--threads-per-rank', type=int, default=TPR)
    run_group.add_argument('-G', '--gpus-per-rank',type=int,default=GPR)

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

    assert args.nodes_per_calculation >= 1
    assert args.ranks_per_node >= 1
    assert args.threads_per_rank >= 1
    assert args.gpus_per_rank >= 1

def main(args, system):

    site_name = args.site_name
    site = Site.objects.get(name=site_name)
    
    COMMON_VASP_JOB_PARAMS = dict(
        num_nodes = args.nodes_per_calculation,
        ranks_per_node = args.ranks_per_node,
        threads_per_rank = args.threads_per_rank,
        threads_per_core=16,
        gpus_per_rank= args.gpus_per_rank,
        launch_params = {"cpu_bind": "depth"},
    )
    
    parents = Job.objects.filter(tags={"system":system})    

    tags = {"system":system,"type":"aniso"}

    ptags = parents.first().tags
    if "group" in ptags:
        tags["group"]=ptags["group"]

    if any([p.tags["name"]=="aniso_scf" for p in parents]):
        print(f"aniso scf job for {system} already created!")
        aniso_scf = parents.filter(tags={"name":"aniso_scf"}).first()
    else:
        parent_ids = [p.id for p in parents if p.tags["name"] in args.parents]
        print(f"Creating aniso scf job for {system} with parent ids {parent_ids}")
    
        job_dat = dict(
            incar_overwrite = SCF_incar_overwrite,
            kpoints_source = args.kpoints,
            scf = True
        )

        aniso_scf = Job(**COMMON_VASP_JOB_PARAMS,
                    app_id = "ElectronicRelaxVasp", 
                    site_name = site_name,
                    workdir=f"{system}/aniso_scf",
                    tags={**tags, "name":"aniso_scf"},
                    parent_ids=parent_ids,
                    data=job_dat,
                    parameters={"ncl":False})
        
        aniso_scf.save()

    print(f"Creating aniso nscf jobs for {system} with parent scf job {aniso_scf.id}")
    nscf_jobs = []
    for ax,saxis in zip(['x','y','z'],[[1,0,0],[0,1,0],[0,0,1]]):
        job_dat = dict(
            incar_overwrite = dict(**NSCF_incar_overwrite,SAXIS=saxis),
            kpoints_source = ".PARENT.",
            scf = False
        )
        nscf_job = Job(**COMMON_VASP_JOB_PARAMS,
                       app_id = "ElectronicRelaxVasp",
                       site_name = site_name,
                       workdir=f"{system}/aniso_nscf_{ax}",
                       tags={**tags,"name":f"aniso_nscf_{ax}"},
                       parent_ids=[aniso_scf.id],
                       data=job_dat,
                       parameters={"ncl":True})
        
        nscf_jobs.append(nscf_job)

    nscf_jobs = Job.objects.bulk_create(nscf_jobs)

if __name__ == "__main__":
    from argparse import Namespace
    #print(' call parser \n')
    parser = get_parser()
    args = parser.parse_args()
    validate_args(args)

    for system in args.system:
        main(args,system=system)