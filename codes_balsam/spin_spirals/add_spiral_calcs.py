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

SCF_incar_overwrite = dict(
    ISTART = 0,
    ICHARG = 2,
    LWAVE = True,
    LCHARG = True,
    LORBIT = 11,
    PREC = "Accurate",
    ENCUT = 500,
    EDIFF = 1.0E-6,
    NELM = 100,
    NSW = 0,
    ISYM = -1,
    IBRION = -1,
    LNONCOLLINEAR = True,
    LSPIRAL = True,
    ENMAX = 600,
    I_CONSTRAINED_M = 1,
    LAMBDA = 0.,
)

def get_rwigs(atoms: Atoms):
    atomic_nums = list(dict.fromkeys(atoms.symbols.numbers))
    rwigs = [covalent_radii[num] for num in atomic_nums]
    return rwigs

def get_qspiral_mag_direction(atoms: Atoms,pos,ref_pos,q: np.ndarray):
    "Give magnitization direction at r=pos for a spin spiral with qspiral=q and m(ref_pos) is rotated like"
    "the y-axis with reference to the spiral."
    "pos, ref_pos: positions in local cell coordinates"
    
    if all(q==0):
        return np.array([0,0,1])
    # z = np.array([0,0,1])
    elif q[2]==0:
        return np.array([0,0,1])
    else:

    # qhat = atoms.cell.cartesian_positions(q)
    # qhat = qhat/np.linalg.norm(qhat)
    # qtheta = np.arccos(qhat[2])
    # qrot = R.from_rotvec(np.cross(z,qhat)*qtheta) # rotate z onto q axis
    # alpha = 2*np.pi*(pos-ref_pos).dot(q)
    # mhat = qrot.apply(np.array([np.cos(alpha),np.sin(alpha),0])) # m = R(qhat)*{cos(q.Ri),sin(q.Ri),0}

        alpha = 2*np.pi*(pos-ref_pos).dot(q)

        mhat = np.array([np.cos(alpha),np.sin(alpha),0]) # m = R(qhat)*{cos(2pi*q.Ri),sin(2*pi*q.Ri),0}
        return mhat
    

def get_magmoms(atoms: Atoms,q,mag_elems=["Mn"],magm=5,nomagm=0.5):
    q = np.array(q)
    
    mag_atoms = [atom for atom in atoms if atom.symbol in mag_elems]
    ref_pos = mag_atoms[0].scaled_position
    magmoms = []
    for at in atoms:
        if at.symbol in mag_elems:
            mm = magm*get_qspiral_mag_direction(atoms,at.scaled_position,ref_pos,q)
        else:
            mm = np.array([0,0,nomagm])
        magmoms.append(mm)
    
    return magmoms

def read_qlist(qlist_file):
    # read a input file containing arguments for qspiral
    with open(qlist_file,'r') as f:
        lines = f.readlines()
    qlist = [[float(num) for num in line.strip().split()] for line in lines]
    return qlist

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--site_name', type=str, required=True)
    parser.add_argument('-i','--input_directory', type=str, nargs='*', required=True, help='<Required> Path to directory containing input files')
    parser.add_argument('-g','--group',type=str)
    parser.add_argument('-q','--qlist', type=str, help='file containing list of qspiral entries')
    parser.add_argument('-n','--num_qpoints', type=int, help='number of q values')
    parser.add_argument('-b','--q_bandpath', type=str, help='bandpath to take q points along')

    parser.add_argument('-k','--kpoints', type=str, required=True, help='KPOINTS file for SCF calculations')
    parser.add_argument('-m','--mag_elems', type=str, nargs='*', required=True, help="List of elements that are magnetic")
    
    run_group = parser.add_argument_group('running arguments')
    run_group.add_argument('-N', '--nodes-per-calculation', type=int, default=NNODES)
    run_group.add_argument('-R', '--ranks-per-node', type=int, default=RPN)
    run_group.add_argument('-T', '--threads-per-rank', type=int, default=TPR)
    run_group.add_argument('-G', '--gpus-per-rank',type=int,default=GPR)

    return parser

def validate_args(args):
    print(args)
    
    assert exists(args.kpoints), f"{args.kpoints} does not exist"
    args.kpoints = abspath(expanduser(args.kpoints))

    start_dir = abspath(expanduser(args.input_dir))
    assert isdir(start_dir), f'Calculation input dir {start_dir} does not exist'
    args.start_dir = start_dir

    args.poscar = join(start_dir, 'POSCAR')
    args.potcar = join(start_dir, 'POTCAR')
    args.incar = join(start_dir, 'INCAR')

    assert exists(args.poscar), f"{args.poscar} does not exist"
    assert exists(args.incar), f"{args.incar} does not exist"
    assert exists(args.potcar), f"{args.potcar} does not exist"
    
    if args.q_bandpath:
        assert args.num_qpoints, "if --q_bandpath is set, -n must also be set!"
    else:
        assert args.qlist, "either -q or both -b and -n must be set!"
        assert exists(args.qlist), f"{args.qlist} does not exist"
        args.qlist = read_qlist(args.qlist)

    assert args.nodes_per_calculation >= 1
    assert args.ranks_per_node >= 1
    assert args.threads_per_rank >= 1
    assert args.gpus_per_rank >= 1

def main(args):

    site_name = args.site_name
    site = Site.objects.get(name=site_name)
    
    with open(args.poscar) as fp:
        poscar_comment = fp.readline()
    system = ''.join(poscar_comment.split())

    if Job.objects.filter(tags={"system":system}).count()>0:
        raise RuntimeError(
            f"workflow for {system} already exists! Use a different script to start from existing relaxation jobs"
        )

    COMMON_JOB_PARAMS = dict(
        num_nodes = args.nodes_per_calculation,
        ranks_per_node = args.ranks_per_node,
        threads_per_rank = args.threads_per_rank,
        threads_per_core=16,
        gpus_per_rank= args.gpus_per_rank,
        launch_params = {"cpu_bind": "depth"},
    )

    tags = {"system":system,"type":"spin_spiral"}

    if args.group:
        tags["group"]=args.group

    jobs = []
    atoms = ase.io.read(args.poscar)

    if args.q_bandpath:
        print(args.q_bandpath)
        qlist = atoms.cell.bandpath(args.q_bandpath,npoints=args.num_qpoints,eps=0.05).kpts[:-1]
    else:
        qlist = args.qlist

    for i,q in enumerate(qlist):

        incar_overwrite = dict(**SCF_incar_overwrite)
        incar_overwrite["RWIGS"] = " ".join([f"{x:.2f}" for x in get_rwigs(atoms)])
        incar_overwrite["QSPIRAL"] = " ".join([f"{x:.4f}" for x in q])
        magmoms = get_magmoms(atoms,q,mag_elems=args.mag_elems)
        print(magmoms)
        incar_overwrite["MAGMOM"] = "  ".join([" ".join([f"{x:.4f}" for x in magmom]) for magmom in magmoms])

        print(f"Creating a spin spiral VASP scf job with q = {q} for {system}")
        q_name = f"q_{i}"

        job_dat = dict(
            source_poscar = args.poscar,
            source_potcar = args.potcar,
            source_incar = args.incar,
            kpoints_source = args.kpoints,
            incar_overwrite = incar_overwrite,
            qspiral = list(q)
        )
        print(q_name)
        job  = Job(app_id="ScfVasp",
                   site_name=site_name,
                   workdir=f"{system}/{q_name}",
                   tags={**tags, "name":q_name},
                   data=job_dat,
                   **COMMON_JOB_PARAMS
                   )
    
        jobs.append(job)

    jobs = Job.objects.bulk_create(jobs)

    print(f"Created spin spiral jobs for {system}")
    

if __name__ == "__main__":
    from argparse import Namespace
    #print(' call parser \n')
    parser = get_parser()
    args = parser.parse_args()

    arg_dict = dict(**vars(args))
    arg_dict.pop("input_directory")

    for input_dir in args.input_directory:
        
        n_args = Namespace(**arg_dict,input_dir=input_dir)
        validate_args(n_args)
        main(n_args)