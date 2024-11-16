import ase.io.vasp
from balsam.api import ApplicationDefinition, Site, Job
import json
import numpy as np
from ase.io.vasp import read_vasp_out, read_vasp
from ase import Atoms
from pymatgen.io.vasp import Incar

from os.path import exists, join
import os
import shutil

from ase.calculators.siesta import Siesta
from ase.units import Ry


site_name = "janus_db"

site = Site.objects.get(site_name)

SIESTA_BIN = "/home/minchp2/.conda/envs/siesta_env/bin/siesta"

def get_mag_ions(atoms: Atoms):
    """returns a list of species of magnetic ions in atoms"""
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
    els = atoms.symbols
    mag_ions = list(set(els).intersection(set(all_mag_elements)))
    return mag_ions

def check_contcar(contcar_path):
    if exists(contcar_path):
        if sum(1 for line in open(contcar_path) if line.strip()) > 1:
            return True
    return False

def generate_ldau_blocks(atoms: Atoms, ldauu, ldaul):
    ldau_blocks=[]

    # get list of unique species following ldauu and ldaul order
    els = list(dict.fromkeys(atoms.symbols))

    for i,(el,u,l) in enumerate(zip(els,ldauu,ldaul)):
        if l<0:
            continue
        else:
            ldau_blocks.append(f"""{el}.{i+1}   1
    n=3      2
    {u:0.3f}    0.000
    0.000    0.000
    """)
    
    return ldau_blocks

SIESTA_BIN = "/home/minchp2/.conda/envs/siesta2/bin/siesta"

class Siesta2J(ApplicationDefinition):
    site = "janus_db"

    def preprocess(self):
        job = self.job
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
        
        ax = job.data["axis"]

        print(f"Reading Poscar and Magmoms from {contcar}")
        parent_atoms = read_vasp_out(join(parent_workdir,"OUTCAR"))
        parent_incar = Incar.from_file(join(parent_workdir,"INCAR"))
        print(f"Creating rotated structure!")
        
        cwd = os.getcwd()
        os.chdir("..")
        parent_atoms.write("relaxed_structure.vasp",format="vasp")
        cmd = "TB2J_rotate.py relaxed_structure.vasp --ftype poscar"
        os.system(cmd)
        ax_id = ['z','y','x'].index(ax)
        
        print(f"copying over atoms_{ax_id}.poscar to {join(cwd,'POSCAR')}")
        shutil.copy(f"atoms_{ax_id}.poscar",join(cwd,'POSCAR'))
        os.chdir(cwd)

        atoms = read_vasp("POSCAR")
        magmom = parent_atoms.get_magnetic_moments()
        mag_ions = get_mag_ions(parent_atoms)
        ldauu = parent_incar["LDAUU"]
        ldaul = parent_incar["LDAUL"]
        
        dat = job.data
        dat.update({"mag_ions":mag_ions})
        job.data={**dat}

        at = atoms
        print("setting magmoms!")
        print(mag_ions)
        at.set_initial_magnetic_moments([(m if ion in mag_ions else 0) for m,ion in zip(magmom,atoms.symbols)])

        print("Defining siesta inputs")
        ldau_blocks = generate_ldau_blocks(at,ldauu,ldaul)

        calc = Siesta(label='s2tb',
                xc='PBE',
                mesh_cutoff=200 * Ry,
                energy_shift=0.01 * Ry,
                basis_set='DZP',
                kpts=job.data["kpoints"],
                spin='spin-orbit',
                fdf_arguments={'CDF.Compress': '9',
                                'CDF.Save': 'True',
                                'MaxSCFIteration': '60',
                                'SCF.DM.Tolerance': '0.0001',
                                'SCF.EDM.Tolerance': '1e-2 eV',
                                'SCF.H.Tolerance': '1e-3 eV',
                                'SCF.Mixer.History': '16',
                                'SCF.Mixer.Method': 'Pulay',
                                'SCF.Mixer.Spin': 'spinor',
                                'SCF.Mixer.Weight': '0.4',
                                'SCF.Mixer.Kick': '50',
                                'SCF.Spin.Fix': 'True',
                                'SaveHS': 'True',
                                'Write.DMHS.Netcdf': 'True',
                                'SCFMustConverge': 'True',
                                'DM.UseSaveDM': 'True',
                                'LDAU.proj':ldau_blocks
                },
                pseudo_path="/grand/MI2Dmaterials/minchp2/siesta_pseudos/nc-fr-04_pbe_standard/"
        )

        at.calc = calc
        calc.write_input(at,properties='potential_energy')

        print(job.data)
        job.state = "PREPROCESSED"
        job.save()

    command_template = f"{SIESTA_BIN} < s2tb.fdf > s2tb.out"

    def postprocess(self):
        self.job.state = "POSTPROCESSED"
        self.job.save()

Siesta2J.sync()

class TB2J(ApplicationDefinition):
    site = "janus_db"

    python_exe = "/home/minchp2/.conda/envs/balsam_env/bin/python"

    def shell_preamble(self):
        return '''
        module use /soft/modulefiles
        module load conda
        conda activate balsam_env

        module load PrgEnv-nvhpc
        '''

    def preprocess(self):
        job = self.job
        parents = job.parent_query()

        mag_ions = parents.first().data["mag_ions"]

        dat = job.data
        dat.update(dict(mag_ions = mag_ions))
        job.data = {**dat}
        job.state = "PREPROCESSED"
        job.save()

    def run(self):
        kmesh = self.job.data["kmesh"]
        mag_ions = self.job.data["mag_ions"]

        cwd = os.getcwd()
        for ax in ['x','y','z']:
            os.chdir(ax)
            cmd = f"siesta2J.py --fdf_fname s2tb.fdf --elements {' '.join(mag_ions)} --kmesh {' '.join([str(k) for k in kmesh])} > tb2j_run.log"
            print(cmd)
            run_code = os.system(cmd)
            os.chdir(cwd)
            if run_code>0:
                raise RuntimeError(f"axis {ax} failed")
                return
            
            print(f"finished TB2J for axis {ax}")
        print("combining tb2j runs!")

        cmd = "TB2J_merge.py x/ y/ z/ --type structure"
        run_code = os.system(cmd)
        if run_code>0:
            raise RuntimeError("TB2J_merge failed!")
            return
        
        print("merge results!")
        return 0
    
TB2J.sync()
            


print("bootstrapped Siesta2J and TB2J apps")
