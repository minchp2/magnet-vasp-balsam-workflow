from pymatgen.core import Structure, Element
import os
import shutil
import copy
import yaml
import argparse
import itertools

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--out-dir', type=str, help='directory to place generated job inputs', required=True)
    parser.add_argument('-t','--template', type=str, help='templte directory', required=True)
    parser.add_argument('-s', '--spec', type=str, help='specification .yaml file', required=True)
    parser.add_argument('-p', '--pseudos', type=str, help='location of pseudopotential directory',default='/home/minchp2/vasp_potcars')
    return parser

def make_potcar(symbols,potcar_folder):
    potcar_text=""
    for s in symbols:
        print(s,potcar_folder)
        with open(os.path.join(potcar_folder,s,"POTCAR")) as pot_file:
            for line in pot_file:
                potcar_text+=line
    return potcar_text

def ordered(elems,group):
    out = {k:[] for k in set(group)}
    for el, group in zip(elems, group):
        out[group].append(el)
    for k in out:
        if sorted(out[k])!=out[k]:
            return False    
    return True

parser = get_parser()
args = parser.parse_args()

with open(args.spec) as f:
    spec = yaml.safe_load(f)

base_elems = spec["base"]['elems']
if "mats" in spec:
    mats = spec["mats"]
elif "subs" in spec:
    sites = spec["base"]["sites"]
    subs = spec["subs"]
    prod = itertools.product(*[subs[site] for site in sites])
    mats = [list(p) for p in prod if ordered(p,sites)]

pseudos = spec["pseudos"]

base_struct = Structure.from_file(os.path.join(args.template,"POSCAR"))

out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)
template_dir = args.template
potcar_folder = args.pseudos

for sub_elems in mats:
    new_struct = copy.deepcopy(base_struct)
    new_struct.replace_species({Element(base_el): Element(sub_el) for base_el,sub_el in zip(base_elems,sub_elems)})
    new_struct.sort(lambda x: sub_elems.index(str(x.specie)))
    unique_species = list(dict.fromkeys([str(x.specie) for x in new_struct]))
    pps=[pseudos[x] for x in unique_species]
    print(pps)
    new_potcar = make_potcar(pps,potcar_folder)
    form = str(new_struct.composition).split(' ')
    new_material=''.join(form)
    print(str(new_struct.composition))
    mat_dir = out_dir + "/" +new_material
    print(new_material)
    shutil.copytree(template_dir, mat_dir)
    with open(mat_dir+"/POTCAR", mode='w') as f:
        print(new_potcar, file=f)
    with open(mat_dir+"/POSCAR", mode='w') as f:
        new_poscar=new_struct.to(fmt="poscar")
        old_name=new_poscar.split("\n")[0]
        print(new_poscar.replace(old_name,new_material), file=f)
