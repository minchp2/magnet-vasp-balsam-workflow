from balsam.api import Site, Job
import argparse

NNODES = 1
RPN = 1
TPR = 64
GPR = 1

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--site_name', type=str, required=True)
    parser.add_argument('--system', nargs='*',type=str,required=True, help='list of systems to add siesta_2j calcs for')
    parser.add_argument('--parents', nargs='*',type=str,required=True, help='name of parents for siesta_2j calcs')
    parser.add_argument('--kmesh',nargs=3,type=int, help='Kmesh for TB2J calculation')
    parser.add_argument('--kpoints',nargs=3,type=int, help='Kpoints for siesta')

    run_group = parser.add_argument_group('running arguments')
    run_group.add_argument('-N', '--nodes-per-calculation', type=int, default=NNODES)
    run_group.add_argument('-R', '--ranks-per-node', type=int, default=RPN)
    run_group.add_argument('-T', '--threads-per-rank', type=int, default=TPR)
    run_group.add_argument('-G', '--gpus-per-rank',type=int,default=GPR)

    return parser

def validate_args(args):
    assert args.nodes_per_calculation == 1, "only serial mode currently"
    assert args.ranks_per_node == 1, "only serial mode currently implemented!!!"
    assert args.threads_per_rank >= 1
    assert args.gpus_per_rank >= 1

def main(args,system):

    site_name = args.site_name
    site = Site.objects.get(name=site_name)
    
    COMMON_JOB_PARAMS = dict(
        num_nodes = args.nodes_per_calculation,
        ranks_per_node = args.ranks_per_node,
        threads_per_rank = args.threads_per_rank,
        threads_per_core=16,
        gpus_per_rank= args.gpus_per_rank,
        launch_params = {"cpu_bind": "depth"},
    )

    parents = Job.objects.filter(site_id=site.id,tags={"system":system})
    group = parents.first().tags["group"]
   

    parent_ids = [p.id for p in parents if p.tags["name"] in args.parents]
    
    jobs = []

    ax_job_ids = []

    for ax in ['x','y','z']:
        ax_id = ['z','y','x'].index(ax)
        ax_job = None
        for p in parents:
            if f"siesta_2j_{ax}" in p.tags["name"]:
                print(f"siesta_2j_{ax} for {system} already created!")
                ax_job = p
                break
        if ax_job:
            ax_job_ids.append(ax_job.id)
            continue
        
        print(f"Creating siesta_2j_{ax} jobs for {system} with parent ids {parent_ids}")

        kpoints = [x for x in args.kpoints]
        kmesh = [x for x in args.kmesh]
        print(ax,kpoints,kmesh)

        ax_job = Job(app_id="Siesta2J",
            site_name=site_name,
            workdir=f"{system}/siesta_2j/{ax}",
            tags={"system":system,"group":group,'name':f'siesta_2j_{ax}'},
            data={'axis':ax,'tb2j_kmesh':kmesh,'kpoints':kpoints},
            parent_ids=parent_ids,
            **COMMON_JOB_PARAMS)
        

        jobs.append(ax_job)

    new_ax_jobs = Job.objects.bulk_create(jobs)
    ax_job_ids+=[j.id for j in new_ax_jobs]

    print(ax_job_ids)

    if any(["tb2j" in p.tags["name"] for p in parents]):
        print(f"TB2J job already created for {system}")
        return

    kmesh = [x for x in args.kmesh]

    tb2j = Job(
            app_id="TB2J",
            site_name = site_name,
            workdir = f"{system}/siesta_2j",
            tags={"system":system,"group":group,'name':"tb2j"},
            data={'kmesh':kmesh},
            parent_ids=ax_job_ids,
            num_nodes = 1,
            ranks_per_node = 1,
            threads_per_rank = 16,
            threads_per_core=16,
            gpus_per_rank= 1,
            launch_params = {"cpu_bind": "depth"}
    )

    print(f"Creating TB2J job for {system} with parent ids {ax_job_ids}")
    Job.objects.bulk_create([tb2j])

if __name__ == "__main__":
    from argparse import Namespace
    #print(' call parser \n')
    parser = get_parser()
    args = parser.parse_args()
    validate_args(args)
    #print('\n ******validate_args******* \n')
    for system in args.system:
        main(args,system=system)