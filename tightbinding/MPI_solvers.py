"""
This module should have some methods which make use of the mpi4py package to parallelize the computations usually done with my tightbinding package.
"""

from mpi4py import MPI
import numpy as np

def getConfig(indices, param_dict):
    config = { param : param[indices[i]] for i , (param, param_val) in enumerate(param_dict.items())}
    return config

def getDimensions(param_dict):
    dims = []
    for key, val in param_dict.items():
        dims += [len(val)]
    return np.array(dims)

# Calculate which processes should do which jobs
def split_work(dims, n_workers):
    # total number of jobs:
    n_jobs = np.prod(dims)
    # number of jobs per worker, rounded up:
    jobs_per_worker = - ( - n_jobs // n_workers)
    workloads = []
    
    dummy_array = np.empty(dims)
    #dummy_array.shape is dims
    
    # Create n_workers number of empty lists
    workloads = [[]] * n_workers
    
    for indices, _ in np.ndenumerate(dummy_array):
        indices = np.array(indices)

        linear_index  = np.ravel_multi_index(indices,dims)
        worker_index = linear_index // n_jobs
        
        # Append the index to the right worker
        workloads[worker_index] += [indices]
        
    return workloads
    
# Get a job that needs to be done; create multiple processes; split the job; ask each process to do the job until it's all done
def getJobs(param_dict, comm):
    # Get the dimensions of the parameters
    dims = getDimensions(param_dict)
    
    print("dims: ", dims)
    
    # get number of workers, ignore root
    n_workers = comm.Get_size() - 1 
    
    # create list of indices of jobs that each worker should do
    workloads = split_work(dims, n_workers)
    
    return workloads

# Dumb way: give each worker a list of jobs (configs) and wait all of them to finish. Each worker does this:
def work_cycle(config_indices, param_dict):    
    
    # Create list to keep the eigenvalues:
    eigs = []
    
    # While the list is not empty, keep calculating
    for indices in config_indices:
        
        # Get the configs from the parameter set
        config = getConfig(indices)
        
        # Initialize a nanowire to that config
        nanowire = nanowires.Nanowire(config)
        
        # Solve
        eigs += [solver.solve_single(nanowire)]
    
    return eigs

if __name__ == '__main__':
    
    test_params = {'n_wire':[100,200], 'mu':[-0.5, 0.0, 0.5]}
    
    param_dict  = test_params
    
    comm = MPI.COMM_WORLD
    if comm.Get_rank() == 0:
        #job = 'Yay!'
        #comm.send(job,dest=1)
        print("Going to enter job")
        workloads = getJobs(param_dict, comm)
        
        for worker_index, job in enumerate(workloads):
            # Give them the appropriate jobs
            print("Sending")
            comm.send(job, dest=worker_index)
        
    else:
        print("Receiving")
        #job = comm.recv(source=0)
        #print(job)
        #eigs = work_cycle(job,param_dict)
        pass