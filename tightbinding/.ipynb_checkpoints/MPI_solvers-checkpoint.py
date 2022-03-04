"""
This module should have some methods which make use of the mpi4py package to parallelize the computations usually done with my tightbinding package.
"""

from mpi4py import MPI
import numpy as np
import nanowires
import solvers
import graphics

def getConfig(indices, param_dict):
    config = { param : param_val[indices[i]] for i , (param, param_val) in enumerate(param_dict.items())}
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
        #indices = np.array(indices)

        linear_index  = np.ravel_multi_index(indices,dims)
        worker_index = linear_index // n_jobs
        
        # Append the index to the right worker
        workloads[worker_index] += [indices]
        
    return workloads


# Dumb way: give each worker a list of jobs (configs) and wait all of them to finish. Each worker does this:
def work_cycle(config_indices, param_dict):    
    
    # Create list to keep the eigenvalues:
    eigs = []
    
    # While the list is not empty, keep calculating
    for indices in config_indices:
        
        # Get the configs from the parameter set
        config = getConfig(indices, param_dict)
        
        # Initialize a nanowire to that config
        nanowire = nanowires.Nanowire(**config)
        
        # Solve
        eigs += [solvers.solve_single(nanowire)]
    
    return np.array(eigs)

def save_results(data, jobs, eigs_list):
    for indices, eig in zip(jobs,eigs_list):
        data[indices] = eig
    return data

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    test_params = {'n_wire':[100], 'mu':[0.0, 0.5], 'delta':[0.25], 'b':np.linspace(0,2)}
    
    param_dict  = test_params
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    if rank == 0:
        print("Going to enter job")
        
        # Get the dimensions of the parameters
        dims = getDimensions(param_dict)

        # get number of workers, ignore root
        n_workers = comm.Get_size() - 1 

        # create list of indices of jobs that each worker should do
        workloads = split_work(dims, n_workers)
        
        # Send each worker its workload
        for worker_index, job in enumerate(workloads):
            print("Sending ", worker_index+1)
            comm.send(job, dest=worker_index+1)
        
        results = []
        spectra = np.zeros(dims, dtype='object')
        # Collect the results
        for worker_index, _ in enumerate(workloads):
            print("Receiving ", worker_index+1)
            results += [comm.recv(source=worker_index+1)]
            
        # Process the results and rearange them into the data ndarray
        print("Saving Results")
        for worker_index, (jobs, eigs_list) in enumerate(zip(workloads, results)):
            print("Processing ", worker_index)
            spectra = save_results(spectra, jobs, eigs_list)
        
        data = {'spectra':spectra, 'params':param_dict}
        graphics.plot_data(data, 'spectra')
        plt.show()
        
    else:
        job = comm.recv(source=0)
        eigs_list = work_cycle(job,param_dict)
        comm.send(eigs_list, dest=0)