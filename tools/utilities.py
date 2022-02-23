import numpy as np
import time
import nanowires
import solvers
from tkinter import Tk, filedialog
from tqdm import tqdm
import itertools


global b_p
b_p = {'n_sc':[100], 'n_normal':[4], 'a':[15], 'mu':[0.5]}

def start_file(**kwargs):
    base_params = kwargs.get('base_params', b_p)
    root = Tk()
    root.withdraw() # Hide the main window.
    root.call('wm', 'attributes', '.', '-topmost', True) # Raise the root to the top of all windows.
    file = filedialog.askopenfilename()
    if file: # File was chosen
        data = np.load(file,allow_pickle=True).item()
    else: # No file is chosen
        file = filedialog.asksaveasfilename() # Create a file instead
        data = {'params': base_params}
        print('Starting from scratch')
    return file, data

def gen_params(params_old, **params_new):   # Breaks down if more than one parameter is passed:
                                            # If (B=b1, mu=m1) was calculated and we ask for
                                            # (B=b2): (B=b2, mu=m1) is calculated
                                            # (mu=m2): (B=b1, mu=m2) is calculated
                                            # (B=b2, mu=m2): (B=b2, mu=m2 is calculated), not
                                            # all permutations of (B=(b1,b2), mu=(mu1,mu2)) as expected
    params = dict(params_old)
    for key in params_new:
        params[key] = np.asarray([item for item in params_new[key] if item not in params_old[key]])
    return params

def remove_duplicates(params, specs):
    for key in params:
        params[key], indices = np.unique(params[key], return_index=True)
        specs = np.swapaxes(specs, 0, list(params).index(key)) # Puts axis of interest first
        specs = specs[indices]
        specs = np.swapaxes(specs, 0, list(params).index(key)) # Swaps back

def merge_observable(params_old, params_new, observable_old, observable_new):
    params = dict(params_old)
    observable = np.ndarray.copy(observable_old)
    for key in params:
        if not np.array_equal(params_old[key],params_new[key]):
            params[key] = np.concatenate((params_old[key], params_new[key])) # Join the old and new parameter space
            axis = list(params_old).index(key) # Find the axis which the current parameter represents
            observable = np.concatenate((observable_old, observable_new), axis=axis) # Join the old and new data correctly

            # Sort and remove duplicates
            
            params[key], indices = np.unique(params[key], return_index=True)
            observable = np.swapaxes(observable, 0, list(params).index(key)) # Puts axis of interest first
            observable = observable[indices]
            observable = np.swapaxes(observable, 0, list(params).index(key)) # Swaps back

    return params, observable

def merge_and_save_data(data, new_data, file_path):
    
    old_params = dict(data['params'])
    for key in data:
        if key!= 'params':
            if data[key] is None: data['params'], data[key] = new_data['params'], new_data[key]
            else: 
                data['params'], data[key] = merge_observable(old_params, new_data['params'], data[key], new_data[key])

    for key in data['params']:
        data['params'][key] = np.asarray(data['params'][key])

    np.save(file_path, data)

# Takes in parameters, a solver, a parameter name to be randomized
# Returns solved values for those parameters, with the random_param picked from a gaussian distribution
def generate_data(params, solver, **kwargs):
    
    random_param = kwargs.get('random_param')
    randomizer = kwargs.get('randomizer')
    random = kwargs.get('random',False)
    samples = kwargs.get('samples',1)

    junction = nanowires.SNS_Junction() # Returns an SNS_Junction object
    keys = list(params)
    observables = np.ndarray([len(p) for p in params.values()], dtype=kwargs.get('dtype', object))
    
    # Initialize the junction to prevent a weird bug
    for index, key in enumerate(keys):
            setattr(junction, key, params[key][0])
            #print(key, params[key][0], getattr(junction, key))

    for indices, observable in tqdm(np.ndenumerate(observables), total=np.prod(observables.shape)):
            
        for index, key in enumerate(keys):
            setattr(junction, key, params[key][indices[index]])
            #print(key, params[key][indices[index]], getattr(junction, key))
        
        if random:
            observables[indices] = np.empty(samples) # Makes the data entry an empty array
         
            for i in range(samples):

                randomizer(junction, params_temp, random_param)
                observables[indices][i] = (solver(junction, **kwargs)) # Adds the result for this parameter value to list
        
        else: 
            observables[indices] = solver(junction, **kwargs)
            #print(np.sort(abs(observables[indices])))
   
    return  observables