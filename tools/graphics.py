import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import constants
import os


matplotlib.rcParams['figure.figsize'] = [8.0, 4.0]
matplotlib.rcParams['figure.dpi'] = 200
matplotlib.rcParams['savefig.dpi'] = 200
#matplotlib.rcParams['ytick.labelsize'] = 10
#matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.direction'] = 'in'
    
scatter_params = {'marker':'o', 'markersize':3, 'linewidth':1}
plot_params_random_average = {'marker':'o', 'markersize':3}
plot_params_random_all = {'linestyle':'None', 'marker':'o', 'markersize':3, 'alpha':0.3}
plot_params_normal = {'marker':'o', 'markersize':3, 'linewidth':1}
plot_params_spectrum = {'linewidth':1}

# Change axes_info['keys'] -> axes_keys; axes_info['values'] -> axes_values;?

def getLims(array):
    return array[0],array[-1]

def find_axes(params, observable, **kwargs):
    
    if kwargs.get('multi_level_observable', False):
        lengths = [len(params[key]) for key in params]
        order = np.argsort(lengths)
        observable = np.transpose(observable, axes = order)
        
        # I think that the fact that we are rearranging observable according to order, but getting the longest lengths in another way is making the bug
        # Solution: use order to get the longest params
        
        keys = list(params.keys())
        ordered_keys = [keys[i] for i in order]
        xaxis = ordered_keys[-1]
        columns = ordered_keys[-2]
        rows = ordered_keys[-3]
        colorbar = ordered_keys[-4]

        axes_info = {'keys': {'xaxis': xaxis, 'colorbar': colorbar, 'columns': columns, 'rows': rows},
            'values':{'xaxis': params[xaxis], 'colorbar': colorbar, 'columns': params[columns], 'rows': params[rows]}}
        axes_keys = {'xaxis': xaxis, 'colorbar': colorbar, 'columns': columns, 'rows': rows}
        axes_values = {'xaxis': params[xaxis], 'colorbar': colorbar, 'columns': params[columns], 'rows': params[rows]}
    
    else:
        lengths = [len(params[key]) for key in params]        
        order = np.argsort(lengths)
        observable = np.transpose(observable, axes = order)
        
        keys = list(params.keys())
        ordered_keys = [keys[i] for i in order]
        xaxis = ordered_keys[-1]
        colorbar = ordered_keys[-2]
        columns = ordered_keys[-3]
        rows = ordered_keys[-4]

        axes_info = {'keys': {'xaxis': xaxis, 'colorbar': colorbar, 'columns': columns, 'rows': rows},
                'values':{'xaxis': params[xaxis], 'colorbar': params[colorbar], 'columns': params[columns], 'rows': params[rows]}}
        axes_keys = {'xaxis': xaxis, 'colorbar': colorbar, 'columns': columns, 'rows': rows}
        axes_values = {'xaxis': params[xaxis], 'colorbar': params[colorbar], 'columns': params[columns], 'rows': params[rows]}
        
    #print(axes_info)

        
    return axes_info, observable

def make_grid(values, axes_info, **kwargs):
    
    if kwargs.get('multi_level_observable', False):
        #print('hmm')
        ncols = len(axes_info['values']['columns'])
        nrows = len(axes_info['values']['rows'])

        spec = matplotlib.gridspec.GridSpec(nrows=nrows, 
                                            ncols=ncols, # Not leaving space for colorbar
                                            )

        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex=True, 
                                 sharey=True, 
                                 squeeze=False)#,
                                 #gridspec_kw={'width_ratios': [1]*ncols})


        cax = None

        fig.subplots_adjust(wspace=0.05)
        
        
    else:

        ncols = len(axes_info['values']['columns'])

        spec = matplotlib.gridspec.GridSpec(nrows=1, 
                                            ncols=ncols+1, # Leaving space for colorbar
                                            width_ratios=[1]*ncols + [0.05])

        fig, axes = plt.subplots(nrows=1,
                                 ncols=ncols+1,
                                 sharex=True, 
                                 sharey=True, 
                                 squeeze=False,
                                 gridspec_kw={'width_ratios': [1]*ncols + [0.05]})


        axes[0,-1].tick_params(bottom=False, left=False, labelbottom=False)
        cax = fig.add_subplot(spec[-1])

        fig.subplots_adjust(wspace=0.05)
        
    
    return fig, axes, cax

def make_plots(values, axes_info, cmap, axes, **kwargs):

    if kwargs.get('random', False):
        for indices in np.ndindex(values.shape[:-1]):
            arr = np.stack(values[indices]) # In case there are several points per param set
            if kwargs.get('average', False): 
                mean, std = np.mean(arr,axis=1), np.std(arr,axis=1)
                axes[0,indices[-2]].errorbar(axes_info['values']['xaxis'], 
                                             mean, 
                                             yerr=std, 
                                             c=cmap(indices[-1]), 
                                             **plot_params_random_average)
            else: axes[0,indices[-2]].plot(axes_info['values']['xaxis'], 
                                           arr, 
                                           c=cmap(indices[-1]), 
                                           **plot_params_random_all)
                
    else:

        if kwargs.get('multi_level_observable', False):
            for indices in np.ndindex(values.shape[:-1]):
                vals = np.stack(values[indices])
                axes[indices[-2],indices[-1]].plot(axes_info['values']['xaxis'], 
                                         vals, 
                                         c='k',
                                         **kwargs.get('plot_params', plot_params_normal))

        else:
            for indices in np.ndindex(values.shape[:-1]):
                vals = np.stack(values[indices])

                axes[indices[-3],indices[-2]].plot(axes_info['values']['xaxis'], 
                                         vals, 
                                         #values[indices], 
                                         c=cmap(indices[-1]),
                                         **kwargs.get('plot_params', plot_params_normal),
                                                  zorder = kwargs.get('zordering',1)*indices[-1])

        axes[0,0].set_xlim(*kwargs.get('xlim', getLims(axes_info['values']['xaxis'])))
        if 'ylim' in kwargs:
            axes[0,0].set_ylim(*kwargs.get('ylim'))
        if 'xticks' in kwargs:
            axes[0,0].set_xticks(kwargs.get('xticks'))
        
def make_colorbar(values, axes_info, cmap, fig, cax): 
    
    cax.tick_params(bottom=False, left=False, labelbottom=False, right=True, direction='inout')
    
    ticks = range(len(axes_info['values']['colorbar']))
    norm = matplotlib.colors.Normalize(vmin=ticks[0]-0.5,vmax=ticks[-1]+0.5)
    mappable = matplotlib.cm.ScalarMappable(cmap=cmap,norm=norm)

    cbar = fig.colorbar(mappable,
                        label=constants.params_labels[axes_info['keys']['colorbar']],
                        cax = cax,
                        ticks=ticks)
                        #format='%.3f')
        
    cbar.ax.set_yticklabels(np.around(axes_info['values']['colorbar'][ticks],decimals=2))
        
def make_labels(observable_name, axes_info, fig):
    ax = fig.add_subplot(111, frameon=False, xticks=[], yticks=[])

    ax.set_xlabel(constants.params_labels[axes_info['keys']['xaxis']], labelpad=20)
    ax.set_ylabel(constants.observable_labels.get(observable_name), labelpad=40)

def make_titles(params, axes_info, fig, axes, **kwargs):
    
    if kwargs.get('multi_level_observable', False):
        
        columns_header = constants.params_labels[axes_info['keys']['columns']]
        if kwargs.get('ParamsTitle',False):fig.text(0.1, 0.9, columns_header)
        
        rows_header = constants.params_labels[axes_info['keys']['rows']]
        if kwargs.get('ParamsTitle',False):fig.text(0.0, 0.9, rows_header)
        
        if len(axes)>1:
            for i, ax in enumerate(axes[:,0]):
                label = str(round(axes_info['values']['rows'][i],2))
                if kwargs.get('ParamsTitle',False):ax.set_ylabel(label)
                #fig.text(0.0, 0.6 - i*(0.9-0.1)/len(axes_info['values']['rows']), label, rotation='vertical')

        for i, ax in enumerate(axes[0,:]):
            if kwargs.get('ParamsTitle',False):ax.set_title(round(axes_info['values']['columns'][i],4))

        suptitle=''
        for key in params:
            #if (key not in axes_info['keys'].values()) | (key is axes_info['keys']['colorbar']) | ((key is axes_info['keys']['rows']) and len(params[key])==1):
            if len(params[key])==1:
                suptitle += constants.params_labels[key] + ': ' + str(np.around(params[key],2)) + '; '
        if kwargs.get('ParamsTitle',False):plt.suptitle(suptitle)
        
        if kwargs.get('subplot_label'):
            for i, ax in enumerate((axes).flatten()):
                label = constants.subplot_labels[i]
                ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', color='k', bbox=dict(facecolor='white'))
                #ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=16, fontweight='bold')

        
    else:

        titles_header = constants.params_labels[axes_info['keys']['columns']]
        if kwargs.get('ParamsTitle',False):fig.text(0.1, 0.9, titles_header)

        for i, ax in enumerate(axes[0,:-1]):
            if kwargs.get('ParamsTitle',False):ax.set_title(axes_info['values']['columns'][i])

        suptitle=''
        for key in params:
            if (key not in axes_info['keys'].values()) | ( key==axes_info['keys']['rows'] and len(axes_info['values']['rows'])==1):
                suptitle += constants.params_labels[key] + ': ' + str(np.around(params[key],2)) + '; '
        if kwargs.get('ParamsTitle',False):plt.suptitle(suptitle)
        
        if kwargs.get('subplot_label', True):
            for i, ax in enumerate(axes[0,:-1]):
                label = constants.subplot_labels[i]
                ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', color='k', bbox=dict(facecolor='white'))
    
def save_fig(observable_name, fig, **kwargs):
    
    plot_path = kwargs.get('plot_path')
    if plot_path:
        filename = plot_path + '_' + observable_name + '.png'
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        fig.savefig(filename, facecolor='w', edgecolor='none')        
    
def plot_data(data, observable_name, **kwargs):
    
    main_observable = observable_name
    #if hasattr(observable_name, '__len__'):
    if isinstance(observable_name, list):
        main_observable = observable_name[0]
    
    
    if main_observable=='spectra':
        spread = 0.5
        kwargs['ylim'] = kwargs.get('ylim', (-spread, spread))
        kwargs['plot_params'] = kwargs.get('plot_params', plot_params_spectrum)
        # Unfortunately I cannot make the bands axis (n_level) become part of the observable ndarray
        # Because the length of this axis (=number of bands in the spectrum) changes with n_wire
        kwargs['multi_level_observable'] = True
        
    params = data['params']
    axes_info, observable = find_axes(params, data[main_observable], **kwargs)
            
    fig, axes, cax = make_grid(observable, axes_info, **kwargs)

    cmap = plt.get_cmap('jet',lut=len(axes_info['values']['colorbar']))
    
    make_plots(observable, axes_info, cmap, axes, **kwargs)
    
    if not kwargs.get('multi_level_observable', False):
        make_colorbar(observable, axes_info, cmap, fig, cax=cax)

    make_labels(main_observable, axes_info, fig)
    make_titles(params, axes_info, fig, axes, **kwargs)
    
    save_fig(main_observable, fig, **kwargs)

    return fig
