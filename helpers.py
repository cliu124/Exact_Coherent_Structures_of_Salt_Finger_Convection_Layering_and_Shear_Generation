########### APS Publication-Ready Matplotlib Configuration ###########
import matplotlib.pyplot as plt
import numpy as np

width_single = 3.38  # APS single-column width in inches
width_double = 6.75   # APS double-column width in inches

# Use Matplotlib's built-in colorblind-friendly palette
# plt.style.use('tableau-colorblind10')

# Configure clean, publication-ready settings
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.size': 8,              # Standard 8 pt text for APS single-column
    # 'mathtext.fontset': 'stix',
    # 'axes.spines.top': False,     # Removes unnecessary top border line
    # 'axes.spines.right': False,   # Removes unnecessary right border line
    'pdf.fonttype': 42            # Keeps text editable in vector PDFs
})


import os
import h5py
import pandas as pd
import matplotlib as mpl
import matplotlib.colors as mcolors
import netCDF4 as nc
import pickle
import glob
from pathlib import Path

# unstable criteria
unstable_criteria = 8e-3 # growth rate threshold to consider a mode unstable, thermal diffusivity scale


def load_file(path):
    if os.path.exists(path):
        with open(path, "r") as f:
            val = float(f.read().strip())
            return val
    return None

def return_coordinates(file_path):
    if file_path.endswith('.nc'):
        with nc.Dataset(file_path, 'r') as ncf:
            x = ncf.variables['x'][:]  # grid coordinates
            z = ncf.variables['z'][:]
            return x.ravel(), z.ravel()
    elif file_path.endswith('.h5'):
        with h5py.File(file_path, 'r') as h5f:
            if 'x' in h5f: x = h5f['x'][:]
            elif 'xg' in h5f: x = h5f['xg'][:] # grid coordinates
            if 'z' in h5f: z = h5f['z'][:]
            elif 'zg' in h5f: z = h5f['zg'][:]
            return x.ravel(),z.ravel()
        
def return_fields(file_path):
    if file_path.endswith('.nc'):
        with nc.Dataset(file_path, 'r') as ncf:
            x = ncf.variables['x'][:]  # grid coordinates
            z = ncf.variables['z'][:]
            te = ncf.variables['te'][:].T  # temperature
            sa = ncf.variables['sa'][:].T  # salinity
            u = ncf.variables['u_0'][:].T  # x-velocity
            w = ncf.variables['u_1'][:].T  # z-velocity
            return x.ravel(), z.ravel(), te, sa, u, w
    elif file_path.endswith('.h5'):
        with h5py.File(file_path, 'r') as h5f:
            # print(h5f.keys())
            if 'x' in h5f: x = h5f['x'][:]
            elif 'xg' in h5f: x = h5f['xg'][:] # grid coordinates
            if 'z' in h5f: z = h5f['z'][:]
            elif 'zg' in h5f: z = h5f['zg'][:]
            te = h5f['t'][:] # temperature
            sa = h5f['s'][:] # salinity
            u = h5f['u'][:] # x-velocity
            w = h5f['w'][:] # z-velocity
            return x.ravel(),z.ravel(), te, sa, u, w

def plot_sample(path, name=None, color_surface = 'PuOr_r',color_S = 'r', 
                plot_mean_profiles=True, symmetric=False):
    fig, ax = plt.subplots(figsize = (3,3))
    [x,z,_,s,_,_]=return_fields(path)
    minvalue = np.min(s)
    maxvalue = np.max(s)
    X, Z = np.meshgrid(x, z) 
    ax.pcolormesh(X, Z, s.T, vmin=minvalue if not symmetric else -maxvalue, vmax=maxvalue, cmap = color_surface, shading='auto')
    ax.set_aspect('equal')
    ax.set_yticks([])
    ax.set_xticks([])
    
    if plot_mean_profiles:
        axs = ax.inset_axes([0, 0, 1, 1], transform=ax.transAxes)
        axs.patch.set_alpha(0)
        total_s = s.mean(axis=0) + z
        axs.plot(total_s, z, color=color_S, lw=5)
        axs.set_ylim(z.min(), z.max())
        axs.set_xticks([])
        axs.set_yticks([])
        
    fig.tight_layout()
    if name is not None: fig.savefig(name, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

#### Load data
def load_energy_eq(file_path,xlabel='D',ylabel='I'):
    df = pd.read_csv(file_path, skipinitialspace=True)
    D, I = df[xlabel].values[-1], df[ylabel].values[-1]
    return D, I

def load_energy_po(file_pattern,xlabel='D',ylabel='I'):
    D = []
    I = []
    for fname in sorted(glob.glob(file_pattern)):
        with h5py.File(fname, 'r') as f:
            D.append(f['tasks'][xlabel][:])
            I.append(f['tasks'][ylabel][:])
    D = np.concatenate(D, axis=0).ravel()
    I = np.concatenate(I, axis=0).ravel()
    return D, I

def load_energy_dns(file_pattern,xlabel='D',ylabel='I'):
    files = glob.glob(file_pattern)
    
    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {file_pattern}")

    t_list, D_list, I_list = [], [], []

    for fname in files:
        with h5py.File(fname, 'r') as f:
            t_list.append(f['scales/sim_time'][:])
            D_list.append(f[f'tasks/{xlabel}'][:])
            I_list.append(f[f'tasks/{ylabel}'][:])
    t = np.concatenate(t_list, axis=0).ravel()
    D = np.concatenate(D_list, axis=0).ravel().real
    I = np.concatenate(I_list, axis=0).ravel().real
    
    # Sort arrays by simulation time
    sort_order = np.argsort(t)
    return D[sort_order], I[sort_order]

def load_snapshot_dns(file_pattern,field='sa'):
    files = glob.glob(file_pattern)
    
    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {file_pattern}")
        
    t_list, snapshot_list = [], []

    loaded_coordinates = False
    for fname in files:
        with h5py.File(fname, 'r') as f:
            if not loaded_coordinates:
                x_key = [k for k in f['scales'].keys() if k.startswith('x_hash')][0]
                x = np.squeeze(f[f'scales/{x_key}'][:])
                z_key = [k for k in f['scales'].keys() if k.startswith('z_hash')][0]
                z = np.squeeze(f[f'scales/{z_key}'][:]) 
                loaded_coordinates = True

            t_list.append(np.atleast_1d(f['scales/sim_time'][:]))
            snapshot_list.append(f[f'tasks/{field}'][:])

        
    t_list = np.concatenate(t_list, axis=0).ravel()
    snapshot_list = np.concatenate(snapshot_list, axis=0)
    # print(f"shape of t_list: {np.shape(t_list)}")  
    # print(f"shape of field_list: {np.shape(snapshot_list)}")  
    
    # Sort arrays by simulation time
    sort_order = np.argsort(t_list)
    
    return x, z, t_list[sort_order], snapshot_list[sort_order]


##########################################
###### Load bifurcation information ######

def get_bifurcation(folder_path, mu = "mu", response="Sh", relativeX=False, relativeZ=False, periodicOrbit=False):
    if os.path.isfile(os.path.join(folder_path, "flow_properties.txt")):
        df0 = pd.read_csv(os.path.join(folder_path, "flow_properties.txt"), skipinitialspace=True)
    elif os.path.isfile(os.path.join(folder_path, "flow_properties.csv")):
        df0 = pd.read_csv(os.path.join(folder_path, "flow_properties.csv"), skipinitialspace=True)
    elif os.path.isfile(os.path.join(folder_path, "continuation.txt")):
        df0 = pd.read_csv(os.path.join(folder_path, "continuation.txt"), skipinitialspace=True)
    else:
        print(f"Warning: {folder_path}flow_properties.txt or flow_properties.csv or continuation.txt not found.")
        return [], [], [], []
    
    # print(df0.columns)
    mu_list = df0[mu] if mu in df0.columns else df0[df0.columns[0]]
    if periodicOrbit:
        response_list = np.zeros((len(mu_list),2))
    else:
        response_list = df0[response] if response in df0.columns else np.zeros(len(mu_list))
    
    dir_list = df0["dir"]
    dir_list = [s.strip() for s in dir_list]

    stability = np.zeros(len(mu_list), dtype=bool)
    growthrate = np.zeros(len(mu_list), dtype=np.complex128)

    # Initialize all optional outputs as None
    Tp = None
    ax = None
    az = None
    # print(f"df0.columns = {df0.columns}")
    if relativeX or relativeZ or periodicOrbit:
        Tp = df0["Tp"] if "Tp" in df0.columns else np.full(len(mu_list), 0.2)
    if relativeX:
        ax = df0["ax"] if "ax" in df0.columns else None
    if relativeZ:
        az = df0["az"] if "az" in df0.columns else None
        
    file_path = "/time-dependent/full_solution/full_solution_s1.h5" if os.path.isdir(folder_path+"search-0/time-dependent/") else "/flowproperties/flowproperties_s1.h5"
    
    for element, dir_element in enumerate(dir_list):
        if periodicOrbit:
            with h5py.File(folder_path+dir_element+file_path, mode='r') as file:
                response_data = file['tasks'][response][:]
            response_list[element,0] = min(np.array(response_data[:,0,0].real))
            response_list[element,1] = max(np.array(response_data[:,0,0].real))
        stability[element], growthrate[element] = check_stability(folder_path+dir_element+"stability/growthrate.mtx")
        # print(f"        Processed element {element+1} in {dir_element}: {mu}={mu_list[element]}, {response}={response_list[element]}, stability={'unstable' if stability[element] else 'stable'}, growthrate={growthrate[element]}")

    return mu_list,response_list,stability,growthrate,Tp,ax,az

def check_stability(growth_rates_path):
    if not os.path.exists(growth_rates_path):
        print(f"Warning: {growth_rates_path} not found.")
        return True, 0
    from scipy.io import mmread
    growth_rates = mmread(growth_rates_path).astype(np.complex128)
    fastest_growth_rate = growth_rates[0,0] # only take the first mode
    if fastest_growth_rate is None:
        return True, None
    if fastest_growth_rate.real > unstable_criteria:
        return True, fastest_growth_rate
    else:
        return False, fastest_growth_rate

def load_bifurcation_info(path_dict, mu='Lx', response='Sh', relativeX=False, relativeZ=False, periodicOrbit=False, save_csv=None):
    mu_list, response_list, is_stable_list, growthrate_list = [], [], [], []
     # use these for traveling waves and relative periodic orbits
    Tp_list = [] if relativeX or relativeZ or periodicOrbit else None
    ax_list = [] if relativeX else None
    az_list = [] if relativeZ else None

    for dir, start, end, upward in path_dict:
        print(f"    Loading continuation data from {dir}...")
        mu_slice,response_slice,is_stable_slice,growthrate_slice,Tp_slice,ax_slice,az_slice = get_bifurcation(folder_path=dir, 
                                            mu=mu, 
                                            response=response, 
                                            relativeX=relativeX, relativeZ=relativeZ, 
                                            periodicOrbit=periodicOrbit)
        mu_slice = mu_slice[start:end] if upward else mu_slice[start:end][::-1]
        response_slice = response_slice[start:end] if upward else response_slice[start:end][::-1]
        is_stable_slice = is_stable_slice[start:end] if upward else is_stable_slice[start:end][::-1]
        growthrate_slice = growthrate_slice[start:end] if upward else growthrate_slice[start:end][::-1]
        if relativeX or relativeZ or periodicOrbit: Tp_slice = Tp_slice[start:end] if upward else Tp_slice[start:end][::-1]
        if relativeX: ax_slice = ax_slice[start:end] if upward else ax_slice[start:end][::-1]
        if relativeZ: az_slice = az_slice[start:end] if upward else az_slice[start:end][::-1]
        
        mu_list.extend(mu_slice)
        response_list.extend(response_slice)
        is_stable_list.extend(is_stable_slice)
        growthrate_list.extend(growthrate_slice)
        if relativeX or relativeZ or periodicOrbit: Tp_list.extend(Tp_slice)
        if relativeX: ax_list.extend(ax_slice)
        if relativeZ: az_list.extend(az_slice)
    is_stable_list = ~np.array(is_stable_list)
    print("")
    # save data of bifurcation to csv file
    if save_csv is not None:
        data2csv(save_csv,mu_list,response_list,is_stable_list,growthrate_list,
                 Tp_list,ax_list,az_list,
                 mu=mu,response=response)
    return mu_list,response_list,is_stable_list,growthrate_list,Tp_list,ax_list,az_list

def merge_branch(mu_list,response_list,is_stable_list,growthrate_list,Tp_list,ax_list,az_list,
                 type='first'):
    """ type = 'first'/'last'/'both' """
    response_array = np.array(response_list)
    if response_array.ndim == 2:  # periodic orbit case
        # add one more point to connect two min. and max. branch of PO
        if type == 'first':
            mu_list = np.concatenate([[mu_list[0]-(mu_list[1]-mu_list[0])], mu_list])
            val = 0.5 * (response_list[0, 0] + response_list[0, 1])
            response_list = np.concatenate([[[val, val]],response_list], axis=0)
            is_stable_list = np.concatenate([[is_stable_list[0]], is_stable_list])
            growthrate_list = np.concatenate([[growthrate_list[0]], growthrate_list])
            if Tp_list is not None:
                Tp_list = np.concatenate([[Tp_list[0]], Tp_list])
            if ax_list is not None:
                ax_list = np.concatenate([[ax_list[0]], ax_list])
            if az_list is not None:
                az_list = np.concatenate([[az_list[0]], az_list])
        elif type == 'last':
            mu_list = np.concatenate([mu_list, [mu_list[-1]-(mu_list[-2]-mu_list[-1])]])
            val = 0.5 * (response_list[-1, 0] + response_list[-1, 1])
            response_list = np.concatenate([response_list, [[val, val]]],axis=0)
            is_stable_list = np.concatenate([is_stable_list, [is_stable_list[-1]]])
            growthrate_list = np.concatenate([growthrate_list, [growthrate_list[-1]]])
            if Tp_list is not None:
                Tp_list = np.concatenate([Tp_list, [Tp_list[-1]]])
            if ax_list is not None:
                ax_list = np.concatenate([ax_list, [ax_list[-1]]])
            if az_list is not None:
                az_list = np.concatenate([az_list, [az_list[-1]]])
        elif type == 'both':
            mu_list, response_list, is_stable_list, growthrate_list, Tp_list, ax_list, az_list = \
                merge_branch(mu_list, response_list, is_stable_list, growthrate_list, Tp_list, ax_list, az_list, type='first')
            mu_list, response_list, is_stable_list, growthrate_list, Tp_list, ax_list, az_list = \
                merge_branch(mu_list, response_list, is_stable_list, growthrate_list, Tp_list, ax_list, az_list, type='last')
    return mu_list,response_list,is_stable_list,growthrate_list,Tp_list,ax_list,az_list

def reduce_branch(mu_list,response_list,is_stable_list,growthrate_list,Tp_list,ax_list,az_list,
           first_idx=0, last_idx=None):
    return mu_list[first_idx:last_idx], response_list[first_idx:last_idx],\
            is_stable_list[first_idx:last_idx], growthrate_list[first_idx:last_idx],\
            Tp_list[first_idx:last_idx] if Tp_list is not None else None,\
            ax_list[first_idx:last_idx] if ax_list is not None else None,\
            az_list[first_idx:last_idx] if az_list is not None else None

def remove_base_cond_flux(mu_list,response_list,is_stable_list,growthrate_list,Tp_list,ax_list,az_list):
    return mu_list,np.array(response_list)-1.0,is_stable_list,growthrate_list,Tp_list,ax_list,az_list


def data2csv(filename,mu_list,response_list,is_stable_list,growthrate_list,
             Tp_list=None, ax_list=None, az_list=None,
             mu='mu',response='response'):
    dir = os.path.dirname(filename)
    if dir:  # only call makedirs if there's actually a directory part
        os.makedirs(dir, exist_ok=True)
    response_array = np.array(response_list)
    if response_array.ndim == 2:  # periodic orbit case
        df = pd.DataFrame({
            f'{mu}': mu_list,
            f'{response}_min': response_array[:, 0],
            f'{response}_max': response_array[:, 1],
            'is_stable': is_stable_list,
            'growthrate_real': np.real(growthrate_list),
            'growthrate_imag': np.imag(growthrate_list)
        })
    else:
        df = pd.DataFrame({
            f'{mu}': mu_list,
            f'{response}': response_list,
            'is_stable': is_stable_list,
            'growthrate_real': np.real(growthrate_list),
            'growthrate_imag': np.imag(growthrate_list)
        })
    if Tp_list is not None and len(Tp_list) > 0:
        df['Tp'] = Tp_list
    if ax_list is not None and len(ax_list) > 0:
        df['ax'] = ax_list
    if az_list is not None and len(az_list) > 0:
        df['az'] = az_list
    df.to_csv(filename, index=False)

def csv2data(filename):
    df = pd.read_csv(filename)
    mu_col, response_col = df.columns[0], df.columns[1]
    
    mu_list = df[mu_col].tolist()
    is_stable_list = df['is_stable'].tolist()
    growthrate_list = (df['growthrate_real'] + 1j * df['growthrate_imag']).tolist()

    # detect periodic orbit case
    if response_col.endswith('_min'):
        response_name = response_col[:-4]  # strip '_min'
        response_list = np.column_stack([
            df[f'{response_name}_min'].to_numpy(),
            df[f'{response_name}_max'].to_numpy()
        ])
    else:
        response_name = response_col
        response_list = df[response_col].tolist()
    
    Tp_list = df['Tp'].tolist() if 'Tp' in df.columns else None
    ax_list = df['ax'].tolist() if 'ax' in df.columns else None
    az_list = df['az'].tolist() if 'az' in df.columns else None
    
    return mu_list, response_list, is_stable_list, growthrate_list, Tp_list, ax_list, az_list

#############################################
######### Plot bifurcation diagram ##########

def plot_branch(figax, bifurcation_info, label=None,
                lw=1., color='k', unstableSolStyle="dashed", alpha=1.):
    mu, response, is_stable, growthrate, Tp, ax, az = bifurcation_info
    response = np.array(response)
    if response.ndim==2: # solution is periodic orbit, plotting min. and max. lines
        plot_bifurcation(figax, mu=mu, response=response[:,0], is_stable=is_stable, label=label,
                        color=color, unstableSolStyle=unstableSolStyle, lw=lw, alpha=alpha)
        plot_bifurcation(figax, mu=mu, response=response[:,1], is_stable=is_stable, label=None,
                        color=color, unstableSolStyle=unstableSolStyle, lw=lw, alpha=alpha)
    else:
        plot_bifurcation(figax, mu=mu, response=response, is_stable=is_stable, label=label,
                        color=color, unstableSolStyle=unstableSolStyle, lw=lw, alpha=alpha)

def plot_bifurcation(ax,mu,response,is_stable,label=None,
                     color='k',unstableSolStyle="dashed",lw=2,alpha=1):
    segments = separate_into_segments(mu,response,is_stable)
    # print(segments)
    if len(segments)>1:
        # this bifurcation has both stable and unstable bifurcations
        # add label to first stable bifurcation
        addedLabel = False
        for i, segment in enumerate(segments):
            if segment['is_stable']:
                if addedLabel:
                    ax.plot(segment['mu'], segment['response'],lw=lw,alpha=alpha,label=None, linestyle='-', color=color, zorder=1)
                else:
                    ax.plot(segment['mu'], segment['response'],lw=lw,alpha=alpha,label=label, linestyle='-', color=color, zorder=1)
                    addedLabel = True
            else:
                ax.plot(segment['mu'], segment['response'],lw=lw,alpha=alpha,label=None, linestyle=unstableSolStyle, color=color, zorder=1)
    else:
        # this bifurcation has only one stable or unstable on the whole bifurcation
        # add label to anywhere 
        for segment in segments:
            if segment['is_stable']:
                ax.plot(segment['mu'], segment['response'],lw=lw,alpha=alpha,label=label, linestyle='-', color=color, zorder=1)
            else:
                ax.plot(segment['mu'], segment['response'],lw=lw,alpha=alpha,label=label, linestyle=unstableSolStyle, color=color, zorder=1)

def separate_into_segments(mu,response,is_stable):
    # Initialize lists for segments
    segments = []
    current_segment = {'mu': [], 'response': [], 'is_stable': is_stable[0]}
    # Separate the data into segments based on the condition
    for i in range(len(mu)):
        if is_stable[i] == current_segment['is_stable']:
            current_segment['mu'].append(mu[i])
            current_segment['response'].append(response[i])
        else:
            segments.append(current_segment)
            current_segment = {'mu': [mu[i-1], mu[i]], 'response': [response[i-1], response[i]], 'is_stable': is_stable[i]}
    segments.append(current_segment)
    return segments

##########################################################
#####  Plot scatters representing bifurcation points #####

def add_scatter(figax, bifurcation_info, pointType='stabilityChange',
                facecolors='None', marker='o', edgecolors='k', 
                markerHopf='d', hopf_edgecolors='k', s=20., lw=1.):
    """ pointType: 'first'/'last'/'stabilityChange'/'Hopf'/'Direct' """
    mu, response, is_stable, growthrate, Tp, ax, az = bifurcation_info
    response = np.array(response)

    def is_Hopf(a:complex, b:complex):
        return abs(a.imag) > 1e-6 or abs(b.imag) > 1e-6

    if response.ndim==2: # solution is periodic orbit
        if pointType == 'first':
            figax.scatter(mu[0], response[0,0],
                        marker=marker, s=s, linewidths=lw, 
                        facecolors=facecolors, edgecolors=edgecolors, 
                        zorder=2)
        elif pointType == 'last':
            figax.scatter(mu[-1], response[-1,0],
                        marker=marker, s=s, linewidths=lw, 
                        facecolors=facecolors, edgecolors=edgecolors, 
                        zorder=2)
        else:
            pass
    else:
        if pointType == 'first':
            figax.scatter(mu[0], response[0],
                        marker=marker, s=s, linewidths=lw, 
                        facecolors=facecolors, edgecolors=edgecolors, 
                        zorder=2)
        elif pointType == 'last':
            figax.scatter(mu[-1], response[-1],
                        marker=marker, s=s, linewidths=lw, 
                        facecolors=facecolors, edgecolors=edgecolors, 
                        zorder=2)
        else: # plot the direct and Hopf bifurcation points
            for i in range(len(mu)-1):
                if is_stable[i] != is_stable[i+1]: # break stability: Stable <-> Unstable
                    if (pointType == 'Hopf' or pointType == 'stabilityChange') and is_Hopf(growthrate[i],growthrate[i+1]):
                        figax.scatter(mu[i],response[i],marker=markerHopf,s=s,linewidths=lw, facecolors=facecolors, edgecolors=hopf_edgecolors, zorder=2)
                    elif pointType == 'Direct' or pointType == 'stabilityChange':
                        figax.scatter(mu[i],response[i],marker=marker,s=s,linewidths=lw, facecolors=facecolors, edgecolors=edgecolors, zorder=2)

