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
  with nc.Dataset(file_path, 'r') as ncf:
    x = ncf.variables['x'][:]  # grid coordinates
    z = ncf.variables['z'][:]
    return x.ravel(), z.ravel()

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

def plot_sample(path, name,color_surface = 'PuOr_r',color_S = 'r', plot_mean_profiles=True, symmetric=False):
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
    fig.savefig(name, dpi=300, bbox_inches='tight')
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

################# 
#### Loading bifurcation data from HDF5 files
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
    
def get_bifurcation(folder_path, mu = "mu", response="Sh", is_po=False):
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
    if is_po:
        response_list = np.zeros((len(mu_list),2))
    else:
        response_list = df0[response] if response in df0.columns else np.zeros(len(mu_list))
    
    dir_list = df0["dir"]
    dir_list = [s.strip() for s in dir_list]

    stability = np.zeros(len(mu_list), dtype=bool)
    growthrate = np.zeros(len(mu_list), dtype=np.complex128)
    file_path = "/time-dependent/full_solution/full_solution_s1.h5" if os.path.isdir(folder_path+"search-0/time-dependent/") else "/flowproperties/flowproperties_s1.h5"
    
    for element, dir_element in enumerate(dir_list):
        if is_po:
            with h5py.File(folder_path+dir_element+file_path, mode='r') as file:
                response_data = file['tasks'][response][:]
            response_list[element,0] = min(np.array(response_data[:,0,0].real))
            response_list[element,1] = max(np.array(response_data[:,0,0].real))
        
        stability[element], growthrate[element] = check_stability(folder_path+dir_element+"stability/growthrate.mtx")
        
        print(f"        Processed element {element+1} in {dir_element}: {mu}={mu_list[element]}, {response}={response_list[element]}, stability={'unstable' if stability[element] else 'stable'}, growthrate={growthrate[element]}")

    return mu_list, response_list, stability, growthrate

def separate_into_segments(mu,refes,isunstable):
    # Initialize lists for segments
    segments = []
    current_segment = {'mu': [], 'refes': [], 'isunstable': isunstable[0]}
    # Separate the data into segments based on the condition
    for i in range(len(mu)):
        if isunstable[i] == current_segment['isunstable']:
            current_segment['mu'].append(mu[i])
            current_segment['refes'].append(refes[i])
        else:
            segments.append(current_segment)
            current_segment = {'mu': [mu[i-1], mu[i]], 'refes': [refes[i-1], refes[i]], 'isunstable': isunstable[i]}
    segments.append(current_segment)
    return segments

def plot_bifurcation(ax,mu,refes,isunstable,label=None,color='k',unstableline="dashed",lw=2,alpha=1):
    segments = separate_into_segments(mu,refes,isunstable)
    # print(segments)
    if len(segments)>1:
        # this bifurcation has both stable and unstable bifurcations
        # add label to first stable bifurcation
        addedLabel = False
        for i, segment in enumerate(segments):
            if segment['isunstable']:
                ax.plot(segment['mu'], segment['refes'],lw=lw,alpha=alpha,label=None, linestyle=unstableline, color=color, zorder=1)
            else:
                if addedLabel:
                    ax.plot(segment['mu'], segment['refes'],lw=lw,alpha=alpha,label=None, linestyle='-', color=color, zorder=1)
                else:
                    ax.plot(segment['mu'], segment['refes'],lw=lw,alpha=alpha,label=label, linestyle='-', color=color, zorder=1)
                    addedLabel = True
    else:
        # this bifurcation has only one stable or unstable on the whole bifurcation
        # add label to anywhere 
        for segment in segments:
            if segment['isunstable']:
                ax.plot(segment['mu'], segment['refes'],lw=lw,alpha=alpha,label=label, linestyle=unstableline, color=color, zorder=1)
            else:
                ax.plot(segment['mu'], segment['refes'],lw=lw,alpha=alpha,label=label, linestyle='-', color=color, zorder=1)

def bifurcation_diagram(data_path, mu = "mu", response="Sh"):
    data = {}
    for key, branch_info in data_path.items():
        mu_fulllist, response_fulllist, stability_fulllist, growthrate_fulllist = [], [], [], []
        # Extract color metadata safely
        color = branch_info.get('color', '#000000')
        alpha = branch_info.get('alpha', 1.0)
        isPO = branch_info.get('isPO', False)
        showLabel = branch_info.get('showLabel', True)
        merge_start = branch_info.get('merge_start', False)
        merge_end = branch_info.get('merge_end', False)
        scatter_first_point = branch_info.get('scatter_first_point', False)
        scatter_end_point = branch_info.get('scatter_end_point', False)
        scatter_stability_change = branch_info.get('scatter_stability_change', False)
        for subdir, subsubdir in branch_info.items():
            # Skip metadata keys during directory parsing
            if subdir == 'color' \
                or subdir == 'alpha' \
                or subdir == 'isPO' \
                or subdir == 'showLabel' \
                or subdir == 'merge_start' or subdir == 'merge_end' \
                or subdir == 'scatter_first_point' or subdir == 'scatter_end_point' or subdir == 'scatter_stability_change':
                continue
            print(f"Loading {key} data from {subdir}...")
            # mu_list, response_list = [], []
            for dir, start, end, upward in subsubdir:
                full_path = os.path.join(subdir, dir)
                print(f"    Loading continuation data from {full_path}...")
                mu_list, response_list, stability_list, growthrate_list = get_bifurcation(full_path, mu=mu, response=response, is_po=isPO)

                mu_slice = mu_list[start:end]
                response_slice = response_list[start:end]
                stability_slice = stability_list[start:end]
                growthrate_slice = growthrate_list[start:end]

                if not upward:
                    mu_slice = mu_list[start:end][::-1]
                    response_slice = response_list[start:end][::-1]
                    stability_slice = stability_list[start:end][::-1]
                    growthrate_slice = growthrate_list[start:end][::-1]

                mu_fulllist.extend(mu_slice)
                response_fulllist.extend(response_slice)
                stability_fulllist.extend(stability_slice)
                growthrate_fulllist.extend(growthrate_slice)
            # print(f"loaded mu_fulllist={mu_fulllist}")
            # ...
        mu_fulllist = np.array(mu_fulllist)
        response_fulllist = np.array(response_fulllist)
        stability_fulllist = np.array(stability_fulllist)
        growthrate_fulllist = np.array(growthrate_fulllist)
        
        # add one more point to connect two min. and max. branch of PO
        if branch_info.get('isPO', False) and branch_info.get('merge_start', False):
            mu_fulllist = np.concatenate([[mu_fulllist[0]-(mu_fulllist[1]-mu_fulllist[0])], mu_fulllist])
            val = 0.5 * (response_fulllist[0, 0] + response_fulllist[0, 1])
            response_fulllist = np.concatenate([[[val, val]],response_fulllist], axis=0)
            stability_fulllist = np.concatenate([[stability_fulllist[0]], stability_fulllist])
            growthrate_fulllist = np.concatenate([[growthrate_fulllist[0]], growthrate_fulllist])
        if branch_info.get('isPO', False) and branch_info.get('merge_end', False):
            mu_fulllist = np.concatenate([mu_fulllist, [mu_fulllist[-1]-(mu_fulllist[-2]-mu_fulllist[-1])]])
            val = 0.5 * (response_fulllist[-1, 0] + response_fulllist[-1, 1])
            response_fulllist = np.concatenate([response_fulllist, [[val, val]]],axis=0)
            stability_fulllist = np.concatenate([stability_fulllist, [stability_fulllist[-1]]])
            growthrate_fulllist = np.concatenate([growthrate_fulllist, [growthrate_fulllist[-1]]])
                
        # collect dataset
        data[key] = {
            'mu': mu_fulllist,
            'response': response_fulllist,
            'Stability': stability_fulllist,
            'Growthrate': growthrate_fulllist,
            'color': color,
            'alpha': alpha,
            'isPO': isPO,
            'showLabel': showLabel,
            'merge_start': merge_start,
            'merge_end': merge_end,
            'scatter_first_point': scatter_first_point,
            'scatter_end_point': scatter_end_point,
            'scatter_stability_change': scatter_stability_change
        }
    return data

def plot_bifurcation_diagram(ax, data, lw=1):
    for key, branch_info in data.items():
        # print(f"Plotting bifurcation for {key}...")
        if branch_info.get('isPO', False):
            response_array = np.array(data[key]['response'])
            plot_bifurcation(ax,
                            mu=data[key]['mu'], 
                            refes=response_array[:,0], # min value
                            isunstable=data[key]['Stability'],
                            label=f'{key}' if branch_info.get('showLabel', False) else None,
                            color=data[key]['color'],
                            unstableline="dashed",
                            lw=lw,
                            alpha=data[key]['alpha'])
            plot_bifurcation(ax,
                            mu=data[key]['mu'], 
                            refes=response_array[:,1], # max value
                            isunstable=data[key]['Stability'],
                            label= None,
                            color=data[key]['color'],
                            unstableline="dashed",
                            lw=lw,
                            alpha=data[key]['alpha'])
        else:    
            plot_bifurcation(ax,
                            mu=data[key]['mu'], 
                            refes=data[key]['response'], 
                            isunstable=data[key]['Stability'],
                            label=f'{key}' if branch_info.get('showLabel', False) else None,
                            color=data[key]['color'],
                            unstableline="dashed",
                            lw=lw,
                            alpha=data[key]['alpha'])


def add_scatter_bifurcation_points(ax,
                                   mu,
                                   refes,
                                   stability,
                                   leading_eigen,
                                   marker_hopf='d',
                                   marker='*',
                                   label=None,
                                   color='none',
                                   hopf_edgecolors='k',
                                   edgecolors='k',
                                   s=20,
                                   lw=1.):
    def break_stability_eigen(a, b):
        if a.real*b.real < 0:
            return True
        return False
    def break_stability(a, b):
        if (a==True and b==False) or (a==False and b==True):
            # print('break stability from ',a, ' to', b)
            return True
        return False
    def Hopf_bifurcation(a, b):
        if (a.imag>0 and abs(a.imag)>1e-6) or (b.imag>0 and abs(b.imag)>1e-6):
            return True
        return False
    for i in range(len(mu)-1):
        if break_stability(stability[i],stability[i+1]):
            if Hopf_bifurcation(leading_eigen[i],leading_eigen[i+1]):
                ax.scatter(mu[i+1],refes[i+1],marker=marker_hopf,s=s,linewidths=lw, color=color, edgecolors=hopf_edgecolors, zorder=2)
            else:
                ax.scatter(mu[i+1],refes[i+1],marker=marker,s=s,linewidths=lw, color=color, edgecolors=edgecolors, zorder=2)

def add_scatter_first_point(ax,
                            mu,
                            refes,
                            marker='*',
                            label=None,
                            color='none', 
                            edgecolors='k',
                            s=20,
                            lw=1.):
    ax.scatter(mu[0],
               refes[0],
               marker=marker,
               s=s, 
               linewidths=lw, 
               facecolors=color, 
               edgecolors=edgecolors, 
               zorder=2)

def plot_scatter_for_stability_change(ax, data, marker_hopf='d', color_hopf='k', marker='o', color='k', s=20, lw=1.):
    for key, branch_info in data.items():
        if branch_info.get('scatter_stability_change', False):
            add_scatter_bifurcation_points(ax, 
                                        mu=data[key]['mu'], 
                                        refes=data[key]['response'],
                                        stability=data[key]['Stability'], 
                                        leading_eigen=data[key]['Growthrate'],
                                        marker_hopf=marker_hopf, 
                                        marker=marker,
                                        color='none',
                                        hopf_edgecolors=color_hopf, 
                                        edgecolors='k', 
                                        s=s,
                                            lw=lw)

def plot_scatter_for_first_point(ax, data, marker='o', color='k', s=20, lw=1.):
    for key, branch_info in data.items():
        if branch_info.get('scatter_first_point', False):
            add_scatter_first_point(ax, 
                                    mu=data[key]['mu'], 
                                    refes=data[key]['response'],
                                    marker=marker, 
                                    color='none', 
                                    edgecolors=color, 
                                    s=s,
                                    lw=lw)
        if branch_info.get('scatter_end_point', False):
            add_scatter_first_point(ax, 
                                    mu=data[key]['mu'][::-1], 
                                    refes=data[key]['response'][::-1],
                                    marker=marker, 
                                    color='none', 
                                    edgecolors=color, 
                                    s=s,
                                    lw=lw)
        