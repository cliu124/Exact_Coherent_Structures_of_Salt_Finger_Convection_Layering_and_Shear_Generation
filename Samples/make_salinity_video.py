"""
Plot 2D cartesian snapshots.

Usage:
    plot_snapshots.py <files>... [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./frames]

"""

import h5py
import numpy as np
from mpi4py import MPI
import matplotlib
from matplotlib.gridspec import GridSpec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dedalus.extras import plot_tools
from matplotlib.animation import FuncAnimation

# Enable LaTeX for all text
import matplotlib as mpl
mpl.rcParams.update({
    'text.usetex': True,           # Enable LaTeX rendering
    'font.family': 'serif',        # Use serif fonts
    'font.size': 14,               # Base font size
    'mathtext.fontset': 'stix',    # STIX font for math symbols
})

def main(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""
    # Plot writes
    with h5py.File(filename, mode='r') as file:
        x_key = [k for k in file['scales'].keys() if k.startswith('x_hash')][0]
        x = file[f'scales/{x_key}'][:]        # usually shape (Nz, 1, 1)
        x = np.squeeze(x)
        z_key = [k for k in file['scales'].keys() if k.startswith('z_hash')][0]
        z = file[f'scales/{z_key}'][:]        # usually shape (Nz, 1, 1)
        z = np.squeeze(z) 
        for index in range(start, start+count):
            time = file['scales/sim_time']
            u = file['tasks']['u'][0]
            w = file['tasks']['u'][1]
            t = file['tasks']['te']
            s = file['tasks']['sa']
        u = np.array(u).real
        w = np.array(w).real
        t = np.array(t).real
        s = np.array(s).real
        time = np.array(time)

        nframes = np.shape(time)[0]

        if MPI.COMM_WORLD.rank == 0:
            print("Data loaded successfully.")
            print(f"(Nx,Nz) = ({np.shape(x)[0]},{np.shape(z)[0]})")
            print("nframes = ", nframes)

            N = 5  # Number of periods to repeat

            #####
            fig = plt.figure(figsize=(4, 3))
            gs = GridSpec(1, 1, figure=fig)
            ax1 = fig.add_subplot(gs[0,0])

            # Initial title
            # title = ax1.set_title(r'$t = 0.00 T_p$')

            ### plot salinity
            smin = np.min(s)
            smax = np.max(s)
            max_abs = max(abs(smin), abs(smax))
            print("min s = ", smin, ", max s = ",smax)
            # ax1.set_title(f'time={time[0]:.2f}')
            X_ax2, Z_ax2 = np.meshgrid(x, z)
            p2 = ax1.pcolormesh(X_ax2,Z_ax2,s[0].T,vmin=-max_abs,vmax=max_abs, shading='auto', cmap = 'PuOr_r')
            fig.colorbar(p2, ax=ax1, label='$s$', ticks=np.linspace(-max_abs, max_abs, 5))
            ax1.set_xlabel('$x$')
            ax1.set_ylabel('$z$')
            ax1.set_ylim([0,1])
            ax1.set_yticks([0,0.5,1],['0','0.5','1'])
            ax1.set_aspect('equal')
            plt.tight_layout()
            
            def update(frame):
                # p2.set_array(s[frame].T)
                p2.set_array(s[frame % nframes].T)

                # Calculate time in multiples of period Tp
                t_in_Tp = frame / nframes
                # ax1.set_title(rf'$t = {t_in_Tp:.2f} T_p$')

            ani = FuncAnimation(fig, 
                                update, 
                                frames=nframes * N,
                                interval=40, 
                                blit=False, 
                                repeat=True)
            ani.save(output / 'Salinity.mp4', 
                     writer='ffmpeg', 
                     fps=24 if nframes==100 else 48, 
                     dpi=300, 
                     bitrate=2000)

if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path)