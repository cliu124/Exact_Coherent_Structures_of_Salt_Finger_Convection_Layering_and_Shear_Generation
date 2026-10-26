This is data and post-processing code for the paper "Exact Coherent Structures of Salt-Finger Convection: Layering and Shear Generation".

Due to the storage limit, we do not provide data for making supplementary video S5 (`Fig3a_video.ipynb`). To avoid the Windows path-length limit (Error 0x80010135) when extracting the zip file, we highly recommend downloading the data via GitHub clone.

```bash
# clone data [run this in WSL terminal, via VSCode]
git clone https://github.com/duc-nguyen1410/Exact_Coherent_Structures_of_Salt_Finger_Convection_Layering_and_Shear_Generation.git
```
Alternatively, extract the zip file using third-party tools (e.g., 7-Zip, WinRAR, ...).

Running Jupyter notebooks requires some libraries:
```bash
# create conda environment with libraries within anaconda
conda create -n test-env python=3.12 numpy matplotlib scipy jupyter ipykernel pandas netCDF4 h5py
conda activate test-env

# or install libraries via pip
pip install numpy matplotlib scipy pandas jupyter ipykernel netCDF4 h5py

## please select exact python kernel when running notebook. 
# for example, here select kernel "test-env", or use your available python kernel
```

We save information of each bifurcation in `csv` files with the following format:
```
Lx,Sh,is_stable,growthrate_real,growthrate_imag,Tp,ax,az
...
```
Here, columns from left to right represent the continuous parameter `Lx`, the response in the bifurcation diagram `Sh`, the stability of the solution (`True`: stable, `False`: unstable), the real and imaginary parts of the fastest growth rate (`growthrate_real`, `growthrate_imag`), and additional parameters of traveling waves and relative periodic orbits: time period `Tp`, horizontal translation `ax` normalized by domain width $L_x$, vertical translation `az` normalized by domain height $L_z$.