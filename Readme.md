This is data and post-processing code for the paper "Exact Coherent Structures of Salt-Finger Convection: Layering and Shear Generation".

Because the storage limit, we do not provide data for making supplementary video S5 (`Fig3a_video.ipynb`). To avoid Windows path-length limit (Error 0x80010135) when extracting zip file, we highly recommend downloading data via Github clone.

```bash
# clone data [run this in WSL terminal, via VScode]
git clone https://github.com/duc-nguyen1410/Exact_Coherent_Structures_of_Salt_Finger_Convection_Layering_and_Shear_Generation.git
```

<!-- Also, you need to install Git LFS (used to save large files on Github) to load full large files from github to your computer after cloning.
```bash
# step 1: clone data, large files will be saved by pointers instead of real files. [run this in WSL terminal, via VScode]
git clone https://github.com/duc-nguyen1410/Exact_Coherent_Structures_of_Salt_Finger_Convection_Layering_and_Shear_Generation.git

# step 2: get full data of large files [run this in Git bash, via VScode]
git lfs install # install if needed
cd Exact_Coherent_Structures_of_Salt_Finger_Convection_Layering_and_Shear_Generation
git lfs pull  # this downloads the actual large files
``` -->

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
Here, columns from left to right present the continuous parameter `Lx`, the response in bifurcation diagram `Sh`, the stability of solution (`True`: unstable, `False`: stable), real and imaginary parts of fastest growth rate (`growthrate_real`,`growthrate_imag`) and additional parameters of traveling waves and relative periodic orbits: time period `Tp`, horizontal translation `ax` normalized by domain width $L_x$, vertical translation `az` normalized by domain height $L_z$.