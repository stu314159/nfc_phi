## validate.py
"""
Python script that will compute the point-wise difference
between a simulation result computed from MATLAB with that
computed using turbineChannel.  

Required: same geometry, same discretization, same simulation parameters,
and same number of time steps for the comparison (obviously...)

The gold_standard.mat file must be in the working directory - this file contains the velocity components and pressure
as computed in the MATLAB-based code that is taken to be correct.

The params.lbm parameter file must be in the folder so that simulation geometric and scaling parameters
associated with the simulation can be gathered and applied to the raw output data.

The user must input the "dump number" for the turbineChannel data so that it can extract the binary simulation data.

The script output is a vtk-file for velocity and pressure component differences.

"""

import numpy as np
import math
import scipy.io
import argparse
from vtkHelper import saveVelocityAndPressureVTK_binary as writeVTK

# the turbineChannel dump number is a required input - use argparse to capture.

parser = argparse.ArgumentParser(prog='validate.py',description='validate turbineChannel output')
parser.add_argument('dump_num',type=int)
# parse input arguments
args = parser.parse_args()
# assign to required variables
dump_num = args.dump_num


# get the data from MATLAB into the Python environment
mat_data = scipy.io.loadmat('gold_standard.mat')
ux_m = np.asarray(list((mat_data['ux_h']).flatten()))
uy_m = np.asarray(list((mat_data['uy_h']).flatten()))
uz_m = np.asarray(list((mat_data['uz_h']).flatten()))
pressure_m = np.asarray(list((mat_data['pressure_h']).flatten()))

# Information about the LBM run that produced the data - I should get this from params.lbm

# Read data from params.lbm
input_file_name = 'params.lbm'
input_data = open(input_file_name,'r')
latticeType = int(input_data.readline())
Num_ts = int(input_data.readline())
ts_rep_freq = int(input_data.readline())
Warmup_ts = int(input_data.readline())
plot_freq = int(input_data.readline())
Cs = float(input_data.readline())
rho_lbm = float(input_data.readline())
u_lbm = float(input_data.readline())
omega = float(input_data.readline())
Nx = int(input_data.readline())
Ny = int(input_data.readline())
Nz = int(input_data.readline())
Lx_p = float(input_data.readline())
Ly_p = float(input_data.readline())
Lz_p = float(input_data.readline())
t_conv_fact = float(input_data.readline())
l_conv_fact = float(input_data.readline())
p_conv_fact = float(input_data.readline())

input_data.close()

u_conv_fact = l_conv_fact/t_conv_fact;

nnodes = Nx*Ny*Nz

# compute geometric data only once
x = np.linspace(0.,Lx_p,Nx).astype(np.float32);
y = np.linspace(0.,Ly_p,Ny).astype(np.float32);
z = np.linspace(0.,Lz_p,Nz).astype(np.float32);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);

XX = np.reshape(X,numEl)
YY = np.reshape(Y,numEl)
ZZ = np.reshape(Z,numEl)


# load the turbineChannel data from *.b_dat files
rho_fn = 'density'+str(dump_num)+'.b_dat'
ux_fn = 'ux'+str(dump_num)+'.b_dat'
uy_fn = 'uy'+str(dump_num)+'.b_dat'
uz_fn = 'uz'+str(dump_num)+'.b_dat'
ux = np.fromfile(ux_fn,dtype=np.float32)
uy = np.fromfile(uy_fn,dtype=np.float32)
uz = np.fromfile(uz_fn,dtype=np.float32)
pressure = np.fromfile(rho_fn,dtype=np.float32)
    
# convert velocity to physical units
ux *= u_conv_fact
uy *= u_conv_fact
uz *= u_conv_fact
pressure *= p_conv_fact

# compute differences
ux_d = ux - ux_m
uy_d = uy - uy_m
uz_d = uz - uz_m
pressure_d = pressure - pressure_m

# write the results to a vtk file
outfilename = 'velocityAndPressureDifferences.vtk'
dims = (Nx,Ny,Nz)
writeVTK(pressure_d,ux_d,uy_d,uz_d,XX,YY,ZZ,outfilename,dims)

# pray that things work better next time.