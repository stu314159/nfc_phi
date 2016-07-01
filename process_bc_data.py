#!/home/sblair/anaconda/bin/python
###!/p/home/sblair/anaconda/bin/python
"""
Data processing script for LBM binary data files holding BC data.

Call this script from the same directory where the LBM input file
resides as well as all of the *.b_dat files from the simulation.

No arguments are needed.

Usage:

>>python ./process_bc_data.py

When done, you should have a set of *.vtk files containing
the snl, inl and onl boundary data files.

launch with 3 MPI processes; one for each boundary condition array

"""


from mpi4py import MPI
import numpy as np
import math
#from vtkHelper import saveVelocityAndPressureVTK_binary as writeVTK

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import sys
sys.path.insert(1,'.')
from vtkHelper import saveScalarStructuredGridVTK_binary as writeVTK

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

#u_conv_fact = l_conv_fact/t_conv_fact;
u_conv_fact = t_conv_fact/l_conv_fact;
if(rank == 0):
   print("u_conv_fact = %15.15f \n"%u_conv_fact)
nnodes = Nx*Ny*Nz

# compute geometric data only once
x = np.linspace(0.,Lx_p,Nx).astype(np.float64);
y = np.linspace(0.,Ly_p,Ny).astype(np.float64);
z = np.linspace(0.,Lz_p,Nz).astype(np.float64);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);

XX = np.reshape(X,numEl)
YY = np.reshape(Y,numEl)
ZZ = np.reshape(Z,numEl)

if rank==0:
  dat_name = 'snl';
elif rank==1:
  dat_name = 'inl';
elif rank==2:
  dat_name = 'onl';

in_fn = dat_name + '.b_dat'
out_fn = dat_name + '.vtk'
my_dat = np.fromfile(in_fn,dtype=np.int32)

dims = (Nx,Ny,Nz)
writeVTK(my_dat,dat_name,XX,YY,ZZ,out_fn,dims)


