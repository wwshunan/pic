import numpy
from scipy import stats
from mayavi import mlab

def fill_up_array(output_arr,source_arr,nx,posn):
    for i in range(nx):
        output_arr[i] = source_arr[i][posn]

filename = raw_input("Please enter the filename: ")
print "you entered", filename

#now read in values from stream.out file to these arrays
full = numpy.loadtxt(filename)
print full.shape[0]

# Find number of lines in data-file
NX = full.shape[0] 

NDIM = 4  # number of columns in data file, 4 for current example

#initialise arrays:
xx = numpy.zeros((NX))
yy = numpy.zeros((NX))
zz = numpy.zeros((NX))
VV = numpy.zeros((NX))

fy = numpy.reshape(full, (NX,NDIM))

fill_up_array(xx,fy,NX,0)
fill_up_array(yy,fy,NX,1)
fill_up_array(zz,fy,NX,2)
fill_up_array(VV,fy,NX,3)

# Plot scatter with mayavi
figure = mlab.figure('Fix-Line')
pts = mlab.points3d(xx, yy, zz, VV,scale_factor=0.25,resolution=5,scale_mode='none')
mlab.axes(xlabel="X",ylabel="Y",zlabel="Z")
mlab.show()

