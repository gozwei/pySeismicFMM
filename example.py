import numpy
import datetime
from pySeismicFMM3D import SeismicFMM3D
from geopy.distance import great_circle


# Init SeismicFMM3D module:
myFMM = SeismicFMM3D()

# Specify used model size (lat, lon, z):
myFMM.SetModelSize(631, 536, 6261)

# Set grid cell size (dlat, dlon, dz):
myFMM.SetGridSize(numpy.deg2rad(0.01), numpy.deg2rad(0.02), 10)

# Read binary files with velocity model, geological model, lat, lon and z (H)
myFMM.ReadVelocityModel("/Users/marcinpolkowski//FMM/MODELVPF.bin", numpy.single)
myFMM.ReadModel("/Users/marcinpolkowski//FMM/MODEL.bin", numpy.uint8)
myFMM.ReadLatVector("/Users/marcinpolkowski//FMM/lat.bin", numpy.double)
myFMM.ReadLonVector("/Users/marcinpolkowski//FMM/lon.bin", numpy.double)
myFMM.ReadHVector("/Users/marcinpolkowski//FMM/H.bin", numpy.double)

# Allocate memory:
myFMM.CreateCalculationVariables()

# Set seismic source:
myFMM.SetSourceGeo(50.0717,19.3396,-1000)

# Perform calculations (number of time steps)
myFMM.Do(400000000)

# Get travel time at any point within the model:
print(myFMM.GetTimeGeo(51.232,21.0068,110))