# pySeismicFMM

## Installation

`3DFMM` needs to be compiled. This process may require installation of additional modules:

```
python3 3Dsetup.py build_ext --inplace
```

## Example usage (example.py)
```
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
```

# Notes
1. Binary files used in example are form 3D model of Poland [Grad, M., Polkowski, M., Ostaficzuk, S.R., High-resolution 3D seismic model of the crustal and uppermost mantle structure in Poland, Tectonophysics, Vol. 666, pp. 188 - 210]
2. Instead of reading binary files models can by set using `SetVelocityModel()` and other similar functions
3. Seismic source can by set at any geographic coordinate within the model
4. Calculated travel time can be red at any geographic coordinate within the model
5. Calculating requires setting number of time steps (grid cells). This will be fixed to provide more friendly solution

# More information
1. This code was described in my PhD thesis (in polish): [http://marcinpolkowski.com/files/dyplom_phd.pdf]
2. It was presented on few posters in Vienna (EGU) and San Francisco (AGU). Poster PDF are available on my website: [http://marcinpolkowski.com/]
