import numpy
import FMM3D
import matplotlib.pyplot as plt
import time
from geopy.distance import great_circle

class SeismicFMM3D():
	def __init__(self):
		print('FMM 3D initiated')

	def SetModelSize(self, lat, lon, z):        
		self.size_lat = lat
		self.size_lon = lon
		self.size_z = z
		
	def SetGridSize(self, dlat, dlon, dz):
		self.dlat = dlat;
		self.dlon = dlon;
		self.dz = dz;

	def CreateCalculationVariables(self):
		self.time = numpy.zeros(self.model_velocity.shape, dtype=numpy.single)
		self.accepted = numpy.zeros(self.model_velocity.shape, dtype=numpy.bool)
		self.trace = numpy.zeros(self.model_velocity.shape, dtype=numpy.int32)

	def __ReadBinaryFile(file, size, dtype):
		return numpy.fromfile(file, dtype=dtype).reshape(size,  order='F')
		   
	def ReadVelocityModel(self, file, dtype=numpy.single):
		print('Reading velocity model from {0:s}..'.format(file), end='\t', flush=True)
		self.model_velocity = SeismicFMM3D.__ReadBinaryFile(file, (self.size_lat, self.size_lon, self.size_z), dtype)
		print('done.', flush=True)        
	def SetVelocityModel(self, M):
		self.model_velocity = M
		
	def ReadModel(self, file, dtype=numpy.single):
		print('Reading model from {0:s}..'.format(file), end='\t', flush=True)
		self.model = SeismicFMM3D.__ReadBinaryFile(file, (self.size_lat, self.size_lon, self.size_z), dtype)
		print('done.', flush=True)        
	def SetModel(self, M):
		self.model = M

	def ReadTime(self, file, dtype=numpy.single):
		print('Reading time from {0:s}..'.format(file), end='\t', flush=True)
		self.time = SeismicFMM3D.__ReadBinaryFile(file, (self.size_lat, self.size_lon, self.size_z), dtype)
		print('done.', flush=True)        
	def SetTime(self, T):
		self.time = T
	
	def ReadLatVector(self, file, dtype=numpy.single):
		print('Reading latitude vector from {0:s}..'.format(file), end='\t', flush=True)
		self.lat = SeismicFMM3D.__ReadBinaryFile(file, (self.size_lat), dtype)
		print('done.', flush=True)      
		
	def SetLatVector(self, V):
		self.lat = V
		
	def ReadLonVector(self, file, dtype=numpy.single):
		print('Reading longitude vector from {0:s}..'.format(file), end='\t', flush=True)
		self.lon = SeismicFMM3D.__ReadBinaryFile(file, (self.size_lon), dtype)
		print('done.', flush=True)    
		
	def SetLonVector(self, V):
		self.lon = V
		
	def ReadHVector(self, file, dtype=numpy.single):
		print('Reading depth vector from {0:s}..'.format(file), end='\t', flush=True)
		self.H = SeismicFMM3D.__ReadBinaryFile(file, (self.size_z), dtype)
		print('done.', flush=True)  
	
	def SetHVector(self, V):
		self.H = V
		
	def GetVelocityAtXYZ(self, lat, lon, z):
		return self.model_velocity[lon,lat,z]
		
	def SetSource(self, lat, lon, z):
		self.time[lon,lat,z] = 10e-6
		
	def SetSourceGeo(self, source_lat, source_lon, source_ele):
		#source_lat = 54.07505;
		#source_lon = 17.63717;
		#source_ele = 195;
		
		La, Lo = numpy.meshgrid(self.lat, self.lon)
		
		D = numpy.zeros_like(La)
		for a in range(self.size_lat):
			for b in range(self.size_lon):
				D[b,a] = great_circle((source_lat, source_lon), (self.lat[a], self.lon[b])).meters

		D2 = D.copy()
		D2 = D2.reshape(self.size_lat * self.size_lon)
		D2.sort()
		idx1 = D<D2[4]
		idx2 = numpy.abs(self.H - source_ele) <= self.dz
		
		selected_lat = La[idx1]
		selected_lon = Lo[idx1]
		selected_ele = self.H[idx2]
		
		for a in range(len(selected_lat)):
			for b in range(len(selected_ele)):
				distance = ((great_circle((source_lat, source_lon), (selected_lat[a], selected_lon[a])).meters)**2+(source_ele-selected_ele[b])**2)**.5
				v = numpy.squeeze(self.model_velocity[self.lat==selected_lat[a], self.lon==selected_lon[a], self.H == selected_ele[b]]);
				print(selected_lat[a], selected_lon[a], selected_ele[b], distance, v, sep='\t')
				if v > 0:
					self.time[self.lat==selected_lat[a], self.lon==selected_lon[a], self.H == selected_ele[b]] = distance / v
		#print(self.time[self.time>0])
		
		
	def SaveTime(self, file):
		print('Saving result to file {0:s}..'.format(file), end='\t', flush=True)
		self.time.reshape((self.size_lon * self.size_lat * self.size_z), order='F').tofile(file)
		print('done.', flush=True) 

	def ReadTime(self, file, dtype=numpy.single):
		print('Reading time from {0:s}..'.format(file), end='\t', flush=True)
		self.time = SeismicFMM3D.__ReadBinaryFile(file, (self.size_lat, self.size_lon, self.size_z), dtype)
		print('done.', flush=True) 
	   
	def SaveTrace(self, file):
		print('Saving result to file {0:s}..'.format(file), end='\t', flush=True)
		self.trace.reshape((self.size_lon * self.size_lat * self.size_z), order='F').tofile(file)
		print('done.', flush=True)  
		
	def Do(self, N, order='F'):
		
		print('Setting model size and grid size...', end='\t', flush=True)
		FMM3D.SetModelSize(self.size_lat, self.size_lon, self.size_z, self.dlat, self.dlon, self.dz)
		print('done.', flush=True)   

		print('Preparing for FMM calculation...', end='\t', flush=True)
		self.model_velocity = self.model_velocity.reshape((self.size_lon * self.size_lat * self.size_z), order=order)
		self.time = self.time.reshape((self.size_lon * self.size_lat * self.size_z), order=order)
		self.accepted = self.accepted.reshape((self.size_lon * self.size_lat * self.size_z), order=order)
		self.trace = self.trace.reshape((self.size_lon * self.size_lat * self.size_z), order=order)
		print('done.', flush=True)        
		
		print('++ FMM START ++', flush=True)
		start_time = time.time()
		FMM3D.FMM3D(self.model_velocity, self.time, self.accepted, self.lat, self.lon, self.H, self.trace, int(N))
		end_time = time.time()
		print('++ FMM FINISH in {0:10.6f}s ++'.format(end_time - start_time), flush=True)
		
		print('Preparing data after FMM calculation...', end='\t', flush=True)
		self.model_velocity = self.model_velocity.reshape((self.size_lat, self.size_lon, self.size_z),order=order)
		self.time = self.time.reshape((self.size_lat, self.size_lon, self.size_z), order=order)
		self.accepted = self.accepted.reshape((self.size_lat, self.size_lon, self.size_z), order=order)
		self.trace = self.trace.reshape((self.size_lat, self.size_lon, self.size_z), order=order)
		print('done.', flush=True)   
		
	def Plot(self):
		plt.imshow(self.time[:,:,300])
		plt.show()

	def GetTimeGeo(self, lat, lon, ele):
		xi = numpy.argmin(abs(self.lat-lat))
		yi = numpy.argmin(abs(self.lon-lon))
		zi = numpy.argmin(abs(self.H-ele))
		#print(xi,yi,zi,self.time[xi,yi,zi])
		min_time_val = 10e9#self.time[xi,yi,zi]
		xm = xi
		ym = yi
		zm = zi

		for x in range(xi-1,xi+2):
			for y in range(yi-1,yi+2):
				for z in range(zi-1,zi+2):
					#print("\t",x,y,z, self.time[x,y,z])
					if self.time[x,y,z] > 0:
						if self.time[x,y,z] < min_time_val:
							min_time_val = self.time[x,y,z]
							xm = x
							ym = y
							zm = z
		#print("\t\t",self.time[xm,ym,zm],self.model_velocity[xm,ym,zm])
		distancev = great_circle((self.lat[xm], self.lon[ym]), (lat, lon)).meters
		distanceh = abs(self.H[zm] - ele)
		distance = numpy.sqrt(distanceh**2+distancev**2)
		time = distance/self.model_velocity[xm,ym,zm]
		#print(distance,time,self.model_velocity[xi,yi,zi])
		return self.time[xm,ym,zm] + time
			

		
