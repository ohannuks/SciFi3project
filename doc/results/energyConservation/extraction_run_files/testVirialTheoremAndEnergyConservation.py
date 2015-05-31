import numpy as np
import pylab as pl

# Calculate virial theorem:
phi_time_mean = 0
v_magnitude_squared_time_mean = 0

for t in np.arange(10)+140:
  phi1 = np.loadtxt('phi' + str(t))[:,3]
  vx1 = np.loadtxt('vx' + str(t))[:,3]
  vy1 = np.loadtxt('vy' + str(t))[:,3]
  vz1 = np.loadtxt('vz' + str(t))[:,3]
  phi_sum1 = np.sum(phi1)
  v_magnitude_squared_sum1 = np.sum(vx1**2 + vy1**2 + vz1**2)
  phi_time_mean = phi_time_mean + phi_sum1
  v_magnitude_squared_time_mean = v_magnitude_squared_time_mean + v_magnitude_squared_sum1


phi_time_mean = phi_time_mean / (float)(len(np.arange(10)))
v_magnitude_squared_time_mean = v_magnitude_squared_time_mean / (float)(len(np.arange(10)))

#print phi_time_mean
#print v_magnitude_squared_time_mean

C = v_magnitude_squared_time_mean / phi_time_mean
print "Scaling factor, calculated from virial theorem: " + str(C)

# Make sure energy is conserved:

for t in np.arange(10)+140:
  phi1 = np.loadtxt('phi' + str(t))[:,3]
  vx1 = np.loadtxt('vx' + str(t))[:,3]
  vy1 = np.loadtxt('vy' + str(t))[:,3]
  vz1 = np.loadtxt('vz' + str(t))[:,3]
  phi2 = np.loadtxt('phi' + str(t+1))[:,3]
  vx2 = np.loadtxt('vx' + str(t+1))[:,3]
  vy2 = np.loadtxt('vy' + str(t+1))[:,3]
  vz2 = np.loadtxt('vz' + str(t+1))[:,3]
  phi_sum1 = C*np.sum(phi1)
  phi_sum2 = C*np.sum(phi2)
  v_magnitude_squared_sum1 = np.sum(vx1**2 + vy1**2 + vz1**2)
  v_magnitude_squared_sum2 = np.sum(vx2**2 + vy2**2 + vz2**2)
  print "Energy conservation at timestep " + str(t) + "; potential: " + str(phi_sum1) + ", v_squared: " + str(v_magnitude_squared_sum1) + "\nTotal Energy: " + str(phi_sum1 + v_magnitude_squared_sum1)




