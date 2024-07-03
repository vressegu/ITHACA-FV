import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import sys
sys.path.insert(0, "./")

#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

time = np.load("./ITHACAoutput/direct/timeVector_mat.npy")
probe_true = np.load("./ITHACAoutput/direct/probe_true_mat.npy")
posteriorMean = np.load("./ITHACAoutput/reconstuction/probe_mean_mat.npy")
minConfidence = np.load("./ITHACAoutput/reconstuction/probe_minConfidence_mat.npy")
maxConfidence = np.load("./ITHACAoutput/reconstuction/probe_MaxConfidence_mat.npy")



fig = plt.figure(1,figsize=(8,6))
plt.plot(time, probe_true,"b--", linewidth = 2, label="T true")

#plt.plot(time, posteriorMean[0,:], "o", linewidth = 2, label="p1")
#plt.plot(time, posteriorMean[1,:], "o", linewidth = 2, label="p2")
#plt.plot(time, posteriorMean[2,:], "o", linewidth = 2, label="p3")
#plt.plot(time, posteriorMean[3,:], "o", linewidth = 2, label="p4")

plt.fill_between(time, minConfidence, maxConfidence, color='b', alpha=.1)
plt.plot(time,posteriorMean, linewidth = 2, color='b', label="T rec" )

#plt.legend()
plt.xlabel('Time [s]', fontsize=25)
plt.grid()

plt.legend()


plt.show()
