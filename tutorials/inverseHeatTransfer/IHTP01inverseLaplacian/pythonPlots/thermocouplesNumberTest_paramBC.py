import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relError_L2norm = np.load("../caseDir/ITHACAoutput/thermocouplesNumberTest_paramBC/relError_L2norm_mat.npy")
relError_LinfNorm = np.load("../caseDir/ITHACAoutput/thermocouplesNumberTest_paramBC/relError_LinfNorm_mat.npy")
TCplane_Y = np.load("../caseDir/ITHACAoutput/thermocouplesNumberTest_paramBC/numberTCperAxis_mat.npy")


f = plt.figure(5,figsize=(12,8))
plt.semilogy(TCplane_Y, relError_L2norm, "bo", markersize=15,label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.semilogy(TCplane_Y, relError_LinfNorm, "kv", markersize=15,label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel(r'Number of thermocouples per axis', fontsize=25)
plt.xlim(0, TCplane_Y[-1])
plt.grid()
plt.title(r"Parameterized BC (LU)", fontsize=25)
plt.legend(fontsize=25)
plt.yscale('log')

plt.show()
