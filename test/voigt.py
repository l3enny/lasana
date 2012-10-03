import lineshapes
import numpy as np
import math
from matplotlib import pyplot as plt

wavelengths = np.linspace(-4, 4, 100)
torr = 3.5
torr2hz = 25.6e6
M = 4.002602 * 1.660538921e-27
c = 299792458
f_0 = c/1082.9091140e-9
kB = 1.3806503e-23
temperatures = np.linspace(300, 500, 3)
wavelengths = wavelengths*1e9
fwhm_p = torr*torr2hz
fwhm_d = np.sqrt((8*math.log(2))*kB*temperatures/(M*c**2)) * f_0
synthetic_profiles = np.zeros((len(temperatures), len(wavelengths)))
for i in range(len(temperatures)):
    voigt = lineshapes.voigt(fwhm_d[i], fwhm_p)
    synthetic_profiles[i, :] = voigt(wavelengths)
norms = np.array(np.max(synthetic_profiles, axis=1))
synthetic_profiles = synthetic_profiles/norms[:, np.newaxis]
plt.hold(True)
for i in range(synthetic_profiles.shape[0]):
    plt.plot(wavelengths, synthetic_profiles[i, :])
plt.show()