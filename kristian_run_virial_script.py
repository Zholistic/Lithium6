##### from scipy import constants as const
from imp import reload
from uncertainties import unumpy as unp
reload(eos_analysis)
mass_li = 9.988346e-27
pi = const.pi
hbar = const.hbar
kb = const.Boltzmann
profile = np.mean(profiles_865[:, :], 1)
prof_offset = np.mean(profile[60:65])
profile -= prof_offset
profile *= 1.20
density = 2 * profile / 7.265 * 10**12
omega_z = unp.uarray((5154) * 2 * pi, 140 * 2 * pi)
a_z = unp.sqrt(hbar / (mass_li * omega_z))
omega_r = eos_analysis.radial_trap(864) * 2 * pi
pixel_area = (13e-6*83.0/400)**2
r = np.linspace(0, len(profile), len(profile))
r *= np.sqrt(pixel_area)
potential = 0.5 * mass_li * omega_r**2 * r**2
potential_nk = (potential/kb)* (10**9)

a2d = eos_analysis.scattering2d(865, a_z)
eb = (hbar**2) / (a2d**2 * mass_li)
cutoff =43
end = 70
dens_fit, fit_result = eos_analysis.virial_fit_residuals2(0.2, 0.6, eb, density[cutoff:end], potential[cutoff:end])
