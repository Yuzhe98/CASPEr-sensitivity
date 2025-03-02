"""
D. F. Jackson Kimball, S. Afach, D. Aybas, J. W. Blanchard, D. Budker, G. Centers, M. Engler, N. L. Figueroa, A. Garcon, P. W. Graham, H. Luo, S. Rajendran, M. G. Sendra, A. O. Sushkov, T. Wang, A. Wicken-brock, A. Wilzewski, and T. Wu, Overview of the cosmic axion spin precession experiment (casper), in Microwave Cavities and Detectors for Axion Research, edited by G. Carosi and G. Rybka (Springer International Publishing, Cham, 2020) pp. 105–121.

.bib citation
@InProceedings{10.1007/978-3-030-43761-9_13,
author="Jackson Kimball, Derek F.
and Afach, S.
and Aybas, D.
and Blanchard, J. W.
and Budker, D.
and Centers, G.
and Engler, M.
and Figueroa, N. L.
and Garcon, A.
and Graham, P. W.
and Luo, H.
and Rajendran, S.
and Sendra, M. G.
and Sushkov, A. O.
and Wang, T.
and Wickenbrock, A.
and Wilzewski, A.
and Wu, T.",
editor="Carosi, Gianpaolo
and Rybka, Gray",
title="Overview of the Cosmic Axion Spin Precession Experiment (CASPEr)",
booktitle="Microwave Cavities and Detectors for Axion Research",
year="2020",
publisher="Springer International Publishing",
address="Cham",
pages="105--121",
abstract="An overview of our experimental program to search for axion and axion-like-particle (ALP) dark matter using nuclear magnetic resonance (NMR) techniques is presented. An oscillating axion field can exert a time-varying torque on nuclear spins either directly or via generation of an oscillating nuclear electric dipole moment (EDM). Magnetic resonance techniques can be used to detect such an effect. The first-generation experiments explore many decades of ALP parameter space beyond the current astrophysical and laboratory bounds. It is anticipated that future versions of the experiments will be sensitive to the axions associated with quantum chromodynamics (QCD) having masses ≲10−9eV∕c2{\$}{\$}{\{}{\backslash}lesssim {\}}10^{\{}-9{\}}{\backslash},{\backslash}mathrm {\{}eV{\}}/c^2{\$}{\$}.",
isbn="978-3-030-43761-9"
}

"""
############################################################
############################################################
import os
import sys
os.chdir("..")  # if you want to go to parent folder
os.chdir("..")  # if you want to go to parent folder
print(os.path.abspath(os.curdir))
sys.path.insert(0, os.path.abspath(os.curdir))

from math import pi, sin, cos, sqrt
from Envelope import *
from functioncache import check

# Different samples
# Methanol properites
rho_M_Methanol = PhysicalQuantity(0.792, "g / cm**3 ")
molmassMethanol = PhysicalQuantity(32.04, "g / mol")
rho_N_MethanolProton = PhysicalQuantity(4.0, "") * rho_M_Methanol / molmassMethanol * NA

# Xe properites
rho_M_LXe = PhysicalQuantity(3.1, "g / cm**3 ")
molmassLXe = PhysicalQuantity(131.29, "g / mol")
rho_N_LXe129 = PhysicalQuantity(0.264 * 1, "") * rho_M_LXe / molmassLXe * NA
# print(rhoN_MethanolProton.convert_to("1 / cm ** 3"))

# sample magnetization
B0 = PhysicalQuantity(0.0317, "T")
Temp = PhysicalQuantity(273.15 - 90, "K")

omega_L = B0 * gamma_p
nu_L = omega_L / (2 * pi)

T2 = PhysicalQuantity(0.71, 's')

NMR_lw = PhysicalQuantity(10, 'ppm')
T2star = (1 / (pi * nu_L * NMR_lw)).convert_to('s')
check(T2star)
p = hbar * gamma_p * B0 / ( 2 * k * Temp)
p = p.convert_to('')
# spin density
ns = rho_N_MethanolProton
ns_SPN = PhysicalQuantity(sqrt(rho_N_MethanolProton.value), " 1 / cm**3 ")
# check(ns)
M0 = (mu_p * p * ns).convert_to("A/m")
# check(M0)

M0_SPN = (mu_p * ns_SPN).convert_to("A/m")



# pickup coil characteristics
gV = PhysicalQuantity(37.0, '1/m') #  For cylindrical sample R=4 mm, H=22.53 mm
Vol = PhysicalQuantity(pi * 4.**2 * 22.53, 'mm**3')
Phi_Gradio = gV * mu_0 * M0 * Vol
Phi_Gradio = Phi_Gradio.convert_to('Phi_0')
# check(PhiGradio90)
# check(PhiGradio90.convert_to("Phi_0"))
Phi_Gradio_SPN = gV * mu_0 * M0_SPN * Vol
Phi_Gradio_SPN = Phi_Gradio_SPN.convert_to('Phi_0')


# SQUID C649_G12 for DM measurement 2022.12.14 or 12.23
Lin = PhysicalQuantity(400, 'nH')
Lgrad = PhysicalQuantity(553, 'nH')
Rf = PhysicalQuantity(3, 'kiloohm')

Mf =  PhysicalQuantity(1 / 43.803, 'Phi_0 / microA')
Min =  PhysicalQuantity(1/ 0.5194 , 'Phi_0/microA')
check(Mf.convert_to('Phi_0 / mA'))
check(Min.convert_to('nH'))
# input coupling between Phi_Gradio and Phi_in
input_coupling = Min / (Lgrad + Lin)
input_coupling = input_coupling.convert_to('')
check(input_coupling)

transferfunction = Rf / Mf * Min / (Lgrad + Lin)
transferfunction = transferfunction.convert_to('V/Phi_0')

Phi_in = input_coupling * Phi_Gradio  # flux in the SQUID
Phi_in = Phi_in.convert_to('Phi_0')

Phi_in_SPN = input_coupling * Phi_Gradio_SPN  # flux in the SQUID
Phi_in_SPN = Phi_in_SPN.convert_to('Phi_0')
# check((mu_0 * M0).convert_to('T'))
# check((Min / (Lgrad + Lin)).convert_to('microV/Phi_0'))
# check(transferfunction.convert_to('microV/Phi_0'))
check(transferfunction )
# check(Phi_in)
# check(Phi_in.convert_to('Wb'))
check(Phi_in_SPN)
SPN_PSDnoise = (Phi_in_SPN)**2 / PhysicalQuantity(10, 'Hz')
check(SPN_PSDnoise**0.5)
check((Rf/Mf).convert_to('V*Phi_0**(-1)'))
# check((transferfunction).convert_to('V*Phi_0**(-1)'))



# Dark matter density in different units
rho_M_DM = PhysicalQuantity(0.4, "GeV / cm**3 / c**2") # 
rho_E_DM = PhysicalQuantity(0.4, "GeV / cm**3") # 

gaNN = PhysicalQuantity(1e-11, 'GeV**-1')

va = PhysicalQuantity(220, 'km / s')

alpha = PhysicalQuantity(2.1, 'rad')

tau_a =  PhysicalQuantity(1 / (pi * 1.34), 's')

Omega_a = PhysicalQuantity(1/2, '') * gaNN * ( PhysicalQuantity(2, '') * hbar * c * rho_E_DM) ** (1/2) * va * sin(alpha.value)
check(gaNN)
check(Omega_a.convert_to('Hz') / (2 * pi))


# tipping angle
T2=PhysicalQuantity(0.71, 's')
# T2star=PhysicalQuantity(1000, 's')
phi_wind = Omega_a * ((1./T2 + 1/tau_a)/(1/T2star + 1/tau_a)) * T2 * (tau_a / (tau_a + T2))**0.5
phi_wind = phi_wind.convert_to('rad')
# check(phi_wind.convert_to(''))
check(phi_wind.convert_to('rad'))



axionPower = 1/2 * Phi_in ** 2 * Omega_a ** 2 * ((1./T2 + 1/tau_a)/(1/T2star + 1/tau_a)) ** 2 * T2 ** 2 * (tau_a / (tau_a + T2))
print(axionPower.convert_to('Phi_0**2'))

SQUID_PSDnoise =  PhysicalQuantity(5e-11, 'Phi_0**2 / Hz')
check((SQUID_PSDnoise)**0.5)
# ofPSDnoise = 0.287195 * PhysicalQuantity(5e-11, 'Phi_0**2 / Hz')
ofPSDnoise = 0.287195 * SQUID_PSDnoise * PhysicalQuantity(1, 'Hz')
# ofPSDnoise = PhysicalQuantity(0.2e-10, 'Phi_0**2')
# print((5. * ofPSDnoise) / (axionPower / (gaNN)**2))
glim_sq = (5. * ofPSDnoise) / (axionPower / gaNN ** 2)
glim = glim_sq**(.5)
glim = glim.convert_to('GeV**-1')
check(glim)


# freq_start = PhysicalQuantity(1, 'kHz')
# freq_stop = PhysicalQuantity(4.3, 'MHz')

# mass_start = h_Planck * freq_start / c**2
# mass_stop = h_Planck * freq_stop / c**2
# check(mass_start.convert_to('eV*c**(-2)'))
# check(mass_stop.convert_to('eV*c**(-2)'))


# freq_start = PhysicalQuantity(1, 'kHz')
# freq_stop = PhysicalQuantity(1, 'GHz')

# mass_start = h_Planck * freq_start / c**2
# mass_stop = h_Planck * freq_stop / c**2
# check(mass_start.convert_to('eV*c**(-2)'))
# check(mass_stop.convert_to('eV*c**(-2)'))