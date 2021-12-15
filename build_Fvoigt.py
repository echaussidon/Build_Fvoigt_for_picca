#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Edmond CHAUSSIDON

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from scipy.special import wofz
import constants

##############
## TO CHANGE :

## it takes 3/5 min for my laptop

path_dla = 'data/master_DLA_4.4.3.fits'
path_qso = 'data/zcat_desi_drq_4.4.3.fits'
path_weight_lambda = 'data/weight_lambda.txt'
path_output = 'output/Fvoigt_example.txt'

################
## DON'T TOUCH :

data = fits.open(path_dla)[1].data
qso = fits.open(path_qso)[1].data

# keep only DLA which are front of a QSO.
data = data[:][np.in1d(data['MOCKID'], qso['THING_ID'])]

nb_dla = data['Z_DLA_RSD'].size
nb_qso = qso['Z'].size # number of line of sight

weight_lambda = np.loadtxt(path_weight_lambda)
lamb_w = weight_lambda[:,0]
weight = weight_lambda[:,1]

N_HI_min = np.min(data['N_HI_DLA'])
N_HI_2 = 20
N_HI_max = np.max(data['N_HI_DLA'])
zdla = np.mean(data['Z_DLA_RSD'])

def build_cddf() :
    y_norm, bins = np.histogram(data['N_HI_DLA'], bins=50, density=True)
    z = bins[:-1] + (bins[1]/2-bins[0]/2)
    c = np.polyfit(z, np.log(y_norm),3)
    pol = np.poly1d(c)
    fit = np.exp(pol(z))
    norma = integrate.trapz(fit, z)
    fit = fit/norma
    return z, fit

def voigt(x, sigma=1, gamma=1):
    return np.real(wofz((x + 1j*gamma)/(sigma*np.sqrt(2))))

def tau(lamb, z, N_hi): # lamb = lambda in A and N_HI in log10 and 10**N_hi in cm^-2
    lamb_rf = lamb/(1+z)
    e = 1.6021e-19 # C
    epsilon0 = 8.8541e-12 # C^2.s^2.kg^-1.m^-3
    f = 0.4164
    mp = 1.6726e-27 # kg
    me = 9.109e-31 # kg
    c = 2.9979e8 # m.s^-1
    k = 1.3806e-23 # m^2.kg.s^-2.K-1
    T = 1e4 # K
    gamma = 6.265e8 # s^-1
    lamb_alpha = constants.absorber_IGM["LYA"] # A
    Deltat_lamb = lamb_alpha/c*np.sqrt(2*k*T/mp) # A

    a = gamma/(4*np.pi*Deltat_lamb)*lamb_alpha**2/c*1e-10
    u = (lamb_rf - lamb_alpha)/Deltat_lamb
    H = voigt(u, np.sqrt(1/2), a)

    absorb = np.sqrt(np.pi)*e**2*f*lamb_alpha**2*1e-10/(4*np.pi*epsilon0*me*c**2*Deltat_lamb)*H
    # 10^N_hi in cm^-2 and absorb in m^2
    return 10**N_hi*1e4*absorb

def profile_voigt_lambda(x, z, N_hi):
    t = tau(x, z, N_hi).astype(float)
    return np.exp(-t)

def profile_lambda_to_r(lamb, profile_lambda, fidcosmo): # for lyman-alpha otherwise use an other emission line
    z = lamb/constants.absorber_IGM["LYA"] - 1
    r = fidcosmo.r_comoving(z)
    rr = np.linspace(r[0], r[-1], r.size)
    profile_r = np.interp(rr, r, profile_lambda) # to have a linear sample
    return rr, profile_r

def fft_profile(profile, dx): # not normalized
    n = profile.size
    tmp = (1-profile)
    ft_profile = dx*np.fft.fftshift(np.fft.fft(tmp))
    k = np.fft.fftshift(np.fft.fftfreq(n, dx))*(2*np.pi)
    return ft_profile, k

def build_mean_density(N, fidcosmo): # compute the mean 'density' using a weighted mean and tranform in r space
    def lambda_to_r(lamb, profile_lambda, fidcosmo): # f(lambda)dlambda = f(r)dr
        z = lamb/constants.absorber_IGM["LYA"] - 1
        r = fidcosmo.r_comoving(z)
        rr = np.linspace(r[0], r[-1], r.size)
        profile_lambda = profile_lambda*fidcosmo.hubble(z)*constants.absorber_IGM["LYA"]/3e5
        profile_r = np.interp(rr,r,profile_lambda)
        return rr, profile_r
    def density(n):
        lamb_alpha = constants.absorber_IGM["LYA"] # A
        dn = 0.1 # atom cm-2
        dlambda = 5 # A
        lamb = np.arange(3600, 5501, dlambda)
        f = np.zeros(lamb.size)
        data_lamb = (1 + data['Z_DLA_RSD'])*lamb_alpha
        qso_lamb = (1 + qso['Z'])*lamb_alpha
        for i in range(lamb.size):
            where_dla = np.where((data_lamb > lamb[i] - dlambda/2) & (data_lamb < lamb[i] + dlambda/2) & (data['N_HI_DLA']>n-dn/2) & (data['N_HI_DLA']<n+dn/2))[0]
            where_qso = np.where(qso_lamb > lamb[i])[0]
            f[i]  = where_dla.size/where_qso.size
        f = f/(dn*dlambda) # probablity to have a dla by n and by lambda bin
        return lamb, f
    mean_density = np.zeros(N.size)
    for i in range(len(N)):
         lamb, f = density(N[i])
         r, f_r = lambda_to_r(lamb, f, fidcosmo)
         r_w, weight_r = profile_lambda_to_r(lamb_w, weight, fidcosmo)
         weight_interp = np.interp(r, r_w, weight_r, left=0, right=0)
         mean_r_weight = np.average(f_r, weights=weight_interp)
         mean_density[i] = mean_r_weight
    return mean_density

def save_function():
    fidcosmo = constants.cosmo(Om=0.3)
    lamb = np.arange(2000, 8000, 1)
    dN, cddf = build_cddf()
    r_density_mean = build_mean_density(dN, fidcosmo) # mpc-1 n-1
    for i in range(dN.size):
        profile_lambda = profile_voigt_lambda(lamb, zdla, dN[i])
        r, profile_r = profile_lambda_to_r(lamb, profile_lambda, fidcosmo) # r is in Mpc h^-1 --> k (from tf) will be in (Mpc h^-1)^-1 = h Mpc^-1 :)
        ft_profile, k = fft_profile(profile_r, np.abs(r[1]-r[0]))
        ft_profile = np.abs(ft_profile)
        if i == 0:
            df = np.array([ft_profile*r_density_mean[i]])
        else:
            df = np.concatenate((df, np.array([ft_profile*r_density_mean[i]])))
    Fvoigt = np.zeros(k.size)
    for i in range(k.size):
        Fvoigt[i] = integrate.trapz(df[:,i], dN)
    Fvoigt = Fvoigt/Fvoigt[k.size//2] #normalization

    save = np.transpose(np.concatenate((np.array([k]), np.array([Fvoigt]))))
    np.savetxt(path_output, save)

##############
## Execution :

if __name__ == "__main__":
    save_function()
