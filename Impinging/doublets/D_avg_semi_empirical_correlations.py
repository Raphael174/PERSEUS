"""
Initially from Dombrowki and Johns 1963 and extended in Ryan 1995

Ryan 1995 measured drop sizes at x=16mm from the impingement point, along the sheet centerline
Therefore, r = 16mm and phi=0°
dj = do (assumption) as done in Ryan 1995

variables:

    theta   : impingement angle for doublets
    sheet_ang     : angular location around doublet sheet, phi=0 along jets' plane
    r       : radial distance from impingement point on doublet sheet
    h       : liquid sheet thickness from doublets
    km      : most unstable wave number [Ryan 1995]
    
    rho_g   : chamber/atmospheric gas density
    rho_L   : propellant density
    sigma   : surfce tension

    d_j     : jet diameter, assumed = d_o 


"""

import numpy as np


def h_like_doublet(d_j, theta, r, sheet_ang):

    theta_half = theta/2
    num = d_j**2 * np.sin(theta_half*np.pi/180)**2
    den = 4*r*(1 - np.cos(sheet_ang*np.pi/180)*np.cos(theta_half*np.pi/180))**2
    return num/den

def km_unstable_wave(rho_g, U_s, sigma):
    return rho_g * U_s**2/(2*sigma)

def D_ligament(d_j, theta, U_s, sigma, rho_g, r=16e-3, sheet_ang=0):

    h  = h_like_doublet(d_j, theta, r, sheet_ang)
    km = km_unstable_wave(rho_g, U_s, sigma)

    return np.sqrt(4*h/km)

def D_drop_mean_LinStab(d_j, theta, U_s, sigma, rho_g, rho_L):
    d_L = D_ligament(d_j, theta, U_s, sigma, rho_g)
    cst_term = (3*np.pi/np.sqrt(2))**(1/3) 
    fluid_term = (1 + (3*U_s*rho_L/(np.sqrt(rho_L*sigma*d_j))))**(1/6)
    return cst_term*d_L*fluid_term

def D_drop_mean_John1995(d_j, theta, U_s, sigma, rho_g, rho_L):
    We = rho_L*U_s**2*d_j/sigma
    f_theta = (1 - np.cos(theta*np.pi/180))**2 / (np.sin(theta*np.pi/180))**3
    rho_g__rho_L = rho_g/rho_L

    cst_term = 2.62/64**(1/3)
    return d_j * cst_term * rho_g__rho_L**(-1/6) * (We * f_theta)**(-1/3)



