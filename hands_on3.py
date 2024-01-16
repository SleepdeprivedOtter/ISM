# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
#import seaborn as sns
import warnings
from numpy.polynomial.polynomial import polyfit
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, CartesianDifferential
from mpl_toolkits import mplot3d

#The parameter name is always _small, not indicative of FoV size

# Ignore warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({'font.size': 25})

# Open the files and extract the data, decide for one file
fits_smaller = fits.open('1702468013097O-result.fits') #smaller FoV
#fits_smaller = fits.open('1702468189469O-result.fits') #larger FoV
data_smaller = fits_smaller[1].data
cols_smaller = fits_smaller[1].columns
fits_smaller.close()

#define working dataframe for the dataset
frame_small = pd.DataFrame({'RA': data_smaller.field(5),'DEC': data_smaller.field(7),
                           'pmra': data_smaller.field(13), 'pmdec': data_smaller.field(15),
                           'parallax': data_smaller.field(9), 'Gmag_mean': data_smaller.field(69) - 0.2478,
                           'bp-g': data_smaller.field(74) - 0.3206, 'g-rp': data_smaller.field(79) - 0.1850,
                           'radial velocity': data_smaller.field(89)},
                         index=np.linspace(1,len(data_smaller.field(0)),len(data_smaller.field(0))))

#plot the positions in RA/DEC coordinates
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['RA'], frame_small['DEC'], color='red')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

#plot the proper motions in RA/DEC coordinates
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['pmra'], frame_small['pmdec'], color='red')
plt.xlabel('pmra in mas/yr', fontsize=25)
plt.ylabel('pmdec in mas/yr', fontsize=25)
plt.xlim((-25, 10))
plt.ylim(-10, 25)
plt.show()

#plot the proper motions in RA vs the parallax
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['pmra'], frame_small['parallax'], color='red')
plt.xlabel('pmra in mas/yr', fontsize=25)
plt.ylabel('parallax in mas', fontsize=25)
plt.xlim((-25, -10))
plt.show()

#plot the proper motions in DEC vs the parallax
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['pmdec'], frame_small['parallax'], color='red')
plt.xlabel('pmdec in mas/yr', fontsize=25)
plt.ylabel('parallax in mas', fontsize=25)
plt.xlim((5, 15))
plt.show()

frame_small_filtered = pd.DataFrame(columns = ['RA', 'DEC', 'pmra', 'pmdec', 
                                               'parallax', 'Gmag_mean', 'bp-g',
                                               'g-rp', 'radial velocity'])
for i in range(len(data_smaller.field(0))):
    if frame_small['pmra'][i+1] >= -21.0 and frame_small['pmra'][i+1] <= -15.0:
        filter_pmra = True
    else:
        filter_pmra = False
        
    if frame_small['pmdec'][i+1] >= 9.0 and frame_small['pmdec'][i+1] <= 13.0:
        filter_pmdec = True
    else:
        filter_pmdec = False
        
    if frame_small['parallax'][i+1] >= 6.0 and frame_small['parallax'][i+1] <= 7.5:
        filter_parallax = True
    else:
        filter_parallax = False
        
    if filter_pmra*filter_pmdec*filter_parallax == True:
        frame_small_filtered = frame_small_filtered._append(frame_small.iloc[i])
    else:
        pass
        



#plot the proper motions in RA/DEC coordinates
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['pmra'], frame_small['pmdec'], color='red')
plt.scatter(frame_small_filtered['pmra'], frame_small_filtered['pmdec'], color='blue')
plt.xlabel('pmra in mas/yr', fontsize=25)
plt.ylabel('pmdec in mas/yr', fontsize=25)
plt.xlim((-25, 10))
plt.ylim(-10, 25)
plt.show()

#plot the proper motions in RA vs the parallax
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['pmra'], frame_small['parallax'], color='red')
plt.scatter(frame_small_filtered['pmra'], frame_small_filtered['parallax'], color='blue')
plt.xlabel('pmra in mas/yr', fontsize=25)
plt.ylabel('parallax in mas', fontsize=25)
plt.xlim((-25, -10))
#plt.ylim((6, 7))
plt.show()

#plot the proper motions in DEC vs the parallax
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['pmdec'], frame_small['parallax'], color='red')
plt.scatter(frame_small_filtered['pmdec'], frame_small_filtered['parallax'], color='blue')
plt.xlabel('pmdec in mas/yr', fontsize=25)
plt.ylabel('parallax in mas', fontsize=25)
plt.xlim((5, 15))
plt.show()

#plot the subset in the RA/DEC positions
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['RA'], frame_small['DEC'], color='red')
plt.scatter(frame_small_filtered['RA'], frame_small_filtered['DEC'], color='blue')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.show()

#get the isochrone data
isochrones = pd.read_table("Isochrones2.dat", sep="\s+", header=0, comment='#')

#split into separate isochrones
logAge_list = np.linspace(6.0, 10.1, 42)
header_isochrones = list(isochrones)
isochrone_list = []
for item in logAge_list: 
    fileName = 'isochrone_%s'%(int(np.round(item*10, 0)),)
    globals()[fileName] = isochrones.loc[np.round(isochrones['logAge'], 1) == item]
    isochrone_list.append(globals()[fileName])

#plot CMD of data
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small['Gmag_mean'] - frame_small['g-rp'], frame_small['Gmag_mean'], color='red')
plt.scatter(frame_small_filtered['Gmag_mean'] - frame_small_filtered['g-rp'], frame_small_filtered['Gmag_mean'], color='blue')
plt.gca().invert_yaxis()
plt.xlabel('G-rp [mag]', fontsize=25)
plt.ylabel('G [mag]', fontsize=25)
plt.show()

#plot best fitting isochrone
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small_filtered['Gmag_mean'] - frame_small_filtered['g-rp'], frame_small_filtered['Gmag_mean']-6, color='blue')
mmm = 15
plt.scatter(isochrone_list[mmm]['Gmag'] - isochrone_list[mmm]['G_RPmag'],
            isochrone_list[mmm]['Gmag'], color='green', label='log(age/yr)=6.15')
plt.xlim((-0.5, 2.5))
plt.ylim(-6, 16)
plt.gca().invert_yaxis()
plt.xlabel('G-rp [mag]', fontsize=25)
plt.ylabel('G [mag]', fontsize=25)
plt.legend()
plt.show()

#plot selection of isochrones to visualize variations
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small_filtered['Gmag_mean'] - frame_small_filtered['g-rp'], frame_small_filtered['Gmag_mean']-6, color='blue')
mmm = 15
plt.scatter(isochrone_list[mmm]['Gmag'] - isochrone_list[mmm]['G_RPmag'],
            isochrone_list[mmm]['Gmag'], color='green', label='log(age/yr)=6.15')
plt.scatter(isochrone_list[mmm+5]['Gmag'] - isochrone_list[mmm+5]['G_RPmag'],
            isochrone_list[mmm+5]['Gmag'], color='red', label='log(age/yr)=6.20')
plt.scatter(isochrone_list[mmm-5]['Gmag'] - isochrone_list[mmm-5]['G_RPmag'],
            isochrone_list[mmm-5]['Gmag'], color='violet', label='log(age/yr)=6.10')
plt.xlim((-0.5, 2.5))
plt.ylim(-6, 16)
plt.gca().invert_yaxis()
plt.xlabel('G-rp [mag]', fontsize=25)
plt.ylabel('G [mag]', fontsize=25)
plt.legend()
plt.show()

#create empty frame for conversion 
frame_small_cartesian = pd.DataFrame(columns = ['X', 'Y', 'Z', 'U', 'V', 'W'])

#actual conversion
for m in range(len(frame_small_filtered['RA'])):
    icrs = ICRS(ra=frame_small_filtered.iloc[m, 0]*u.degree,
                dec=frame_small_filtered.iloc[m, 1]*u.degree,
                distance=(1000/frame_small_filtered.iloc[m, 4])*u.pc,
                pm_ra_cosdec=frame_small_filtered.iloc[m, 2]*u.mas/u.yr,
                pm_dec=frame_small_filtered.iloc[m, 3]*u.mas/u.yr,
                radial_velocity=frame_small_filtered.iloc[m, 8]*u.km/u.s)
    icrs.representation_type = 'cartesian'
    frame_small_cartesian = frame_small_cartesian._append({'X':icrs.x.value,
                                                           'Y':icrs.y.value,
                                                           'Z':icrs.z.value,
                                                           'U':icrs.v_x.value,
                                                           'V':icrs.v_y.value,
                                                           'W':icrs.v_z.value},
                                                          ignore_index=True)

#3D plot of the positions
fig = plt.figure(figsize=(25,10))
ax = plt.axes(projection='3d')
ax.scatter(frame_small_cartesian['X'],
            frame_small_cartesian['Y'],
            frame_small_cartesian['Z'], color='blue', s = 50)
ax.set_xlabel('X', labelpad=20)
ax.set_ylabel('Y', labelpad=20)
ax.set_zlabel('Z', labelpad=20)
ax.set_xlim(-70, -54)
ax.set_ylim(16, 27)
ax.set_zlim(-155, -120)
ax.set_xticks([-70, -66, -62, -58])
ax.set_zticks([-120, -130, -140, -150])
ax.dist = 11.5
plt.show()

#3D plot of the velocities
fig = plt.figure(figsize=(25,10))
ax = fig.add_subplot(projection='3d')
ax.scatter(frame_small_cartesian['U'],
            frame_small_cartesian['V'],
            frame_small_cartesian['W'], color='blue', s = 50)
ax.set_xlabel('u', labelpad=20)
ax.set_ylabel('v', labelpad=20)
ax.set_zlabel('w', labelpad=20)
ax.set_xlim(-30, 10)
ax.set_ylim(10, 25)
ax.set_zlim(-40, 20)
ax.set_xticks([-25, -15, -5, 5])
ax.set_yticks([10, 15, 20, 25])
ax.set_zticks([-40, -20, 0, 20])
ax.view_init(elev=30, azim=30, roll=0)
plt.show()

#plot U vs V
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small_cartesian['U'], frame_small_cartesian['V'], color='blue')
plt.xlabel('U', fontsize=25)
plt.ylabel('V', fontsize=25)
plt.close()

#plot U vs W
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small_cartesian['U'], frame_small_cartesian['W'], color='blue')
plt.xlabel('U', fontsize=25)
plt.ylabel('W', fontsize=25)
plt.close()

#plot V vs W
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frame_small_cartesian['Y'], frame_small_cartesian['Z'], color='blue')
plt.xlabel('Y', fontsize=25)
plt.ylabel('Z', fontsize=25)
plt.show()

#how many stars are in our subset
print(len(frame_small_cartesian['X']), len(frame_small_cartesian['U']))

#plot X vs U and linear fit of it
frameXU = pd.DataFrame({'X':frame_small_cartesian['X'],
                        'U':frame_small_cartesian['U']})
frameXU= frameXU.dropna(axis= 0, how='any')
print(len(frameXU)) #how much velocity data do we have
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frameXU['X'],frameXU['U'], color='blue')
fitxu = polyfit(frameXU['X'],frameXU['U'], deg=1)
xurange = np.linspace(-70,-50, 100)
plt.plot(xurange,fitxu[1]*xurange+fitxu[0], color='red')
plt.xlabel('X in pc', fontsize=25)
plt.ylabel('U in km/s', fontsize=25)
plt.show()

#plot Y vs V and linear fit of it
frameYV = pd.DataFrame({'Y':frame_small_cartesian['Y'],
                        'V':frame_small_cartesian['V']})
frameYV= frameYV.dropna(axis= 0, how='any')
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frameYV['Y'],frameYV['V'], color='blue')
fityv = polyfit(frameYV['Y'],frameYV['V'], deg=1)
yvrange = np.linspace(10,35, 100)
plt.plot(yvrange,fityv[1]*yvrange+fityv[0], color='red')
plt.xlabel('Y in pc', fontsize=25)
plt.ylabel('V in km/s', fontsize=25)
plt.show()

#plot Z vs W and linear fit of it
frameZW = pd.DataFrame({'Z':frame_small_cartesian['Z'],
                        'W':frame_small_cartesian['W']})
frameZW= frameZW.dropna(axis= 0, how='any')
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(frameZW['Z'],frameZW['W'], color='blue')
fitzw = polyfit(frameZW['Z'],frameZW['W'], deg=1)
zwrange = np.linspace(-155,-120, 100)
plt.plot(zwrange,fitzw[1]*zwrange+fitzw[0], color='red')
plt.xlabel('Z in pc', fontsize=25)
plt.ylabel('W in km/s', fontsize=25)
plt.ylim(-50, 50)
plt.show()

