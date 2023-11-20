# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 10:16:59 2023

@author: 96fan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sns
import warnings
import matplotlib as mpl



# Ignore warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({'font.size': 25})

# Open the files and extract the data
fits_2mass = fits.open('2MASS.fits')
#fits_2mass.info()
#print(repr(fits_2mass[0].data)) #NOTE:HEADER EMPTY
data_2mass = fits_2mass[1].data
cols_2mass = fits_2mass[1].columns
fits_2mass.close()

fits_gaia = fits.open('GAIA_DR3.fits')
#fits_gaia.info()
data_gaia = fits_gaia[1].data
cols_gaia = fits_gaia[1].columns
fits_gaia.close()

fits_apass = fits.open('APASS.fits')
#fits_apass.info()
data_apass = fits_apass[1].data
cols_apass = fits_apass[1].columns
fits_apass.close()

# SHOW COLUMN NAMES IF NEEDED
#print(cols_gaia.names)

# PART 2
print('PART 2')

# Plot a 2D-space diagram of the data
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(data_gaia.field(0), data_gaia.field(1) ,s=1 , c='green') 
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

# Calculate the histogram of the luminosities in the G-band
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
mag_range = np.linspace(8, 22, 15)
n, bins, _ = plt.hist(data_gaia.field('Gmag'), bins=mag_range, color='blue', density=False)
plt.close()

# Calculate the slope of the histogram and plot both together
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
# Completeness until 18mg, used for fit: 12-17mag
linfit_gaia_gmag, cov_gaia = np.polyfit(np.linspace(12.5, 16.5, 5), np.log10(n[4:9]), 1, cov=True)
plt.stairs(np.log10(n), mag_range, color='blue')
plt.plot(mag_range, (linfit_gaia_gmag[0] *mag_range + linfit_gaia_gmag[1]), color='red')
#plt.yscale('log')
plt.xlabel('mag', fontsize=25)
plt.ylabel('logN', fontsize=25)
plt.close()
print('The b-value is:', np.round(linfit_gaia_gmag[0], 3))
print('Uncertainties fitting parameters:', np.sqrt(np.diag(cov_gaia)))
b_gaia = np.round(linfit_gaia_gmag[0], 3)

# Plot the 2D density plot with seaborn
sns.set_style('whitegrid')
Gmag_array = np.array(data_gaia.field('Gmag'))
pd_RA_DEC = pd.DataFrame({'RA': data_gaia.field(0),'DEC': data_gaia.field(1)},
                         index=np.linspace(1,len(data_gaia.field(0)),len(data_gaia.field(0))))
pd_RA_DEC['Gmag'] = Gmag_array.byteswap().newbyteorder()
np_RA = np.array(data_gaia.field(0))
new_np_RA = np_RA.byteswap().newbyteorder()
np_DEC = np.array(data_gaia.field(1))
new_np_DEC = np_DEC.byteswap().newbyteorder()
#sns.kdeplot(x=new_np_RA, y=new_np_DEC, cmap="Reds", shade=True, bw_adjust=.5)
plt.close()

def total_extinction(b, N_control, N_science):
    A = (1/b)*np.log10(N_control/N_science)
    return A

def create_subset_gaia(mid_RA, mid_DEC, side_RA, side_DEC):
    right_RA = mid_RA + (side_RA/2)
    left_RA = mid_RA - (side_RA/2)
    upper_DEC = mid_DEC + (side_DEC/2)
    lower_DEC = mid_DEC - (side_DEC/2)
    filter_RA = pd_RA_DEC['RA'].between(left_RA, right_RA, inclusive = 'both')
    pd_RA_DEC['filter_RA'] = filter_RA
    filter_DEC = pd_RA_DEC['DEC'].between(lower_DEC, upper_DEC, inclusive = 'both')
    pd_RA_DEC['filter_DEC'] = filter_DEC
    np__RA_DEC = pd_RA_DEC.to_numpy()
    #print(np.shape(np__RA_DEC))
    #print(pd_RA_DEC[1])
    #print(filter_RA)
    
    subset_RA = []
    subset_DEC = []
    subset_Gmag = []
    subset_testRA = []
    subset_testDEC = []
    for i in range(len(data_gaia.field(0))):
        if np__RA_DEC[i][3]==1 and np__RA_DEC[i][4]==1:
            row_got = np__RA_DEC[i]
            subset_RA.append(row_got[0])
            subset_DEC.append(row_got[1])
            subset_Gmag.append(row_got[2])
            subset_testRA.append(row_got[3])
            subset_testDEC.append(row_got[4])
        else:
            pass
    
    subset = pd.DataFrame({'RA': subset_RA, 'DEC': subset_DEC,'Gmag': subset_Gmag,
                           'boolRA': subset_testRA,'boolDEC': subset_testDEC},
                          index=np.linspace(1,len(subset_RA),len(subset_RA)))
    
    return subset

# Define the fields to check
control_left = create_subset_gaia(192, -77.75, 0.5, 0.25)
control_right = create_subset_gaia(196, -76.65, 0.5, 0.25)
subset1 = create_subset_gaia(193.8, -76.7, 0.5, 0.25)
subset2 = create_subset_gaia(194.9, -77.20, 0.5, 0.25)
subset3 = create_subset_gaia(193.3, -77.15, 0.5, 0.25)
subset4 = create_subset_gaia(193.75, -77.15, 0.5, 0.25)
subset5 = create_subset_gaia(194.7, -77.5, 0.5, 0.25)
subset6 = create_subset_gaia(196, -77.47, 0.5, 0.25)

# Plot the fields 
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(data_gaia.field(0), data_gaia.field(1) ,s=1 , c='green') 
plt.scatter(subset1['RA'], subset1['DEC'],s=1 , c='red')
plt.scatter(subset2['RA'], subset2['DEC'],s=1 , c='black')
plt.scatter(subset3['RA'], subset3['DEC'],s=1 , c='gold')
plt.scatter(subset4['RA'], subset4['DEC'],s=1 , c='hotpink')
plt.scatter(subset5['RA'], subset5['DEC'],s=1 , c='darkviolet')
plt.scatter(subset6['RA'], subset6['DEC'],s=1 , c='chocolate')
plt.scatter(control_left['RA'], control_left['DEC'],s=1 , c='deepskyblue')
plt.scatter(control_right['RA'], control_right['DEC'],s=1 , c='navy')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.show()

# Calculate the A_G and print them out
A_1_left = total_extinction(b_gaia, len(control_left), len(subset1))/0.83627
A_1_right = total_extinction(b_gaia, len(control_right), len(subset1))/0.83627
A_2_left = total_extinction(b_gaia, len(control_left), len(subset2))/0.83627
A_2_right = total_extinction(b_gaia, len(control_right), len(subset2))/0.83627
A_3_left = total_extinction(b_gaia, len(control_left), len(subset3))/0.83627
A_3_right = total_extinction(b_gaia, len(control_right), len(subset3))/0.83627
A_4_left = total_extinction(b_gaia, len(control_left), len(subset4))/0.83627
A_4_right = total_extinction(b_gaia, len(control_right), len(subset4))/0.83627
A_5_left = total_extinction(b_gaia, len(control_left), len(subset5))/0.83627
A_5_right = total_extinction(b_gaia, len(control_right), len(subset5))/0.83627
A_6_left = total_extinction(b_gaia, len(control_left), len(subset6))/0.83627
A_6_right = total_extinction(b_gaia, len(control_right), len(subset6))/0.83627

print('N values:', len(control_left), len(control_right), len(subset1), len(subset2),
      len(subset3), len(subset4), len(subset5), len(subset6))
print('A values for left control field:', np.round(A_1_left, 3), np.round(A_2_left, 3),
      np.round(A_3_left, 3), np.round(A_4_left, 3), np.round(A_5_left, 3), np.round(A_6_left, 3))
print('A values for right control field:', np.round(A_1_right, 3), np.round(A_2_right, 3),
      np.round(A_3_right, 3), np.round(A_4_right, 3), np.round(A_5_right, 3), np.round(A_6_right, 3))

# 2D distribution plot with A_v color code (need more fields for that,
#       normal denisty plot for now)
#print(cols_2mass.names)
sns.set_style('whitegrid')
pd_RA_DEC_2mass = pd.DataFrame({'RA': data_2mass.field(3),'DEC': data_2mass.field(4)},
                         index=np.linspace(1,len(data_2mass.field(3)),len(data_2mass.field(3))))
np_RA_2mass = np.array(data_2mass.field(3))
new_np_RA_2mass = np_RA_2mass.byteswap().newbyteorder()
np_DEC_2mass = np.array(data_2mass.field(4))
new_np_DEC_2mass = np_DEC_2mass.byteswap().newbyteorder()
#sns.kdeplot(x=new_np_RA_2mass, y=new_np_DEC_2mass, cmap="Reds", shade=True, bw_adjust=.5)
plt.close()

# PART 3
print('PART 3')

# Plot a 2D-space diagram of the data
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(data_apass.field(1), data_apass.field(2) ,s=1 , c='green') 
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

# Calculate the histogram of the luminosities in the V-band
#print(cols_apass.names)
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
mag_range = np.linspace(10, 19, 10)
n_apass, bins, _ = plt.hist(data_apass.field('Vmag'), bins=mag_range, color='blue', density=False)
plt.close()

# Calculate the slope of the histogram and plot both together
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
# Completeness until 15mg, used for fit: 10-15mag
linfit_apass_vmag, cov_apass = np.polyfit(np.linspace(11.5, 14.5, 4), np.log10(n_apass[1:5]), 1, cov=True)
plt.stairs(np.log10(n_apass), mag_range, color='blue')
plt.plot(mag_range, (linfit_apass_vmag[0] *mag_range + linfit_apass_vmag[1]), color='red')
plt.xlabel('mag', fontsize=25)
plt.ylabel('logN', fontsize=25)
plt.close()
print('The b-value is:', np.round(linfit_apass_vmag[0], 3))
print('Uncertainties fitting parameters:', np.sqrt(np.diag(cov_apass)))
b_apass = np.round(linfit_apass_vmag[0], 3)

# Plot the 2D density plot with seaborn and color code by B-V
BVmag_array_apass = np.array(data_apass.field('B-V'))
pd_RA_DEC_apass = pd.DataFrame({'RA': data_apass.field(1),'DEC': data_apass.field(2)},
                         index=np.linspace(1,len(data_apass.field(1)),len(data_apass.field(1))))
pd_RA_DEC_apass['BVmag'] = BVmag_array_apass.byteswap().newbyteorder()
np_RA_apass = np.array(data_apass.field(1))
new_np_RA_apass = np_RA_apass.byteswap().newbyteorder()
np_DEC_apass = np.array(data_apass.field(2))
new_np_DEC_apass = np_DEC_apass.byteswap().newbyteorder()
new_BVmag_array_apass = BVmag_array_apass.byteswap().newbyteorder()

fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_apass['RA'], pd_RA_DEC_apass['DEC'], s=100,
            c=pd_RA_DEC_apass['BVmag'], cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax, label='E(B-V)')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

# Define function for the apass to create subsets
def create_subset_apass(mid_RA, mid_DEC, side_RA, side_DEC):
    right_RA = mid_RA + (side_RA/2)
    left_RA = mid_RA - (side_RA/2)
    upper_DEC = mid_DEC + (side_DEC/2)
    lower_DEC = mid_DEC - (side_DEC/2)
    filter_RA = pd_RA_DEC_apass['RA'].between(left_RA, right_RA, inclusive = 'both')
    pd_RA_DEC_apass['filter_RA'] = filter_RA
    filter_DEC = pd_RA_DEC_apass['DEC'].between(lower_DEC, upper_DEC, inclusive = 'both')
    pd_RA_DEC_apass['filter_DEC'] = filter_DEC
    np__RA_DEC = pd_RA_DEC_apass.to_numpy()
    #print(np.shape(np__RA_DEC))
    #print(pd_RA_DEC[1])
    #print(filter_RA)
    
    subset_RA = []
    subset_DEC = []
    subset_BVmag = []
    subset_testRA = []
    subset_testDEC = []
    for i in range(len(data_apass.field(0))):
        if np__RA_DEC[i][3]==1 and np__RA_DEC[i][4]==1:
            row_got = np__RA_DEC[i]
            subset_RA.append(row_got[0])
            subset_DEC.append(row_got[1])
            subset_BVmag.append(row_got[2])
            subset_testRA.append(row_got[3])
            subset_testDEC.append(row_got[4])
        else:
            pass
    
    subset = pd.DataFrame({'RA': subset_RA, 'DEC': subset_DEC,'BVmag': subset_BVmag,
                           'boolRA': subset_testRA,'boolDEC': subset_testDEC},
                          index=np.linspace(1,len(subset_RA),len(subset_RA)))
    
    return subset

# Create the same subsets as before, just for apass
# NOTE: To get better results, the subsets should be larger
control_left_apass = create_subset_apass(192, -77.75, 0.5, 0.25)
control_right_apass = create_subset_apass(196, -76.65, 0.5, 0.25)
subset1_apass = create_subset_apass(193.8, -76.7, 0.5, 0.25)
subset2_apass = create_subset_apass(194.9, -77.20, 0.5, 0.25)
subset3_apass = create_subset_apass(193.3, -77.15, 0.5, 0.25)
subset4_apass = create_subset_apass(193.75, -77.15, 0.5, 0.25)
subset5_apass = create_subset_apass(194.7, -77.5, 0.5, 0.25)
subset6_apass = create_subset_apass(196, -77.47, 0.5, 0.25)

# Define function to calculate mean and std of the Mag of a given field
def mean_std_fields(field):
    mean = np.mean(field.iloc[:,2])
    std = np.std(field.iloc[:,2])
    return np.round(mean, 3), np.round(std, 3)

# Get the mean and std for all fields
BV_mean_CL_apass, BV_std_CL_apass = mean_std_fields(control_left_apass)
BV_mean_CR_apass, BV_std_CR_apass = mean_std_fields(control_right_apass)
BV_mean_s1_apass, BV_std_s1_apass = mean_std_fields(subset1_apass)
BV_mean_s2_apass, BV_std_s2_apass = mean_std_fields(subset2_apass)
BV_mean_s3_apass, BV_std_s3_apass = mean_std_fields(subset3_apass)
BV_mean_s4_apass, BV_std_s4_apass = mean_std_fields(subset4_apass)
BV_mean_s5_apass, BV_std_s5_apass = mean_std_fields(subset5_apass)
BV_mean_s6_apass, BV_std_s6_apass = mean_std_fields(subset6_apass)

print('Means:', BV_mean_CL_apass, BV_mean_CR_apass, BV_mean_s1_apass, BV_mean_s2_apass, 
      BV_mean_s3_apass, BV_mean_s4_apass, BV_mean_s5_apass, BV_mean_s6_apass)
print('Std:', BV_std_CL_apass, BV_std_CR_apass, BV_std_s1_apass, BV_std_s2_apass, 
      BV_std_s3_apass, BV_std_s4_apass, BV_std_s5_apass, BV_std_s6_apass)

# Calculate the extinction from it
A_L_s1_apass = np.round(3.1*(BV_mean_s1_apass - BV_mean_CL_apass), 3)
A_R_s1_apass = np.round(3.1*(BV_mean_s1_apass - BV_mean_CR_apass), 3)
A_L_s2_apass = np.round(3.1*(BV_mean_s2_apass - BV_mean_CL_apass), 3)
A_R_s2_apass = np.round(3.1*(BV_mean_s2_apass - BV_mean_CR_apass), 3)
A_L_s3_apass = np.round(3.1*(BV_mean_s3_apass - BV_mean_CL_apass), 3)
A_R_s3_apass = np.round(3.1*(BV_mean_s3_apass - BV_mean_CR_apass), 3)
A_L_s4_apass = np.round(3.1*(BV_mean_s4_apass - BV_mean_CL_apass), 3)
A_R_s4_apass = np.round(3.1*(BV_mean_s4_apass - BV_mean_CR_apass), 3)
A_L_s5_apass = np.round(3.1*(BV_mean_s5_apass - BV_mean_CL_apass), 3)
A_R_s5_apass = np.round(3.1*(BV_mean_s5_apass - BV_mean_CR_apass), 3)
A_L_s6_apass = np.round(3.1*(BV_mean_s6_apass - BV_mean_CL_apass), 3)
A_R_s6_apass = np.round(3.1*(BV_mean_s6_apass - BV_mean_CR_apass), 3)

# Give out the results
print('A values for left control field:', A_L_s1_apass, A_L_s2_apass, A_L_s3_apass,
      A_L_s4_apass, A_L_s5_apass, A_L_s6_apass)
print('A values for right control field:', A_R_s1_apass, A_R_s2_apass, A_R_s3_apass,
      A_R_s4_apass, A_R_s5_apass, A_R_s6_apass)

# Write color excess into the dataframe
pd_RA_DEC_apass['BVmag_correctedL'] = pd_RA_DEC_apass['BVmag'] - BV_mean_CL_apass
pd_RA_DEC_apass['BVmag_correctedR'] = pd_RA_DEC_apass['BVmag'] - BV_mean_CR_apass

# Calculate extinction via dataframe
A_CL = 3.1*pd_RA_DEC_apass['BVmag_correctedL']
A_CR = 3.1*pd_RA_DEC_apass['BVmag_correctedR']

# Plot extinction map for left control field
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_apass['RA'], pd_RA_DEC_apass['DEC'], s=100,
            c=A_CL, cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax)
plt.close()

# Plot extinction map for right control field
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_apass['RA'], pd_RA_DEC_apass['DEC'], s=100,
            c=A_CR, cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax, label='A_BV')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

# PART 4
print('PART 4')

# Calculate the histogram of the luminosities in the V-band
#print(cols_2mass.names)

# Create subset with only Qflg = AAA
subset_2mass_RA = []
subset_2mass_DEC = []
subset_2mass_Jmag = []
subset_2mass_Jmag_True = []
for i in range(len(data_2mass['Qflg'])):
    if data_2mass[i]['Qflg']=='AAA':
        subset_2mass_RA.append(data_2mass[i]['RAJ2000'])
        subset_2mass_DEC.append(data_2mass[i]['DEJ2000'])
        Hmag = data_2mass[i]['Hmag']
        Kmag = data_2mass[i]['Kmag']
        HKmag = Hmag - Kmag
        subset_2mass_Jmag.append(HKmag)
        subset_2mass_Jmag_True.append(data_2mass[i]['Jmag'])
    else:
        pass

# Put data(subset) into a dataframe
pd_RA_DEC_2mass = pd.DataFrame({'RA': subset_2mass_RA, 'DEC': subset_2mass_DEC,
                      'Jmag': subset_2mass_Jmag}, 
                      index=np.linspace(1,len(subset_2mass_RA),len(subset_2mass_RA)))

# Plot a 2D-space diagram of the data
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(subset_2mass_RA, subset_2mass_DEC ,s=1 , c='green') 
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

# Calculate the number of stars per luminosity in Jmag
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
mag_range = np.linspace(8, 17, 10)
n_2mass, bins, _ = plt.hist(subset_2mass_Jmag_True, bins=mag_range, color='blue', density=False)
plt.close()


# Calculate the slope of the histogram and plot both together
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
# Completeness until 15mg, used for fit: 10-15mag
linfit_2mass_jmag, cov_2mass = np.polyfit(np.linspace(10.5, 13.5, 4), np.log10(n_2mass[2:6]), 1, cov=True)
plt.stairs(np.log10(n_2mass), mag_range, color='blue')
plt.plot(mag_range, (linfit_2mass_jmag[0] *mag_range + linfit_2mass_jmag[1]), color='red')
plt.xlabel('mag', fontsize=25)
plt.ylabel('logN', fontsize=25)
plt.close()
print('The b-value is:', np.round(linfit_2mass_jmag[0], 3))
print('Uncertainties fitting parameters:', np.sqrt(np.diag(cov_2mass)))
b_2mass = np.round(linfit_2mass_jmag[0], 3)

# Plot the 2D density plot with seaborn and color code by J
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_2mass['RA'], pd_RA_DEC_2mass['DEC'], s=20,
            c=pd_RA_DEC_2mass['Jmag'], cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax, label= 'H-K')
plt.close()

# Create fields for the 2mass data
def create_subset_2mass(mid_RA, mid_DEC, side_RA, side_DEC):
    right_RA = mid_RA + (side_RA/2)
    left_RA = mid_RA - (side_RA/2)
    upper_DEC = mid_DEC + (side_DEC/2)
    lower_DEC = mid_DEC - (side_DEC/2)
    filter_RA = pd_RA_DEC_2mass['RA'].between(left_RA, right_RA, inclusive = 'both')
    pd_RA_DEC_2mass['filter_RA'] = filter_RA
    filter_DEC = pd_RA_DEC_2mass['DEC'].between(lower_DEC, upper_DEC, inclusive = 'both')
    pd_RA_DEC_2mass['filter_DEC'] = filter_DEC
    np__RA_DEC = pd_RA_DEC_2mass.to_numpy()
    #print(np.shape(np__RA_DEC))
    #print(pd_RA_DEC[1])
    #print(filter_RA)
    
    subset_RA = []
    subset_DEC = []
    subset_Jmag = []
    subset_testRA = []
    subset_testDEC = []
    for i in range(len(pd_RA_DEC_2mass['RA'])):
        if np__RA_DEC[i][3]==1 and np__RA_DEC[i][4]==1:
            row_got = np__RA_DEC[i]
            subset_RA.append(row_got[0])
            subset_DEC.append(row_got[1])
            subset_Jmag.append(row_got[2])
            subset_testRA.append(row_got[3])
            subset_testDEC.append(row_got[4])
        else:
            pass
    
    subset = pd.DataFrame({'RA': subset_RA, 'DEC': subset_DEC,'Jmag': subset_Jmag,
                           'boolRA': subset_testRA,'boolDEC': subset_testDEC},
                          index=np.linspace(1,len(subset_RA),len(subset_RA)))
    
    return subset

# Create the fields 
control_left_2mass = create_subset_2mass(192, -77.75, 0.5, 0.25)
control_right_2mass = create_subset_2mass(196, -76.65, 0.5, 0.25)
subset1_2mass = create_subset_2mass(193.8, -76.7, 0.5, 0.25)
subset2_2mass = create_subset_2mass(194.9, -77.20, 0.5, 0.25)
subset3_2mass = create_subset_2mass(193.3, -77.15, 0.5, 0.25)
subset4_2mass = create_subset_2mass(193.75, -77.15, 0.5, 0.25)
subset5_2mass = create_subset_2mass(194.7, -77.5, 0.5, 0.25)
subset6_2mass = create_subset_2mass(196, -77.47, 0.5, 0.25)

# Calculate the mean and std for the Jmag
J_mean_CL_2mass, J_std_CL_2mass = mean_std_fields(control_left_2mass)
J_mean_CR_2mass, J_std_CR_2mass = mean_std_fields(control_right_2mass)
J_mean_s1_2mass, J_std_s1_2mass = mean_std_fields(subset1_2mass)
J_mean_s2_2mass, J_std_s2_2mass = mean_std_fields(subset2_2mass)
J_mean_s3_2mass, J_std_s3_2mass = mean_std_fields(subset3_2mass)
J_mean_s4_2mass, J_std_s4_2mass = mean_std_fields(subset4_2mass)
J_mean_s5_2mass, J_std_s5_2mass = mean_std_fields(subset5_2mass)
J_mean_s6_2mass, J_std_s6_2mass = mean_std_fields(subset6_2mass)

print('Means:', J_mean_CL_2mass, J_mean_CR_2mass, J_mean_s1_2mass, J_mean_s2_2mass, 
      J_mean_s3_2mass, J_mean_s4_2mass, J_mean_s5_2mass, J_mean_s6_2mass)
print('Std:', J_std_CL_2mass, J_std_CR_2mass, J_std_s1_2mass, J_std_s2_2mass, 
      J_std_s3_2mass, J_std_s4_2mass, J_std_s5_2mass, J_std_s6_2mass)

# Calculate extinction from it
# NOTE: Not sure if this is right as its the excess of J and not J-V, but there
#           is no data from the satellite
A_L_s1_2mass = np.round(15.87*(J_mean_s1_2mass - J_mean_CL_2mass), 3)
A_R_s1_2mass = np.round(15.87*(J_mean_s1_2mass - J_mean_CR_2mass), 3)
A_L_s2_2mass = np.round(15.87*(J_mean_s2_2mass - J_mean_CL_2mass), 3)
A_R_s2_2mass = np.round(15.87*(J_mean_s2_2mass - J_mean_CR_2mass), 3)
A_L_s3_2mass = np.round(15.87*(J_mean_s3_2mass - J_mean_CL_2mass), 3)
A_R_s3_2mass = np.round(15.87*(J_mean_s3_2mass - J_mean_CR_2mass), 3)
A_L_s4_2mass = np.round(15.87*(J_mean_s4_2mass - J_mean_CL_2mass), 3)
A_R_s4_2mass = np.round(15.87*(J_mean_s4_2mass - J_mean_CR_2mass), 3)
A_L_s5_2mass = np.round(15.87*(J_mean_s5_2mass - J_mean_CL_2mass), 3)
A_R_s5_2mass = np.round(15.87*(J_mean_s5_2mass - J_mean_CR_2mass), 3)
A_L_s6_2mass = np.round(15.87*(J_mean_s6_2mass - J_mean_CL_2mass), 3)
A_R_s6_2mass = np.round(15.87*(J_mean_s6_2mass - J_mean_CR_2mass), 3)

# Give out the results
print('A values for left control field:', A_L_s1_2mass, A_L_s2_2mass, A_L_s3_2mass,
      A_L_s4_2mass, A_L_s5_2mass, A_L_s6_2mass)
print('A values for right control field:', A_R_s1_2mass, A_R_s2_2mass, A_R_s3_2mass,
      A_R_s4_2mass, A_R_s5_2mass, A_R_s6_2mass)

# Write color excess into dataframe
pd_RA_DEC_2mass['Jmag_correctedL'] = pd_RA_DEC_2mass['Jmag'] - J_mean_CL_2mass
pd_RA_DEC_2mass['Jmag_correctedR'] = pd_RA_DEC_2mass['Jmag'] - J_mean_CR_2mass

# Calculate extinction via dataframe
A_CL = 15.87*pd_RA_DEC_2mass['Jmag_correctedL']
A_CR = 15.87*pd_RA_DEC_2mass['Jmag_correctedR']

# Plot extinction map for left control field 
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_2mass['RA'], pd_RA_DEC_2mass['DEC'], s=20,
            c=A_CL, cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax)
plt.close()
# Calculate extinction map for the right control field
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_2mass['RA'], pd_RA_DEC_2mass['DEC'], s=20,
            c=A_CR, cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax, label='A_J')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.close()

'''
# Test for high resolution data in gaia, this part takes about 2h to execute, only
#    do so if really needed
print('TEST')
gaia_control = len(create_subset_gaia(196, -76.65, 0.5, 0.25)['RA'])

# Iterate over all stars, create one field centered at each star, get extinction 
#    for it, put that into list
Av_gaia = []
empty_area = []
np_pd_RA_DEC = pd_RA_DEC.to_numpy()
for k in range(len(pd_RA_DEC['RA'])):
    single_RA = np_pd_RA_DEC[k][0]
    single_DEC = np_pd_RA_DEC[k][1]
    single_area = create_subset_gaia(single_RA, single_DEC, 0.5, 0.25)
    #counts_gaia.append(len(single_area['RA']))
    if len(single_area) != 0.0:
        Av = total_extinction(b_gaia, gaia_control, len(single_area))/0.83627
        Av_gaia.append(Av)
    elif len(single_area) == 0.0:
        empty_area.append(k)
        Av_gaia.append(np.nan)
    else:
        pass
        
# Plot the extinction per star as a colormap
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC['RA'], pd_RA_DEC['DEC'], s=1,
            c=Av_gaia, cmap='viridis')
#color_lims = mpl.cm.ScalarMappable(norm=None, cmap='viridis')
plt.colorbar(None, ax=ax, label='A_V')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.show()


print('TEST_APASS')
apass_control = create_subset_apass(196, -76.65, 2, 1)
apass_mean, apass_std = mean_std_fields(apass_control)

# Iterate over all stars, create one field centered at each star, get extinction 
#    for it, put that into list
Av_apass = []
empty_area = []
np_pd_RA_DEC_apass = pd_RA_DEC_apass.to_numpy()
for k in range(len(pd_RA_DEC_apass['RA'])):
    single_RA = np_pd_RA_DEC_apass[k][0]
    single_DEC = np_pd_RA_DEC_apass[k][1]
    single_area = create_subset_apass(single_RA, single_DEC, 2, 1)
    if len(single_area) != 0.0:
        mean, sdt = mean_std_fields(single_area)
        Av = np.round(3.1*(mean - apass_mean), 3)
        Av_apass.append(Av)
    elif len(single_area) == 0.0:
        empty_area.append(k)
        Av_apass.append(np.nan)
    else:
        pass
        
# Plot the extinction per star as a colormap
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_apass['RA'], pd_RA_DEC_apass['DEC'], s=100,
            c=Av_apass, cmap='viridis')
plt.colorbar(None, ax=ax, label='A_V')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.show()

'''

print('TEST_2MASS')
mass_control = create_subset_2mass(196, -76.65, 1.5, 0.75)
mass_mean, mass_std = mean_std_fields(mass_control)

# Iterate over all stars, create one field centered at each star, get extinction 
#    for it, put that into list
Av_2mass = []
empty_area = []
np_pd_RA_DEC_2mass = pd_RA_DEC_2mass.to_numpy()
for k in range(len(pd_RA_DEC_2mass['RA'])):
    single_RA = np_pd_RA_DEC_2mass[k][0]
    single_DEC = np_pd_RA_DEC_2mass[k][1]
    single_area = create_subset_2mass(single_RA, single_DEC, 1.5, 0.75)
    single_mag = np_pd_RA_DEC_2mass[k][2]
    if len(single_area) != 0.0:
        mean, sdt = mean_std_fields(single_area)
        Av = np.round(15.87*(mean - mass_mean), 3)
        Av_2mass.append(Av)
    elif len(single_area) == 0.0:
        empty_area.append(k)
        Av_2mass.append(np.nan)
    else:
        pass
        
# Plot the extinction per star as a colormap
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
plt.scatter(pd_RA_DEC_2mass['RA'], pd_RA_DEC_2mass['DEC'], s=20,
            c=Av_2mass, cmap='viridis')
plt.colorbar(None, ax=ax, label='A_V')
plt.xlabel('RA', fontsize=25)
plt.ylabel('DEC', fontsize=25)
plt.show()


#print(cols_2mass.names)
sns.set_style('whitegrid')
#sns.kdeplot(x=pd_RA_DEC_2mass['RA'], y=pd_RA_DEC_2mass['DEC'], cmap="Reds", shade=True, bw_adjust=.5)
plt.close()

#print(cols_2mass.names)
sns.set_style('whitegrid')
np_RA_apass = pd_RA_DEC_apass['RA'].to_numpy()
new_np_RA_apass = np_RA_apass.byteswap().newbyteorder()
np_DEC_apass = pd_RA_DEC_apass['DEC'].to_numpy()
new_np_DEC_apass = np_DEC_apass.byteswap().newbyteorder()
#sns.kdeplot(x=new_np_RA_apass, y=new_np_DEC_apass, cmap="Reds", shade=True, bw_adjust=.5)
plt.close()











