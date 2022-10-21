#!/usr/bin/env python3.6
"""
getcztmkfpar.py: Computes average values of various parameters of 
AstroSat CZTI orbit like lat/long/elv for time bins specified by the 
user, using values available in the MKF file. It also computes 
additional parameters like geomagnetic lat/long and cut-off rigidity 

Mithun N P S (06/12/2020)

Revisions:

v0.1 06/12/2020  Initial version
v0.2 15/12/2020  Time output is now T-Tstart; Added computation of mag lat/lon
                 by dipole approximation making AACGM lat/lon optional output. 
                 AACGM lat/lon now computed for altitude of zero.
                 COR computation done for DP maglat/lon and optionally for 
                 AACGM lat/lon.
v0.3 18/01/2021  Modified handling of tstart to accept it if 5000<tstart<mkfstart.
                 In that case output file not to have entries for time bins without 
                 any or partial MKF data.
                 Included computation of SAA crossing times based on lat/lon, optional 
                 output to be generated with -s / -saatimes option.
                 Included computation of Zenith distance of Sun, seem to have some 
                 discrepancy with SUNELV available in the MKF file, to be verified. 

Additions:

Vipul Prasad M & Rahul G (12/08/2022)

           --    Binned noise cleaned veto spectrum.
           --    Functionality for selecting desired interval between SAA crosssing times.

"""


from astropy.table import Table, hstack
from astropy.table import Column, MaskedColumn
import numpy as np
import argparse, sys
from datetime import datetime,timedelta
import aacgmv2
from scipy.spatial.transform import Rotation as Rot
import matplotlib.pyplot as plt
import pandas as pd
import Astrosat_Time as astime

from astropy.time import Time
from astropy.coordinates import EarthLocation,SkyCoord,get_sun,AltAz
from astropy import units as u
from astropy.io import fits
# To avoid autodownload of earth rotation files for altaz transform
from astropy.utils import iers
iers.conf.auto_download = False


# Function to compute mean of longitudes by getting the mean of
# corresponding complex numbers and getting back to angle

def longitude_mean(long_deg):
    mean_Z = np.sum(np.exp(1j * long_deg * np.pi/180.0))/float(len(long_deg))
    return np.angle(mean_Z, deg=True)

# Function to convert geographic lat/lon to geomagnetic lat/lon considering
# dipole approximation. Convert geo lat/lon to cartesian coordinates and then
# transform to magnetic coordinates by two rotations. First about z-axis from
# 0 deg lon to Magnetic North Pole (NP) meridian and then about the new y-axis
# from Geographic NP latitude to Magnetic NP latitude. Convert the transformed
# coordinates back to spherical system (lat/lon)

def geo2dipolemag(glat,glon):

    # Lat/lon of magnetic north pole
    NPlon=287.32  # 72.68 W or 287.32
    NPlat=80.65

    glat=np.array(glat)*np.pi/180.0
    glon=np.array(glon)*np.pi/180.0

    mlat=[]
    mlon=[]

    for i,templat in enumerate(glat):
        xyz=[np.cos(glat[i])*np.cos(glon[i]),np.cos(glat[i])*np.sin(glon[i]),np.sin(glat[i])] # To cartesian
        r1 = Rot.from_euler('z', -NPlon, degrees=True)  # First rotation
        r2 = Rot.from_euler('y', -(90.0-NPlat), degrees=True)   # Second rotation
        transxyz= r2.apply(r1.apply(xyz))   # Apply both rotations
        # Get mag lat/lon from transformed cartesian coordinates
        mlat.append(np.arctan2(transxyz[2],np.sqrt(transxyz[1]**2+transxyz[0]**2))*180.0/np.pi)
        mlon.append(np.arctan2(transxyz[1],transxyz[0])*180.0/np.pi)

    return([mlat,mlon])

def veto_noiseclean(bandlist, time_list):

    #function will clean the data.
    #bandlist = list of data of all the energy bands

    new_bandlist = np.copy(bandlist)

    for i in range(5,len(time_list)-5):

        avg_list = []
        flag_list = []
        # For a particluar time, calculate the average count of former 5 sec and later 5 sec.
        # lets say it is ri,
        # if it is higher than 10 sigma (can be changed) as compared to avi (that is (ri-avi) > 10*sqrt(avi)) then flag it.
        # If all three rates (for the 3 energy ranges) are flagged high, that is genuine. If the flag is for one or
        # two energy channels, deem to be anomalous and replace every ri with the corresponding avi. 

        for b in range(len(bandlist)):
            temp1 = bandlist[b][i-5:i]
            temp2 = bandlist[b][i+1:i+6]
            temp = [*temp1,*temp2]

            avg = sum(temp)/len(temp)
            avg_list.append(avg)

            if ((bandlist[b][i] - avg) > sigma*np.sqrt(avg)):
                flag_list.append(1)

        if ((len(flag_list) != len(bandlist))&(len(flag_list) != 0)):
            for b in range(len(bandlist)):
                new_bandlist[b][i] = avg_list[b]

    return new_bandlist


def get_evt_bandwise_data(evtfname, start_time, end_time):
    # Gives the binned noise-cleaned veto counts in different bands.
    #
    # Declaring table.
    veto_full_data = Table()

    # Read evt file 
    vetoSpec = fits.open(evtfname)[5].data

    # selecting veto data 
    vetoSpec = vetoSpec[(vetoSpec['time'] >= start_time) & (vetoSpec['time'] <= end_time)]
    veto_start = vetoSpec['Time'][0]#int(vetoSpec['Time'][0])   
    veto_end = vetoSpec['Time'][-1]#np.floor(vetoSpec['Time'][-1]+1)

    tbin_edges = np.array([start_time+(i*tbin) for i in range(0, int((end_time - start_time)/tbin)+1)],)
    tmid = 0.5*(tbin_edges[1:] + tbin_edges[:-1])

    veto_full_data.add_column(np.array(tmid), name= 'Time(veto)')

    for quad in quad_list:

        vetoSpec_quad = vetoSpec['VetoSpec'][vetoSpec['QuadID'] == quad]
        vetoSpec_quad_time = vetoSpec['Time'][vetoSpec['QuadID'] == quad]


        # Creating channnel edges corresponding to energy bin edges.
        # ph = ( energy - offset)/gain 
        channel_edges_quad = [int((float(val) - quad_offset[quad])/quad_gain[quad]) for val in ebin_edges]

        vetoSpect_band_list = []

        for i in range(len(channel_edges_quad)-1):
            # creating total counts in diffent bands.
            vetoSpec_band = np.sum(vetoSpec_quad[:, channel_edges_quad[i]:channel_edges_quad[i+1]], axis = 1)
            vetoSpect_band_list.append(vetoSpec_band)


        # noise cleaning 
        vetoSpect_band_list_noise_clean = veto_noiseclean(vetoSpect_band_list, vetoSpec_quad_time)
        
        
        for b,band in enumerate(vetoSpect_band_list_noise_clean):
            binned_band = []
            binned_band_time = []

            for i in range(len(tbin_edges)-1):
                
                data = pd.DataFrame({'time':vetoSpec_quad_time, 'counts':band})
                condition = (data['time'] >= tbin_edges[i]) & (data['time'] < tbin_edges[i+1])
                interval = len(data['time'][condition])
                counts_in_bin = sum(data['counts'][condition])
                if counts_in_bin != 0:
                    binned_band.append(counts_in_bin/interval)
        
                else:
                    binned_band.append(0) # null values to keep track of time bins with no events.


            veto_full_data.add_column(binned_band, name = 'Eband_{}_Quad_{}'.format(b+1, quad))

    return veto_full_data, veto_start, veto_end


def plot_veto(tabledata):

    bands = ['{} - {} KeV'.format(ebin_edges[i], ebin_edges[i+1]) for i in range(len(ebin_edges) - 1)]

    plt.figure(figsize = (10, 12))
    plt.subplots_adjust(hspace=0.4)
    
    for b in range(len(ebin_edges) - 1):
        plt.subplot(3,1,b+1)
        plt.title('Energy band {}'.format(bands[b]))
        for quad in quad_list:
            plt.plot(tabledata['Time(veto)'], tabledata['Eband_{}_Quad_{}'.format(b+1, quad)], label = 'Q_{}'.format(quad))
            plt.ylabel('Counts rate')
            plt.legend(loc = 'upper right')
    plt.xlabel('Time')
    plt.savefig('{}_veto_plots.pdf'.format(args.outfile))
    plt.show()

def get_mkf_params(tstart, tstop, tbin):

    #tbin_edges = np.linspace(start_time, end_time, int((end_time -start_time)/tbin))
    
    # Select the MKF table within tstart to tstop
    ind=((mkfdata['TIME'] >= tstart) & (mkfdata['TIME'] <= tstop))
    selmkfdata=mkfdata[ind]

    # Mask obviously junk CPM Rates above threshold (set to 10,000)
    selmkfdata['CPM_Rate']= MaskedColumn(selmkfdata['CPM_Rate'], mask=(selmkfdata['CPM_Rate'] > 10000.0))

    # Get time bin number for each row of selected table
    # and get the average value table (groups rows having same
    # time bin and then average the values for each group)

    mkftbin=np.trunc((selmkfdata['TIME']-tstart)/tbin)
    # print(mkftbin)
    grpselmkfdata=selmkfdata.group_by(mkftbin)
    binmkfdata = grpselmkfdata.groups.aggregate(np.mean)


    # Earth longitude needs special treatment, arithmetic mean will give incorrect
    # values while dateline crossing. Compute the mean for long correctly
    # and replace in the binned table
    # MAY NEED SIMILAR TREATMENT FOR Roll_RA/ELV/SUN_ANGLE FOR OTHER DATA SETS,
    # NEED TO CHECK.

    binearthlong = grpselmkfdata['EARTHLON'].groups.aggregate(longitude_mean)
    binmkfdata['EARTHLON']=binearthlong

    # Compute Magnetic lat/lon and COR by dipole approximation

    [mlat, mlon] = geo2dipolemag(np.array(binmkfdata['EARTHLAT']),np.array(binmkfdata['EARTHLON']))
    cor_dp=14.5*((np.cos(np.array(np.array(mlat)*np.pi/180.0)))**2.0)/(1.1**2.0)
    #denom=(1+np.sqrt(1+(np.cos(np.array(np.array(mlat)*np.pi/180.0)))**3.0))**2.0
    #cor_dp2=(60.0/1.1**2.0)*((np.cos(np.array(np.array(mlat)*np.pi/180.0)))**4.0)/denom

    binmkfdata.add_column(np.array(mlat),name='DPMLAT')
    binmkfdata.add_column(np.array(mlon),name='DPMLON')
    binmkfdata.add_column(cor_dp,name='COR_DP')

    # Compute AACGM lat/lon from geo lat/lon as well as COR and add them to the
    # MKF table if input argument --aacgm is set

    if args.aacgm:

        dtimeStart=datetime(2010,1,1)+timedelta(seconds=tstart)     ## Start time as epoch

        [aacgmlat,aacgmlon,mlt]=aacgmv2.get_aacgm_coord_arr(binmkfdata['EARTHLAT'], \
                        binmkfdata['EARTHLON'], np.zeros(len(binmkfdata['EARTHLON'])), dtimeStart)

        binmkfdata.add_column(aacgmlat,name='AACGMLAT')
        binmkfdata.add_column(aacgmlon,name='AACGMLON')

        # Compute COR as per Eqn A1 Tawa et al 2008, Suzaku
        # 1.1 is approximate radius (600 km+6371 km) in units of R_earth

        cor_aacgm=14.5*((np.cos(np.array(binmkfdata['AACGMLAT']*np.pi/180.0)))**2.0)/(1.1**2.0)
        binmkfdata.add_column(cor_aacgm,name='COR_AACGM')

    # Compute Zenith distance of Sun: Get RA,DEC of Sun at each time bin and then
    # transform this to altaz coordinates for the lat/lon/alt of satellite, zenith distance
    # is 90-altitude

    refdtime=tobs=datetime(2010,1,1)

    sun_zdist=[]
    sang=[]

    for tbin,tobs in enumerate(binmkfdata['TIME']):
        obsTime=Time(refdtime+timedelta(seconds=tobs))
        coord=get_sun(obsTime)
        obsloc=EarthLocation(lat=binmkfdata['EARTHLAT'][tbin],lon=binmkfdata['EARTHLON'][tbin],\
                height=binmkfdata['ALTITUDE'][tbin]*u.km)
        azcoord=coord.transform_to(AltAz(location=obsloc, obstime=obsTime))
        sun_zdist.append(90.0-azcoord.alt.value)
        #srccoord=SkyCoord(ra=binmkfdata['Roll_RA'][tbin]*u.degree,dec=binmkfdata['Roll_DEC'][tbin]*u.degree,frame='icrs')
        #sang.append(coord.separation(srccoord))
        #print(coord.ra.value,coord.dec.value) #,azcoord.az.value,azcoord.alt.value,srccoord.separation(coord).value)

    binmkfdata.add_column(np.array(sun_zdist),name='SUN_ZDIST')

    # Change Time to Time-Tstart

    #binmkfdata['TIME']-=tstart

    # Now keep only columns that are needed in the final table
    # and write to output text file. Add/remove column names
    # to be written out in the list below

    if args.aacgm:
        binmkfdata = binmkfdata['TIME', 'Roll_RA', 'Roll_DEC', 'ROLL_ROT','EARTHLAT','EARTHLON',\
                'SUN_ANGLE','CPM_Rate','ELV','DPMLAT','DPMLON','COR_DP','AACGMLAT','AACGMLON','COR_AACGM','SUN_ZDIST','SUNELV']

    else:
        binmkfdata = binmkfdata['TIME', 'Roll_RA', 'Roll_DEC', 'ROLL_ROT','EARTHLAT','EARTHLON',\
                'SUN_ANGLE','CPM_Rate','ELV','DPMLAT','DPMLON','COR_DP','SUN_ZDIST','SUNELV']

    #binmkfdata.write(args.outfile,format='ascii.fixed_width',overwrite=True,delimiter=' ')
    return binmkfdata#, np.array(list(set(mkftbin)))


# function to enable clicking to select interval from plot.

# def on_click(event):
#     global xpoint
#     xpoint = event.xdata
#     temp = []
#     temp.append(xpoint)
#     if len(temp) == 1:
#         fig.canvas.mpl_disconnect(cid)
#         plt.close(1)
#     return xpoint




# Parse user input parameters



parser = argparse.ArgumentParser(epilog="""
    Function: Get average values of MKF parameters including mag lat/lon & COR for specified
    time duration and binning. Output parameters are written to text output file""")
parser.add_argument("mkffile", help="Input MKF FITS file", type=str)
parser.add_argument("evtfile", help="Input EVT FITS file", type=str)
parser.add_argument("outfile", help="Output text file name", type=str)
parser.add_argument("tbin", help="Time bin size (s)", type=float)
parser.add_argument("--sigma", help="standard deviation for noise clean(sigma)", type=float, default = 10)
parser.add_argument("--quad", help="Specify needed quadrants [eg and default value: 0,1,2,3]", type=str, default = '0,1,2,3')
parser.add_argument("-a", "--aacgm", action="store_true", help="Include AACGM lat/lon values")
parser.add_argument('--ebin_edges',type = str, help = "eg and default value: 100.0,250.0,350.0,500.0", default = '100.0,250.0,350.0,500.0')
parser.add_argument('-np', '--no_plot', help = 'Specify this to not plot the veto spectrum', action='store_false')

args = parser.parse_args()


# parameters for functions get_evt_bandwise_data and veto_noiseclean

quad_gain = [5.591, 5.594, 5.943, 5.222]
quad_offset = [-56.741, -41.239, -41.682, -26.528]
evtfname = args.evtfile
ebin_edges = args.ebin_edges.split(',')
num_ebands = len(ebin_edges)-1
sigma = args.sigma
tbin = int(args.tbin)
quad_list = [int(i) for i in args.quad.split(',')]



## Compute SAA crossing times

# SAA crossing is defined as the time when the spacecraft crosses
# the line lon=(((lat+6)*slope)-50.0) with slope=20.0/12.0 in the
# latitude, longitude space. Crossing is identified by evaluating
# the line equation with two successive lat,lon and when it changes
# sign from negative to positive.

mkfdata = Table.read(args.mkffile)

slope=20./12.0
lsign=np.sign(mkfdata['EARTHLON']-(((mkfdata['EARTHLAT']+6)*slope)-50.0))

saa_index=[]
for n in range(0,len(lsign)-1):
    if((lsign[n] < 0) & (lsign[n+1] > 0)):
        saa_index.append(n)

saaCrossTimes=np.array(mkfdata['TIME'][saa_index])
tdiff=np.ediff1d(saaCrossTimes)

# Prepend and append SAA cross times before tstart and after tstop; assumes orbital period of 6300 s
# if only one time is available, otherwise the different of first two and last two SAA times are used as
# orbital periods.

if (len(saaCrossTimes)==1):
    saaCrossTimes=np.concatenate(([saaCrossTimes[0]-6300.0],saaCrossTimes,[saaCrossTimes[0]+6300.0]))
else:
    saaCrossTimes=np.concatenate(([saaCrossTimes[0]-tdiff[0]],saaCrossTimes,[saaCrossTimes[-1]+tdiff[-1]]))

#Save output to a text file named MKFBasename_saaCrossTimes.txt
#saatimefile=args.mkffile[:-4]+'_saaCrossTimes.txt'
#np.savetxt(saatimefile,saaCrossTimes-tstart,fmt="%f")

x = mkfdata['TIME']
y = mkfdata['CPM_Rate']


int_no = 0
for i in range(1, len(saaCrossTimes)):
    utctime_start = astime.convertAStoAll(float(saaCrossTimes[i-1]))[2]
    utctime_end = astime.convertAStoAll(float(saaCrossTimes[i]))[2]
    print('\nInterval: {} = {}  - {}'.format(int_no, utctime_start, utctime_end))
    int_no += 1

fig = plt.figure()
plt.plot(x,y)

# plot with saa crossing times and cpm rate.
for i in range(0, len(saaCrossTimes)-1):
    centre = (saaCrossTimes[i+1] - saaCrossTimes[i])/2.0
    #print(saaCrossTimes[i])
    plt.axvline(saaCrossTimes[i], color = 'black')
    plt.text(saaCrossTimes[i]+centre , (np.max(mkfdata['CPM_Rate']) - (1/10.0)*np.max(mkfdata['CPM_Rate'])), s = '{}'.format(i))


plt.axvline(saaCrossTimes[-1], color = 'black')
# cid = fig.canvas.mpl_connect('button_press_event', on_click)
plt.title('Not down required interval index')
plt.xlabel('TIME')
plt.ylabel('CPM_Rate')
plt.show()


saa_start = int(input("\nSelect interval: "))

start_time = saaCrossTimes[saa_start]
end_time = saaCrossTimes[saa_start + 1]

print("\nSelected start time: ", astime.convertAStoAll(start_time)[2])
print("Selected end time: ", astime.convertAStoAll(end_time)[2])


# if saa selected interval is first or last.
# make start time and end_time the earliest data available time.
if start_time < mkfdata['TIME'][0]:
    start_time = mkfdata['TIME'][0]

if end_time > mkfdata['TIME'][-1]:
    end_time = mkfdata['TIME'][-1]



veto_bandwise_data, veto_start, veto_end = get_evt_bandwise_data(args.evtfile, start_time, end_time)
mkf_param_data = get_mkf_params(start_time,end_time, tbin)
#mkf_param_data = mkf_param_data[(mkf_param_data['TIME'] >= veto_start) & (mkf_param_data['TIME'] <= veto_end)]
combined_data = hstack([veto_bandwise_data, mkf_param_data])
combined_data.write(args.outfile,format='ascii.fixed_width',overwrite=True,delimiter=' ')

if args.no_plot is True:
    plot_veto(veto_bandwise_data)

print(combined_data['Time(veto)', 'TIME'][0:5])

#----------------------- END OF CODE ------------------------------------------------#
