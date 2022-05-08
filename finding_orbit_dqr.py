#! /usr/bin/env python3
#Vipul

# Finds orbit number corresponding that covers a time using the data quality report page.

import argparse
import pandas as pd
import numpy as np
import datetime as dt
from astropy.time import Time
import astropy.units as u
from time import time

def convert_time(time):
    time = time[0:18]
    try: 
        time_formatted = dt.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')    
    except:
        time_formatted = dt.datetime.strptime(time, '%Y-%m-%d %H-%M-%S')

    return time_formatted

parser	=	argparse.ArgumentParser()
#parser.add_argument("triggertime", help="Triggertime in Astrosat seconds, eg: 360158420.0", type=float)
parser.add_argument("--tast",type=float,help='Trigger time in astrosat seconds')
parser.add_argument("--tutc",type=str,help='Trigger time in UTC')
args = parser.parse_args()

run_start = time()
if args.tast == None:
    trigtime_utc  = convert_time(args.tutc)
else:
    #trigtime_utc = Time('2010-01-01') + (args.tast+3)*u.second
    total_seconds = (dt.datetime.strptime('2010-01-01', '%Y-%m-%d').timestamp() + args.tast + 3)
    trigtime_utc = dt.datetime.fromtimestamp(total_seconds)#.strftime('%Y-%m-%d %H-%M-%S')

data = pd.read_html("https://www.iucaa.in/~astrosat/czti_dqr/", header = 0)[0]
#data.to_csv("dqr.csv")
#data = pd.read_csv("dqr.csv")
#data = pd.read_csv("orbitinfo.csv", names = ['Folder', 'OBSID', 'Observer', 'Object', 'RA', 'Dec', 'Exposure time', 'Date/time start', 'Date/time end'], skiprows = 1)
data['Date/time start'] = data['Date/time start'].apply(convert_time)
data['Date/time end'] = data['Date/time end'].apply(convert_time)

orbit_data = data[data["Folder"].str.len() == 43]
orbit_data = orbit_data.dropna()

observation_data = data[data['Folder'].str.len() != 43]
observation_data = observation_data.dropna()

prev_index = 0
obs_index = 0
for index, row in observation_data.iterrows():
    if (trigtime_utc >= row['Date/time start']) and (trigtime_utc <= row['Date/time end']):
        #print(index, row['Folder'])
        obs_index = prev_index
        #print(obs_index)
        break
    prev_index = index

num = 1
for orbit_index, row in orbit_data.loc[obs_index:].iterrows():
    #print(orbit_index, row['Folder'])
    if (trigtime_utc >= row['Date/time start']) and (trigtime_utc <= row['Date/time end']):
        #print('Folder name: ',row['Folder'][0:-6])
        #print('Orbit: ', row['Folder'][-5:])
        print('Path {}: '.format(num), '/data2/czti/level2/'+row['Folder'][0:-6]+'/czti/orbit/'+row['Folder'][-5:], '({} - {})'.format(row['Date/time start'], row['Date/time end']))
        num+=1

run_end = time()
print(run_end - run_start)
        


