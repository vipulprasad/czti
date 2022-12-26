import pandas as pd
import numpy as np
import datetime as dt

num_orbits = int(input("Enter the number of orbits to load:"))

#start_folder = input("Enter starting folder folder name(eg:20210726_T04_030T01_9000004588_level2_31503): ")
#end_folder = input("Enter the ending folder name(eg:20210726_T04_030T01_9000004588_level2_31503) or type a for all remaining orbits: ")


data = pd.read_html("https://www.iucaa.in/~astrosat/czti_dqr/")[0].head(num_orbits)
#data.to_csv("dqr.csv")

"""start_ind = dqr.index[dqr["Folder"] == start_folder].tolist()[0]
if end_folder == 'a':
	end_ind = 0
else:
	end_ind = dqr.index[dqr["Folder"] == end_folder].tolist()[0]"""

#data = dqr[end_ind: start_ind]
data = data[data["Folder"].str.len() == 43]
data = data.dropna()
#data.to_csv("dqr.csv")

noise_ratios = []

print('ObsID     OrbitID    Quadrant       DetID  PixID')

for folder in data['Folder']:

    orbit = folder[38:43]
    obs = folder[26:30] #folder[9:12]+'_'+
    year = folder[0:4]
    url = 'https://www.iucaa.in/~astrosat/czti_dqr/{}/{}/index.html'.format(year, folder)
    orbit_data = pd.read_html(url)

    if len(orbit_data) == 10:
        n = 5
    else:
        n = 4

    noise_tab = orbit_data[n][[orbit_data[n].columns[0], orbit_data[n].columns[6]]]
    pix_tab = orbit_data[n+1]['Pixel properties']



    for i in range(4):

    	if float(noise_tab['Noise dominated (detector-on time)'][i].strip('%')) >= 5.0:

    		noise_ratios.append([obs, orbit, noise_tab['Quadrant'][i]+'='+noise_tab['Noise dominated (detector-on time)'][i], pix_tab['DetID'][i*3], pix_tab['PixID'][i*3]])
    		#print('{}\t{}\t\t{}\t\t{}\t{}'.format(obs, orbit, noise_tab['Quadrant'][i]+'='+noise_tab['Noise dominated (detector-on time)'][i], pix_tab['DetID'][i*3], pix_tab['PixID'][i*3]))

noise_ratios = pd.DataFrame(noise_ratios, columns = ['ObsID', 'OrbitID', 'Quadrant', 'DetID', 'PixID'])
sorted_data = noise_ratios.sort_values(['OrbitID'], ascending = [True])
sorted_data['spc'] = '    '


print(sorted_data[['ObsID', 'spc', 'OrbitID', 'spc', 'Quadrant', 'spc', 'DetID', 'spc', 'PixID']].to_string(index=False, header = False))


"""
data['Date/time start'] = pd.to_datetime(data['Date/time start'], format = '%Y-%m-%d %H:%M:%S')
data['Date/time end'] = pd.to_datetime(data['Date/time end'], format = '%Y-%m-%d %H:%M:%S')
data['orbits'] = data['Folder'].str.slice(38,43).astype(int)
data['obs'] = data['Folder'].str.slice(26,30).astype(int)

sorted_data = data[['obs', 'orbits', 'Date/time start', 'Date/time end']].sort_values(['orbits', 'obs'], ascending = [False, False])

data_gap = []

for i in range(len(sort_data)):

	if sort_data['Date/time start'][i] > sort_data['Date/time end'][i]:"""
