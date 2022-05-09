import os
import subprocess
import argparse
import datetime as dt
import pandas as pd


parser	=	argparse.ArgumentParser()
parser.add_argument("grbname", help="Name of the GRB(eg:GRB211213A)", type=str)
parser.add_argument("trigtime", help="Trigger time in astrosat seconds", type=float)
#parser.add_argument("orbit_folder", help="Folder name of the orbit", type=str)
args = parser.parse_args()

GRB_name = args.grbname
trigtime = args.trigtime

def convert_time(time):
    try: 
        time_formatted = dt.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')    
    except:
        time_formatted = dt.datetime.strptime(time, '%Y-%m-%d %H-%M-%S')

    return time_formatted


def find_orbit(time):

	total_seconds = (dt.datetime.strptime('2010-01-01', '%Y-%m-%d').timestamp() + time + 3)
	trigtime_utc = dt.datetime.fromtimestamp(total_seconds)
    #data = pd.read_html("https://www.iucaa.in/~astrosat/czti_dqr/", header = 0)[0]
	#data.to_csv("dqr.csv")
	data = pd.read_csv("dqr.csv")
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

	for orbit_index, row in orbit_data.loc[obs_index:].iterrows():
	    #print(orbit_index, row['Folder'])
	    if (trigtime_utc >= row['Date/time start']) and (trigtime_utc <= row['Date/time end']):
	        print('Folder name: ',row['Folder'])
	        path = '/data2/czti/level2/'+row['Folder'][0:-6]+'/czti/orbit/'+row['Folder'][-5:]
	        break
	return path




#foldername = input("Enter the folder path (eg:/data2/czti/level2/20180313_T02_013T01_9000001976_level2/czti/orbit/13290): ")

folder_path = find_orbit(args.trigtime).replace(level2, level1)
orbit_folder = folder_path.split('/')[4]
#orbit_no = foldername.split('/')[-1]
rootname = 'AS1'+orbit_folder[9:30]

print('Folder path: ',folder_path)
print('Root name: ',rootname)

os.system('mkdir {}'.format(GRB_name))
os.system('cp -r {}_V1.0 {}/level1'.format(folder_path, GRB_name))
os.system('mkdir {}/level2'.format(GRB_name))
os.system('mkdir {}/level2/aux'.format(GRB_name))
os.system('mkdir {}/level2/moseM0'.format(GRB_name))
os.system('mkdir {}/level2/moseSS'.format(GRB_name))
os.system('mkdir {}/level2/modeM9'.format(GRB_name))

level1_dir = GRB_name+'/level1'
level2_dir = GRB_name+'/level2'

print('########## cztscience2event############\n')

infile=level1_dir+'/modeM0/'+rootname+'cztM0_level1.fits'
TCTfile=level1_dir+'/aux/'+rootname+'czt_level1.tct'
outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2.fits'
hdrInfoFile=level2_dir+'/modeM0/'+rootname+'cztM0_level2.hdr'
bunchfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bunch.fits'

cztscience2event_cmd_M0 = 'cztscience2event infile={} TCTfile={} outfile={} hdrInfoFile={} bunchfile={} nPackets=100000 BigEndian=y clobber=y history=y debug=y'.format(infile, TCTfile, outfile, hdrInfoFile, bunchfile)

cztscience2event_cmd_SS = 'cztscience2event infile={} TCTfile={} outfile={} hdrInfoFile={} bunchfile={} nPackets=100000 BigEndian=y clobber=y history=y debug=y'.format(infile.replace('M0', 'SS'), TCTfile, outfile.replace('M0', 'SS'), hdrInfoFile.replace('M0', 'SS'), bunchfile.replace('M0', 'SS'))

os.system(cztscience2event_cmd_M0)
os.system(cztscience2event_cmd_SS)



print('########## cztbunchclean ############\n')

infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2.fits'
bunchfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bunch.fits'
outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc.fits'
livetimefile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc_livetime.fits'

cztbunchclean_cmd = 'cztbunchclean {} {} {} {} 30 20 1 clobber=y history=y'.format(infile, outfile, bunchfile, livetimefile)

os.system(cztbunchclean_cmd)




print('########## cztpha2energy ##########\n')

inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc.fits'
outevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc.evt'
tempfile=level2_dir+'/modeSS/'+rootname+'cztSS_level2.fits'

cztpha2energy_cmd = 'cztpha2energy {} _ _ {} {} y y'.format(inevtfile, outevtfile, tempfile)

os.system(cztpha2energy_cmd)



print('########## cztgtigen ##########\n')

inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc.evt'
mkffile=level2_dir+'/'+rootname+'czt_level2.mkf'
mkfthreshold='mkfThresholds.txt' # str(os.getcwd())+'/mkfThresholds.txt'
outgtifile=level2_dir+'/modeM0/'+rootname+'cztM0_level2.gti'

cztgtigen_cmd = 'cztgtigen {} {} {} {} - y y'.format(inevtfile, mkffile, mkfthreshold, outgtifile) 

call(cztgtigen_cmd)



print('######## cztdatasel #########\n')

inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc.evt'
gitfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2.gti'
gtitype='QUAD' #'COMMON'

if gtitype == 'QUAD':
	outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds.evt'
else:
	outfile=outfile_quad.replace(quad, common)

cztdatasel_cmd = 'cztdatasel {} {} {} {} y y'.format(inevtfile, gitfile, gtitype, outfile)

os.system(cztdatasel_cmd)


print('\n\n ## Running the T90 script #### \n\n')

print('########### cztpixclean ###########')

inlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc_livetime.fits'

caldbbadpix='/home/cztipoc/CALDB/data/as1/czti/bcf/AS1cztbadpix20160908v01.fits'


inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds.evt'

outevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc.evt'
outlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_livetime.fits'
outbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix.fits'
#outdblevt=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad.dblevt'

cztpixclean_cmd_t90 = 'cztpixclean par_infile={} par_inlivetimefile={} par_outfile1={} par_outlivetimefile={} par_badpixfile={} par_nsigma=5 par_det_tbinsize=1 par_pix_tbinsize=1 par_det_count_thresh=1000 par_pix_count_thresh=100'.format(inevtfile, inlivetime, outevtfile, outlivetime, outbadpix, det_th_1, pix_th_1)
os.system(cztpixclean_cmd_t90)

print('####### cztevtclean ########\n')

infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc.evt'
outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.evt'
#infile_dbl=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad.dblevt'
#outfile_dbl=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.dblevt'

cztevtclean_cmd_t90='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"'.format(infile, outfile)

t90_output=GRB_name+'/{}_T90'.format(GRB_name)

t90_cmd = 'python calculate_T90_test_Gfit.py {} {} {} {} {} --tmark {}'.format(outfile, caldbbadpix, outlivetime, t90_output, GRB_name, trigtime)

os.system(t90_cmd)



print('\n\n ## Running the Compton counts script ### \n\n')

print('########### cztpixclean ###########')

inlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc_livetime.fits'
det_th_1, det_th_2 = 200, 400
pix_th_1, pix_th_2 = 4, 10


inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds.evt'

outevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc_200_4.evt'
outlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_livetime_200_4.fits'
outbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix_200_4.fits'
outdblevt=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_200_4.dblevt'

outevtfile_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc_400_10.evt'
outlivetime_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_livetime_400_10_400_10.fits'
outbadpix_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix_400_10.fits'
outdblevt_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_400_10.dblevt'


cztpixclean_cmd_1 = 'cztpixclean par_infile={} par_inlivetimefile={} par_writedblevt=y par_outfile1={} par_outfile2={} par_outlivetimefile={} par_badpixfile={} par_nsigma=5 par_det_tbinsize=1 par_pix_tbinsize=1 par_det_count_thresh={} par_pix_count_thresh={}'.format(inevtfile, inlivetime, outevtfile, outdblevt, outlivetime, outbadpix, det_th_1, pix_th_1)

cztpixclean_cmd_2 = 'cztpixclean par_infile={} par_inlivetimefile={} par_writedblevt=y par_outfile1={} par_outfile2={} par_outlivetimefile={} par_badpixfile={} par_nsigma=5 par_det_tbinsize=1 par_pix_tbinsize=1 par_det_count_thresh={} par_pix_count_thresh={}'.format(inevtfile, inlivetime, outevtfile_400_10, outdblevt_400_10, outlivetime_400_10, outbadpix_400_10, det_th_2, pix_th_2)

os.system(cztpixclean_cmd_1)
os.system(cztpixclean_cmd_2)


print('####### cztevtclean ########\n')


infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc_200_4.evt'
outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean_200_4.evt'
infile_dbl=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_200_4.dblevt'
outfile_dbl=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean_200_4.dblevt'


infile_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc_400_10.evt'
outfile_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean_400_10.evt'
infile_dbl_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_400_10.dblevt'
outfile_dbl_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean_400_10.dblevt'



cztevtclean_cmd_1='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"'.format(infile, outfile)
cztevtclean_cmd_2='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="y" history="y"'.format(infile_dbl, outfile_dbl)
cztevtclean_cmd_3='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"'.format(infile_400_10, outfile_400_10)
cztevtclean_cmd_4='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="y" history="y"'.format(infile_dbl_400_10, outfile_dbl_400_10)

os.system(cztevtclean_cmd_1)
os.system(cztevtclean_cmd_2)
os.system(cztevtclean_cmd_3)
os.system(cztevtclean_cmd_4)


compton_output_1=GRB_name+'/{}_{}_{}'.format(GRB_name, det_th_1, pix_th_1)
compton_output_2=GRB_name+'/{}_{}_{}'.format(GRB_name, det_th_2, pix_th_2)


compton_cmd_1 = 'python latest_compton_events.py {} {} {} --tmark {}'.format(outfile_dbl, caldbbadpix, compton_output_1, trigtime)
compton_cmd_2 = 'python latest_compton_events.py {} {} {} --tmark {}'.format(outfile_dbl_400_10, caldbbadpix, compton_output_2, trigtime)


os.system(compton_cmd_1)
os.system(compton_cmd_2)
