import os
import argparse

parser	=	argparse.ArgumentParser()
parser.add_argument("grbname", help="Name of the GRB(eg:GRB211213A)", type=str)
parser.add_argument("orbit_folder", help="Folder name of the orbit", type=str)
args = parser.parse_args()

GRB_name = args.grbname
caldbbadpix='/home/cztipoc/CALDB/data/as1/czti/bcf/AS1cztbadpix20160908v01.fits'



#foldername = input("Enter the folder path (eg:/data2/czti/level2/20180313_T02_013T01_9000001976_level2/czti/orbit/13290): ")

observation_folder = args.orbit_folder  #foldername.split('/')[4]
#orbit_no = foldername.split('/')[-1]
rootname = 'AS1'+observation_folder[9:30]


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
mkfthreshold=mkfThresholds.txt # str(os.getcwd())+'/mkfThresholds.txt'
outgtifile=level2_dir+'/modeM0/'+rootname+'cztM0_level2.gti'

cztgitgen_cmd = 'cztgitgen inevtfile mkffile mkfThreshold outgtifile - y y'

os.system(cztgitgen_cmd)


print('########## cztgaas ##########\n')

evtdatafile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc.evt'
mkffile=level2_dir+'/'+rootname+'czt_level2.mkf'
outaspect=level2_dir+'/modeM0/'+rootname+'cztM0_level2.aspect'

cztgaas_cmd = 'cztgaas {} _ _ _ _ {} {} y y'.format(evtdatafile, mkffile, outaspect)

os.system(cztgaas_cmd)


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


print('########### cztpixclean ###########')

inlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_bc_livetime.fits'
det_th_1, det_th_2 = 200, 400
pix_th_1, pix_th_2 = 4, 10

if gtitype == 'QUAD':
	inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds.evt'

	outevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc_200_4.evt'
	outlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_livetime_200_4.fits'
	outbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix_200_4.fits'
	outdblevt=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_200_4.dblevt'

	outevtfile_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc_400_10.evt'
	outlivetime_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_livetime_400_10_400_10.fits'
	outbadpix_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix_400_10.fits'
	outdblevt_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_400_10.dblevt'

else:
	inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_bc_ds.evt'

	outevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_bc_ds_pc_200_4.evt'
	outlivetime=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_livetime_200_4.fits'
	outbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_badpix_200_4.fits'
	outdblevt=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_200_4.dblevt'

	outevtfile_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_bc_ds_pc_400_10.evt'
	outlivetime_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_livetime_400_10.fits'
	outbadpix_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_badpix_400_10.fits'
	outdblevt_400_10=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_400_10.dblevt'



cztpixclean_cmd_1 = 'cztpixclean par_infile={} par_inlivetimefile={} par_writedblevt=y par_outfile1={} par_outfile2={} par_outlivetimefile={} par_badpixfile={} par_nsigma=5 par_det_tbinsize=1 par_pix_tbinsize=1 par_det_count_thresh={} par_pix_count_thresh={}'.format(inevtfile, inlivetime, outevtfile, outdblevt, outlivetime, outbadpix, det_th_1, pix_th_1)

cztpixclean_cmd_2 = 'cztpixclean par_infile={} par_inlivetimefile={} par_writedblevt=y par_outfile1={} par_outfile2={} par_outlivetimefile={} par_badpixfile={} par_nsigma=5 par_det_tbinsize=1 par_pix_tbinsize=1 par_det_count_thresh={} par_pix_count_thresh={}'.format(inevtfile, inlivetime, outevtfile_400_10, outdblevt_400_10, outlivetime_400_10, outbadpix_400_10, det_th_2, pix_th_2)

os.system(cztpixclean_cmd_1)
os.system(cztpixclean_cmd_2)

print('######### cztflagbadpix ################')

if gtitype == 'QUAD':
	inbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix.fits'
	outbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix.fits_'

else:
	inbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_badpix.fits'
	outbadpix=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_badpix.fits_'

cztflagbadpix_cmd = 'nbadpixFiles=1 {} {} clobber=y history=y'.format(inbadpix, outbadpix)

os.systme(cztflagbadpix_cmd)

print('####### cztevtclean ########\n')

if gtitype == 'QUAD':

	infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_bc_ds_pc.evt'
	outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.evt'
	infile_dbl=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad.dblevt'
	outfile_dbl=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.dblevt'
else:
	infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_bc_ds_pc.evt'
	outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.evt'

cztevtclean_cmd_1='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"'.format(infile, outfile)
cztevtclean_cmd_2='cztevtclean infile={} outfile={} alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"'.format(infile_dbl, outfile_dbl)

os.system(cztevtclean_cmd)


print('####### cztdpigen #########\n')

if gtitype == 'QUAD':

	inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.evt'
	inbadpixfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_badpix.fits'
	outdphfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.dph'
	outdpifile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.dpi'

else:

	inevtfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.evt'
	inbadpixfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_badpix.fits'
	outdphfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.dph'
	outdpifile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.dpi'


cztdpigen_cmd = 'cztdpigen par_infile={} par_badpixFile={} par_badpixThreshold=0 par_outDPHfile={} par_outDPIfile={} par_quadsToProcess=0,1,2,3 par_ebins=- par_timerange=- par_clobber=y par_history=y'.format(inevtfile, inbadpixfile, outdphfile, outdpifile)

os.system(cztdpigen_cmd)


print('####### cztimage #########\n')

intype='dpi'

aspectfile_q0=level2_dir+'/modeM0/'+rootname+'cztM0_level2.aspect_Q0'
aspectfile_q1=level2_dir+'/modeM0/'+rootname+'cztM0_level2.aspect_Q1'
aspectfile_q2=level2_dir+'/modeM0/'+rootname+'cztM0_level2.aspect_Q2'
aspectfile_q3=level2_dir+'/modeM0/'+rootname+'cztM0_level2.aspect_Q3'

if intype == 'dpi':
	if gtitype == 'QUAD':

		infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.dpi'
		outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.img'
	else:
		infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.dpi'
		outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.img'
else:
	if gtitype == 'QUAD':

		infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.dph'
		outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_quad_clean.img'
	else:
		infile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.dph'
		outfile=level2_dir+'/modeM0/'+rootname+'cztM0_level2_common_clean.img'



