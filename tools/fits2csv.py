#!/usr/bin/python2
#
# to convert a table fits file to a csv file.
# usage: fits2csv.py fitsname.fit [-c col1,col2,col3] outname.csv
# arguments:
# 	-c 	specify columns in the fits file
# 		default = all columns
#

from astropy.table import Table
import argparse


################## specify the arguments ###################

parser = argparse.ArgumentParser()
parser.add_argument("input", help="input fits file, e.g. table.fits")
parser.add_argument("output", help="output csv file, e.g. table.csv")
parser.add_argument("-c", "--cols", default = 'all columns',
		help="[-c ra,dec,mag,etc] \
		specify colmuns in the input file.\
		Defaul = all columns")
args = parser.parse_args()

################## getting the arguments ###################

fitname = args.input
outname = args.output
columns = tuple(col for col in args.cols.split(',') if col.strip())
print 'extracting',columns,'from', fitname, 'to', outname

################ convert specified columns to csv #################  

t = Table.read(fitname)
if columns != ('all columns',):
	t = t[columns]

t.write(outname,format='ascii.csv')

exit(0)

