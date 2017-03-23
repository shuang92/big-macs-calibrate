#!/usr/bin/python2
#
# to convert a csv file to a table fits file.
# usage: fits2csv.py cvsname.cvs [-c col1,col2,col3] fitsname.fits
# arguments:
# 	-c 	specify columns in the cvs file
# 		default = all columns
#

from astropy.table import Table
import argparse


################## specify the arguments ###################

parser = argparse.ArgumentParser()
parser.add_argument("input", help="input csv file, e.g. table.csv")
parser.add_argument("output", help="output fits file, e.g. table.fits")
parser.add_argument("-c", "--cols", default = 'all columns',
		help="[-c ra,dec,mag,etc] \
		specify colmuns in the input file.\
		Defaul = all columns")
args = parser.parse_args()

################## getting the arguments ###################

csvname = args.input
outname = args.output
columns = tuple(col for col in args.cols.split(',') if col.strip())
print 'extracting',columns,'from', csvname, 'to', outname

################ convert specified columns to csv #################  

t = Table.read(csvname, format='ascii.csv', guess=False)
if columns != ('all columns',):
	t = t[columns]

t.write(outname,format='fits')

exit(0)

