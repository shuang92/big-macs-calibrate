#!/usr/local/bin/python
""">> sqlcl << command line query tool by Tamas Budavari <budavari@jhu.edu>
Usage: sqlcl [options] sqlfile(s)

Options:
        -s url	   : URL with the ASP interface (default: pha)
        -f fmt     : set output format (html,xml,csv - default: csv)
        -q query   : specify query on the command line
        -l         : skip first line of output with column names
        -v	   : verbose mode dumps settings in header
        -h	   : show this message"""

formats = ['csv','xml','html']

astro_url='http://cas.sdss.org/astro/en/tools/search/x_sql.asp'
public_url='http://cas.sdss.org/public/en/tools/search/x_sql.asp'

astro_url='http://skyserver.sdss.org/dr13/en/tools/search/x_results.aspx?searchtool=SQL&'
public_url='http://skyserver.sdss.org/dr13/en/tools/search/x_results.aspx?searchtool=SQL&'


default_url=astro_url
default_fmt='csv'

def usage(status, msg=''):
    "Error message and usage"
    print __doc__
    if msg:
        print '-- ERROR: %s' % msg
    sys.exit(status)

def filtercomment(sql):
    "Get rid of comments starting with --"
    import os
    fsql = ''
    for line in sql.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep;
    return fsql

def query(sql,url=default_url,fmt=default_fmt):
    "Run query and return file object"
    import urllib
    fsql = filtercomment(sql)
    params = urllib.urlencode({'cmd': fsql, 'format': fmt})
    return urllib.urlopen(url+params)    

def panstarrs_ebv(lon, lat, coordsys='equ', mode='full'):
    import json, requests
    '''
    Send a line-of-sight reddening query to the Argonaut web server.
    
    Inputs:
      lon, lat: longitude and latitude, in degrees.
      coordsys: 'gal' for Galactic, 'equ' for Equatorial (J2000).
      mode: 'full', 'lite' or 'sfd'
    
    In 'full' mode, outputs a dictionary containing, among other things:
      'distmod':    The distance moduli that define the distance bins.
      'best':       The best-fit (maximum proability density)
                    line-of-sight reddening, in units of SFD-equivalent
                    E(B-V), to each distance modulus in 'distmod.' See
                    Schlafly & Finkbeiner (2011) for a definition of the
                    reddening vector (use R_V = 3.1).
      'samples':    Samples of the line-of-sight reddening, drawn from
                    the probability density on reddening profiles.
      'success':    1 if the query succeeded, and 0 otherwise.
      'converged':  1 if the line-of-sight reddening fit converged, and
                    0 otherwise.
      'n_stars':    # of stars used to fit the line-of-sight reddening.
      'DM_reliable_min':  Minimum reliable distance modulus in pixel.
      'DM_reliable_max':  Maximum reliable distance modulus in pixel.
    
    Less information is returned in 'lite' mode, while in 'sfd' mode,
    the Schlegel, Finkbeiner & Davis (1998) E(B-V) is returned.
    '''
    
    url = 'http://argonaut.skymaps.info/gal-lb-query-light'
    
    payload = {'mode': mode}
    
    if coordsys.lower() in ['gal', 'g']:
        payload['l'] = lon
        payload['b'] = lat
    elif coordsys.lower() in ['equ', 'e']:
        payload['ra'] = lon
        payload['dec'] = lat
    else:
        raise ValueError("coordsys '{0}' not understood.".format(coordsys))
    
    headers = {'content-type': 'application/json'}
    
    r = requests.post(url, data=json.dumps(payload), headers=headers)
    
    try:
        r.raise_for_status()
    except requests.exceptions.HTTPError as e:
        print('Response received from Argonaut:')
        print(r.text)
        raise e
    
    ebv = json.loads(r.text)
    return ebv['EBV_SFD']

def pan_catalog_cut(cat_raw_name, RA, DEC):
    "Apply several cuts and extinction correction to panstarrs catalog"
    from astropy.table import Table
    import itertools, numpy as np
        
    catalog_raw = Table.read(cat_raw_name, format='ascii.csv', guess=False)
    N = len(catalog_raw) 
    
    #colors = ['g','r','i','z','y']
    colors = ['r']

    psfMags = [c +'FPSFMag' for c in colors]
    KronMags = [c + 'FKronMag' for c in colors] 
    

    ## separate stars from galaxies
    ## https://confluence.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies
    flag = np.ones(N, dtype=bool)
    for psfMag, KronMag in zip(psfMags, KronMags):
        	flag *= (catalog_raw[psfMag] - catalog_raw[KronMag] < 0.05)
        	flag *= (catalog_raw[psfMag] - catalog_raw[KronMag] > -0.2)
    
    index = np.where(flag==False)[0]
    catalog_raw.remove_rows(index)    ## remove galaxies
    
    for KronMag in KronMags:
    	catalog_raw.remove_column(KronMag)  ## remove KronMags

    ## dust extinction correction: http://argonaut.skymaps.info/
    ## coefficients: Schlafly & Finkbeiner, 2011
    EBV = panstarrs_ebv(RA,DEC,mode='sfd')
    coeffs = {'g':3.172, 'r':2.271, 'i':1.682, 'z':1.322, 'y':1.087}
    for psfMag, color in zip(psfMags, colors):
        catalog_raw[psfMag] -= EBV * coeffs[color]
        print 'dust extinction for PanSTARRS band ' + color + ':', EBV*coeffs[color]

    ## use psfMags as Pogson Mags
    #for c in colors:
    #    catalog_raw['psfPogCorr_'+c] = catalog_raw['psfMagCorr_'+c]
    #    catalog_raw['psfPogErr_'+c] = catalog_raw['psfMagErr_'+c]

    
    catalog_raw.write("cat_pan.csv", format='ascii.csv') 
    return "cat_pan.csv"

def pan_query(query, RA, DEC):
    "Run panstarrs query via Casjobs"
    import os, glob

    bashCommand = []
    bashCommand.append('rm *.csv')
    bashCommand.append('java -jar ${CasJobs} run ' + "\'" + query + "\'")
    bashCommand.append('java -jar ${CasJobs} extract -b big_macs_db -type csv -d -F')
    bashCommand.append('java -jar ${CasJobs} execute -t mydb/1 \"drop table big_macs_db\" ')

    for c in bashCommand:
        os.system(c)

    cat_raw_name = glob.glob("./*.csv")[0]
    cat_pan_name = pan_catalog_cut(cat_raw_name, RA, DEC)

    return cat_pan_name


def write_header(ofp,pre,url,qry):
    import  time
    ofp.write('%s SOURCE: %s\n' % (pre,url))
    ofp.write('%s TIME: %s\n' % (pre,time.asctime()))    
    ofp.write('%s QUERY:\n' % pre)
    for l in qry.split('\n'):
        ofp.write('%s   %s\n' % (pre,l))
    
def main(argv):
    "Parse command line and do it..."
    import os, getopt, string
    
    queries = []
    url = os.getenv("SQLCLURL",default_url)
    fmt = default_fmt
    writefirst = 1
    verbose = 0
    
    # Parse command line
    try:
        optlist, args = getopt.getopt(argv[1:],'s:f:q:vlh?')
    except getopt.error, e:
        usage(1,e)
        
    for o,a in optlist:
        if   o=='-s': url = a
        elif o=='-f': fmt = a
        elif o=='-q': queries.append(a)
        elif o=='-l': writefirst = 0
        elif o=='-v': verbose += 1
        else: usage(0)
        
    if fmt not in formats:
        usage(1,'Wrong format!')

    # Enqueue queries in files
    for fname in args:
        try:
            queries.append(open(fname).read())
        except IOError, e:
            usage(1,e)

    # Run all queries sequentially
    for qry in queries:
        ofp = sys.stdout
        if verbose:
            write_header(ofp,'#',url,qry)
        file = query(qry,url,fmt)
        # Output line by line (in case it's big)
        line = file.readline()
        if line.startswith("ERROR"): # SQL Statement Error -> stderr
            ofp = sys.stderr
        if writefirst:
            ofp.write(string.rstrip(line)+os.linesep)
        line = file.readline()
        while line:
            ofp.write(string.rstrip(line)+os.linesep)
            line = file.readline()


if __name__=='__main__':
    import sys
    main(sys.argv)
