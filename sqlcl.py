
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
    print(__doc__)
    if msg:
        print('-- ERROR: %s' % msg)
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
    print(url+params)
    return urllib.urlopen(url+params)    

def gaia_query(file, query, EBV, DR):
    from astroquery.gaia import Gaia
    from astropy.table import Table
    if (DR == 2):
        Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"  # Reselect Data Release 2
    elif (DR == 3):
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
    import numpy as np
    ''' Que Gaia SQL server '''
    job = Gaia.launch_job_async(query,verbose=True)
    gaia_data = job.get_results() #Accessing data from query might be different for DR3 (different units??)
    print("obtained gaia data")

    colors = ['g','bp','rp']

    Av = 3.1 * EBV #holden# is this still correct for DR3? Looks like it, both extinction at 550nm
    bp_rp = gaia_data['bp_rp']
	# calculate the extinction (Gaia Data Release 2:Observational Hertzsprung-Russell diagrams)
    if (DR == 2):
        coeffs = {  'kg':[0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099], \
                    'kbp':[1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043], \
                    'krp':[0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006] }
        c_terms = [np.ones(bp_rp.shape), bp_rp, bp_rp**2, bp_rp**3, Av, Av**2, bp_rp*Av] 
        #zps_ab = { 'g':25.7934, 'bp':25.3806, 'rp':25.1161} #DR2 - Don't Use
        zps_ab = { 'g':25.7916, 'bp':25.3862, 'rp':25.1162} #DR2 Revised https://www.cosmos.esa.int/web/gaia/iow_20180316
    elif (DR == 3):
        #DR3 Coeffs https://www.cosmos.esa.int/web/gaia/edr3-extinction-law - Main Sequence - BPRP XName 
        coeffs = {  'kg':[0.99596972154, -0.15972646030, 0.01223807382, 0.00090726555, -0.03771602639, 0.00151347495, -0.00002523645, 0.01145226581, -0.00093691499, -0.00026029677], \
                    'kbp':[1.15363197483, -0.08140129917, -0.03601302398, 0.01921435856, -0.02239754824, 0.00084056268, -0.00001310180, 0.00660124080, -0.00088224750, -0.00011121576], \
                    'krp':[0.66320787941, -0.01798471649, 0.00049376945, -0.00267994406, -0.00651422147, 0.00003301799, 0.00000157894, -0.00007980090, 0.00025567981, 0.00001104766] }
        c_terms = [np.ones(bp_rp.shape), bp_rp, bp_rp**2, bp_rp**3, Av, Av**2, Av**3, Av*bp_rp, Av*(bp_rp**2), (Av**2)*bp_rp] #DR3 https://www.cosmos.esa.int/web/gaia/edr3-extinction-law
        zps_ab = { 'g':25.8010 , 'bp':25.3540, 'rp':25.1040} #DR3 https://www.cosmos.esa.int/web/gaia/edr3-passbands AB ZPs



    k_g, k_bp, k_rp = 0.0, 0.0, 0.0
    for i in range(len(c_terms)):
        k_g += coeffs['kg'][i] * c_terms[i]
        k_bp += coeffs['kbp'][i] * c_terms[i]
        k_rp += coeffs['krp'][i] * c_terms[i]
    
    a_g = Table.Column( name = 'a_g', data = k_g * Av)
    a_bp = Table.Column( name = 'a_bp', data = k_bp * Av)
    a_rp = Table.Column( name = 'a_rp', data = k_rp * Av)
    
    gaia_data.add_column(a_g)
    gaia_data.add_column(a_bp)
    gaia_data.add_column(a_rp)

	# calculate magnitude  and err
    for c in colors:
        ab_mag = Table.Column( name='ab_' + c, data = -2.5*np.log10( gaia_data['phot_' + c + '_mean_flux'] ) + zps_ab[c]  - gaia_data['a_' + c]) ## zp and ext
        ab_mag = Table.Column( name='ab_' + c, data = -2.5*np.log10( gaia_data['phot_' + c + '_mean_flux'] ) + zps_ab[c] ) ## zp but no ext
        mag_err = Table.Column( name = 'phot_'+ c + '_mean_mag_error', data = 2.5 * gaia_data['phot_'+ c +'_mean_flux_error'] / gaia_data['phot_' + c +'_mean_flux'] )
        
        gaia_data.add_column(ab_mag)
        gaia_data.add_column(mag_err)
        
        gaia_data.remove_column('phot_' + c +'_mean_flux')
        gaia_data.remove_column('phot_' + c +'_mean_flux_error')
        
    gaia_data.write(file + '.cut.csv', format='ascii.csv', overwrite=True)
    #gaia_data.write(file + '.csv', format='ascii.csv', overwrite=True)


def panstarrs_ebv(lon, lat, coordsys='equ', mode='full'):
    '''
    import json, requests
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
        print(r.text) #requests.exceptions.HTTPError: 500 Server Error: INTERNAL SERVER ERROR for url: http://argonaut.skymaps.info/api/v2/sfd/query
        raise e
    
    ebv = json.loads(r.text)
    return ebv['EBV_SFD']
    '''
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as units
    from dustmaps.bayestar import BayestarQuery
    bayestar = BayestarQuery(map_fname="/fs/ddn/sdf/group/kipac/u/awright/bayestar2019.h5")
    coords = SkyCoord(ra=lon*units.deg, dec=lat*units.deg,
                    frame='icrs')

    reddening = bayestar(coords, mode='median')
    print("REDDENING")
    print(reddening)
    print(lon)
    print(lat)
    #print(reddening[best])
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as units
    from dustmaps.sfd import SFDQuery
    from dustmaps.sfd import fetch
    fetch() #get the sfd map
    sfd = SFDQuery()
    if coordsys.lower() in ['gal', 'g']:
        coords = SkyCoord(l=lon*units.deg, b=lat*units.deg,
                    frame='galactic')
    elif coordsys.lower() in ['equ', 'e']:
        coords = SkyCoord(ra=lon*units.deg, dec=lat*units.deg,
                    frame='icrs')

    ebv = sfd(coords)
    return ebv

def pan_catalog_cut(file, cat_raw_name, RA, DEC):
    "Apply several cuts and extinction correction to panstarrs catalog"
    from astropy.table import Table
    import itertools, numpy as np
        
    print(cat_raw_name)
    catalog_raw = Table.read(cat_raw_name, format='ascii.csv', guess=False)
    catalog_raw.rename_column('raMean', 'ra')
    catalog_raw.rename_column('decMean', 'dec')

    N = len(catalog_raw) 
    
    #colors = ['g','r','i','z','y']
    colors = ['r']

    psfMags = [c +'PSFMag' for c in colors]
    KronMags = [c + 'KronMag' for c in colors] 
    
    ## delete none-detections
    flag = np.ones(N, dtype=bool)
    for psfMag in psfMags:
        flag *= (catalog_raw[psfMag] > 0)
        flag *= (catalog_raw[psfMag] < 30)

    index = np.where(flag==False)
    catalog_raw.remove_rows(index[0])

    ## separate stars from galaxies
    ## https://confluence.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies

    N = len(catalog_raw) 
    flag = np.ones(N, dtype=bool)
    for psfMag, KronMag in zip(psfMags, KronMags):
        flag *= (catalog_raw[psfMag] - catalog_raw[KronMag] < 0.05)
        flag *= (catalog_raw[psfMag] - catalog_raw[KronMag] > -0.2)
    
    index = np.where(flag==False)
    catalog_raw.remove_rows(index[0])    ## remove galaxies
    
    for KronMag in KronMags:
        catalog_raw.remove_column(KronMag)  ## remove KronMags
    
    ## dust extinction correction: http://argonaut.skymaps.info/
    ## coefficients: Schlafly & Finkbeiner, 2011
    EBV = panstarrs_ebv(RA,DEC,mode='sfd') 

    coeffs = {'g':3.172, 'r':2.271, 'i':1.682, 'z':1.322, 'y':1.087}
    for psfMag, color in zip(psfMags, colors):
        catalog_raw[psfMag] -= EBV * coeffs[color]
        print('dust extinction for PanSTARRS band ' + color + ':', EBV*coeffs[color])
    catalog_raw.write(file + ".csv", format='ascii.csv', overwrite=True)
    return file + ".csv"

def pan_query(file, cmd, RA, DEC):

    import os, glob
    if not os.path.exists(file +'.pan_raw.csv'):
        os.system(cmd)

    cat_raw_name = file + '.pan_raw.csv'
    cat_pan_name = pan_catalog_cut(file, cat_raw_name, RA, DEC)

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
    except getopt.error as e:
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
        except IOError as e:
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
