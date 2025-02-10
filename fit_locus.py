
#def import_lib():

if __name__ != '__main__':
    print('importing modules')
    import os, re, string
    import random, scipy #dbm #commands #anydbm
    from astropy.io import fits
    import numpy as np
    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib import pyplot as plt
    from scipy import linalg
    from scipy import optimize
    from glob import glob
    from copy import copy
    from functools import reduce
    import utilities
    
    print('finished importing modules')
    #import_lib()

global itr
itr = 0


def fix_kpno():
    
    p = fits.open('./EXAMPLES/kpno.fits')

    for color in ['g','r','i','z']: 
        mask = p[1].data['FLAGS_reg1_' + color] != 0  
        p[1].data['MAG_APERCORR_reg1_' + color][mask] = 99

        mask = p[1].data['IMAFLAGS_ISO_reg1_' + color] != 0  
        print(mask)
        p[1].data['MAG_APERCORR_reg1_' + color][mask] = 99

    p.writeto('./EXAMPLES/kpno_fixed.fits')
        

def join_cats(cs,outputfile):
    from astropy.io import fits
    tables = {}
    i = 0
    cols = []
    seqnr = 0 
    for c in cs:
        if len(c) == 2:
            TAB = c[1]
            c = c[0]
        else: TAB = 'STDTAB'
        i += 1
        print(c)
        tables[str(i)] = fits.open(c)
        for column in  tables[str(i)][TAB].columns:           
            if column.name == 'SeqNr':
                if not seqnr:
                    seqnr += 1
                else:                                            
                    column.name = column.name + '_' + str(seqnr)
                    seqnr += 1

            cols.append(column)

    print(cols)
    print(len(cols))
    hdu = fits.PrimaryHDU()
    #hduSTDTAB = fits.new_table(cols) 
    hduSTDTAB = fits.BinTableHDU.from_columns(cols) 
    hdulist = fits.HDUList([hdu])
    hdulist.append(hduSTDTAB)
    hdulist[1].header.update('EXTNAME','STDTAB')
    import os
    os.system('rm ' + outputfile)
    print(outputfile)
    hdulist.writeto(outputfile)

def get_survey_stars(file, inputcat, racol, deccol, necessary_columns, EBV, survey='SDSS', sdssUnit=False): 

    import scipy, math
    from astropy.io import fits
    from functools import reduce

    RA, DEC, RADIUS = get_catalog_parameters(inputcat, racol, deccol)
    #print("Radius: " + RADIUS)

    print('WILL SEARCH FOR STARS WITHIN ' + str(RADIUS) + 'min OF ' + str(RA) + ' ' + str(DEC))

    if survey == 'SDSS':
	    #colors = ['u','g','r','i','z']
	    #color_AB = [['u',-0.04],['g',0],['r',0],['i',0],['z',0.02]]
        colors = ['r']
        color_AB = [['r',0]]

        ''' includes conversion to Pogson magnitudes from luptitudes '''
        keys = ['ra','dec']
        keys += ['psfMag_%(color)s - extinction_%(color)s + %(AB).2f as psfMagCorr_%(color)s' % {'color':color,'AB':AB} for color, AB in color_AB ] 
        keys += ['psfMagErr_%(color)s' % {'color':color,'AB':AB} for color, AB in color_AB ] 
        keys += ['-2.5*LOG10(psfFlux_%(color)s) + 22.5 - extinction_%(color)s + %(AB).2f as psfPogCorr_%(color)s' % {'color':color,'AB':AB} for color, AB in color_AB ] 
        keys += ['1.085/SQRT(psfFluxIvar_%(color)s)/psfFlux_%(color)s as psfPogErr_%(color)s' % {'color':color,'AB':AB} for color, AB in color_AB ] 

        wherekeys = ['psfFlux_%(color)s > 0 ' % {'color':color,'AB':AB} for color, AB in color_AB ] 
        import sqlcl
        ''' includes AB correction and extinction correction , require g mag (luptitude) error less than 0.1  '''   
        #query = 'select ra, dec, s.psfMag_u - extinction_u - 0.04, s.psfMag_g - extinction_g, s.psfMag_r - extinction_r, s.psfMag_i - extinction_i, s.psfMag_z - extinction_z + 0.02, s.psfMagErr_u, s.psfMagErr_g, s.psfMagErr_r, s.psfMagErr_i, s.psfMagErr_z from star as s JOIN dbo.fGetNearbyObjEq(' + str(RA) + ',' + str(DEC) + ',' + str(RADIUS) + ' ) AS GN ON s.objID = GN.objID where s.clean=1 and s.psfMagErr_g < 0.1' # and s.psfMagErr_z < 0.1'
        query = 'select ' + reduce(lambda x,y: x + ',' + y, keys) + ' from star as s JOIN dbo.fGetNearbyObjEq(' + str(RA) + ',' + str(DEC) + ',' + str(RADIUS) + ' ) AS GN ON s.objID = GN.objID where s.clean=1 and ' + reduce(lambda x,y: x + ' and ' + y,wherekeys) 
        #query = 'select ra, dec FROM star'
        ''' cannot query SDSS database more than once per second '''
        print(query)
        
        lines = sqlcl.query(query).readlines()
        print(len(lines) - 2, 'STAR(S) FOUND')
        #print(lines[1])

        #returned_keys = re.split(r'\,',lines[1][:-1])
        returned_keys = re.split(r'\,', lines[1][:-1].decode('utf-8'))
        print(returned_keys)
        saveKeys = returned_keys[2:]

        ''' make a array with empty list with an entry for each key '''
        catalogStars = dict(zip(returned_keys,list([[] for x in returned_keys])))
        print(catalogStars.keys())

        if lines[1] == 'N' or len(lines) -2  < 5:
            print('NO USABLE SDSS DATA FOUND, PROCEEDING') 
            matched = False
            returnCat = inputcat
        else:
            matched = True
            for line in lines[2:]:
                line = line.decode('utf-8').replace('\n', '')
                res = re.split(',',line)
                for i in range(len(res)): 
                    catalogStars[returned_keys[i]].append(float(res[i]))  

    elif survey == 'PanSTARRS':
        #colors = ['g','r','i','z','y']   
        colors = ['r']   

        keys = ['raMean','decMean']
        keys += ['%(color)sPSFMag' % {'color':color} for color in colors ] 
        keys += ['%(color)sKronMag' % {'color':color} for color in colors ]
        keys += ['%(color)sPSFMagErr' % {'color':color} for color in colors ] 
        
        columns = '['
        for key in keys:
            columns += key + ','
        columns = columns[:-1] + ']'
        
        cmd = "curl -g \'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.csv"
        cmd += "?ra=" + str(RA) + "&dec=" + str(DEC) + "&radius=" + str(RADIUS/60) + "&columns=" + columns   ## RADIUS in arcmin -> degree
        cmd += "&nDetections>1"
        for c in colors:
            cmd += "&" + c + "KronMag>0"
            cmd += "&n" + c + ">0"
        cmd += "\' > " + file + ".pan_raw.csv"
        print(cmd)
        
        pan_bands = ''
        for c in colors:
            pan_bands += c
        print("Query PanSTARRS " + pan_bands + " for reference")
        print(cmd)
        
        import sqlcl

        ref_cat_name = sqlcl.pan_query(file, cmd, RA, DEC)
        
        with open(ref_cat_name) as ref_cat:
            lines = ref_cat.readlines()
            print(lines)
        print(len(lines) - 1, 'STAR(S) FOUND')
        
        print(lines[0]) #not actually getting all the columns
        print("Re=", lines[0][:-1])
        returned_keys = re.split(r'\,',lines[0][:-1])
        saveKeys = returned_keys[2:]

        ''' make a array with empty list with an entry for each key '''
        catalogStars = dict(zip(returned_keys,list([[] for x in returned_keys])))
        print(catalogStars.keys())

        if lines[0] == 'N' or len(lines) -1  < 5:
            print('NO USABLE PanSTARRS DATA FOUND, PROCEEDING')
            matched = False
            returnCat = inputcat
        else:
            matched = True
            for line in lines[1:]:
                line = line.replace('\n','')
                res = re.split(',',line)
                for i in range(len(res)): 
                    catalogStars[returned_keys[i]].append(float(res[i]))

    elif survey == 'Gaia':
        import sqlcl
        ''' Gaia ADQL, Radius in degrees. Color excess cut:
        https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_qa/ssec_cu5pho_excessflux.html '''
        DR = 3 #Release of Gaia to be Used, make sure to change res file too #cDR
        if DR == 3:
            color_range = "AND bp_rp >  -0.06 AND bp_rp < 2.5"
        else:
            #color_range = "AND bp_rp >  0.6 AND bp_rp < 1.6"
            color_range = ""
        
        RAD = RADIUS / 60
        if DR == 2:
            query = "SELECT ra, dec, bp_rp, \
                phot_g_mean_flux, phot_g_mean_flux_error,  \
                            phot_bp_mean_flux, phot_bp_mean_flux_error, \
                            phot_rp_mean_flux, phot_rp_mean_flux_error \
                            FROM gaiadr" + str(DR) + ".gaia_source \
                            WHERE 1=CONTAINS( POINT('ICRS',ra,dec), BOX('ICRS'," + str(RA) + "," + str(DEC) + "," + str(RAD) + ", " + str(RAD) + ")) \
                            AND phot_g_mean_mag<=19 AND phot_bp_mean_mag>=5 AND phot_rp_mean_mag>=5 \
                            AND phot_bp_rp_excess_factor > (1.0 + 0.015*bp_rp*bp_rp) AND phot_bp_rp_excess_factor < (1.3 + 0.06*bp_rp*bp_rp) " \
                            + str(color_range)
            #AND phot_g_mean_mag<=19
        elif DR==3: 
            #holden# Using C* to filter off of bp_rp_excess_factor https://www.aanda.org/articles/aa/full_html/2023/06/aa43680-22/aa43680-22.html#R27
            #Looks like that is the way to filter of off BP/RP flux excess, rather than equations similar to DR2("do not take in consideration the uncertainties on the flux excess factor")
            #https://gea.esac.esa.int/archive/documentation/GEDR3/Data_processing/chap_cu5pho/cu5pho_sec_photProc/cu5pho_ssec_photVal.html
    
            #c_star mean is 0 (the ABS stuff is c_star, the excess factor - parens)
            query = f"""
            SELECT dr3.ra, dr3.dec, dr3.bp_rp,
                dr3.phot_g_mean_flux, dr3.phot_g_mean_flux_error,
                dr3.phot_bp_mean_flux, dr3.phot_bp_mean_flux_error,
                dr3.phot_rp_mean_flux, dr3.phot_rp_mean_flux_error
            FROM gaiadr3.gaia_source AS dr3
            WHERE 1 = CONTAINS(
                    POINT('ICRS', ra, dec),
                    BOX('ICRS', {RA}, {DEC}, {RAD}, {RAD})
                )
                AND phot_g_mean_mag <= 21
                AND phot_bp_mean_mag >= 5
                AND phot_rp_mean_mag >= 5
                AND (
                   ( bp_rp >= 4.0 AND ABS(dr3.phot_bp_rp_excess_factor-(1.057572+0.1405537*bp_rp))<(0.0059898 + 8.817481e-12 * POWER(dr3.phot_g_mean_mag,7.618399)) ) OR
                   ( bp_rp < 0.5 AND ABS(dr3.phot_bp_rp_excess_factor-(1.154360+0.033772*bp_rp+0.32277*bp_rp*bp_rp))<(0.0059898 + 8.817481e-12 * POWER(dr3.phot_g_mean_mag,7.618399)) ) OR
                   ( bp_rp >= 0.5 AND bp_rp < 4.0 AND ABS(dr3.phot_bp_rp_excess_factor-(1.162004+.011464*bp_rp+0.049255*bp_rp*bp_rp-0.005879*bp_rp*bp_rp*bp_rp))<(0.0059898 + 8.817481e-12 * POWER(dr3.phot_g_mean_mag,7.618399))  )
                )
                {color_range}
            """
            #see https://scholar.google.com/scholar_lookup?title=Gaia+Early+Data+Release+3+-+Photometric+content+and+validation&author=Riello+M.+De+Angeli+F.+Evans+D.+W.&journal=A%26A&volume=649&pages=A3&publication_year=2021&issn=0004-6361%2C1432-0746&doi=10.1051%2F0004-6361%2F202039587
            #sction 9.4 for information about filtering with c_star, stddev is estimated using a power-law



# Uncomment the following line if additional filters are needed
# query += "AND c_star < 20000 AND bp_rp > 0.6 AND bp_rp < 1.6"
        print(query)
        
        EBV, gallong, gallat = galactic_extinction_and_coordinates(RA,DEC)
        print("FILE=",file)
        sqlcl.gaia_query(file, query, EBV, DR=DR)	

        with open(file + '.cut.csv') as ref_cat:
            lines = ref_cat.readlines()
        print(len(lines) - 1, 'STAR(S) FOUND')
        
        returned_keys = re.split(r'\,',lines[0][:-1])
        saveKeys = returned_keys[2:]

        ''' make a array with empty list with an entry for each key '''
        catalogStars = dict(zip(returned_keys,list([[] for x in returned_keys])))
        print(catalogStars.keys())

        if lines[0] == 'N' or len(lines) -1  < 5:
            print('NO USABLE PanSTARRS DATA FOUND, PROCEEDING')
            matched = False
            returnCat = inputcat
        else:
            matched = True
            for line in lines[1:]:
                line = line.replace('\n','')
                res = re.split(',',line)
                for i in range(len(res)): 
                    catalogStars[returned_keys[i]].append(float(res[i]))
               
    elif survey == '2MASS':
        if True or RADIUS > 59:
            RADIUS = 59
            coordinate = str(RA) + '+' + str(DEC)
            catalog = '2MASS_stars.cat'

            ''' NOTE 2MASS MAGS NOT CORRECTED FOR DUST -- SHOULD BE CORRECTED '''

            ''' select 2MASS stars with ph_qual=A for J band (includes a S/N cut) and use_src=1 '''
            command = "wget \"http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?outfmt=1&objstr=" + coordinate + "&spatial=Cone&radius=" + str(RADIUS) + "&radunits=arcmin&catalog=fp_psc&selcols=ph_qual,ra,dec,j_m,j_cmsig&constraints=ph_qual+like+%27A__%27+and+use_src%3D1\" -O "  + catalog
            print(command)

            import os
            os.system(command)

            lines = open(catalog,'r').readlines()
        
            keyDict = {} 
            saveKeys = ['ra','dec','j_m','j_cmsig']
            for line in lines:
                if line[0] == '|' and keyDict == {}:
                    returned_keys_full = re.split(r'\|',line)[1:]
                    returned_keys = [r.replace(' ','') for r in returned_keys_full]
                    index = 1
                    for key_full in returned_keys_full:
                        indexStart = index 
                        index += len(key_full) + 1
                        indexEnd = index  
                    
                        keyDict[key_full.replace(' ','')] = {'indexStart': indexStart, 'indexEnd': indexEnd}

                    catalogStars = dict(zip(returned_keys,list([[] for x in returned_keys])))

                elif line[0] != '#' and line[0] != '|' and line[0] != '\\':
                    for key in saveKeys:
                        value = float(line[keyDict[key]['indexStart']:keyDict[key]['indexEnd']])
                        if key == 'j_m':
                            ''' correct J magnitude for extinction '''
                            value = value - EBV * 0.709
                        catalogStars[key].append(value)
            
            if catalogStars.values()[0]:
                print('NO USABLE 2MASS DATA FOUND, PROCEEDING')
                matched = False 
                returnCat = inputcat
            else: 
                matched = True

    if matched:
        print('making KDTrees')                                                                                                                                                                                                                                
        if survey == 'SDSS' and sdssUnit:
            ''' make a catalog of all SDSS stars (i.e., not just those matched against catalog stars) ''' #LOOK# is SDSS handled differently?                                                       
            cols = []
            for column_name in returned_keys[2:]: 
                cols.append(fits.Column(name=column_name,format='1E',array=np.array(catalogStars[column_name])))

            coldefs = fits.ColDefs(cols)
	    #hdu_new = fits.new_table(coldefs)
            hdu_new = fits.BinTableHDU.from_columns(coldefs)

            returnCat = hdu_new

            matched = True

        else:
            from scipy import spatial 
            
            data_catalog = list(zip(catalogStars['ra'],catalogStars['dec']))

            data_inputcat = list(zip(inputcat.data.field(racol),inputcat.data.field(deccol)))

            kdtree_catalog = spatial.KDTree(data_catalog)
            kdtree_inputcat = spatial.KDTree(data_inputcat)
            match = kdtree_catalog.query_ball_tree(kdtree_inputcat,2./3600.) #LOOK# finding matches with nearest neighbors within 2-arcsec
            print(match)

            ''' make catalog with same number of row as inputcat and columns for catalog mags  '''
            rows = len(inputcat.data)
            #for column_name in returned_keys[2:]: 
            #    print(column_name)
            
            cols = []
            for column_name in necessary_columns: #inputcat.columns:
                #cols.append(column)
                cols.append(fits.Column(name=column_name,format='1E',array=inputcat.data.field(column_name)))


            necessary_columns += saveKeys

            for column_name in saveKeys: #Adding new hold/ext columns
                #print(column_name)
                array = np.ones(rows) * -99
                cols.append(fits.Column(name=column_name,format='1E',array=array))
    

            coldefs = fits.ColDefs(cols)
            hdu_new = fits.TableHDU.from_columns(coldefs)

            matchedStars = 0

            for i in range(len(match)): #LOOK# probably need to change something here for double matches
                if len(match[i]) == 1:
                    matchedStars += 1
                    for column_name in saveKeys: 
                        #print(column_name)
                        hdu_new.data.field(column_name)[match[i][0]] = catalogStars[column_name][i] #adding rPSF mags to matched stars
                        '''
                        if column_name == 'rPSFMag':
                            print(match[i][0]) # In current form, you can get into the scenario where match[i][0] can be the same for multiple i
                            # this means that an input cat star is the NN for multiple ref cat stars
                            # could be fixed by not incrementing matched stars if hdu_new.data.field(column_name)[match[i][0]] != -99 or something
                            print(catalogStars[column_name][i])
                        #catalog stars is ref cat
                        '''
            #print(len(match))
            #print(len(saveKeys))
            #print(matchedStars)
            #print("rPSF Mag ", hdu_new.data.field('rPSFMag')[24])
            #print("W-J-V Mag ", hdu_new.data.field('W-J-V')[24])

            ''' require at least five matched stars '''
            if matchedStars > 3:
                if False: #turn to true if you want to remove nonmatches
                    temp_table = Table(hdu_new.data)

                    # Filter out rows where rPSFMag == 99
                    filtered_table = temp_table[temp_table['rPSFMag'] != -99]

                    # Convert back to a NumPy structured array
                    filtered_data = filtered_table.as_array()

                    # Create a new BinTableHDU with the filtered data
                    hdu_new = fits.BinTableHDU(data=filtered_data)

                matched = matchedStars 
                hdu = fits.PrimaryHDU()               
                hdulist = fits.HDUList([hdu,hdu_new])
                print(len(match))
                import os
                os.system('rm merge.fits')
                hdulist.writeto('merge.fits')
                returnCat = hdu_new #hdu_new has data from all stars (both ref cat and normal cat, INCLUDING UNMATCHED)

            else: 
                print(str(matchedStars) + ' MATCHES WITH ' + survey  + ' CATALOG')
                returnCat = inputcat
                matched = 0 

    
    print(returnCat)
    
    ''' #print all, INCLUDING UNMATCHED, input mags in hdu_new
    for column_name in hdu_new.columns.names:  # Iterating over all column names
        column_data = hdu_new.data.field(column_name)  # Access the data for the current column
        print(f"Column: {column_name}")
        print(column_data)
        print("-" * 50)  # Separator for readability
    '''
    
    return returnCat, matched, necessary_columns


''' retrieve SFD dust extinction and galactic coordinates from NED '''
def galactic_extinction_and_coordinates(RA,DEC): 
    
        print('RETRIEVING DUST EXTINCTION AT RA=' + str(RA) + ' DEC=' + str(DEC) + ' FROM NED')
        import urllib, os, re, string, time #anydbm dbm 
        from urllib.request import urlopen
        from functools import reduce

        form = [
            "in_csys=Equatorial",
            "in_equinox=J2000.0",
            "obs_epoch=2000.0",    
            "lon=%(ra).7fd" % {'ra':float(RA)},
            "lat=%(dec).7fd" % {'dec':float(DEC)},
            "pa=0.0",
            "out_csys=Galactic",
            "out_equinox=J2000.0",
        ]

        response = urlopen('http://nedwww.ipac.caltech.edu/cgi-bin/nph-calc?' + reduce(lambda x,y: str(x) + '&' + str(y),form) + '"')  
        text = response.read()

        ''' scan for Galactic coordinates '''
        found = False 

        text = text.decode('utf-8')
        for l in text.split('\n'): 
            if found:
                res = re.split(r'\s+',l)
                gallong = float(res[0])
                gallat = float(res[1])
                break                    
            if l.find('Galactic') != -1 and l.find('Output:') != -1:
                found = True
            
        ''' find extinction in each band '''
        dict = {}
        for q in ['U','B','V','R','J','H','K']:
            for m in text.split('\n'):
                if m[0:11] == 'Landolt ' + q + ' (' :
                    line = re.split(r'\s+', m)       
                    dict[q] = line[3]         

        ebv = float(dict['B']) - float(dict['V'])

        print('EBV=', ebv)
        print('GAL LONG', gallong)
        print('GAL LAT', gallat)

        return ebv, gallong, gallat


''' sort each set by center wavelength '''
def sort_wavelength(x,y):
    if x['center wavelength']>y['center wavelength']:
        return 1
    else: return -1

def assign_zp(filt,pars,zps,zps_hold):
    if filt in zps:
        out = pars[zps[filt]]
    else: 
        out = zps_hold[filt] 
    return out

def get_kit(): 
    f = open(os.environ['kpno'] + '/process_kpno/locuskit','r')
    m = pickle.Unpickler(f)
    locus = m.load()
    return locus


def get_catalog_parameters(fulltable, racol, deccol):
    ''' calculate field center '''
    import numpy as np

    import scipy
    #DEC = np.median(fulltable.data.field(deccol))
    DEC = (fulltable.data.field(deccol).min() + fulltable.data.field(deccol).max())/2.
    DEC_DIFF_SQ = ((fulltable.data.field(deccol) - DEC) * 60.)**2.

    #RA = np.median(fulltable.data.field(racol))
    RA = (fulltable.data.field(racol).min() + fulltable.data.field(racol).max())/2.
    RA_DIFF_SQ = ((fulltable.data.field(racol) - RA) * 60. * np.cos(DEC))**2.

    RADII = (DEC_DIFF_SQ + RA_DIFF_SQ)**0.5

    return RA, DEC, RADII.max() 


def run(file,columns_description,output_directory=None,plots_directory=None,extension='OBJECTS',racol=None,deccol=None,end_of_locus_reject=1,plot_iteration_increment=50, min_err=0.02, bootstrap_num=0, snpath=None, night=None, run=None, prefix='',data_from_sdss=False, addSDSS=False, addPanSTARRS=False, addGaia=False, number_of_plots=10, add2MASS=False, sdssUnit=False, twoStep=False):
    from importlib import reload
    print(add2MASS, addGaia, addPanSTARRS, addSDSS)
    try: 
        extension = int(extension)
    except: pass

    print('trying to open file', file)
    fulltable = fits.open(file)[extension]

    input_info = utilities.parse_columns(columns_description)
    necessary_columns = [racol, deccol] + [x['mag'] for x in input_info] + [x['mag_err'] for x in input_info]

    print(necessary_columns)

    RA, DEC, RADIUS = get_catalog_parameters(fulltable, racol, deccol) 

    if RA is not None and DEC is not None:
        EBV, gallong, gallat = galactic_extinction_and_coordinates(RA,DEC)

    print('EBV=',EBV)
    #add in projection
    #inputcat.data.field(racol) - RA)**2. + (inputcat.data.field(deccol) - DEC)**2.)**0.5
    
    if twoStep:
        pass #do nothing, get relative ZPs later

    if output_directory is None:
        fs = file.split('/')
        if len(fs) > 1:
            output_directory = '/'.join(fs[:-1])
        else:
            output_directory = './'

        offsets_file = output_directory + '/' + file.split('/')[-1]  + '.offsets.list'


    fitSDSS = False
    foundSDSS = 0 
    if addSDSS:
        fulltable, foundSDSS, necessary_columns = get_survey_stars(file, fulltable, racol, deccol, necessary_columns, EBV, survey='SDSS', sdssUnit=sdssUnit)
        if foundSDSS: fitSDSS = True
            
    foundPanSTARRS = 0 
    if addPanSTARRS:
        fulltable, foundPanSTARRS, necessary_columns = get_survey_stars(file, fulltable, racol, deccol, necessary_columns, EBV, survey='PanSTARRS')

    foundGaia = 0
    if addGaia:
        fulltable, foundGaia, necessary_columns = get_survey_stars(file, fulltable, racol, deccol, necessary_columns, EBV, survey='Gaia')
            
    found2MASS = 0 
    if add2MASS:
        fulltable, found2MASS, necessary_columns = get_survey_stars(file, fulltable, racol, deccol, necessary_columns, EBV, survey='2MASS')
        if found2MASS: fit2MASS = True

    if output_directory is None:
        fs = file.split('/')
        if len(fs) > 1:
            output_directory = '/'.join(fs[:-1])
        else:
            output_directory = './'

    if plots_directory is None: 
        plots_directory = output_directory + '/PLOTS/'

    reload(utilities)

    ''' if SDSS stars, hold no other filter ZPs fixed (avoid tension) '''
    if addSDSS and foundSDSS:
        #input_info = utilities.parse_columns(columns_description,fitSDSS=False,noHoldExceptSDSS=True)
        for i in range(len(input_info)):
            input_info[i]['HOLD_VARY'] = 'VARY'

        if sdssUnit: 
            ''' include only SDSS magnitudes in unit test '''
            input_info = [] 

	    #sdss_info = [{'mag':'psfPogCorr_' + c, 'plotName':'SDSS ' + c, 'filter': 'SDSS-' + c + '.res', 'mag_err': 'psfPogErr_' + c, 'HOLD_VARY':'HOLD', 'ZP':0.} for c in ['g','r','i','z'] ]
        sdss_info = [{'mag':'psfPogCorr_' + c, 'plotName':'SDSS ' + c, 'filter': 'SDSS-' + c + '.res', 'mag_err': 'psfPogErr_' + c, 'HOLD_VARY':'HOLD', 'ZP':0.} for c in ['r'] ]

        for filt_dict in sdss_info:
            ''' avoid duplicate filters -- will override '''
            if filt_dict['mag'] not in [f['mag'] for f in input_info]:
                input_info += [filt_dict]

        ''' if SDSS unit test, hold only z-band zeropoint constant '''            
        if sdssUnit:
            for i in range(len(input_info)):                
                if input_info[i]['mag'] != 'psfPogCorr_z': 
                    input_info[i]['HOLD_VARY'] = 'VARY'
                    
    if addPanSTARRS and foundPanSTARRS:
        for i in range(len(input_info)):
            input_info[i]['HOLD_VARY'] = 'VARY'

        #panstarrs_info = [{'mag':c + 'PSFMag', 'plotName':'PanSTARRS ' + c, 'filter': 'PAN-STARRS.PS1.' + c + '.res', 'mag_err': c + 'PSFMagErr', 'HOLD_VARY':'HOLD', 'ZP':0.} for c in ['r','i','z','y'] ]
        panstarrs_info = [{'mag':c + 'PSFMag', 'plotName':'PanSTARRS ' + c, 'filter': 'PAN-STARRS.PS1.' + c + '.res', 'mag_err': c + 'PSFMagErr', 'HOLD_VARY':'HOLD', 'ZP':0.} for c in ['r'] ]
	    #panstarrs_info += [{'mag':c + 'PSFMag', 'plotName':'PanSTARRS ' + c, 'filter': 'PAN-STARRS.PS1.' + c + '.res', 'mag_err': c + 'PSFMagErr', 'HOLD_VARY':'VARY', 'ZP':0.} for c in ['z'] ]

        for filt_dict in panstarrs_info:
            ''' avoid duplicate filters -- will override '''
            if filt_dict['mag'] not in [f['mag'] for f in input_info]:
                input_info += [filt_dict]

    if addGaia and foundGaia:
        for i in range(len(input_info)):
            input_info[i]['HOLD_VARY'] = 'VARY'

        #gaia_info = [{'mag':'ab_g', 'plotName':'Gaia G' , 'filter': 'Gaia_dr2_revised.g.res', 'mag_err': 'phot_g_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.}]
        #gaia_info = [ {'mag':'ab_bp', 'plotName':'Gaia Gbp' , 'filter': 'Gaia_dr2_revised.bp.res', 'mag_err': 'phot_bp_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.}]
        #gaia_info = [{'mag':'ab_rp', 'plotName':'Gaia Grp' , 'filter': 'Gaia_dr2_revised.rp.res', 'mag_err': 'phot_rp_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.}]
        
        #gaia_info = [{'mag':'ab_g', 'plotName':'Gaia G' , 'filter': 'Gaia_dr2_revised.g.res', 'mag_err': 'phot_g_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.},\
        #        {'mag':'ab_bp', 'plotName':'Gaia Gbp' , 'filter': 'Gaia_dr2_revised.bp.res', 'mag_err': 'phot_bp_mean_mag_error', 'HOLD_VARY':'VARY', 'ZP':0.},\
        #        {'mag':'ab_rp', 'plotName':'Gaia Grp' , 'filter': 'Gaia_dr2_revised.rp.res', 'mag_err': 'phot_rp_mean_mag_error', 'HOLD_VARY':'VARY', 'ZP':0.} ]

        '''
        gaia_info = [{'mag':'ab_g', 'plotName':'Gaia G' , 'filter': 'Gaia_dr3.g.res', 'mag_err': 'phot_g_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.},\
                {'mag':'ab_bp', 'plotName':'Gaia Gbp' , 'filter': 'Gaia_dr3.bp.res', 'mag_err': 'phot_bp_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.},\
                {'mag':'ab_rp', 'plotName':'Gaia Grp' , 'filter': 'Gaia_dr3.rp.res', 'mag_err': 'phot_rp_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.} ]
        '''
        
        DR = 3 #cDR
        if DR == 3:
            gaia_info = [{'mag':'ab_g', 'plotName':'Gaia G' , 'filter': 'Gaia_dr3.g.res', 'mag_err': 'phot_g_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.}]
        elif DR == 2:
            gaia_info = [{'mag':'ab_g', 'plotName':'Gaia G' , 'filter': 'Gaia_dr2_revised.g.res', 'mag_err': 'phot_g_mean_mag_error', 'HOLD_VARY':'HOLD', 'ZP':0.}]

        for filt_dict in gaia_info:
            ''' avoid duplicate filters -- will override '''
            if filt_dict['mag'] not in [f['mag'] for f in input_info]:
                    input_info += [filt_dict]
                    
    if add2MASS and found2MASS:
        ''' if no SDSS, see if there are 2MASS matches '''
        #input_info = utilities.parse_columns(columns_description,fitSDSS=False,noHoldExcept2MASS=True)
        for i in range(len(input_info)):
            if input_info[i]['filter'].find('SDSS') == -1:
                input_info[i]['HOLD_VARY'] = 'VARY'

        sdss_info = [{'mag':'j_m', 'plotName':'2MASS J', 'filter': 'J2MASS.res', 'mag_err': 'j_cmsig', 'HOLD_VARY':'HOLD', 'ZP':0.} ]
        for filt_dict in sdss_info:
            ''' avoid duplicate filters -- will override '''
            if filt_dict['mag'] not in [f['mag'] for f in input_info]:
                input_info += [filt_dict]


    
    ''' check to see if at least one but not all filter is held constant '''
    if not list(filter(lambda x: x['HOLD_VARY'] == 'HOLD', input_info)): 
        raise Exception('None of your magnitudes is held fixed (i.e., HOLD_VARY HOLD)')
    #if not list(filter(lambda x: x['HOLD_VARY'] == 'VARY', input_info)): 
	    #raise Exception('All of your magnitudes are held fixed (i.e., HOLD_VARY VARY)')

    filters = utilities.get_filters([[a['mag'], a['filter']] for a in input_info])

    ''' update input_info with filter functions '''
    for i in range(len(filters)):
        input_info[i].update(filters[i])


    ''' separate into mag ZPs to be held fixed and varied '''
    info_hold = list(filter(lambda x: x['HOLD_VARY'] == 'HOLD',input_info))     
    info_vary = list(filter(lambda x: x['HOLD_VARY'] == 'VARY',input_info))        

    #info_hold.sort(sort_wavelength) 
    #info_vary.sort(sort_wavelength)

    info_hold.sort(key=lambda x: x['center wavelength'])
    info_vary.sort(key=lambda x: x['center wavelength'])


    ''' recombine '''
    input_info = info_hold + info_vary


    if RA is not None and DEC is not None:
        #EBV, gallong, gallat = galactic_extinction_and_coordinates(RA,DEC)
                                                                           
        for i in range(len(input_info)):
            print(input_info[i]['mag'])
            coeff = utilities.compute_ext(input_info[i])
            extinction = coeff * EBV 
            input_info[i]['extinction'] = extinction
            input_info[i]['gallong'] = gallong 
            input_info[i]['gallat'] = gallat 
            print(input_info[i]['mag'], extinction, ' (mag) in field', coeff)
    

    print(input_info)
    mag_locus = utilities.synthesize_expected_locus_for_observations(input_info)

    print(mag_locus)

    #fulltable = fits.open(file)[extension]

    #mask = ((fulltable.data.field('Xpos-SUBARU-W-J-V')- 5000)**2. +  (fulltable.data.field('Ypos-SUBARU-W-J-V') - 5000.)**2.)**0.5 < 2000
    #fulltable.data = fulltable.data[mask]
    #mask = fulltable.data.field('MaxVal-SUBARU-W-S-Z+') < 10 
    #fulltable.data = fulltable.data[mask]

    #fulltable.data = fulltable.data[:100]

    ''' if not SeqNr column, add one '''
    #if not list(filter(lambda x: x.name=='SeqNr', fulltable.columns)): 
    if not any(x.name == 'SeqNr' for x in fulltable.columns):
        cols = []
        for col in fulltable.columns:
            cols.append(col)
        cols.append(fits.Column(name='SeqNr',format='J',array=np.arange(len(fulltable.data))))
        hdu = fits.PrimaryHDU()
        hdulist = fits.HDUList([hdu])
        fulltable = fits.BinTableHDU.from_columns(cols)

    table = fulltable.data 

    print('INPUT CATALOG', file, 'EXTENSION', extension)

    red_input_info = []
    blue_input_info = []
    for mag in input_info: 
        if mag['center wavelength'] > 4000:
            mag['blue/red'] = 'REDDER'
            red_input_info.append(mag)
        else: 
            #holden# is there any band that is appropriate for bigmacs that would fit this? 
            #u-band does, but it shouldn't be used (also isn't working, says RC is BLUER for some reason), need to test two-step for this
            mag['blue/red'] = 'BLUER/RESTRICTED' 
            blue_input_info.append(mag)

    ''' designate which filter zeropoints to be held constant when matching bands '''
    zps_dict_all = {} 
    zps_dict_all_err = {} 
    cal_type = {}

    def update_zps(zps_dict_all,zps_dict_all_err, cal_type, results, red_or_blue):
        #if not combo['hold'] in zps_dict_all:
        #    zps_dict_all[combo['hold']] = 0.
        for key in results['full'].keys(): #list(filter(lambda x: x['HOLD_VARY']=='VARY', input_info)): 
            if results['hold_vary'][key] == 'VARY':
                zps_dict_all[key] = results['full'][key]
                zps_dict_all_err[key] = results['errors'][key]
                cal_type[key] = red_or_blue
        return zps_dict_all, zps_dict_all_err, cal_type

    ''' clear out plotting directory '''
    import os        
    os.system('rm ' + plots_directory + '/qc_*png')                                                    

    #program already ran to get relative ZPs, now get absolute ZPs
    def parse_file(file_path): #function to extract relative ZPs from output file
        red_band_order = [] #bands have to be in right order for fit function to work
        blue_band_order = []
        #band_order = ['W-J-B', 'W-J-V', 'W-C-RC', 'W-C-IC', 'W-S-Z+']
        red_vary_input_info = list(filter(lambda x: x['HOLD_VARY'] == 'VARY', red_input_info))
        blue_vary_input_info = list(filter(lambda x: x['HOLD_VARY'] == 'VARY', blue_input_info))
        for i in range(len(red_vary_input_info)):
            band = red_vary_input_info[i]['mag']
            red_band_order.append(band)

        for i in range(len(blue_vary_input_info)):
            band = blue_vary_input_info[i]['mag']
            blue_band_order.append(band)  
        
        print(red_band_order)
        print(blue_band_order)

        red_bands = [None] * len(red_band_order)
        red_zps = [None] * len(red_band_order)
        red_errors = [None] * len(red_band_order)

        blue_bands = [None] * len(blue_band_order)
        blue_zps = [None] * len(blue_band_order)
        blue_errors = [None] * len(blue_band_order)
            
        with open(file_path, 'r') as results_file:
            for line in results_file:
                # Skip comment lines and metadata
                if line.strip().startswith('#') or "RESULTS" in line:
                    continue

                # Split the line into components
                parts = line.split()
                if len(parts) >= 4:  # Ensure it's a valid data line
                    if ("REDDER" in line):
                        band = parts[0]  # First column is the band
                        zp = float(parts[1])  # Second column is the ZP
                        error = float(parts[3])  # Fourth column is the error on ZP
                        
                        red_bands[red_band_order.index(band)] = band
                        red_zps[red_band_order.index(band)] = zp
                        red_errors[red_band_order.index(band)] = error
                    elif ("BLUER" in line):
                        band = parts[0]  # First column is the band
                        zp = float(parts[1])  # Second column is the ZP
                        error = float(parts[3])  # Fourth column is the error on ZP
                        
                        blue_bands[blue_band_order.index(band)] = band
                        blue_zps[blue_band_order.index(band)] = zp
                        blue_errors[blue_band_order.index(band)] = error

        red_result_dict = {band: (zp, error) for band, zp, error in zip(red_bands, red_zps, red_errors)}
        blue_result_dict = {band: (zp, error) for band, zp, error in zip(blue_bands, blue_zps, blue_errors)}
        return red_result_dict, blue_result_dict

    ''' first calibrate redder filters '''
    # Where red ZPs get calculated, similar for blue_input
    if (twoStep):
        red_relative_zps_info, blue_relative_zps_info = parse_file(offsets_file)
        
        red_relative_zps = [info[0] for info in red_relative_zps_info.values()]
        red_relative_zps = [np.float64(val) for val in red_relative_zps]
        blue_relative_zps = [info[0] for info in blue_relative_zps_info.values()]
        blue_relative_zps = [np.float64(val) for val in blue_relative_zps]
        
        results, ref_mags, SeqNr = fit(table, red_input_info, mag_locus, min_err=min_err, end_of_locus_reject=end_of_locus_reject, plot_iteration_increment=plot_iteration_increment, bootstrap=True, bootstrap_num=bootstrap_num, plotdir=plots_directory, pre_zps=None, number_of_plots=number_of_plots, relative_zps=red_relative_zps)
    else:
        results, ref_mags, SeqNr = fit(table, red_input_info, mag_locus, min_err=min_err, end_of_locus_reject=end_of_locus_reject, plot_iteration_increment=plot_iteration_increment, bootstrap=True, bootstrap_num=bootstrap_num, plotdir=plots_directory, pre_zps=None, number_of_plots=number_of_plots)
    
    zps_dict_all, zps_dict_all_err, cal_type = update_zps(zps_dict_all,zps_dict_all_err,cal_type,results,'REDDER')

    print(len(ref_mags), len(SeqNr))

    ''' calibrate using bright u-band stars: removed'''

    if len(blue_input_info):

        gmr = ref_mags[:,1] - ref_mags[:,2]
        print(gmr)
        table = table[SeqNr]
        #table.field(blue_input_info[0]['mag'])[gmr > 0.5] = 99.

        print(table.field(blue_input_info[0]['mag']))
        print(SeqNr)


        print(len(table))

        for i in range(len(red_input_info)):
            red_input_info[i]['HOLD_VARY'] = 'HOLD'
            if red_input_info[i]['mag'] in zps_dict_all:
                red_input_info[i]['ZP'] = zps_dict_all[red_input_info[i]['mag']]
                red_input_info[i]['ZPERR'] = zps_dict_all_err[red_input_info[i]['mag']]
            else:
                red_input_info[i]['ZP'] = 0. 
                red_input_info[i]['ZPERR'] = 0. 

            print(zps_dict_all)
            print(red_input_info[i]['mag'], red_input_info[i]['ZP'])#, zps_dict_all[red_input_info[i]['mag']]


        print(red_input_info)
        
        if twoStep:
            results, ref_mags, SeqNr = fit(table, red_input_info + blue_input_info, mag_locus, min_err=min_err, end_of_locus_reject=end_of_locus_reject, plot_iteration_increment=plot_iteration_increment, bootstrap=True, bootstrap_num=bootstrap_num, plotdir=plots_directory, pre_zps=None, number_of_plots=number_of_plots, relative_zps=red_relative_zps + blue_relative_zps)
        else:
            results, ref_mags, SeqNr = fit(table, red_input_info + blue_input_info, mag_locus, min_err=min_err, end_of_locus_reject=end_of_locus_reject, plot_iteration_increment=plot_iteration_increment, bootstrap=True, bootstrap_num=bootstrap_num, plotdir=plots_directory, pre_zps=None, number_of_plots=number_of_plots)

        zps_dict_all, zps_dict_all_err, cal_type = update_zps(zps_dict_all,zps_dict_all_err,cal_type, results,'BLUER')

        print(results)
        print(zps_dict_all)


    output_string = '' 

    if twoStep and False:
        #red_relative_zps_info[band][0] ZP
        #red_relative_zps_info[band][1] ZP error
        abs_zp_offset = 0
        for band in red_relative_zps_info.keys():   
            if red_relative_zps_info[band][0] == 0: #HOLD band
                abs_zp_offset =  zps_dict_all[band]
                zps_dict_all_err[band] = 0

        for band in red_relative_zps_info.keys():   
            if red_relative_zps_info[band][0] != 0: #VARY bands
                zps_dict_all[band] = abs_zp_offset + red_relative_zps_info[band][0]
                zps_dict_all_err[band] = red_relative_zps_info[band][1]  

    if foundSDSS: 
        output_string += '#  USED ' + str(foundSDSS) + ' MATCHED SDSS STARS \n'
    if found2MASS:
        output_string += '#  USED ' + str(found2MASS) + ' MATCHED 2MASS STARS \n'
    if foundPanSTARRS:
        output_string += '#  USED ' + str(foundPanSTARRS) + ' MATCHED PanSTARRS STARS \n'
    if foundGaia:
        output_string += '#  USED ' + str(foundGaia) + ' MATCHED Gaia STARS \n'

    output_string += '# ' + str(gallat) + ' ' + str(gallong) + ' galactic latitude longitude \n'
    output_string += '# ' + str(results['redchi']) + ' reduced chi squared value \n'
    output_string += '# ' + str(results['num']) + ' number of stars \n'

    output_string += '# RESULTS: (ADD THESE ZP ADJUSTMENTS TO CATALOG MAGNITUDES) \n'
    for key in zps_dict_all.keys():    
        #print(key + ' ' + str(zps_dict_all[key]) + ' ' + str(zps_dict_all_err[key]) + ' ' + cal_type[key])
        output_string += key + ' ' + str(zps_dict_all[key]) + ' +- ' + str(zps_dict_all_err[key]) + ' ' + cal_type[key] + '\n'              

    ''' write out the magnitude zeropoints that were held constant during the fit '''
    for filt_hold in info_hold:    
        #print(filt_hold['mag'] + ' HELD ' + str(filt_hold['ZP']) ) 
        for band in input_info:
            if band['name'] == filt_hold['mag']:
                hold_cal_type = band['blue/red']
                if hold_cal_type == "BLUER/RESTRICTED":
                    hold_cal_type == "BLUER"
                break;

        output_string += filt_hold['mag'] + ' ' + str(filt_hold['ZP']) + ' +- -99 ' + hold_cal_type + '\n'              
   
    print(output_string)
 
    print('NUMBER OF BOOTSTRAPS:', bootstrap_num)
    print('IF ERROR IS -99, NEED TO HAVE > 1 BOOTSTRAP')
    offset_list = output_directory + '/' + file.split('/')[-1]  + '.offsets.list'
    offset_list_file = open(offset_list,'w')
    offset_list_file.write(output_string)
    offset_list_file.close()

    print('LIST OF ZEROPOINTS WRITTEN TO', offset_list)

    print(snpath)


def fit(table, input_info_unsorted, mag_locus,  
        end_of_locus_reject=3,
        min_err=0.02, 
        min_bands_per_star=3, 
        startingzps=None, 
        plot_iteration_increment=50, 
        max_err=0.3, 
        bootstrap=False, 
        bootstrap_num=0, 
        plotdir='.', 
        save_bootstrap_plots=False, 
        pre_zps=None, #HOLDEN# Look into this
        number_of_plots = 10,
        fast=True,
        publish=True ,
        relative_zps=None           
        ):

    os.system('mkdir -p ' + plotdir)
    import numpy as np
    from functools import reduce
    import time


    params = {'backend' : 'ps',
         'text.usetex' : False,
          'ps.usedistiller' : 'xpdf',
          'ps.distiller.res' : 6000}
    mpl.rcParams.update(params)

    fig_size = [5,5]
    params = {'axes.labelsize' : 14,
              'font.size' : 14,
	      'legend.fontsize' : 12,
	      'xtick.labelsize' : 10,
	      'ytick.labelsize' : 10,
              'scatter.marker': 'o',
	      'figure.figsize' : fig_size}
    mpl.rcParams.update(params)

    vary_input_info = list(filter(lambda x: x['HOLD_VARY'] == 'VARY', input_info_unsorted))
    hold_input_info = list(filter(lambda x: x['HOLD_VARY'] == 'HOLD', input_info_unsorted))

    print([a['filter'] for a in vary_input_info])
    print([a['filter'] for a in hold_input_info])

    input_info = hold_input_info + vary_input_info

    zps ={} 
    for i in range(len(vary_input_info)):
        zps[vary_input_info[i]['mag']] = i

    print(zps)

    number_locus_points = len(mag_locus) 
    number_all_stars = len(table.field(input_info[0]['mag']))

    ''' for each point in locus, make a list of the locus in each color (locus has same number of points in each color) '''
    ''' just a rearrangement '''
    locus_list = [] #Locus that is being fit
    ref_locus_list = [] #locus being fit to (not ref catalog, I think)
    for j in range(number_locus_points):
        o = []
        o_ref = []
        for c in input_info:
            print(mag_locus[j].keys())
            o.append(mag_locus[j][c['mag']])
        for c in ['USDSS','GSDSS','RSDSS','ISDSS','ZSDSS']:
            o_ref.append(mag_locus[j][c])
        locus_list.append(o)
        ref_locus_list.append(o_ref)

    results = {} 

    if bootstrap:
        cycles = ['full'] + ['bootstrap' + str(i) for i in range(bootstrap_num)] 
    else:        
        cycles = ['full']

    for iteration in cycles:       


        zps_hold={} 
        for i in range(len(hold_input_info)):
            zps_hold[hold_input_info[i]['mag']] =  hold_input_info[i]['ZP']

            ''' if bootstrap, sample uncertaintities of HOLD bands in bootstrap '''
            if 'ZPERR' in hold_input_info[i] and iteration.find('bootstrap') != -1:
                    
                    import random as rd 
                    zp_err = float(hold_input_info[i]['ZPERR'])
                    if zp_err > 0: 
                        zps_hold[hold_input_info[i]['mag']] += rd.gauss(0,zp_err)
                    #print(zps_hold[hold_input_info[i]['mag']])
                    #print(iteration)
                
                                                                                                                                                                   
        print(zps, zps_hold)

        ''' make matrix with a full set of locus points for each star '''    
        locus_matrix = np.array(number_all_stars*[locus_list])
        ref_locus_matrix = np.array(number_all_stars*[ref_locus_list])
	#print(locus_matrix.shape)

        ''' assemble matricies to make instrumental measured bands '''
        SeqNr = table.field('SeqNr')
        A_band = np.swapaxes(np.swapaxes(np.array(number_locus_points*[[table.field(a['mag']) for a in input_info]]),0,2),1,2)
        n = len(table.field(input_info[0]['mag']))
        def isitJ(name):
            if name.find('JCAT') != -1:
                return np.ones(n)
            else: 
                return np.zeros(n)                

        A_err = np.swapaxes(np.swapaxes(np.array(number_locus_points*[[table.field(a['mag_err']) for a in input_info]]),0,2),1,2)
        #print(A_err.shape)
        ''' only use stars with errors less than max_err '''            

        print(A_band.shape)
        if True:
            mask = A_err > max_err  
            mask[A_err > 1.5] = 1  
            A_band[mask] = 99

        ''' make matrix specifying good values '''

        good = np.ones(A_band.shape)
        good[abs(A_band) == 99] = 0
        good[abs(A_band) == 0] = 0
        good = good[:,0,:]
        good_bands_per_star = good.sum(axis=1) # sum all of the good bands for any given star
        

        ''' figure out the cut-off '''
        index = np.where(good_bands_per_star>=min_bands_per_star)
        SeqNr = SeqNr[index]
        A_band = A_band[index]
        A_err = A_err[index]
        locus_matrix = locus_matrix[index]
        ref_locus_matrix = ref_locus_matrix[index]

        A_err[A_err<min_err] = min_err 

        ''' if a bootstrap iteration, bootstrap with replacement '''
        if iteration.find('bootstrap') != -1:
            length = len(A_band)
            random_indices = []
	    #unique_indices = {}
            for e in range(length): 
                index = int(random.random()*length - 1)
		#unique_indices[index] = 'yes'
                random_indices.append(index)

            #print(random_indices, len(unique_indices.keys()))

            SeqNr = np.array([SeqNr[i] for i in random_indices])             
            A_band = np.array([A_band[i] for i in random_indices])             
            A_err = np.array([A_err[i] for i in random_indices])
            locus_matrix = np.array([locus_matrix[i] for i in random_indices])
            ref_locus_matrix = np.array([ref_locus_matrix[i] for i in random_indices])
        
        bands = A_band 
        bands_err = A_err

        ''' set errors on bad measurements (value=+-99) equal to 100000. and bands equal to 0 '''
        bands_err[abs(A_band) == 99] = 1000.   
        bands[abs(A_band) == 99] = 0.   
        #print(A_band)

        #print(bands.shape, locus_matrix.shape)
        number_good_stars = len(locus_matrix)
        print(number_good_stars) # still have all stars here for panstarrs and gaia

        ''' update good matrix after masking '''
        good = np.ones(A_band.shape) 
        good[abs(A_band) == 99] = 0
        good[abs(A_band) == 0] = 0

        global itr
        itr = 0

        keep_fitting = True
        fit_num = 0
        outliers = 'no outlier rejection'

        while keep_fitting:
            def optimize_offset_errfunc(offset,pars,residuals=False, savefig=None):
                return errfunc(pars,offset=offset)

            def errfunc(pars,residuals=False,savefig=None, offset=0):
                if not offset == 0:
                    pars = [x + offset for x in pars]
                 
                global itr 
                stat_tot = 0
                zp_bands = np.zeros((number_good_stars,number_locus_points,len(input_info))) 
                for i in range(len(input_info)):
                    a = input_info[i]['mag']
                    zp_bands[:,:,i] = assign_zp(a,pars,zps,zps_hold)
                num_prelim = (bands - locus_matrix + zp_bands) / bands_err**2. #holden# numerator of eq 6
                num_prelim[good == 0] = 0.
                num = (num_prelim.sum(axis=2))
                denom_prelim = 1. / bands_err**2. #holden# denom of eq 6
                denom_prelim[good == 0] = 0.
                denom = (denom_prelim.sum(axis=2))
                mean = num / denom
                mean_array = np.dstack(len(input_info)*[mean]) #holden# eq 6

                ds_prelim = (bands - locus_matrix + zp_bands - mean_array)**2. #/ ds_err**2. #holden# numerator of eq 5
                ds_prelim[good == 0] = 0
                ''' calculate reduced chi squared '''
                ds = ds_prelim.sum(axis=2)**0.5 
                resid_prelim = (bands - locus_matrix + zp_bands - mean_array )**2. / bands_err**2. #holden# eq 5
                plot = (bands -locus_matrix + zp_bands - mean_array ) 
                resid_prelim[good == 0] = 0
                resid = resid_prelim.sum(axis=2) / good.sum(axis=2) 

                resid_sum = resid_prelim.sum(axis=2) #/ good.sum(axis=2) #holden# eq 5

                ''' these two are not necessarily the same star '''
                match_locus_index = resid.argmin(axis=1) ## closest locus to each star
                select_diff = resid[np.arange(len(match_locus_index)),match_locus_index]
                select_sum = resid_sum[np.arange(len(match_locus_index)),match_locus_index] #holden# eq 5

                select_good = good[np.arange(len(match_locus_index)),match_locus_index]

                dist = ds[np.arange(len(match_locus_index)),match_locus_index]
                spectrum_normalization = mean[np.arange(len(match_locus_index)),match_locus_index]

                print('good', good.sum())
                
                                                                                   
                #chi_squared_total = select_sum.sum()
                #print("HERE")
                #print(select_sum)
                #print(len(select_sum))
                
                chi_squared_total = np.nansum(select_sum) #holden# verify that it's okay to skip over nan values #holden# eq 5 #means some stars have nans
                data_points = select_good.sum()
                print('data points', data_points)
                print('stars', len(select_good))
                degrees_of_freedom = data_points - (bands.shape[-1] - 1) - 2*len(select_good) - 1
                # degrees of freedom = datapoints - parameters - 1
                
                ''' two fit parameters for each star: median and choice of closest locus point (I think) '''
                #redchi = stat_tot / float(max(1,len(bands) - 1))
                redchi = chi_squared_total / float(degrees_of_freedom)

                ''' compute reference apparent magnitudes of stars ''' 
                norm = np.swapaxes(np.array([spectrum_normalization.tolist()]*5),0,1)
                ref_locus_mags = ref_locus_matrix[np.arange(len(match_locus_index)),match_locus_index,:]
                ref_mags =  norm + ref_locus_mags 

                stat_tot = chi_squared_total #select_diff.sum()

                print('ZPs', dict(zip([a['mag'] for a in input_info] ,([zps_hold[a['mag']] for a in hold_input_info] + ['%.6f' % a for a in list(pars)]))))
                #print("zps_hold values:", [zps_hold[a['mag']] for a in hold_input_info])
                #print("pars:", list(pars))


                

                print('CURRENT TASK:', iteration)
                print('STARS:', len(bands))

                #print('select sum len:', len(select_sum))
                #print('select sum:', select_sum)

                print('chi^2', '%.5f' % stat_tot, )
                print('degrees of freedom', '%d' % degrees_of_freedom, )
                print('red chi^2', '%.5f' % redchi)
                print('iteration', itr)
                
                #holden# toggle plots
                #if iteration == 'full' and (itr % plot_iteration_increment == 0 or savefig is not None):
                #   plot_progress(pars,stat_tot,savefig)
                itr += 1

                if residuals:
                    #print(end_of_locus_reject)
                    end_of_locus = np.array([reduce(lambda x,y: x*y, [match_locus_index[i] != x for x in range(end_of_locus_reject)]) for i in range(len(match_locus_index))])
                    print(select_diff.shape) 
                    print(dist.shape)
                    print(redchi.shape)
                    print(end_of_locus.shape)
                    print(len(bands)) 
                    print(ref_mags.shape)
                    return select_diff, dist, redchi, end_of_locus, len(bands), ref_mags 
                else: return stat_tot

            def plot_progress(pars,stat_tot=None,savefig=None):
                zp_bands = np.zeros((number_good_stars,number_locus_points,len(input_info))) 
                for i in range(len(input_info)):
                    a = input_info[i]['mag']
                    zp_bands[:,:,i] = assign_zp(a,pars,zps,zps_hold)

                if pre_zps:
                    pre_zp_bands = np.swapaxes(np.swapaxes(np.array(number_locus_points*[number_good_stars*[[assign_zp(a[0],pars,pre_zps,zps_hold) for a in input_info]]]),0,1),0,0)
                    pre_zp_bands = np.zeros((number_good_stars,number_locus_points,len(pre_zpz))) 
                    for i in range(len(pre_zps)):
                        a = pre_zps[i]
                        zp_bands[:,:,i] = assign_zp(a[0][0],pars,zps,zps_hold)-assign_zp(a[1][0],pars,zps,zps_hold)
                                                                                                                                                                                                      
                oa = copy(input_info)
                #oa.sort(sort_wavelength)
                oa.sort(key=lambda x: x['center wavelength'])
                
                oa_no_ref = list(filter(lambda x: x['mag'].find('psfMag') == -1 and x['mag'].find('phot_g_mean_mag') == -1 and x['mag'].find('PSFMag') == -1, oa))


                def plot_combinations(input):
                    list = []
                    index = []
                    N = len(input)
                    for a in range(N):
                        for b in range(a+1,N):
                            for c in range(N):
                                for d in range(c+1,N): 
                                    if (len(index)==0) and ( [a,b] != [c,d] ):
                                        index.append( [[a,b], [c,d]] )
                                    else:
                                        flag = True
                                        flag *= not ([a,b] == [c,d])
                                        for [a0,b0], [c0,d0] in index:
                                            flag *= not ([a0,b0,c0,d0] == [c,d,a,b])
                                        if flag:
                                            index.append( [[a,b], [c,d]] )
                    for [a,b], [c,d] in index:
                        list.append([[input[c],input[d]],[input[a],input[b]]])
                    return list


                if savefig is not None:
                    index_list = plot_combinations(oa_no_ref)
                    index_list += plot_combinations(oa)

                else: 
                    #index_list = [[[oa[1]['mag'],oa[2]['mag']],[oa[0]['mag'],oa[1]['mag']]]]
                    index_list = [[[oa[0],oa[1]],[oa[1],oa[2]]]]

                def ind(filt):
                    for j in range(len(input_info)):
                        if input_info[j]['mag'] == filt:
                            return j

		#for [c1_1, c1_2], [c2_1,c2_2] in index_list[:number_of_plots]: 
        #print(c1_2)
                
                for [c1_1, c1_2], [c2_1,c2_2] in index_list: 
                    c1_band1 = c1_1['mag']
                    c1_band2 = c1_2['mag']
                    c2_band1 = c2_1['mag']
                    c2_band2 = c2_2['mag']

                    #print(input_info)
                    #print(ind(c1_band1), ind(c1_band2))
                    #print(ind(c2_band1), ind(c2_band2))
                    #print(c2_band1, c2_band2)

                    print(c1_band1, c1_band2, c2_band1, c2_band2)

                    if ind(c1_band1) is not None and ind(c1_band2) is not None and ind(c2_band1) is not None and ind(c2_band2) is not None:
                        x_color = np.array(bands + zp_bands)[:,0,ind(c1_band1)] - np.array(bands + zp_bands)[:,0,ind(c1_band2)]

                        y_app_mag = np.array(bands + zp_bands)[:,0,ind(c2_band1)] 
                        #print(ind(c2_band1), ind(c2_band2))

                        y_color = (bands + zp_bands)[:,0,ind(c2_band1)] - (bands + zp_bands)[:,0,ind(c2_band2)]

                        if pre_zps:
                            pre_x_color = np.array((bands + pre_zp_bands)[:,0,color1_index].tolist())
                            pre_y_color = (bands + pre_zp_bands)[:,0,color2_index]

                        x_err_1 = (bands_err)[:,0,ind(c1_band1)]
                        x_err_2 = (bands_err)[:,0,ind(c1_band2)]
                        y_err_1 = (bands_err)[:,0,ind(c2_band1)]
                        y_err_2 = (bands_err)[:,0,ind(c2_band2)] #elimator of nonref catalogs
                        

                        mask = (x_err_1<100)*(x_err_2<100)*(y_err_1<100)*(y_err_2<100)
            
                        print(len(x_color))
                        print(len(y_color))
                        x_color = x_color[mask]
                        y_color = y_color[mask]
                        print(len(x_color))
                        print(len(y_color)) #unmatchedstars withh get masked out if rPSF is included in c1_1, c1_2, ..., but if not they stay
                        #print(c1_band1, c1_band2, c2_band1, c2_band2)
                        
                        y_app_mag = y_app_mag[mask]
                        x_err = (x_err_1**2. + x_err_2**2.)**0.5
                        y_err = (y_err_1**2. + y_err_2**2.)**0.5
                        y_err = y_err[mask]
                        x_err = x_err[mask]
                        
                        if pre_zps:
                            pre_x_color = pre_x_color[mask]
                            pre_y_color = pre_y_color[mask]

                        #print(len(x_color), len(x_color))

                        plt.clf()                                                                            
                        plt.axes([0.15,0.125,0.95-0.15,0.95-0.125])

                        x_a = c1_1['plotName'] 
                        x_b = c1_2['plotName'] 
                        y_a = c2_1['plotName'] 
                        y_b = c2_2['plotName'] 
                                       
                        if 'extinction' in input_info[ind(c1_band1)]:
                            x_extinct = input_info[ind(c1_band1)]['extinction'] - input_info[ind(c1_band2)]['extinction']
                            y_extinct = input_info[ind(c2_band1)]['extinction'] - input_info[ind(c2_band2)]['extinction']

                            gallat = '%.1f' % input_info[ind(c1_band1)]['gallat']
                        else: x_extinct, y_extinct, gallat = 0,0,'NA'

                        y_app_mag_name = y_a + ' (mag)'
        
                        units = " ${\\rm (mag)}$"
                        x_color_name = x_a + '  -  ' + x_b + units
                        y_color_name = y_a + '  -  ' + y_b + units 
                        plt.xlabel(x_color_name)
                        plt.ylabel(y_color_name)

                        if len(x_color):
                            plt.scatter(x_color,y_color,color='#0066ff',s=4,marker='o', zorder=20)
                            plt.errorbar(x_color,y_color,xerr=x_err,yerr=y_err,marker=None,fmt='o',ecolor="#e8e8e8",ms=1, mew=1, zorder=1) #,mc='none')   

                            c1_locus = locus_matrix[0,:,ind(c1_band1)] - locus_matrix[0,:,ind(c1_band2)]
                            c2_locus = locus_matrix[0,:,ind(c2_band1)] - locus_matrix[0,:,ind(c2_band2)]
                            plt.plot(c1_locus,c2_locus,'r-',linewidth=1,zorder=30)
                            plt.scatter(c1_locus,c2_locus,color='red',s=7,marker='o',zorder=30)

                            if pre_zps:
                                plt.errorbar(pre_x_color,pre_y_color,xerr=x_err,yerr=y_err,fmt='o',c='green')
                                plt.scatter(pre_x_color,pre_y_color,c='green')

                            x_diff = (c1_locus[-1] - c1_locus[0])
                            y_diff = (c2_locus[-1] - c2_locus[0])
                            plt.arrow(c1_locus[0]+x_diff*0.1,c2_locus[-1]-y_diff*0.1,x_extinct,y_extinct,width=0.01,color='black')

                            if not publish:
                                plt.text(c1_locus[0]+x_diff*0.1 + x_extinct,c2_locus[-1]-y_diff*0.1 + y_extinct,'  ext. vec.',color='black')
                            if stat_tot is not None and not publish:
                                plt.title('N=' + str(len(x_color)) + ' chi$^{2}$=' + ('%.1f' % stat_tot) + ' ' + iteration + ' ' + outliers + ' GALLAT=' + str(gallat))

                            fit_band_zps = reduce(lambda x,y: x + y, [z[-2:].replace('C','').replace('-','') for z in [a['mag'] for a in input_info]])
                            print('savefig', savefig)

                            ''' only save figure if savefig is not None '''
                            if savefig is not None: 
                                if (iteration.find('bootstrap')==-1 or save_bootstrap_plots):
                                    file = plotdir + '/qc_' + fit_band_zps + '_' + x_color_name.replace(units,'').replace(' ','') + '_' + y_color_name.replace(units,'').replace(' ','') + '_' + savefig.replace(' ','_')
                                    file = file.replace('$','')
                                    print(file)

                                   
                                    print(mpl.rcParams['figure.figsize'])
                                    plt.savefig(file)
                                    
                        


                def order_plots(a):
                    # Define a numerical value based on the type of string for sorting
                    if a.find('psfMag') != -1:
                        return 1  # Prioritize 'psfMag' entries
                    elif a.find('PSFMag') != -1:
                        return 2  # Prioritize 'PSFMag' entries
                    elif a.find('phot_g_mean_mag') != -1:
                        return 3  # Prioritize 'phot_g_mean_mag' entries
                    else:
                        return 4  # For other entries, the least priority

                if savefig is not None:
                                                            
                    fs = glob(plotdir + '/qc_*png')                                                                 
                    print(fs)
                    ''' put non-reference plots first '''
                    fs.sort(key=order_plots)
                    print(fs)
                    html = open(plotdir + '/all.html','w')
                    html.write('<html>\n')
                    for type in ['full_outliers_removed','full_egregious_outliers_removed','no_outlier_rejection']:
                        html.write('<h1>' + type + '</h1>\n')
                        for f in fs:                                                     
                            if f.split('/')[-1].find(type) != -1:
                                html.write('<img src=' + f.split('/')[-1] + '></img><br>\n')
                    html.close()

                    
            ''' starting guess for zeropoint : median hold instrumental magnitude - median hold locus magnitude ''' 
            #print(A_band.shape)


            if iteration == 'full': 
                # Could move this block to the not twoStep section?
                if True:
                    print("input info", input_info)
                    print("hold input info", hold_input_info)
                    pinit = []                                                                                                               
                    for i in range(len(hold_input_info),len(input_info)):
                        key = input_info[i]['mag'] ## varying magnitudes
                        info_hold = list(filter(lambda x: x['HOLD_VARY'] == 'HOLD', input_info)) #[0]['mag']
                        ''' calculate average color for actual and model locus '''
                        diff = A_band[:,0,i] - A_band[:,0,0]
                        print(A_band.shape, good.shape)
                        good_diff = good[:,0,i] + good[:,0,0]
                        print(np.sum(good[:,0,i]), 'number of good measurements in band')
                        #diff = diff[good_diff == 2]
                        ind = np.where(good_diff == 2)
                        diff = diff[ind[0]]

                        print(key, len(diff))
                  
                        if len(diff) == 0:
                            print('no stars have good measurements in relevant bands') #error here
                            raise Exception 
                        median_instrumental = np.nanmedian(diff) #holden# verify that this is okay
                        locus_here = [mag_locus[x][input_info[i]['mag']] - mag_locus[x][info_hold[0]['mag']] for x in range(len(mag_locus))]
                        median_locus = np.median(locus_here)
                        #print('diff:', diff)
                        print('median indtrumental', median_instrumental)
                        print('median locus', median_locus)
                        pinit.append(median_locus - median_instrumental)
                #$pinit = [0 for key in [a['mag'] for a in vary_input_info]]
            else:
                ''' add random offset of 1.0 mag '''

                if not fast:
                    pinit = [results['full'][key] + random.random()*1.0 for key in [a['mag'] for a in vary_input_info]]
                else: 
                    pinit = [results['full'][key] for key in [a['mag'] for a in vary_input_info]]

            print(pinit)
            if relative_zps is not None:
                initial_offset = 0
                new_offset = scipy.optimize.fmin(optimize_offset_errfunc,initial_offset,maxiter=10000,maxfun=100000,ftol=0.00001,xtol=0.00001,args=(relative_zps,)) 
                out = np.array([np.float64(val + new_offset) for val in relative_zps])
                print("rel zps: ", relative_zps)
                print("abs zps: ", out)
                #exit()
            else:
                out = scipy.optimize.fmin(errfunc,pinit,maxiter=10000,maxfun=100000,ftol=0.00001,xtol=0.00001,args=())

            print('STARTINGOUT:', out)
            if iteration == 'full':
                errfunc(out,savefig=(iteration+'_'+outliers+'.png').replace('$',''))
            #print(out)

            #print('starting')      

            print(out, )
            print([zps_hold[a['mag']] for a in hold_input_info] )

            #[zps_hold[a['mag']] for a in hold_input_info] + 
            print('ENDINGOUT:', list(out))
            residuals,dist,redchi,end_of_locus, num, ref_mags = errfunc(pars=list(out),residuals=True)
            #print(dist)
            #print('finished')
            #print('bands' , len(bands))
                                                                                      
            #print(end_of_locus)
            #print(bands.shape)
            #print(dist.shape, residuals.shape)

            if fit_num == 0:                                                                          
                resid_thresh = 30
                print(residuals)
                bands = bands[residuals < resid_thresh]
                bands_err = bands_err[residuals < resid_thresh]
                locus_matrix = locus_matrix[residuals < resid_thresh]
                ref_locus_matrix = ref_locus_matrix[residuals < resid_thresh]
                SeqNr = SeqNr[residuals < resid_thresh]
                good = good[residuals < resid_thresh]
                end_of_locus = end_of_locus[residuals < resid_thresh]
                dist = dist[residuals < resid_thresh]
                ref_mags = ref_mags[residuals < resid_thresh]

            else: 
                resid_thresh = 6
                bands = bands[residuals < resid_thresh]
                bands_err = bands_err[residuals < resid_thresh]
                locus_matrix = locus_matrix[residuals < resid_thresh]
                ref_locus_matrix = ref_locus_matrix[residuals < resid_thresh]
                SeqNr = SeqNr[residuals < resid_thresh]
                good = good[residuals < resid_thresh]
                end_of_locus = end_of_locus[residuals < resid_thresh]
                dist = dist[residuals < resid_thresh]
                ref_mags = ref_mags[residuals < resid_thresh]

                ''' first filter on distance '''
                ind = np.where(dist<3)
                bands = bands[ind]
                bands_err = bands_err[ind]
                locus_matrix = locus_matrix[ind]
                ref_locus_matrix = ref_locus_matrix[ind]
                SeqNr = SeqNr[ind]
                good = good[ind]
                residuals = residuals[ind]
                end_of_locus = end_of_locus[ind]
                ref_mags = ref_mags[ind]

                if True:
                    ''' filter on end of locus '''                      
                    bands = bands[end_of_locus]
                    bands_err = bands_err[end_of_locus]
                    locus_matrix = locus_matrix[end_of_locus]
                    ref_locus_matrix = ref_locus_matrix[end_of_locus]
                    SeqNr = SeqNr[end_of_locus]
                    good = good[end_of_locus]
                    ref_mags = ref_mags[end_of_locus]

            #print(number_good_stars, len(locus_matrix))


            fit_num += 1                                                                                      

            if fit_num == 1:
                print('REFITTING AFTER REMOVING EGREGIOUS OUTLIERS ')
                outliers = 'egregious outliers removed '# + str(resid_thresh)
                number_good_stars = len(locus_matrix)
                print(str(number_good_stars), 'STARS LEFT')
            elif number_good_stars > len(locus_matrix) or len(list(filter(lambda x: x is False,end_of_locus.tolist()))) > 0 :
                print('REFITTING AFTER REMOVING ' + str(number_good_stars - len(locus_matrix) ) + ' OUTLIERS AND STARS MATCHING BLUE END OF LOCUS')
                number_good_stars = len(locus_matrix)

                print(str(number_good_stars), 'STARS LEFT')
                #print('bands' , len(bands))                                     
                #print(bands.shape, locus_matrix.shape)
                pinit = out 
                outliers = 'outliers removed'
    
            else:
                print('NO OUTLYING STARS OR STARS MATCHING BLUE END OF LOCUS, PROCEEDING')
                keep_fitting = False

        print(str(number_good_stars), 'STARS LEFT') #110 stars here for pan, 692 for Gaia

        results[iteration] = dict(zip([a['mag'] for a in input_info],([zps_hold[a['mag']] for a in hold_input_info] + out.tolist())))
        results['ref_mags_' + iteration] = copy(ref_mags)
        results['SeqNr_' + iteration] = copy(SeqNr)


        mask = bands_err < 100 #Filters out nonmatches in ext cat plots (plots w/ those band(s))
        print()

    results['redchi'] = redchi
    results['num'] = num        

    #print(results) 
    errors = {}
    bootstraps = {}
    #print('BOOTSTRAPPING ERRORS:')
    print(input_info)
    
    print([a['mag'] for a in input_info])
    
    for key in [a['mag'] for a in input_info]:
        l = []
        print(results.keys())
        for r in results.keys():
            if r != 'full' and r != 'redchi' and r != 'num' and r.find('ref_mags') == -1 and r.find('SeqNr') == -1:
                print(r, key)
                l.append(results[r][key])
        #print(key+':', np.std(l), 'mag')
       
        if len(l) > 1: 
            errors[key] = '%.4f' % np.std(l)
        else: errors[key] = -99

        #scipy.cov(np.array(l))


        if bootstrap_num > 0 and len(l) > 0:
            bootstraps[key] = reduce(lambda x,y: x + ',' + y, [str(z) for z in l])
        else: bootstraps[key] = 'None'

    results['hold_vary'] = dict(zip([a['mag'] for a in input_info],[a['HOLD_VARY'] for a in input_info]))
    results['bootstraps'] = bootstraps
    results['errors'] = errors
    results['bootstrapnum'] = bootstrap_num 

    print(str(number_good_stars), 'STARS LEFT') #108 stars here for pan, 715 for gaia

    return results, results['ref_mags_full'], results['SeqNr_full']

if __name__ == '__main__':
    #all(subarudir,cluster,DETECT_FILTER,AP_TYPE,magtype)

#if True: # __name__ == '__main__':

    from optparse import OptionParser

    usage = "usage: python fit_locus.py [options] --help \n\nGiven catalog of stellar magnitudes and total (atmosphere+mirrors+optics+filter+CCD) response, \ncomputes expected stellar locus, and fits for zeropoint calibration. \nRequires description of the columns in the input FITS table.\n\nExample: python fit_locus.py -f stars.fits -c stars.columns -e 1 -b 10"

    parser = OptionParser(usage)
    parser.add_option("-f","--file",help="FITS catalog file")
    parser.add_option("-e","--extension",help="extension of FITS file containing stellar magnitudes (number or name) (default: 1)",default=1)
    parser.add_option("-b","--bootstrap",type="int",help="number of bootstraps for error estimation (default: 0)",default=0)
    parser.add_option("-c","--columns",help="column description file")
    parser.add_option("-o","--output",help="output calibration file directory location, if different from directory containing catalog file",default=None)
    parser.add_option("-p","--plots",help="destination directory for plots, if different from /output directory/PLOTS",default=None)
    parser.add_option("-r","--racol",help="name of column in FITS file with object RA in DEGREES (default: X_WORLD)",default='X_WORLD')
    parser.add_option("-d","--deccol",help="name of column in FITS file with object DEC in DEGREES (default: XWORLD)",default='Y_WORLD')
    parser.add_option("-z","--SN",help="snpath",default=None)
    parser.add_option("-t","--run",help="run",default=None)
    parser.add_option("-n","--night",help="night",default=None)
    parser.add_option("-s","--addSDSSgriz",action='store_true',help="automatically search for and add SDSS griz stellar photometry, if available")
    parser.add_option("-j","--add2MASSJ",action='store_true',help="automatically search for and add 2MASS J stellar photometry, if available")
    parser.add_option("-a","--addPanSTARRS",action='store_true',help="automatically search for and add PanSTARRS griz stellar photometry, if available")
    parser.add_option("-g","--addGaia",action='store_true',help="automatically search for and add Gaia dr2 G band stellar photometry")
    parser.add_option("-w","--numberofplots",help="number of plots to make (default: 10)",default=10)
    parser.add_option("-u","--sdssUnit",help="run SDSS unit test (only works if in coverage)",action='store_true')
    parser.add_option("-2","--twoStep",action='store_true',help="run 2 step process with external catalogs. Relative ZPs calculated first, then scaled by an external catalog for abs ZPs. Hold 1 band in columns file and add 1 external catalog.")
    
    import sys

    args = sys.argv     
    
    #args = ['-f','stars.fits','-c','sdss.columns','-e','1','-a']
    #args = ['-f','A383.fits','-c','sdss.columns','-e','1','0','-a']
    #args = ['-f','MACS0717+37.stars.calibrated.cat','-c','sdss.columns','-e','1','-b','0','-a']
    #args = ['-f','HDFN.fits','-c','sdss.columns','-e','1','-b','0','-a']
    #args = ['-f','A2552.stars.calibrated.cat','-c','A2552.columns','-e','1','-b','0']
    #args = ['-f','A2552.stars.calibrated.cat','-c','A2552.columns','-e','1','-b','0']
    #args = ['-f','MACS1347-11.stars.calibrated.cat','-c','MACS1347-11.columns','-e','1','-b','0','-l','False']

    (options, args) = parser.parse_args(args)

    if options.file is None: 
        parser.error('you must specify an input FITS catalog file')
    elif options.columns is None: 
        parser.error('you must specify a file specifying the input magnitudes and corresponding filters')
    elif options.extension is None: 
        parser.error('you must specify the extension of the input FITS catalog containing the stellar magnitudes')

    print('importing libraries')
    import os, re, string
    import random, scipy #dbm #commands,anydm
    from astropy.io import fits
    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib import pyplot as plt
    import numpy as np
    from astroquery.gaia import Gaia
    from astropy.table import Table
    from scipy import linalg
    from scipy import optimize
    from glob import glob
    from copy import copy
    import utilities
    print('finished importing libraries')

    if options.twoStep: #could put in the fit function, like how blue locus is done (blue locus might also validate two step strategy)
        #run without external catalog first for relative ZP's
        run(options.file,options.columns,output_directory=options.output,plots_directory=options.plots,extension=options.extension,racol=options.racol,deccol=options.deccol,bootstrap_num=options.bootstrap, add2MASS=False, addSDSS=False, addPanSTARRS=False, addGaia=False, number_of_plots=options.numberofplots, sdssUnit=False, twoStep=False)  
        
        def write_external_columns_file(input_filename, output_filename):
            try:
                with open(input_filename, 'r') as infile:
                    lines = infile.readlines()

                # find the filter being held, should only be one
                new_lines = []

                for line in lines:
                    if "HOLD" in line:
                        # make the HELD filter variable for the absolute ZP calibration
                        hold_pos = line.find("HOLD")
                        modified_line = line[:hold_pos] + "VARY\n" # Replace everything after and including HOLD with VARY
                        new_lines.append(modified_line)
                    else:
                        new_lines.append(line)

                with open(output_filename, 'w') as outfile:
                    outfile.writelines(new_lines)
            except FileNotFoundError:
                print(f"Error: The file '{input_filename}' was not found.")
            except Exception as e:
                print(f"An error occurred: {e}")

        old_columns_file = options.columns
        options.columns += ".exter_anchor"

        write_external_columns_file(old_columns_file, options.columns)

    run(options.file,options.columns,output_directory=options.output,plots_directory=options.plots,extension=options.extension,racol=options.racol,deccol=options.deccol,bootstrap_num=options.bootstrap, add2MASS=options.add2MASSJ, addSDSS=options.addSDSSgriz, addPanSTARRS=options.addPanSTARRS, addGaia=options.addGaia, number_of_plots=options.numberofplots, sdssUnit=options.sdssUnit, twoStep=options.twoStep)
       
