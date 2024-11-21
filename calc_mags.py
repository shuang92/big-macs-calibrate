import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from glob import glob


###########################################
def get_filter_name(cat):

    mags = []
    flux_cols = []
    flux_errs = []
    
    for col in cat.columns:
        if 'flux' in col and 'err' not in col:
            flux_cols.append(col)
            mags.append(col[5:])
        if 'flux' in col and 'err' in col:
            flux_errs.append(col)

    return mags, flux_cols, flux_errs

def calc_mag(cat, mags, flux_cols, flux_errs):

    for mag, flux, err in zip(mags, flux_cols, flux_errs):

        cat[mag] = -2.5*np.log10(cat[flux]) + 27
        cat[mag + '_err'] = 2.5*cat[err]/cat[flux]
        cat.remove_columns([flux, err])

    return cat

##########################################
file_list = glob('flux_files/*')
for file in file_list:
    print(file)
    file_name = file[len('flux_files/'):]
    print(file_name)
    cat = Table.read(file, format='fits')
    mags, flux_cols, flux_errs = get_filter_name(cat)
    cat = calc_mag(cat, mags, flux_cols, flux_errs)
    for mag in mags:
        plt.scatter(cat[mag], cat[mag+'_err'], s=5, label=mag)
    plt.legend(fontsize=14)
    plt.xlabel('Magnitudes', fontsize=14)
    plt.ylabel('Mag Err', fontsize=14)
    plt.title(file_name)
    plt.savefig(file_name + '.png')
    plt.clf()
    for mag in mags:
        ind, = np.where(cat[mag+'_err'] > 0.1)
        cat.remove_rows(ind)
        plt.scatter(cat[mag], cat[mag+'_err'], s=5, label=mag)
    plt.legend(fontsize=14)
    plt.xlabel('Magnitudes', fontsize=14)
    plt.ylabel('Mag Err', fontsize=14)
    plt.title(file_name + ' Mag Err < 0.1')
    plt.savefig(file_name + '_clipped.png')
    plt.clf()
    cat.write('mag_files/' + file_name, format='fits', overwrite=True)
