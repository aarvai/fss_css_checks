import asciitable
import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate
import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime

from bad_times import bad_times

SAFEMODE_2012150 = '2012:150:03:33:29'

def plot_fss_errors(out, savefigs=False):
    ok = ~out['spm_act_bad']
    times= out['times'][ok]
    spm_act = out['spm_act'][ok]
    sun_prs = out['alpha_sun'][ok] & out['beta_sun'][ok]
    fss_pitch_err = out['beta'][ok] - out['pitch'][ok]
    fss_roll_err = out['alpha'][ok] - out['roll'][ok]
    fss_err = sqrt(fss_pitch_err**2 + fss_roll_err**2)
    pitch = out['pitch'][ok]
    roll = out['roll'][ok]
    print 'Points w/o sun presence:'
    print(sum(~sun_prs))
    print 'Points w/ sun presence:'
    print(sum(sun_prs))
    
    #Plot FSS Errors vs Time
    figure()
    plot_cxctime(times[sun_prs & ~spm_act], fss_err[sun_prs & ~spm_act], 'b.', 
                 mec='b', label='SPM Disabled')
    plot_cxctime(times[sun_prs & spm_act], fss_err[sun_prs & spm_act], 'r.', 
                 mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('FSS Error [deg]')
    title('Total FSS Error (with Sun Presence) \n fed into Sun Position Monitor')    
    if savefigs==True:
        savefig('fss_errors_vs_time.png')
    
    figure()
    subplot(2, 1, 1)
    plot_cxctime(times[sun_prs & ~spm_act], fss_roll_err[sun_prs & ~spm_act], 'b.', 
                 mec='b', label='SPM Disabled')
    plot_cxctime(times[sun_prs & spm_act], fss_roll_err[sun_prs & spm_act], 'r.', 
                 mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('FSS Roll Error [deg]')
    title('FSS Roll and Pitch Errors (with Sun Presence) \n fed into Sun Position Monitor')    
    subplot(2, 1, 2)
    plot_cxctime(times[sun_prs & ~spm_act], fss_pitch_err[sun_prs & ~spm_act], 'b.', 
                 mec='b', label='SPM Disabled')
    plot_cxctime(times[sun_prs & spm_act], fss_pitch_err[sun_prs & spm_act], 'r.', 
                 mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('FSS Pitch Error [deg]')
    if savefigs==True:
        savefig('fss_errors_vs_time2.png')

        
def plot_css_errors(out, savefigs=False):
    ok = ~out['spm_act_bad'] & ~out['eclipse'] & ~out['low_alt']
    times= out['times'][ok]
    spm_act = out['spm_act'][ok]
    eclipse = out['eclipse'][ok]
    low_alt = out['low_alt'][ok]
    css_pitch_err = out['pitch_css'][ok] - out['pitch'][ok]
    css_roll_err = out['roll_css'][ok] - out['roll'][ok]
    css_err = sqrt(css_pitch_err**2 + css_roll_err**2)
    pitch = out['pitch'][ok]
    roll = out['roll'][ok]
    
    #Plot CSS Errors vs Time
    figure()
    plot_cxctime(times[~spm_act], css_err[~spm_act], 
                 'b.', mec='b', label='SPM Disabled')
    plot_cxctime(times[spm_act ], css_err[spm_act], 
                 'r.', mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('CSS Error [deg]')
    title('Total CSS Error (excluding eclipses and low altitudes) \n fed into Sun Position Monitor') 
    if savefigs==True:
        savefig('css_errors_vs_time.png')
    
    figure()
    subplot(2, 1, 1)
    plot_cxctime(times[~spm_act], css_roll_err[~spm_act],  
                 'b.', mec='b', label='SPM Disabled')
    plot_cxctime(times[spm_act], css_roll_err[spm_act],  
                 'r.', mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('CSS Roll Error [deg]')
    title('CSS Roll and Pitch Errors (excluding eclipses and low altitudes) \n fed into Sun Position Monitor')    
    subplot(2, 1, 2)
    plot_cxctime(times[~spm_act], css_pitch_err[~spm_act], 
                 'b.', mec='b', label='SPM Disabled')
    plot_cxctime(times[spm_act], css_pitch_err[spm_act], 
                 'r.', mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('CSS Pitch Error [deg]')
    if savefigs==True:
        savefig('css_errors_vs_time2.png')    
    
    zipvals = zip((css_err, css_roll_err, css_pitch_err),
                  ('Total CSS Error', 'CSS Roll Error', 'CSS Pitch Error'),
                  ('css_errors', 'css_roll_errors', 'css_pitch_errors'))
    
    #Plot CSS Errors vs Attitude
    for var, name, plot_name in zipvals:
        figure()
        scatter(roll, pitch, c=var, edgecolors='none')
        c = colorbar()
        c.set_label(name + ' [deg]')
        xlabel('Roll [deg]')
        ylabel('Pitch [deg]')
        title(name + ' vs Attitude \n (Excludes eclipses and low altitudes)')
        grid()
        if savefigs==True:
            savefig(plot_name + '_vs_att.png')

def plot_css_errors_by_year(out, savefigs=False):
    if min(out['times']) > 63158464:
        print('Warning:  plot_css_errs_by_year assumes a start time of 2000:001')
    ok = ~out['eclipse'] & ~out['low_alt']
    times= out['times'][ok]
    css_pitch_err = out['pitch_css'][ok] - out['pitch'][ok]
    css_roll_err = out['roll_css'][ok] - out['roll'][ok]
    css_err = sqrt(css_pitch_err**2 + css_roll_err**2)
    pitch = out['pitch'][ok]
    roll = out['roll'][ok]
    t = times[0]
    dt = 3600 * 24 * 365
    yr = 0

    while t < times[-1]:
        i = (times > t) & (times < t + dt)
      
        figure()
        scatter(roll[i], pitch[i], c=css_roll_err[i], edgecolors='none')
        title(str(2000 + yr) + ' CSS Roll Errors \n (Excludes eclipses and low altitudes)')
        xlabel('Roll Angle [deg]')
        ylabel('Pitch Angle [deg]')
        c = colorbar()
        clim([floor(min(css_roll_err)), ceil(max(css_roll_err))])
        c.set_label('CSS Roll Error [deg]')
        xlim([-30,30])
        ylim([20,200])
        grid()
        savefig('css_roll_errors_' + str(yr+2000) + '.png')
        close()
        
        figure()
        scatter(roll[i], pitch[i], c=css_pitch_err[i], edgecolors='none')
        title(str(2000 + yr) + ' FSS Pitch Errors \n (Excludes eclipses and low altitudes)')
        xlabel('Roll Angle [deg]')
        ylabel('Pitch Angle [deg]')
        c = colorbar()
        clim([floor(min(css_pitch_err)), ceil(max(css_pitch_err))])
        c.set_label('CSS Pitch Error [deg]')
        xlim([-30,30])
        ylim([20,200])
        grid()
        savefig('css_pitch_errors_' + str(yr+2000) + '.png')
        close()
        
        t = t + dt
        yr = yr + 1


def plot_css_errors_bin(out, savefigs=False, pitch_bin=90, roll_bin=0):
    ok = ~out['eclipse'] & ~out['low_alt']
    times= out['times'][ok]
    css_pitch_err = out['pitch_css'][ok] - out['pitch'][ok]
    css_roll_err = out['roll_css'][ok] - out['roll'][ok]
    css_err = sqrt(css_pitch_err**2 + css_roll_err**2)
    pitch = out['pitch'][ok]
    roll = out['roll'][ok]
    bin = (abs(pitch - pitch_bin) < 1) & (abs(roll - roll_bin) < 1)
    
    figure()    
    subplot(2, 1, 1)
    plot_cxctime(times[bin], abs(css_roll_err[bin]), '.')
    grid()
    ylabel('CSS Roll Error [deg]')
    title('CSS Roll and Pitch Errors (excluding eclipses and low altitudes) \n' + 
          'within 1 deg of ' + str(pitch_bin) + ' deg pitch and ' + str(roll_bin) + ' deg roll')    
    subplot(2, 1, 2)
    plot_cxctime(times[bin], abs(css_pitch_err[bin]), '.')
    grid()
    ylabel('CSS Pitch Error [deg]')
    if savefigs==True:
        savefig('css_errors_' + str(pitch_bin) + '_' + str(roll_bin) + '.png')  
    
    
def get_spm_data(start='2000:001', stop=SAFEMODE_2012150, interp=32.8,
             pitch0=45, pitch1=180):
    msids = ('aopssupm', 'aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aoalpang', 'aobetang', 'aoalpsun', 'aobetsun',
             'pitch_css', 'roll_css', 'Dist_SatEarth', 'Sun_EarthCentAng')
    print 'fetching data'
    x = fetch.MSIDset(msids, start, stop)

    # Resample MSIDset (values and bad flags) onto a common time sampling
    print 'starting interpolate'
    x.interpolate(interp, filter_bad=False)

    # Remove data during times of known bad or anomalous data (works as of
    # Ska.engarchive 0.19.1)
    x.filter_bad_times(table=bad_times)

    # Select data only in a limited pitch range
    ok = ((x['pitch'].vals > pitch0) &
          (x['pitch'].vals < pitch1))

    # Determine the logical-or of bad values for all MSIDs and use this
    # to further filter the data sample
    nvals = np.sum(ok)
    bads = np.zeros(nvals, dtype=bool)
    for msid in x.values():
        # Ignore sun position monitor for bad data because it is frequently
        # bad (not available in certain subformats including SSR)
        if msid.MSID == 'AOPSSUPM':
            continue
        print msid.msid, sum(msid.bads[ok])
        bads = bads | msid.bads[ok]
    ok[ok] = ok[ok] & ~bads

    nvals = np.sum(ok)
    colnames = ('times',
                'pitch', 'roll', 'alpha', 'beta', 'pitch_css', 'roll_css',
                'alpha_sun', 'beta_sun', 'spm_act', 
                'spm_act_bad', 'kalman', 'eclipse', 'low_alt')
    dtypes = ('f8',
              'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
              'bool', 'bool', 'bool', 
              'bool', 'bool', 'bool', 'bool')
    out = np.empty(nvals, dtype=zip(colnames, dtypes))
    
    # Define eclipse flag using ephemeris data
    # Add 2 angular degree buffer since ephemeris data only at 5 min intervals
    rad_earth = 6378100 
    ang_rad_earth = arctan(rad_earth / x['Dist_SatEarth'].vals[ok]) * 180 / pi
    ang_rad_sun = .25 
    eclipse = x['Sun_EarthCentAng'].vals[ok] < ang_rad_earth + ang_rad_sun + 2

    # Define low altitude as being below 25,000 km
    low_alt = x['Dist_SatEarth'].vals[ok] < 25000000

    out['times'][:] = x['pitch'].times[ok]
    out['pitch'][:] = x['pitch'].vals[ok]
    out['roll'][:] = x['roll'].vals[ok]
    out['alpha'][:] = -x['aoalpang'].vals[ok]
    out['beta'][:] = 90 - x['aobetang'].vals[ok]
    out['pitch_css'][:] = x['pitch_css'].vals[ok]
    out['roll_css'][:] = x['roll_css'].vals[ok]
    out['alpha_sun'][:] = x['aoalpsun'].vals[ok] == 'SUN '
    out['beta_sun'][:] = x['aobetsun'].vals[ok] == 'SUN '
    out['spm_act'][:] = x['aopssupm'].vals[ok] == 'ACT '
    out['spm_act_bad'][:] = x['aopssupm'].bads[ok]
    out['kalman'][:] = ((x['aoacaseq'].vals[ok] == 'KALM') &
                        (x['aopcadmd'].vals[ok] == 'NPNT'))
    out['eclipse'][:] = eclipse    
    out['low_alt'][:] = low_alt
    return out


def filter_bad_times(msid_self, start=None, stop=None, table=None):
    """Filter out intervals of bad data in the MSID object.

    There are three usage options:

    - Supply no arguments.  This will use the global list of bad times read
      in with fetch.read_bad_times().
    - Supply both ``start`` and ``stop`` values where each is a single
      value in a valid DateTime format.
    - Supply an ``table`` parameter in the form of a 2-column table of
      start and stop dates (space-delimited) or the name of a file with
      data in the same format.

    The ``table`` parameter must be supplied as a table or the name of a
    table file, for example::

      bad_times = ['2008:292:00:00:00 2008:297:00:00:00',
                   '2008:305:00:12:00 2008:305:00:12:03',
                   '2010:101:00:01:12 2010:101:00:01:25']
      msid.filter_bad_times(table=bad_times)
      msid.filter_bad_times(table='msid_bad_times.dat')

    :param start: Start of time interval to exclude (any DateTime format)
    :param stop: End of time interval to exclude (any DateTime format)
    :param table: Two-column table (start, stop) of bad time intervals
    """
    if table is not None:
        bad_times = asciitable.read(table, Reader=asciitable.NoHeader,
                                    names=['start', 'stop'])
    elif start is None and stop is None:
        raise ValueError('filter_times requires 2 args ')
    elif start is None or stop is None:
        raise ValueError('filter_times requires either 2 args '
                         '(start, stop) or no args')
    else:
        bad_times = [(start, stop)]

    ok = np.ones(len(msid_self.times), dtype=bool)
    for start, stop in bad_times:
        tstart = DateTime(start).secs
        tstop = DateTime(stop).secs
        if tstart > tstop:
            raise ValueError("Start time %s must be less than stop time %s"
                             % (start, stop))

        if tstop < msid_self.times[0] or tstart > msid_self.times[-1]:
            continue

        i0, i1 = np.searchsorted(msid_self.times, [tstart, tstop])
        ok[i0:i1 + 1] = False

    colnames = (x for x in msid_self.colnames)
    for colname in colnames:
        attr = getattr(msid_self, colname)
        if isinstance(attr, np.ndarray):
            setattr(msid_self, colname, attr[ok])       
            