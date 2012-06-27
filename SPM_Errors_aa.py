import asciitable
import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate
import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime

from bad_times import bad_times

SAFEMODE_2012150 = '2012:150:03:33:29'

def plot_spm_input_errs(out, savefigs=False):
    ok = ~out['spm_act_bad']
    times= out['times'][ok]
    spm_act = out['spm_act'][ok]
    sun_prs = out['alpha_sun'][ok] & out['beta_sun'][ok] 
    fss_pitch_err = out['beta'][ok] - out['pitch'][ok]
    fss_roll_err = out['alpha'][ok] - out['roll'][ok]
    fss_err = sqrt(fss_pitch_err**2 + fss_roll_err**2)
    css_pitch_err = out['pitch_css'][ok] - out['pitch'][ok]
    css_roll_err = out['roll_css'][ok] - out['roll'][ok]
    css_err = sqrt(css_pitch_err**2 + css_roll_err**2)
    pitch = out['pitch'][ok]
    roll = out['roll'][ok]
    print 'Points w/o sun presence:'
    print(sum(~sun_prs))
    print 'Points w/ sun presence:'
    print(sum(sun_prs))
    
    #Plot FSS Errors vs Time
    figure(1)
    plot_cxctime(times[sun_prs & ~spm_act], fss_err[sun_prs & ~spm_act], 'b.', 
                 mec='b', label='SPM Disabled')
    plot_cxctime(times[sun_prs & spm_act], fss_err[sun_prs & spm_act], 'r.', 
                 mec='r', label='SPM Enabled')
    legend(loc='upperleft')
    grid()
    ylabel('FSS Error [deg]')
    title('FSS Error (with Sun Presence) \n fed into Sun Position Monitor')    
    if savefigs==True:
        savefig('spm_fss_errors_vs_time.png')
        
    #Plot CSS Errors vs Time
    figure(2)
    plot_cxctime(times[sun_prs & ~spm_act], css_err[sun_prs & ~spm_act], 'b.',
                 mec='b', label='SPM Disabled')
    plot_cxctime(times[sun_prs & spm_act], css_err[sun_prs & spm_act], 'r.',
                 mec='r', label='SPM Enabled')
    legend(loc='upper left')
    grid()
    ylabel('CSS Error [deg]')
    title('CSS Error (with FSS Sun Presence) \n fed into Sun Position Monitor') 
    if savefigs==True:
        savefig('spm_css_errors_vs_time.png')
        
    #Plot CSS Errors vs Attitude
    figure(3)
    scatter(roll, pitch, c=css_err, edgecolors='none')
    c = colorbar()
    c.set_label('CSS Error [deg]')
    xlabel('Roll [deg]')
    ylabel('Pitch [deg]')
    title('CSS Errors vs Attitude')
    grid()
    if savefigs==True:
        savefig('spm_css_errors_vs_att.png')
        
    #Plot Attitudes vs Time
    figure(4)
    scatter(css_roll_err, css_pitch_err, c=times, edgecolors='none')
    c = colorbar()
    xlabel('Roll (deg)')
    ylabel('Pitch (deg)')
    title('CSS Errors vs Time')
    grid()
    if savefigs==True:
        savefig('spm_something.png')

def plot_css_by_year(out, savefigs=False):
    print('Warning:  These plots assume a start time of 2000:001')
    ok = ones(len(out)).astype(bool)
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
        scatter(roll, pitch, c=css_roll_err, edgecolors='none')
        title(str(2000 + yr) + ' CSS Roll Errors')
        xlabel('Roll Angle [deg]')
        ylabel('Pitch Angle [deg]')
        colorbar()
        #clim([0.1,0.5])
        xlim([-30,30])
        ylim([20,200])
        savefig('css_roll_errors_' + str(yr+2000) + '.png')
        close()
        
        figure()
        scatter(roll, pitch, c=css_pitch_err, edgecolors='none')
        title(str(2000 + yr) + ' FSS Pitch Errors')
        xlabel('Roll Angle [deg]')
        ylabel('Pitch Angle [deg]')
        colorbar()
        #clim([-.5,2])
        xlim([-25,20])
        ylim([40,160])
        savefig('css_pitch_errors_' + str(yr+2000) + '.png')
        close()
        
        t = t + dt
        yr = yr + 1
    
def get_data(start='2005:001', stop=SAFEMODE_2012150, interp=32.8,
             pitch0=45, pitch1=180):
    msids = ('aopssupm', 'aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aoalpang', 'aobetang', 'aoalpsun', 'aobetsun',
             'pitch_css', 'roll_css')
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
                'alpha_sun', 'beta_sun', 'spm_act', 'spm_act_bad', 'kalman')
    dtypes = ('f8',
              'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
              'bool', 'bool', 'bool', 'bool', 'bool')
    out = np.empty(nvals, dtype=zip(colnames, dtypes))

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
        bad_times = msid_bad_times.get(msid_self.MSID, [])
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

    