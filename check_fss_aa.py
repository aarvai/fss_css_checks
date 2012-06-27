def plot_pitches(out, angle_err_lim=8.0, savefigs=False):
    times= out['times']
    pitch = out['pitch']
    alpha_err = out['alpha'] - out['roll']
    sun = out['alpha_sun'] & out['beta_sun']
    kalman = out['kalman']
    bad = abs(alpha_err) > angle_err_lim

    # Sun Presence
    zipvals = zip((~sun, sun),
                  ('c.', 'r.'),
                  ('c', 'r'),
                  ('No Sun Presence', 'Sun Presence'))
    figure(1)
    for filt, mark, mec, label in zipvals:
        plot_cxctime(times[bad & filt], pitch[bad & filt], mark, 
                     mec=mec, label=label)
        legend(loc='lower left')
        grid('on') 
        title("Pitch Angles when Alpha's Error > " + str(angle_err_lim) + 
              " deg \n Sun Presence vs No Sun Presence")
        ylabel("Pitch Angle [deg]")    
    if savefigs==True:  savefig('pitches_sunprs.png')
    
    # Kalman
    zipvals = zip((~kalman, kalman),
                  ('c.', 'r.'),
                  ('c', 'r'),
                  ('Not Kalman', 'Kalman'))
    figure(2)
    for filt, mark, mec, label in zipvals:
        plot_cxctime(times[bad & filt], pitch[bad & filt], mark, 
                     mec=mec, label=label)
        legend(loc='lower left')
        grid('on') 
        title("Pitch Angles when Alpha's Error > " + str(angle_err_lim) + 
              " deg \n Kalman vs Not Kalman")
        ylabel("Pitch Angle [deg]")  
    if savefigs==True:  savefig('pitches_kalm.png')
    
def plot_temps(out, angle_err_lim=8.0, savefigs=False):
    print('Warning:  These plots assume 138 < pitch < 140 and data from 2011:001 - 2012:150')
    times= out['times']
    bkt_temp1 = out['bkt_temp1']
    bkt_temp2 = out['bkt_temp2']
    tcylaft6 = out['tcylaft6']
    fsse_temp = out['fsse_temp']
    alpha_err = out['alpha'] - out['roll']
    bad = abs(alpha_err) > angle_err_lim

    zipvals = zip((bkt_temp1, bkt_temp2, tcylaft6, fsse_temp),
                  ('Bracket Temp 1', 'Bracket Temp 2', 'TCYLAFT6', 'FSS Electronics Temp'),
                  (3, 4, 5, 6),
                  ('temp_bkt1', 'temp_bkt2', 'temp_tcylaft6', 'temp_fsse'))
    for msid, name, fig, figname in zipvals:
        figure(fig)
        subplot(2,1,1)
        hist(msid, 50, range=(msid.min(), msid.max()))
        title(name + ':\n All Points (2011:001 - 2012:150, 138 < pitch < 140.5)')
        x = xlim()
        subplot(2,1,2)
        hist(msid[bad], 50, range=(msid.min(), msid.max()))
        title('Above + Alpha Error > ' + str(angle_err_lim) + ' deg')
        xlabel('Deg F')
        xlim(x)
        if savefigs==True:  savefig(figname + '.png')
   
def plot_temps_log(out, angle_err_lim=8.0, savefigs=False):
    print('Warning:  These plots assume 138 < pitch < 140 and data from 2011:001 - 2012:150')
    times= out['times']
    bkt_temp1 = out['bkt_temp1']
    bkt_temp2 = out['bkt_temp2']
    tcylaft6 = out['tcylaft6']
    fsse_temp = out['fsse_temp']
    alpha_err = out['alpha'] - out['roll']
    bad = abs(alpha_err) > angle_err_lim

    zipvals = zip((bkt_temp1, bkt_temp2, tcylaft6, fsse_temp),
                  ('Bracket Temp 1', 'Bracket Temp 2', 'TCYLAFT6', 'FSS Electronics Temp'),
                  (3, 4, 5, 6),
                  ('temp_bkt1', 'temp_bkt2', 'temp_tcylaft6', 'temp_fsse'))
    for msid, name, fig, figprefix in zipvals:
        figure(fig, figsize=(8,11))
        subplot(2,1,1)
        hist_all = hist(msid, log=True, bins=50, range=(msid.min(), msid.max()), color='b', 
                   label='All Points (2011:001 - 2012:150, 138 < pitch < 140.5)')
        hist_bad = hist(msid[bad], log=True, bins=50, range=(msid.min(), msid.max()), color='r', 
                   label='Above + Alpha Error > ' + str(angle_err_lim) + ' deg')
        ylim(10**-1, 10**6)
        grid()
        title(name)
        xlabel('Deg F')
        legend()
        x = xlim()
        subplot(2,1,2)
        if all(hist_all[1] == hist_bad[1]):
            ratio = 100 * hist_bad[0].astype(float) / hist_all[0]
            midpts = hist_all[1][:-1] + (hist_all[1][1] - hist_all[1][0]) / 2
            plot(midpts, ratio, 'r-*')
            title('Likelihood of Bad Data (Alpha Error > ' + str(angle_err_lim) + ' deg) \n' + 
                  'vs ' + name)
            xlabel(name + ' (deg F)')
            ylabel('Likelihood of Bad Data (%)')
            ylim(0,100)
            grid()
            tight_layout()
            xlim(x)
        else:  print('Warning:  Histogram bins did not match up.')
        if savefigs==True:  
            savefig(figprefix + '_log.png')

def plot_binary(out, angle_err_lim=8.0, savefigs=False):
    print('Warning:  This plot should only be used when interp=4.1 sec')
    alpha_err = out['alpha'] - out['roll']
    bad = abs(alpha_err) > angle_err_lim
    alpha_raw = out['alpha_raw'].tolist()
    alpha_binary = [bin(raw)[2:].zfill(16) for raw in alpha_raw]
    sun_prs = array([bits[0]=='1' for bits in alpha_binary])
    alpha_array = array([array(list(bits)).astype(int) for bits in alpha_binary])
    # First Figure:  Transitions to Bad Data
    figure()
    just_before_bad = append(~bad[:-1] & bad[1:], False)
    just_went_bad = insert(~bad[:-1] & bad[1:],0,False)
    bits_just_before_bad = alpha_array[just_before_bad]
    bits_just_went_bad = alpha_array[just_went_bad]
    bits_different = bits_just_before_bad != bits_just_went_bad
    bits_flipped = sum(bits_different, 0)
    total_num_events = float(sum(just_went_bad))
    bar(range(15, -1, -1), 100 * bits_flipped / total_num_events, align='center', color='m')
    title('Bit Analysis for AOALPHA \n (AOALPHA = Root MSID for AOALPANG and AOALPSUN)')
    ylabel('% of time bit changed when alpha angle went > ' + str(angle_err_lim) + ' deg')
    xlabel('Bit in AOALPHA \n' + 
           'bit 0 = angle LSB, bit 14 = angle MSB, bit 15 = sun presence')
    grid()
    ylim((0,100))
    if savefigs==True:  
        savefig('bit_analysis_bad.png')
    # Second Figure:  In general
    figure()
    bits_different_gen = alpha_array[1:] != alpha_array[:-1]
    bits_flipped_gen = sum(bits_different_gen, 0)
    total_num_events_gen = float(size(alpha_array, 0) - 1)
    bar(range(15, -1, -1), 100 * bits_flipped_gen / total_num_events_gen, align='center', color='orange')
    title('Bit Analysis for AOALPHA \n (AOALPHA = Root MSID for AOALPANG and AOALPSUN)')
    ylabel('% of time bit changed in general')
    xlabel('Bit in AOALPHA \n' + 
           'bit 0 = angle LSB, bit 14 = angle MSB, bit 15 = sun presence')
    grid()
    ylim((0,100))
    if savefigs==True:  
        savefig('bit_analysis_gen.png')
    
def get_data(start='2005:001', stop=SAFEMODE_2012150, interp=32.8,
             pitch0=100, pitch1=144):
    msids = ('aopssupm', 'aopcadmd', 'aoacaseq', 'pitch', 'roll',
             'aoalpang', 'aobetang', 'aoalpsun', 'aobetsun',
             'tfssbkt1', 'tfssbkt2', 'tcylaft6', 'tpc_fsse', 'aoalpha')
    print 'fetching data'
    x = fetch.MSIDset(msids, start, stop)

    # Resample MSIDset (values and bad flags) onto a common time sampling
    print 'starting interpolate'
    x.interpolate(interp, filter_bad=False)

    # Remove data during times of known bad or anomalous data (works as of
    # Ska.engarchive 0.19.1)
    x.filter_bad_times()

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
                'pitch', 'roll', 'alpha', 'beta',
                'bkt_temp1', 'bkt_temp2', 'tcylaft6', 'fsse_temp',
                'alpha_sun', 'beta_sun', 'spm_act', 'spm_act_bad', 'kalman',
                'alpha_raw')
    dtypes = ('f8',
              'f4', 'f4', 'f4', 'f4',
              'f4', 'f4', 'f4', 'f4',
              'bool', 'bool', 'bool', 'bool', 'bool', 
              'uint16')
    out = np.empty(nvals, dtype=zip(colnames, dtypes))

    out['times'][:] = x['pitch'].times[ok]
    out['pitch'][:] = x['pitch'].vals[ok]
    out['roll'][:] = x['roll'].vals[ok]
    out['alpha'][:] = -x['aoalpang'].vals[ok]
    out['beta'][:] = 90 - x['aobetang'].vals[ok]
    out['bkt_temp1'][:] = x['tfssbkt1'].vals[ok]
    out['bkt_temp2'][:] = x['tfssbkt2'].vals[ok]
    out['tcylaft6'][:] = x['tcylaft6'].vals[ok]
    out['fsse_temp'][:] = x['tpc_fsse'].vals[ok]    
    out['alpha_sun'][:] = x['aoalpsun'].vals[ok] == 'SUN '
    out['beta_sun'][:] = x['aobetsun'].vals[ok] == 'SUN '
    out['spm_act'][:] = x['aopssupm'].vals[ok] == 'ACT '
    out['spm_act_bad'][:] = x['aopssupm'].bads[ok]
    out['kalman'][:] = ((x['aoacaseq'].vals[ok] == 'KALM') &
                        (x['aopcadmd'].vals[ok] == 'NPNT'))
    out['alpha_raw'][:] = x['aoalpha'].vals[ok]                    
    return out