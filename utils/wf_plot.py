from obspy import UTCDateTime
import matplotlib.pyplot as plt
from utils.basic_utils import read_sac_ref_time




def waveform_plot(st,method='dist'):
    if len(st) == 0:
        print("Waveform plot function error, length of Stream is 0!")
    #Inititae parameters
    min_time = st[0].stats.starttime
    max_time = st[0].stats.endtime
    sampling_rate = st[0].stats.sampling_rate
    for tr in st[1:]:
        if tr.stats.starttime < min_time:
            min_time = tr.stats.starttime
        if tr.stats.endtime > max_time:
            max_time = tr.stats.endtime
    
    fig = plt.figure(figsize=(8,10))
    plt.xlim(min_time-min_time,max_time-min_time)
    plt.xlabel("Seconds (s)")
    plt.ylabel("Distance (km)")
    for tr in st:
        origin_time = read_sac_ref_time(tr)
        try:
            dist = tr.stats.sac.dist
        except:
            evla = tr.stats.sac.evla
            evlo = tr.stats.sac.evlo
            stla = tr.stats.sac.stla
            stlo = tr.stats.sac.stlo
            dist = spherical_dist(evlo,evla,stlo,stla)*111
        starttime = tr.stats.starttime
        x_start = starttime - min_time
        tr.data = tr.data*2/(max(tr.data)-min(tr.data))
        plt.plot(np.arange(0,len(tr.data))/sampling_rate+x_start,tr.data+dist,color = 'k',linewidth = 0.5)
        plt.plot([origin_time-min_time,origin_time - min_time],[dist-0.5,dist+0.5],color='k',linewidth = 2)

        try:
            a = tr.stats.sac.a
            rela_a = origin_time - min_time + a
            plt.plot([rela_a,rela_a],[dist-0.5,dist+0.5],color='b',linewidth = 2)
        except:
            pass
        try:
            t0 = tr.stats.sac.t0
            rela_t0 = origin_time - min_time + t0
            plt.plot([rela_t0,rela_t0],[dist-0.5,dist+0.5],color='r',linewidth = 2)
        except:
            pass
