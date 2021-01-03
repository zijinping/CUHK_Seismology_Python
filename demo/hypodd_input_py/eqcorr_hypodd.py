#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create input event.dat, phase pair (dt.ct) and waveform cross-correlation
(dt.cc) files for hypodd 

    Processes:
        0 Subset events
        1 ph2dt
        2.1 Generate event streams (Optional for prebuilt already)
        2.2 Cross-correlation

    Configure the script after __main__

    Input: 
        Catalog file (readable via obspy.read_events)
            - Formats: .pha, .xml, etc...
        Continuous Waveform
            - Store in daily base (year/doy/*.[sac/mseed])

    Output: (In the same directory of the script)
        Event.dat   - Subsetted event with assigned id
        dt.ct       - phase-pair
        dt.cc       - cross-correlation pair

        [ For intermediate step ]
        event_stream.p      - temporal storage of event wf
        event_stream_p.p    - temporal storage of P arrival wf
        event_stream_s.p    - temporal storage of S arrival wf
        event_id_mapper.p   - temporal storage of hypoDD id to catalog id mapper
        hypodd_subset.p     - temporal storage of subseted obspy.catalog object


Created on Fri Nov  6 11:57:19 HKT 2020

@author: jw
"""
from obspy import read, read_events, Stream, UTCDateTime as UTC
from eqcorrscan.utils.catalog_to_dd import write_event,\
                                        write_catalog, \
                                        read_phase,\
                                        write_correlations
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.cross_correlation import correlate, xcorr_max
import os, glob, pickle
from tqdm import tqdm


def _subset_event(event, lat, lon, verbose=True):
    _event_subset = event.filter(
                        "latitude > " +str(lat[0]),
                        "latitude < "+str(lat[1]),
                        "longitude > "+str(lon[0]),
                        "longitude < "+str(lon[1]))

    if verbose:
        print("Spatial subset event: %i > %i" % (len(event), len(_event_subset)))
    return _event_subset


def _subset_event_time(event, starttime, endtime, verbose=True):
    _event_subset = event.filter(
                        "time >= " + str(starttime),
                        "time <= " + str(endtime))
    if verbose:
        print("Time subset event: %i > %i" % (len(event), len(_event_subset)))
    return _event_subset


def get_event_stream(event_catalog, starttime, endtime, waveform_dir, 
    time_before_arrival=2, time_after_arrival=3,
    verbose=False, dump=False, filter_kwargs=None):
    """ Generate dictionary of event waveforms based on event.resource_id.id
    Return:
        Dictionary with keys of event.resoruce_id of 
            1) event waveform
            2) P arrivals
            3) S arrivals
    """
    print("CREATE EVENT WAVEFORM")
    event_stream = dict()
    event_stream_p = dict()
    event_stream_s = dict()

    pbar = tqdm(range(0, int(endtime - starttime), int(24*60*60)))

    for tt in pbar:
        _sttime = starttime + tt
        _edtime = _sttime + 24*60*60
        year = str(_sttime.year)
        doy = str(_sttime.julday).zfill(3)
        month = str(_sttime.month)
        day = str(_sttime.day)


        day_hypodd_subset = _subset_event_time(event_catalog, 
                                                _sttime, _edtime, verbose=False)

        pbar.set_description("Proceesing %s-%s, Events: %i" % (year, doy, len(day_hypodd_subset)))

        if len(day_hypodd_subset)==0:
            continue

        st = Stream()
        for _tr_file in glob.glob(os.path.join(waveform_dir, year, doy, "*seed")):
            try:
                st += read(_tr_file)
            except Exception as ex:
                print(ex)

        st_starttime = st[0].stats.starttime
        st_endtime = st[-1].stats.endtime
        st.detrend('linear')
        if filter_kwargs is not None:
            st.filter(**filter_kwargs)
        
        for event in day_hypodd_subset:
            event_id = event.resource_id.id
            origin_time = event.origins[0].time
            _event_stream = Stream()
            _event_stream_P = Stream()
            _event_stream_S = Stream()
            
            # Faster marching (skip event without trace)
            #if event.origins[0].time > st_endtime or event.origins[0].time < st_starttime:
            #    print("event exceed waveform")
            #    continue
                
            for pick in event.picks:
                pick.waveform_id.network_code = 'SC'
                sta = pick.waveform_id.station_code[-3:]         
                pick_time = pick.time
                _event_stream += st.select(station=sta).slice(origin_time, 
                                                                    origin_time+60)
                if pick.phase_hint == 'P':
                    _event_stream_P += st.select(station=sta).slice(pick_time-time_before_arrival,
                                                                    pick_time+time_after_arrival)
                elif pick.phase_hint == 'S':
                    _event_stream_S += st.select(station=sta).slice(pick_time-time_before_arrival,
                                                                    pick_time+time_after_arrival)
            if verbose:
                print(_event_stream_P)
            _event_stream.merge()
            event_stream[event_id] = _event_stream 
            event_stream_p[event_id] = _event_stream_P
            event_stream_s[event_id] = _event_stream_S 
    
    if dump:
        ### temporal store event_id_mapper and event_stream to save time
        pickle.dump(hypodd_subset, open("hypodd_subset.p", "wb"))
        pickle.dump(event_id_mapper, open("event_id_mapper.p", "wb"))
        pickle.dump(event_stream, open("event_stream.p", "wb"))
        pickle.dump(event_stream_p, open("event_stream_p.p", "wb"))
        pickle.dump(event_stream_s, open("event_stream_s.p", "wb"))


    return event_stream, event_stream_p, event_stream_s


def _gen_event_lat_lon(catalog):
    """ Create dictonary of event latitude and longitude using 
    event.resource_id.id as key"""
    event_lat = dict()
    event_lon = dict()
    for event in catalog.events:
        eve_id = event.resource_id.id
        event_lat[eve_id] = event.origins[0].latitude
        event_lon[eve_id] = event.origins[0].longitude
    
    return event_lat, event_lon


def _correlation_stream(st1, st2, shift=2):
    """ 
    para:   shift - maxshift in seconds
    Return: shifttime, cc value
    """
    sampling_rate = st1.stats.sampling_rate
    cc = correlate(st1.data, st2.data, int(sampling_rate*shift))
    shift, value = xcorr_max(cc)
    return shift/sampling_rate, value


def correlate_event_stream(catalog, event_stream_p, event_stream_s, 
                            event_id_mapper, max_sep=20, min_link=5, min_cc=0.5,
                            max_shift=2,verbose=False, vverbose=False):
    event_lat, event_lon = _gen_event_lat_lon(catalog)
    event_list = list(event_stream.keys())
    event_pair_id = dict()
    event_pair_cc = dict()

    pbar = tqdm(total = len(event_list)*(len(event_list)-1)/2)
    for i, event_id_1 in enumerate(event_list):
        for event_id_2 in event_list[(i+1):]:
            # Check event distance
            dist, _, _ = gps2dist_azimuth(event_lat[event_id_1],
                                    event_lon[event_id_1],
                                    event_lat[event_id_2],
                                    event_lon[event_id_2])
            if dist/1000 > max_sep: 
                continue
             
            # waveform correlation 
            pair_id = event_id_1+event_id_2
            pair_id_mapper = [event_id_1, event_id_2, event_id_mapper[event_id_1], event_id_mapper[event_id_2]]
            event_pair_id[pair_id] = pair_id_mapper

            # Progress bar
            pbar.update(1)
            pbar.set_description("Processing %s " % pair_id)

            link_cc = []
            for tr1 in event_stream_p[event_id_1].select(component="*Z"):
                 for tr2 in event_stream_p[event_id_2].select(component="*Z"):
                    if tr1.stats.station == tr2.stats.station:
                        shift, value = _correlation_stream(tr1, tr2, shift=max_shift)
                        if value < min_cc:
                            continue
                        if vverbose:
                            print([tr1.stats.station, shift, value, 'P'])
                        link_cc.append([tr1.stats.station, shift, value, 'P'])

            for tr1 in event_stream_s[event_id_1].select(component="*Z"):
                 for tr2 in event_stream_s[event_id_2].select(component="*Z"):
                    if tr1.stats.station == tr2.stats.station:
                        shift, value = _correlation_stream(tr1, tr2, shift=max_shift)
                        if value < min_cc:
                            continue
                        link_cc.append([tr1.stats.station, shift, value, 'S'])

            if len(link_cc) < min_link:
                link_cc = []
            if verbose:
                print(pair_id)
                print(link_cc)
            event_pair_cc[pair_id] = link_cc

    pbar.close()
 
    return event_pair_id, event_pair_cc


def _parse_to_dtcc(event_pair_id, event_pair_cc, f_dtcc='dt.cc'):
    with open(f_dtcc, 'w+') as dt_cc_f:
        event_id_list = list(event_pair_id.keys())
        for eve_id in event_id_list:
            [_, _, id1, id2] = event_pair_id[eve_id] 
            link_cc = event_pair_cc[eve_id]
            if len(link_cc)==0:
                continue
            dt_cc_f.write("#   %i  %i   0 \n" % (id1, id2))
            for link in link_cc:
                dt_cc_f.write("%s %.3f %.3f %s \n" % tuple(link))
    print("Created dt.cc")

    return None
     

if __name__ == "__main__":
    ######## Configurations ########
    hypodd_pha = read_events('/NAS2/Sichuan_data/phasenet/event_waveform/HYPODD2/weiyuan_subset.pha')


    ### 0 Subset events with latlon and time ranges
    lat = [29.48, 29.73]
    lon = [104.65, 104.95]

    starttime = UTC(2018,9,1)
    endtime = UTC(2019,3,1)

    ### 1 Generate hypodd output (ph2dt)
    max_sep = 20            # Maximum seperations between event paris in km
    min_link = 8            # minimum links for an event to be paired,
                            # e.g. min. no. of picks with same sta& chan
                            # shared between two paired events


    ### 2 Cross-correlation Configurations
    load_prebuilt = False           # Select true if event waveforms are stored
                                    #   event_stream.p, ...
    
    ## Continuous wf directory
    # Waveform store in waveform_dir/year/doy
    waveform_dir = '/NAS2/Sichuan_data/event_waveform_PN'

    ## Waveform CC configurations
    before_arrival = 0.5
    after_arrival = 1.5

    filter_kwargs = {'type':'bandpass',
                     'freqmin': 2,
                     'freqmax': 15}

    print('Filter settings')
    
    ## Cross-correlation pair configurations
    max_sep = 20
    min_link_cc = 4
    min_cc = 0.5
    max_shift=1.0



    ######## End of Configurations ########
    

    ## 0 Subset events
    hypodd_subset = _subset_event(hypodd_pha, lat, lon)
    hypodd_subset = _subset_event_time(hypodd_subset, starttime, endtime)


    ## 1 ph2dt
    event_id_mapper = write_event(hypodd_subset)
    event_id_mapper = write_catalog(hypodd_subset, \
                                event_id_mapper=event_id_mapper)


    ## 2.1 Generate event_streams with corresponding event_id_mapper
    # read stream and extract event waveform corresponding to the 
    # event_id_mapper

    print(filter_kwargs)

    # Generate event stream
    if load_prebuilt == False:
        event_stream, event_stream_p, event_stream_s = get_event_stream(hypodd_subset,
                                                starttime, endtime, waveform_dir,
                                                time_before_arrival=before_arrival,
                                                time_after_arrival=after_arrival,
                                                filter_kwargs=filter_kwargs, dump=True)
    

    ### load previous event_id_mapper, event_stream
    if load_prebuilt == True:
        print('Load pickle dump')
        with open("hypodd_subset.p", "rb") as f_event_list:
            hypodd_subset = pickle.load(f_event_list)
        with open("event_id_mapper.p", "rb") as f_event_id:
            event_id_mapper = pickle.load(f_event_id)
        with open("event_stream.p", "rb") as f_event_st:
            event_stream = pickle.load(f_event_st)
        with open("event_stream_p.p", "rb") as f_event_st:
            event_stream_p = pickle.load(f_event_st)
        with open("event_stream_s.p", "rb") as f_event_st:
            event_stream_s = pickle.load(f_event_st)


    ### Compute 2.2 Conduct CC over event streams
    # Conduct CC
    event_pair_id, event_pair_cc = correlate_event_stream(hypodd_subset, 
                                                        event_stream_p, 
                                                        event_stream_s,
                                                        event_id_mapper,
                                                        max_sep=max_sep,
                                                        min_cc=min_cc,
                                                        min_link=min_link_cc,
                                                        max_shift=max_shift,
                                                        verbose=False)  
    # Create file (./dt.cc)
    _parse_to_dtcc(event_pair_id, event_pair_cc, f_dtcc='dt.cc')

