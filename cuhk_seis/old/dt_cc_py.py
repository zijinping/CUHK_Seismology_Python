from obspy import read, read_events
from obspy.signal.cross_correlation import correlate, xcorr_max
import numpy as np
from eqcorrscan.utils.catalog_to_dd import write_event,\
                                        write_catalog, read_phase

def correlate(st1, st2, st1_time, st2_time, conf):
    """
    Correlate two trace with their pick time
    returns 
        shifted in seconds, cc_value
    """
    sample_rate = st1[0].stats.sample_rate
    shift = CC_config["max_shift"]*sample_rate
    st1_sliced = st1.slice(
                    st1_time + CC_config['cc_time_before'],
                    st1_time + CC_config["cc_time_after"])
    st2_sliced = st2.slice(
                    st2_time + CC_config['cc_time_before'],
                    st2_time + CC_config["cc_time_after"])

    # correlate
    _cc = correlate(st1_sliced[0], st2_sliced[1], int(shift))
    # Squared _cc
    _cc_shift, _cc_value = xcorr_max(_cc*_cc)
    return _cc_shift/sample_rate, _cc_value
    

def parse_to_dtcc(correlation_list, dt_cc_file):
    """
    dt.cc format:
    Event: #, ID1, ID2, OTC
    CC Pair: sta, dt, wight, pha 
    """
    return None

def _subset_event(event, lat, lon):
    _event_subset = event.filter(
                        "latitude > " +str(lat[0]), 
                        "latitude < "+str(lat[1]), 
                        "longitude > "+str(lon[0]),
                        "longitude < "+str(lon[1]))

    print("Event size: %i > %i" % (len(event), len(_event_subset)))
    return _event_subset

def _subset_event_time(event, starttime, endtime):
    _event_subset = event.filter(
                        "time >= " + str(starttime), 
                        "time <= " + str(endtime))
    print("Event size: %i > %i" % (len(event), len(_event_subset)))
    return _event_subset

if __name__ == "__main__":# 
    hypodd = read_events('../../catalog/pengcheng/weiyuan.pha_201809_201903')

    
    ### Subset events to mulin fault and time ranges
    lat = [29.4, 29.57]
    lon = [104.4, 104.6]

    hypodd_subset = _subset_event(hypodd, lat, lon)

    waveform_dir = ''


    ### Generate hypodd output (ph2dt)
    max_sep = 50            # Maximum seperations between event paris in km
    min_link = 8            # minimum links for an event to be paired, 
                            # e.g. min. no. of picks with same sta& chan
                            # shared between two paired events

    # ph2dt
    event_id_mapper = write_events(hypodd_subset)
    event_id_mapper = write_catalog(hypodd_subseti, \
                                event_id_mapper=event_id_mapper)
    
