from obspy import read,read_events,UTCDateTime,Stream
from obspy.geodetics.base import gps2dist_azimuth
import glob
from obspy.signal.cross_correlation import correlate, xcorr_max
import multiprocessing
import time
import re
import os

def run_cc(i,event_id_mapper,event_id_1,event_id_2,event_p_st,event_s_st,target_folder):
    link_cc = []
    for tr1 in event_p_st[event_id_1].select(component="*Z"):
        for tr2 in event_p_st[event_id_2].select(component="*Z"):
            if tr1.stats.station == tr2.stats.station:
                sampling_rate = tr1.stats.sampling_rate
                cc=correlate(tr1.data,tr2.data,int(2*sampling_rate))
                shift,value = xcorr_max(cc)
                shift = shift/sampling_rate
                if value > min_cc:
                    link_cc.append([tr1.stats.station,shift,value,'P'])
    for tr1 in event_s_st[event_id_1].select(component="*Z"):
        for tr2 in event_s_st[event_id_2].select(component="*Z"):
            if tr1.stats.station == tr2.stats.station:
                sampling_rate = tr1.stats.sampling_rate
                cc=correlate(tr1.data,tr2.data,int(2*sampling_rate))
                shift,value = xcorr_max(cc)
                shift = shift/sampling_rate
                if value > min_cc:
                    link_cc.append([tr1.stats.station,shift,value,'S'])
    if len(link_cc) < min_link:
        link_cc = []
    if link_cc != []:
        with open(target_folder+"dt.cc",'a') as f:
            f.write("#   %i %i   0 \n" %(event_id_mapper[event_id_1],event_id_mapper[event_id_2]))
            for link in link_cc:
                f.write("%s %6.3f %6.3f %s \n" %(link[0][0:2]+link[0][0:5],link[1],link[2],link[3]))
        f.close()

def run_parallel(i,event_id_1,remain_event_list,event_lat,event_lon,event_pair_id,event_id_mapper,event_p_st,event_s_st,target_folder):
    for event_id_2 in remain_event_list:
            dist,_,_ = gps2dist_azimuth(event_lat[event_id_1],event_lon[event_id_1],event_lat[event_id_2],event_lon[event_id_2])
            if dist/1000 > max_sep:
                continue
            pair_id = event_id_1+event_id_2
            pair_id_mapper = [event_id_1,event_id_2,event_id_mapper[event_id_1],event_id_mapper[event_id_2]]
            event_pair_id[pair_id]=pair_id_mapper
            run_cc(i,event_id_mapper,event_id_1,event_id_2,event_p_st,event_s_st,target_folder)

if __name__ == "__main__":
    """
    This demo script modified from Jeremy's dt_cc_py.py script
    The input file is .pha file that is used for further hypoDD relocation process
    run: python phase_cross_correlation.py 
    """
    cores = 4
    target_folder="./"
    cata = read_events(target_folder+"out.pha")

    cata_subset = cata.filter("latitude > 29.305",'latitude < 29.501','longitude > 104.36','longitude < 104.61')
    print("Events qty to be processed: ",len(cata_subset))
    max_sep = 3 #Maximum seperation distance
    min_link = 5
    min_cc = 0.7
    event_id_mapper = {}
    #Keep the same event id as input out.pha file
    for event in cata:
        event_id=event.resource_id.id
        event_id_mapper[event_id]=int(re.split("/",event_id)[2])

    event_list=[]
    event_p_st={}
    event_s_st={}
    event_count=0
    for event in cata_subset:
        event_count+=1
        print(f"Read in event {event_count}         ",end="\r")
        event_id = event.resource_id.id
        st = Stream()
        p_st = Stream()
        s_st = Stream()
    
        str_etime = str(event.origins[0].time) #Event time 
        #Waveform folder was named based on event time
        w_folder = str_etime[0:4]+str_etime[5:7]+str_etime[8:10]+\
            str_etime[11:13]+str_etime[14:16]+str_etime[17:19]+str_etime[20:22]
        #Waveform folder was placed under day_folder
        day_folder = w_folder[0:8]
        try:  #Read in Z component
            st += read("e_waveform/"+day_folder+"/"+w_folder+"/"+"*Z")
        except: #Without corresponding waveform
            continue
        st.detrend("constant")
        st.detrend("linear")
        st.filter('bandpass',freqmin=5,freqmax=14)
        for pick in event.picks:
            sta = pick.waveform_id.station_code[2:]
            pick_time = pick.time
            pick_type = pick.phase_hint
            if pick_type == "P":
                p_st += st.select(station=sta).slice(pick_time-0.5,pick_time+1)
            elif pick_type == "S":
                s_st += st.select(station=sta).slice(pick_time-1,pick_time+3)
    
        event_list.append(event_id)
        event_p_st[event_id] = p_st
        event_s_st[event_id] = s_st
    print("")

    #initiate dt.cc file
    f = open(target_folder+"dt.cc",'w')
    f.close()

    event_lat={}
    event_lon={}
    for event in cata_subset:
        eve_id = event.resource_id.id
        event_lat[eve_id]=event.origins[0].latitude
        event_lon[eve_id]=event.origins[0].longitude
    
    event_pair_id={}
    event_pair_cc={}

    #Below start to run parallel processing
    pool = multiprocessing.Pool(processes = cores)
    tasks = []
    for i,event_id_1 in enumerate(event_list):
        remain_event_list = event_list[(i+1):]
        tasks.append((i,event_id_1,remain_event_list,event_lat,event_lon,event_pair_id,event_id_mapper,event_p_st,event_s_st,target_folder))
    rs = pool.starmap_async(run_parallel,tasks,chunksize = 1)
    pool.close()
    while(True):
        remaining = rs._number_left
        print("Finished: {0}/{1}    ".format(len(tasks) - remaining,len(tasks)),end = '\r')
        if(rs.ready()):
            break
        time.sleep(0.5)
    pool.join()
    print("")
    print("Done!")
    
