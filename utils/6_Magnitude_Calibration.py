#!/usr/bin/env python
# coding: utf-8

# In[1]:


from obspy import UTCDateTime,Stream,read
from utils.basic_utils import _region_subset,_time_subset,load_hypoDD,load_sum,month_day
from utils.basic_utils import spherical_dist,str_2_time,find_nearest,time_2_str,read_sac_ref_time
import os
from math import sqrt
from math import log
import numpy as np
import matplotlib.pyplot as plt
starttime = UTCDateTime("2020-04-22") 
endtime= UTCDateTime("2020-06-13")


# In[2]:


sum_list = load_sum("INV_DD/inv/out.sum") #The waveform data was cut based on the inverse results
ref_file = '/home/zijinping/Desktop/projects/wy_eq/6_wy_eq/data/Molin_Case/hypoDD.reloc' #reference file
ref_dict,ref_df = load_hypoDD(ref_file,shift_hours=8) 


# In[3]:


filt_dict=_time_subset(ref_dict,starttime,endtime)
tmp_dict = filt_dict.copy()
for eve in tmp_dict:
    if filt_dict[eve][3]<2.5: #Accept events with magnitudes larger than 2.2 
        filt_dict.pop(eve)
#Correlate the SC catalog events with hypoinverse results
for ref_eve in filt_dict:
    ref_time = UTCDateTime(ref_eve)
    min_gap = 10000 # set abonral initiate number
    for eve_id in sum_list:
        eve_time = str_2_time(sum_list[eve_id])
        gap = abs(ref_time - eve_time)
        if gap < min_gap:
            min_gap = gap
            min_eve_id = eve_id
    filt_dict[ref_eve].append(min_eve_id)
print(filt_dict)


# In[4]:


with open("Calib_Reference_Events.txt","w") as f:
    for eve in filt_dict:
        lon = filt_dict[eve][0]
        lat = filt_dict[eve][1]
        dep = filt_dict[eve][2]
        mag = filt_dict[eve][3]
        evid = filt_dict[eve][5]
        evfolder = sum_list[evid]
        f.write(eve+" "+format(lon,'7.3f')+" "+format(lat,'6.3f')+" "+format(dep,'5.2f')+" "+format(mag,'4.2f')+" "               +format(evid,'5d')+" "+sum_list[evid]+"\n")
f.close()


# In[5]:


eq_file = "INV_DD/dd/SS_subset_hypoDD_cc_sum.reloc"
eq_dict, eq_df = load_hypoDD(eq_file)


# In[6]:


'''
if not os.path.exists("caldist2SCseq.txt"):
    f = open("caldist2SCseq.txt",'w')
    f.close()
else:
    os.remove("caldist2SCseq.txt")
    f = open("caldist2SCseq.txt",'w')
    f.close()  
for eq in eq_list:
    eq_lon = eq_list[eq][0]
    eq_lat = eq_list[eq][1]
    min_dist = 50 #set a large value
    for ref_eve in filt_dict:
        ref_lon = filt_dict[ref_eve][0]
        ref_lat = filt_dict[ref_eve][1]
        dist = spherical_dist(eq_lon,eq_lat,ref_lon,ref_lat)*111
        if dist < min_dist:
            min_dist=dist
            min_ref_lon=ref_lon
            min_ref_lat=ref_lat
            min_ref_eve=ref_eve
    with open("caldist2SCseq.txt","a") as f:
        f.write(eq+" "+format(eq_lon,"7.3f")+" "+format(eq_lat,'6.3f')+" "+
                min_ref_eve+" "+format(min_ref_lon,'7.3f')+" "+
                format(min_ref_lat,'6.3f')+" "+format(min_dist,'7.3f')+"\n")
    f.close()
'''


# In[7]:


sta_list = ["GS010","GS012","GS020","GS021","GS022","GS023","GS030",               "GS031","GS033","GS040","GS041","GS043","ML02","ML17",               "ML34","ML50","ML66","ML77","GS050","GS051","GS052","GS053"]


# In[11]:


def load_st(date,folder):
    st = read('e_waveform/'+date+"/"+folder+"/*BH*")
    st.detrend('linear')
    st.detrend('constant')
    st.filter('bandpass',freqmin=5,freqmax=14)
    return st
def trace_mapper(st1,st2): #get stream with the same station
    mapper_1 = Stream()
    mapper_2 = Stream()
    #Get unique station list
    sta_list_1 = []
    pair_sta_list = []
    for tr in st1:
        if tr.stats.station not in sta_list_1:
            sta_list_1.append(tr.stats.station)
    for sta in sta_list_1:
        sta_st2 = st2.select(station = sta)
        if len(sta_st2) != 0:
            pair_sta_list.append(sta)
            mapper_1 += st1.select(station = sta)
            mapper_2 += st2.select(station = sta)
    return pair_sta_list,mapper_1, mapper_2

def get_ratio_zone(st,b,e,eve_id,sum_list):
    
    this_event = sum_list[eve_id]
    this_time = str_2_time(this_event)
    next_event = sum_list[eve_id+1]
    next_time = str_2_time(next_event)
    next_origin = next_time - this_time
    
    a_list = [] #p_arrival time list
    t_list = [] #s_arrival time list
    d_list = [] #distance list
    for tr in st.select(component="*Z"):
        try:
            a=tr.stats.sac.a
        except:
            a=0
        try:
            t0=tr.stats.sac.t0
        except:
            t0=0
        evla = tr.stats.sac.evla
        evlo = tr.stats.sac.evlo
        stla = tr.stats.sac.stla
        stlo = tr.stats.sac.stlo
        dist = spherical_dist(evlo,evla,stlo,stla)
        a_list.append(dist)
        t_list.append(a)
        d_list.append(t0)
    a_max = max(a_list)
    t_max = max(t_list)
    d_max = max(d_list)
    #get time region
    if a_max == 0: #No P picks
        Tb = t_max-d_max/10-1
        Te = t_max+5
    else:
        Tb = a_max-1
        Te = a_max+d_max/10+5
    #Check if excess the length of waveform
    if b > Tb:
        Tb = b
    if e < Te:
        Te = e
    #Don't overlap with next event    
    if Te > next_origin:
        Te = next_origin
    if Tb > Te:
        print("A mistake error in confining the region")
    return Tb,Te

def D3_data(st):
    if len(st)!=3:
        print("Error in get combined ENZ data, not enough component")
    N_trace = st.select(component="*N")[0]
    N_data = N_trace.data
    E_trace = st.select(component="*E")[0]
    E_data = E_trace.data
    Z_trace = st.select(component="*Z")[0]
    Z_data = Z_trace.data
    new_data = []
    for i in range(len(N_data)):
        new_data.append(sqrt(N_data[i]**2+E_data[i]**2+Z_data[i]**2))
    return np.array(new_data)


# In[ ]:


#sampling_rate = 100
i=0
for eq in eq_dict:
    eq_lon,eq_lat,eq_dep,_,eq_id = eq_dict[eq]
    eq_folder=sum_list[eq_id]
    eq_date=eq_folder[0:8]
    eq_st = load_st(eq_date,eq_folder)

    min_dist = 50 #set a large value
    for ref_eve in filt_dict:
        ref_lon = filt_dict[ref_eve][0]
        ref_lat = filt_dict[ref_eve][1]
        dist = spherical_dist(eq_lon,eq_lat,ref_lon,ref_lat)*111
        if dist < min_dist:
            min_dist=dist
            min_ref_lon=ref_lon
            min_ref_lat=ref_lat
            min_ref_eve=ref_eve
    ref_mag = filt_dict[min_ref_eve][3]
    ref_id = filt_dict[min_ref_eve][5]
    ref_folder = sum_list[ref_id]
    ref_date = ref_folder[0:8]
    ref_st = load_st(ref_date,ref_folder)
    pair_sta_list,pair_eq_st, pair_ref_st = trace_mapper(eq_st,ref_st)
    r_start = ref_st[0].stats.sac.b
    r_end = ref_st[0].stats.sac.e
    e_start = eq_st[0].stats.sac.b
    e_end = eq_st[0].stats.sac.e
    
    ref_Tb, ref_Te = get_ratio_zone(ref_st,r_start,r_end,ref_id,sum_list)
    eq_Tb, eq_Te = get_ratio_zone(eq_st,e_start,e_end,eq_id,sum_list)

    ratio_list=np.array([])
    eq_plot_waveform = Stream()
    ref_plot_waveform = Stream()
    for sta in pair_sta_list:
        sta_eq_st = pair_eq_st.select(station = sta)
        eq_starttime = sta_eq_st[0].stats.starttime + 20
        sta_eq_st.trim(starttime = eq_starttime + eq_Tb, endtime = eq_starttime+eq_Te)
        eq_plot_waveform += sta_eq_st.select(component="*Z")
        eq_D3_data = np.array(D3_data(sta_eq_st))
        #sta_eq_st.plot()
        sta_ref_st = pair_ref_st.select(station = sta)
        ref_starttime = sta_ref_st[0].stats.starttime+20
        sta_ref_st.trim(starttime=ref_starttime+ref_Tb,endtime =ref_starttime+ref_Te)
        ref_plot_waveform += sta_ref_st.select(component="*Z")

        #sta_ref_st.plot()
        ref_D3_data = np.array(D3_data(sta_ref_st))
        ratio = max(eq_D3_data)/max(ref_D3_data)
        ratio_list = np.append(ratio_list,ratio)
        #print("Data length is ",len(ref_D3_data))
    #plt = waveform_plot(eq_plot_waveform)
    #plt.show()
    #plt2 = waveform_plot(ref_plot_waveform)
    #plt2.show()
    mean_ratio = np.mean(ratio_list)
    calib_mag = ref_mag + log(mean_ratio,10)
    print(f"  Processed {i}/{len(eq_dict)}    ",end = '\r')
    i = i+1
    eq_dict[eq][3]=calib_mag
    


# In[ ]:


print(eq)


# In[ ]:


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
    return plt
    


# In[ ]:


def read_sac_ref_time(tr):
    nzyear = tr.stats.sac.nzyear
    nzjday = tr.stats.sac.nzjday
    nzhour = tr.stats.sac.nzhour
    nzmin = tr.stats.sac.nzmin
    nzsec = tr.stats.sac.nzsec
    nzmsec = tr.stats.sac.nzmsec*0.001
    year,month,day = month_day(nzyear,nzjday)
    sac_ref_time = UTCDateTime(year,month,day,nzhour,nzmin,nzsec)+nzmsec
    return sac_ref_time


# In[ ]:




