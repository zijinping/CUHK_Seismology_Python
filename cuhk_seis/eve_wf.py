import os
import glob
import obspy
from obspy import UTCDateTime,read
from obspy.io.sac.sactrace import SACTrace
import shutil
import re
import matplotlib.pyplot as plt
from cuhk_seis.basic_utils import read_sac_ref_time,spherical_dist,get_st
import numpy as np
from tqdm import tqdm
import pickle
import logging

def load_sta(sta_file):
    sta_dict = {}
    with open(sta_file,'r') as f:
        for line in f:
            line = line.rstrip()
            net,sta,_lon,_lat,_ele = re.split(" +",line)[:5]
            sta_dict[sta] = [float(_lon),float(_lat),int(_ele),net]
    f.close()
    return sta_dict

def load_y2000(y2000_file):
    phs_cont = []
    with open(y2000_file,"r") as f1:
        for line in f1:
            phs_cont.append(line.rstrip())
    f1.close()
    phs_dict = {}
    event_count = 0

    print(">>> Loading phases ... ")
    for line in tqdm(phs_cont):
        f_para = line[0:2]     # first two characters as first parameter(f_para)
        if re.match("\d+",f_para):    # event line
            event_count += 1
            _yr=line[0:4];_mo=line[4:6];_day=line[6:8]
            _hr=line[8:10];_minute=line[10:12];
            yr = int(_yr); mo = int(_mo); day=int(_day); hr=int(_hr);minute=int(_minute);
            _seconds=line[12:14]+"."+line[14:16]
            evla=float(line[16:18])+(float(line[19:21])+float(line[21:23])*0.01)/60
            evlo=float(line[23:26])+(float(line[27:29])+float(line[29:31])*0.01)/60
            evdp=float(line[32:36])/100; evid = int(line[136:146])
           # _,no,year,month,day,o_time,ab_sec,res,evla,evlo,evdp,mag,mag_res,np,ns,nt,sta_gap=re.split(" +",line)
            #e_hr,e_min,e_seconds=re.split(":",o_time)
            e_secs = float(_seconds)
            e_time = UTCDateTime(yr,mo,day,hr,minute,0)+e_secs

            str_time = e_time.strftime('%Y%m%d%H%M%S%f')
            str_time = str_time[:16]
            phs_dict[str_time] = {}
            phs_dict[str_time]["eve_loc"] = [evlo,evla,evdp]
            phs_dict[str_time]["phase"] = []
            phs_dict[str_time]["evid"] = evid

        elif re.match("[A-Z]+",f_para) and f_para != "  ": # phase line
            net = line[5:7]
            sta = re.split(" +",line[0:5])[0]
            year = int(line[17:21])
            month = int(line[21:23])
            day = int(line[23:25])
            hour = int(line[25:27])
            minute = int(line[27:29])
            if line[14]==" ":
                #if sec or msec is 0, it will be "  " in out.arc file
                _sec = line[41:44]; _sec_m = line[44:46]
                if _sec == "   ":
                    _sec = "000"
                if _sec_m == "  ":
                    _sec_m = "00"
                p_type="S"
                phs_time = UTCDateTime(year,month,day,hour,minute,0)+\
                                   (int(float(_sec))+int(_sec_m)*0.01)
                phs_dict[str_time]["phase"].append([net,sta,p_type,phs_time-e_time])
            
            else:
                p_type="P"
                #if sec or msec is 0, it will be "  " in out.arc file
                _sec = line[29:32]; _sec_m = line[32:34]
                if _sec == "   ":
                    _sec = "000"
                if _sec_m == "  ":
                    _sec_m = "00"
                phs_time = UTCDateTime(year,month,day,hour,minute,0)+\
                                   (int(float(_sec))+int(_sec_m)*0.01)
                phs_dict[str_time]["phase"].append([net,sta,p_type,phs_time-e_time])
            
    out_name = y2000_file+".pkl"
    out_file = open(out_name,'wb')
    pickle.dump(phs_dict,out_file)
    out_file.close()
 
    return phs_dict

def cut_eve_wf(phs_dict,sta_dict,src_root,tar_root,tb,te):
    eve_list = phs_dict.keys()
    print(">>> cut event waveform ... ")
    for eve in tqdm(eve_list):
        str_day = eve[:8]
        #create day folder
        if os.path.exists(os.path.join(tar_root,str_day)):
            pass
        else:
            os.makedirs(os.path.join(tar_root,str_day))
        if os.path.exists(os.path.join(tar_root,str_day,eve)):
            shutil.rmtree(os.path.join(tar_root,str_day,eve))
        os.makedirs(os.path.join(tar_root,str_day,eve))
        e_time = UTCDateTime.strptime(eve,"%Y%m%d%H%M%S%f")
        starttime = e_time + tb
        endtime = e_time + te
        phases = phs_dict[eve]["phase"]
        processed_sta = []
        for [net,sta,pha,rel_time] in phases:
            src_folder=os.path.join(src_root,sta)
            if sta not in processed_sta:
                processed_sta.append(sta)
                st = get_st(net,sta,starttime,endtime,src_folder)
                if len(st)==0:
                    continue
                for trace in st:
                    chn = trace.stats.channel
                    out_file = os.path.join(tar_root,str_day,eve,f"{sta}.{chn}")
                    trace.write(out_file,format="SAC")
                    sac = SACTrace.read(out_file,headonly=True)
                    sac.lovrok = 1
                    sac.evlo = phs_dict[eve]["eve_loc"][0]
                    sac.evla = phs_dict[eve]["eve_loc"][1]
                    sac.evdp = phs_dict[eve]["eve_loc"][2]
                    sac.stlo = sta_dict[sta][0]
                    sac.stla = sta_dict[sta][1]
                    sac.stel = sta_dict[sta][2]
                    sac.nzyear = e_time.year
                    sac.nzjday = e_time.julday
                    sac.nzhour = e_time.hour
                    sac.nzmin = e_time.minute
                    sac.nzsec = e_time.second
                    sac.nzmsec = e_time.microsecond/1000
                    sac.b = sac.b+tb
                    sac.o = 0
                    sac.write(out_file)

def write_eve_wf(phs_dict,tar_root,depth=2):
    '''
    depth = 2 means tar_root/str_day/eve/sta.*
    depth = 1 means tar_root/eve/sta.*
    '''
    eve_list = phs_dict.keys()
    print(">>> write arrival time information ...")
    for eve in tqdm(eve_list):
        str_day = eve[:8]
        phases = phs_dict[eve]["phase"]
        for [net,sta,pha,rel_time] in phases:
            if depth == 2:
                sac_files = glob.glob(os.path.join(tar_root,str_day,eve,sta+"*"))
            if depth == 1:
                sac_files = glob.glob(os.path.join(tar_root,eve,sta+"*"))

            if len(sac_files) == 0:
                continue
            if pha == "P":
                for sac_file in sac_files:
                    sac = SACTrace.read(sac_file,headonly=True)
                    sac.a = rel_time
                    logging.info(f">>> write in {eve} {sta} with P arrival {sac.a}")
                    sac.write(sac_file)
            if pha == "S":
                for sac_file in sac_files:
                    sac = SACTrace.read(sac_file,headonly=True)
                    sac.t0 = rel_time
                    logging.info(f">>> write in {eve} {sta} with S arrival {sac.t0}")
                    sac.write(sac_file)
            
def gen_eve_wf(y2000_file,sta_file,src_root="day_data",tar_root="eve_wf",tb = -20,te = 60):
    """
    cut event waveform from continous waveform based on phase file generated by 
    HYPOINVERSE

    Parameters:
    y2000_file: phase file, which is the output from HYPOINVERSE
    sta_file: station file in format: net sta lon lat ele
    src_root: source path for continous waveform
    tar_root: target path of eve based waveform
    tb: waveform begin time in relative to the event occurence time
    te: waveform end time in relative to the event ouccurence time
    """
    if not os.path.exists(tar_root):
        os.mkdir(tar_root)
    base_path = os.path.abspath("./")
    sta_dict = load_sta(sta_file)
    phs_dict = load_y2000(y2000_file)

    cut_eve_wf(phs_dict,sta_dict,src_root,tar_root,tb,te)
    write_eve_wf(phs_dict,tar_root)


def wf_dist_plot(st,length=10,color=None,label_sta=True,out_format="PNG",scaling_factor=2):
    '''    
    Description
    ------------
    Plot event waveform by distance. The start time is event origin time.
    Data should be in sac format. The output is a saved file with the title
    of reference time.

    Parameters
    -----------
                st: obspy Stream object
            length: The time window is defined by length in seconds.
             color: The usage of color is the same as matplotlib.pyplot. 
                    Using default color if not defined.
         label_sta: Whether label station name on the plot.
        out_format: "png","jpg","pdf"..., The same as matplotlib.pyplot.savefig
    scaling_factor: The waveform are normalized, increase scaling_facotr to
                    make the waveform plot more obvious

    Below data information needed:
    |   P arrival: tr.stats.sac.a
    |   S arrival: tr.stats.sac.t0
    |        evla: tr.stats.sac.evla
    |        evlo: tr.stats.sac.evlo
    |        stla: tr.stats.sac.stla
    |        stlo: tr.stats.sac.stlo
    '''
    st.detrend("linear")
    st.detrend("constant")
    try:
        e_mag = st[0].stats.sac.mag
    except:
        e_mag = -9
    starttime=st[0].stats.starttime
    endtime =st[0].stats.endtime
    #Reference time shoule be the same for all traces.
    sac_ref_time = read_sac_ref_time(st[0]) #In UTCDateTime format
    o_value = st[0].stats.sac.o 
    event_time = sac_ref_time + o_value #event origin time
    st.trim(starttime = event_time, endtime = event_time + length)
    if len(st) == 0:
        print("Error: Nothing to plot!")
    #Inititae parameters
    for tr in st:
        evla = tr.stats.sac.evla
        evlo = tr.stats.sac.evlo
        stla = tr.stats.sac.stla
        stlo = tr.stats.sac.stlo
        #It is recommend to set tr.stats.distance in meters by osbpy guideline
        tr.stats.distance = spherical_dist(evlo,evla,stlo,stla)*111*1000
    #Get mininumtime, maximum time, and max distance
    min_time = st[0].stats.starttime
    max_time = st[0].stats.endtime
    max_dist = st[0].stats.distance/1000 #in km
    for tr in st[1:]:
        if tr.stats.starttime < min_time:
            min_time = tr.stats.starttime
        if tr.stats.endtime > max_time:
            max_time = tr.stats.endtime
        if tr.stats.distance/1000 > max_dist:
            max_dist = tr.stats.distance/1000
    sampling_rate = st[0].stats.sampling_rate

    #Initiate plot parameters
    plt.figure(figsize = (8,10))
    plt.xlim(0,max_time-min_time)
    plt.ylim(0,max_dist+3)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel("Time (s)",fontdict={'size':16})
    plt.ylabel("Distance (km)",fontdict={'size':16})
    #Plot trace by trace
    for tr in st:
        sta = tr.stats.station
        tr_ref_time = read_sac_ref_time(tr) 
        tr_o_value = tr.stats.sac.o
        event_time = tr_ref_time + tr_o_value
        x_start = event_time - min_time
        dist = tr.stats.distance/1000
        #Normalize the event
        disp_data = tr.data/(max(tr.data) - min(tr.data))
        disp_data = disp_data*scaling_factor
        if color == None:
            plt.plot(np.arange(0,len(tr.data))/sampling_rate+x_start,
                    disp_data+dist,
                    linewidth = 0.5)
        else:
            plt.plot(np.arange(0,len(tr.data))/sampling_rate+x_start,
                    disp_data+dist,
                    color = color,
                    linewidth = 0.5)
        if label_sta:
            plt.text(0.1,dist+0.2,sta,fontsize=12)
        #Plot P arrival if available
        try:
            a = tr.stats.sac.a
            rela_a = tr_ref_time + a - min_time
            gap = 0.5*max_dist/25
            plt.plot([rela_a,rela_a],[dist-gap,dist+gap],color='b',linewidth=2)
        except:
            pass
        #Plot S arrival if available
        try:
            t0 = tr.stats.sac.t0
            rela_t0 = tr_ref_time + t0 - min_time
            gap = 0.5*max_dist/25
            plt.plot([rela_t0,rela_t0],[dist-gap,dist+gap],color='r',linewidth=2)
        except:
            pass
        if e_mag != -9:
            plt.title(str(tr_ref_time)+f"_M{e_mag}",fontdict={'size':18})
        else:
            plt.title(str(tr_ref_time),fontdict={'size':18})
    plt.savefig(f"{sac_ref_time.strftime('%Y%m%d%H%M%S')}.png",format=out_format)
    plt.close()            
