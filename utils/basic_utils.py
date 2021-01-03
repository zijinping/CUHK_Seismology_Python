import os
import numpy as np
from math import radians,cos,acos,sin,asin,sqrt,ceil,pi
import obspy
from obspy import Stream
import glob
import re
from obspy import UTCDateTime
import pandas as pd

root_path="/NAS1/Sichuan_data/continous_waveform_sac/"

def pha_subset(pha_file,lon_min,lon_max,lat_min,lat_max,filt=8):
    out_file = pha_file+".st"
    f = open(out_file,"w")
    f.close()
    pha_content = []
    with open(pha_file,"r") as f:
        for line in f:
            pha_content.append(line.rstrip())
    f.close()
    i = 0
    j = 0
    record_list=[]
    for line in pha_content:
        if line[0]=="#":
            if i>0 and len(record_list) > (filt+1):
                j=j+1
                with open(out_file,"a") as f:
                    for record in record_list:
                        f.write(record+"\n")
                f.close()
                record_list = []
                record_list.append(line)
            else:
                record_list = []
                record_list.append(line)
            i=i+1
            lat = re.split(" +",line)[7]
            lon = re.split(" +",line)[8]
            if float(lat)>lat_min and float(lat)<lat_max and float(lon)>lon_min and float(lon)<lon_max:
                region_pass = True
            else:
                region_pass = False
        else:
            if region_pass:
                record_list.append(line)
    if i>0 and len(record_list) > (filt+1):
        j=j+1
        with open(out_file,"a") as f:
            for record in record_list:
                f.write(record+"\n")
        f.close()
    print("Event before filtering",i)
    print("Events qty after filtering",j)
     

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


def spherical_dist(lon_1,lat_1,lon_2,lat_2):
    """
    test
    """
    lon_1,lat_1,lon_2,lat_2 = map(radians,[lon_1,lat_1,lon_2,lat_2])
    a=acos(sin(lat_1)*sin(lat_2)+cos(lat_1)*cos(lat_2)*cos(lon_2-lon_1))
    return a*180/pi


def get_st(net,sta,starttime,endtime,f_folder,time_zone_effect=0):
    '''
    Parameters list: net,sta,chn,starttime,endtime,f_path,time_zone_effect=0
    '''
    starttime = starttime-time_zone_effect*60*60
    endtime = endtime-time_zone_effect*60*60
    
    inc_list=[]

    for file in os.listdir(f_folder):
        file_path = os.path.join(f_folder,file)
        try:
            st = obspy.read(file,headonly=True)
        except:
            continue
        t1,t2 = st[0].stats.starttime,st[0].stats.endtime
        if t2 < starttime or t1 > endtime:
            continue
        else:
            inc_list.append(file_path)

    #Read in data
    st = Stream()
    for path in inc_list:
        st += obspy.read(path)
    if len(st) == 0:
        pass
    else:
        st.trim(starttime,endtime)
    return st

def julday(year,month,day):
    ref_time=UTCDateTime(year,1,1)
    tar_time=UTCDateTime(year,month,day)
    julday=(tar_time-ref_time)/(24*60*60)+1
    return int(julday)

def month_day(year,julday):
    #check if year is leap year
    leap=False
    if year%100==0:
        if year%400==0:
            leap=True
    else:
        if year%4==0:
            leap=True
    normal_list=[0,31,59,90,120,151,181,212,243,273,304,334,365]
    leap_list=[0,31,60,91,121,152,182,213,244,274,305,335,366]
    if leap:
        i=0
        while leap_list[i]<julday:
            i=i+1
        month=i
        day=julday-leap_list[i-1]
        return year,month,day
    else:
        i=0
        while normal_list[i]<julday:
            i=i+1
        month=i
        day=julday-normal_list[i-1]
        return year,month,day

def find_nearest(array,value):
    '''
    find the nearest array value and index
    '''
    if type(array) != np.ndarray:
        array=np.array(array)
    idx=np.abs(array-value).argmin()
    diff=array[idx]-value
    return idx,diff

def str2time(name):
    year = name[0:4]
    mo = name[4:6]
    dy = name[6:8]
    hr = name[8:10]
    min = name[10:12]
    sec = name[12:]
    time = UTCDateTime(year+"-"+mo+"-"+dy+"T"+hr+":"+min+":"+sec[0:2]+"."+sec[2:])
    return time

def time2str(time):
    year = time.year
    month = time.month
    day = time.day
    hour = time.hour
    minute = time.minute
    second = time.second
    m_second = time.microsecond
    string = str(year)+str(month).zfill(2)+str(day).zfill(2)+str(hour).zfill(2)+\
            str(minute).zfill(2)+str(second).zfill(2)+str(m_second).zfill(6)
    return string

def load_sum(sum_file):
    '''
    out.sum file after hypoinverse run
    '''
    sum_list = {}
    with open(sum_file,'r') as f:
        for line in f:
            eve_id=int(line[136:146])
            eve_folder = line[0:16]
            sum_list[eve_id] = eve_folder
    return sum_list

def load_hypoDD(reloc_file="hypoDD.reloc",shift_hours=8):
    """
    return eve_list, df
    """
    eve_list={}
    columns = ["ID","LAT","LON","DEPTH","X","Y","Z","EX","EY","EZ",\
           "YR","MO","DY","HR","MI","SC","MAG",\
           "NCCP","NCCS","NCTP","NCTS","RCC","RCT","CID"]
    with open(reloc_file,"r") as f:
        for line in f:
            data=re.split(" +",line.rstrip())[1:]
            try:
                data_arr = np.vstack((data_arr,data))
            except:
                data_arr = np.array(data)
            eve_id = data[0]
            eve_lat = data[1]
            eve_lon = data[2]
            eve_dep = data[3]
            eve_time = UTCDateTime(int(data[10]),int(data[11]),int(data[12]),int(data[13]),int(data[14]),0) +float(data[15])- shift_hours*60*60
            eve_time_str=str(eve_time)[0:-1]
            eve_mag =data[16]
            eve_list[eve_time_str]=[float(eve_lon),float(eve_lat),float(eve_dep),float(eve_mag),int(eve_id)]
    f.close()
    df = pd.DataFrame(data=data_arr,columns=columns)
    return eve_list,df

def hypoDD_mag_mapper(reloc_file,out_sum):
    """
    get magnitude from out.sum file after hyperinverse
    """
    #get the magnitude dictionary
    event_mag_list = {}
    with open(out_sum,"r") as f_obj:
        for line in f_obj:
            event_id = int(line[136:146])
            event_mag = int(line[123:126])*0.01
            event_mag_list[event_id]=event_mag
    f_obj.close()
    #add in the magnitude
    new_dd = []
    with open(reloc_file,"r") as f_obj:
        for line in f_obj:
            dd_event_id = int(line[0:9])
            dd_event_mag = event_mag_list[dd_event_id]
            new_line=line[:128]+format(dd_event_mag,'4.2f')+line[132:]
            new_dd.append(new_line)
    f_obj.close()
    with open(reloc_file,"w") as f_obj:
        for line in new_dd:
            f_obj.write(line)
    f_obj.close()

def hypoDD_ref_days(reloc_file,ref_time,shift_hours=0):
    new_add=[]
    with open(reloc_file,"r") as f:
        for line in f:
            year = int(line[103:107])
            month = int(line[108:110])
            day = int(line[111:113])
            hour = int(line[114:116])
            minute = int(line[117:119])
            seconds = float(line[120:126])
            eve_time = UTCDateTime(year,month,day,hour,minute,0)+seconds - shift_hours*60*60
            days = (eve_time - ref_time)*1.0/(24*60*60)
            new_line=line[:-1]+" "+format(days,'4.2f')
            new_add.append(new_line)
    f.close()
    with open(reloc_file+".add","w") as f:
        for line in new_add:
            f.write(line+"\n")
    f.close()

def _region_subset(eve_dict,filt_condition="-999/-999/-999/-999/-999/-999"):
    """
    eve_dict is a dictionary,depth in km. Return a selected dict
    """
    #copy original dictionary to avoid action on original dictionary
    filt_dict = eve_dict.copy()
    #filt_condition
    lon_min,lon_max,lat_min,lat_max,dep_min,dep_max=re.split("/",filt_condition)
    if lon_min != "-999":
        for eve in eve_dict:
            if filt_dict[eve][0]<float(lon_min) or filt_dict[eve][0]>float(lon_max):
                filt_dict.pop(eve)
    tmp_dict = filt_dict.copy()
    if lat_min != "-999":
        for eve in tmp_dict:
            if filt_dict[eve][1]<float(lat_min) or filt_dict[eve][1]>float(lat_max):
                filt_dict.pop(eve)
    tmp_dict = filt_dict.copy()
    if dep_min != "-999":
        for eve in tmp_dict:
            if eve_dict[eve][2]<float(dep_min) or eve_dict[eve][2]>float(dep_max):
                filt_dict.pop(eve)
    return filt_dict
    
def _time_subset(eve_dict,starttime,endtime):
    """
    eve_dict is a dictionary, return a selected dict
    """
    filt_dict=eve_dict.copy()
    for eve in eve_dict:
        eve_time = UTCDateTime(eve)
        if eve_time < starttime or eve_time > endtime:
            filt_dict.pop(eve)
    return filt_dict
