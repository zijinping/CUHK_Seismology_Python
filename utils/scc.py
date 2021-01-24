import obspy
from obspy import Stream,UTCDateTime
from math import sqrt
import numpy as np

def wf_scc(tmplt_sta,st_sta,ncom):
    """
    Sliding-window cross-correlation between template and target waveform
    reference: Yang et al.,2009, BSSA.
    tmplt_sta: template waveform of one station
    st_sta: target waveform of the same station
    ncom: number of component, n = 3 means 3 components cross-correlation
    return:
        ccmax: maximum cross-correlation coefficient
        aamax: amplitude ratio at ccmax
        i0: the shifting index at ccmax
    """
    temp=Stream()
    stream = Stream()
    temp.append(tmplt_sta.select(component="*N")[0])
    temp.append(tmplt_sta.select(component="*E")[0])
    temp.append(tmplt_sta.select(component="*Z")[0])
    stream.append(st_sta.select(component="*N")[0])
    stream.append(st_sta.select(component="*E")[0])
    stream.append(st_sta.select(component="*Z")[0])
    mm = temp[0].stats.npts
    dt = temp[0].stats.delta
    normMaster = 0.
    ic = 0
    while ic < ncom:
        k = 0
        while k < mm:
            normMaster += temp[ic].data[k]*temp[ic].data[k]
            k += 1
        ic += 1
    normMaster = sqrt(normMaster)
    
    starttime = stream[0].stats.starttime
    npts = stream[0].stats.npts
    norm = 0
    j = 0
    while j < mm-1:
        ic=0
        while ic < ncom:
            norm += stream[ic][j]*stream[ic][j]
            ic+=1
        j=j+1
    ccmax = 0
    j=0
    cc_list=[]
    while j<=npts-mm:
        cc = 0
        ic = 0
        while ic < ncom:
            norm+=stream[ic].data[j+mm-1]*stream[ic].data[j+mm-1]
            k=0
            while k <mm:
                cc += temp[ic].data[k]*stream[ic][j+k]
                k+=1
            ic+=1
        aa = sqrt(norm)/normMaster
        cc = cc*aa/norm #cc = <f|g>/sqrt((f|f)(g|g))
        if(cc>=ccmax):
            ccmax = cc
            aamax = aa
            i0 = j
        cc_list.append(cc)
        ic = 0
        while ic < ncom:
            norm -= stream[ic][j]*stream[ic][j]
            ic+=1
        j=j+1
    return ccmax,aamax,i0

def data_scc(tmplt_data,st_data,ncom):
    """
    Sliding-window cross-correlation between template and target waveform
    reference: Yang et al.,2009, BSSA.
    tmplt_data: template waveform of one station
    st_data: target waveform of the same station
    ncom: number of component, n = 3 means 3 components cross-correlation
    return:
        ccmax: maximum cross-correlation coefficient
        aamax: amplitude ratio at ccmax
        i0: the shifting index at ccmax
    """
    normMaster = 0.
    ic = 0
    mm = len(tmplt_data[0])
    while ic < ncom:
        k = 0
        while k < mm:
            normMaster += tmplt_data[ic][k]*tmplt_data[ic][k]
            k += 1
        ic += 1
    normMaster = sqrt(normMaster)
    
    npts = len(st_data[0])
    norm = 0
    j = 0
    while j < mm-1:
        ic=0
        while ic < ncom:
            norm += st_data[ic][j]*st_data[ic][j]
            ic+=1
        j=j+1
    ccmax = -1
    aamax = -1
    j=0
    cc_list=[]
    while j<=npts-mm:
        cc = 0
        ic = 0
        while ic < ncom:
            norm+=st_data[ic][j+mm-1]*st_data[ic][j+mm-1]
            k=0
            while k <mm:
                cc += tmplt_data[ic][k]*st_data[ic][j+k]
                k+=1
            ic+=1
        aa = sqrt(norm)/normMaster
        cc = cc*aa/norm #cc = <f|g>/sqrt((f|f)(g|g))
        if(cc>=ccmax):
            ccmax = cc
            aamax = aa
            i0 = j
        cc_list.append(cc)
        ic = 0
        while ic < ncom:
            norm -= st_data[ic][j]*st_data[ic][j]
            ic+=1
        j=j+1
    return ccmax,aamax,i0,cc_list

