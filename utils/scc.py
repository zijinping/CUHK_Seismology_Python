# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
#    Author: ZI, Jinping
# Institute: The Seismology Lab, the Chinese University of Hong Kong
# Reference: Yang, H. L. Zhu, and R. Chu(2009), Fault-Plane Determination of
#            the 18 April 2008 Mt. Carmel, Illinois, Earthquake by Detection
#            and Relocating Aftershocks, BSSA, 99(6), 3413-3420
#            doi:10.1785/0120090038
#   History: 
#       2021-01-25 Initial coding
#
#     Usage: python scc.py [-Ccc] [-E] [-Mn] [-O] [-Tlength] [-Wt1/t2[/maxShift]]
#            -C: cross-correlation threshold (default = 0.7)
#-----------------------------------------------------------------------------

import obspy
from obspy import Stream,UTCDateTime
from math import sqrt
import numpy as np
import sys

def wf_scc(tmplt_st,sta_st,ncom):
    """
    Sliding-window cross-correlation between template and target waveform
    reference: Yang et al.,2009, BSSA.

    Parameters
    -----------
    tmplt_st: template waveform of one station
       sta_st: target waveform of the same station
         ncom: number of component, n = 3 means 3 components cross-correlation

    Return
    -----------
        ccmax: maximum cross-correlation coefficient
        aamax: amplitude ratio at ccmax
        i0: the shifting index at ccmax
    """

    tmplt_data = []
    sta_data = []
    if ncom == 3:
        tmplt_data.append([tmplt_st.select(component="*N")[0].data])
        tmplt_data.append([tmplt_st.select(component="*E")[0].data])
        tmplt_data.append([tmplt_st.select(component="*Z")[0].data])
        sta_data.append([sta_st.select(component="*N")[0].data])
        sta_data.append([sta_st.select(component="*E")[0].data])
        sta_data.append([sta_st.select(component="*Z")[0].data])
    elif ncom == 1:
        tmplt_data.append([tmplt_sta[0].data])
    dt = temp[0].stats.delta
    ccmax,aamax,i0,cc_list = data_scc(tmplt_data,st_data,ncom)


def data_scc(tmplt_data,st_data,ncom):
    """
    Sliding-window cross-correlation between template and target waveform
    reference: Yang, H. et al.,2009, BSSA.

    Parameters
    -----------
    tmplt_data: template waveform of one station
       st_data: target waveform of the same station
          ncom: number of component, n = 3 means 3 components cross-correlation

    return
    ----------
         ccmax: maximum cross-correlation coefficient
         aamax: amplitude ratio at ccmax
            i0: the shifting index at ccmax
       cc_list: the list of cross-correlation coefficient in each step
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

if __name__ == "__main__":
    """
    Usage
    """
    para_len = len(sys.args)