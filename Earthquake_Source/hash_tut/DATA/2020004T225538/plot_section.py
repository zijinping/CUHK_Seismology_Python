#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Section of the sac files

@auther: jw
"""
from obspy import read

import matplotlib;matplotlib.use('tkagg')


_sac_files="sac/*HHZ*"

st = read(_sac_files)

print(st)
# convert sac.dist to tr.stats.distance
for tr in st:
    tr.stats.distance = tr.stats.sac['dist']*1000

#import pdb; pdb.set_trace()
st.filter('highpass', freq = 1)
# plot section of the waveform
st.plot(type='section',orientation='horizontal',time_down=True, recordlength=180)
