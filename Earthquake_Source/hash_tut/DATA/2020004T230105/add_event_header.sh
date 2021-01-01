# !/bin/sh
# Update event header in sac and trim the waveform 

# Event information
ev_time="GMT 2020 004 23 01 05 641"
ev_lat=22.086
ev_lon=113.846
ev_depth=10.0
ev_mag=2.53

SAC_file="2020004225500.00*"

echo r $SAC_file
echo evla $ev_lat evlo $ev_lon evdp $ev_depth mag $ev_mag
echo  O $ev_time
sac << EOF
r $SAC_file
ch evla $ev_lat evlo $ev_lon evdp $ev_depth mag $ev_mag
ch o $ev_time
chnhdr allt (0 - &1,o&) IZTYPE IO
wh
cut o -30 180
r $SAC_file
w over
q
EOF
