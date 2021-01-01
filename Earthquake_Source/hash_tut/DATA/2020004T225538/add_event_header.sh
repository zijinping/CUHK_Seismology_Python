# !/bin/sh

# Event information
ev_time="GMT 2020 004 20 55 38 387"
ev_lat=22.072
ev_lon=113.846
ev_depth=12
ev_mag=3.47

SAC_file="tmp0*"

echo r $SAC_file
echo evla $ev_lat evlo $ev_lon evdp $ev_depth mag $ev_mag
echo  O $ev_time
sac << EOF
r $SAC_file
ch evla $ev_lat evlo $ev_lon evdp $ev_depth mag $ev_mag
cd O $ev_time
wh
q
EOF
