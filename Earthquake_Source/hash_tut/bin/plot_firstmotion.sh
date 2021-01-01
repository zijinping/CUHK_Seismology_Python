# !/bin/sh
# Plot seismic waves on a regionalmap
# Required P arrival header in the sac file

### CONFIGURE THE SAC FILES DIRECTORY AND FILE NAMES
sacfolder='../DATA/2020004T225538'
#sacfolder='../../tmp/2020/20200105225538'
sac_file='*HZ*'



# Map range
# Automatic from station metadata
lon_min=`saclst stlo f $sacfolder/$sac_file|  awk '{print $2}' | sort -n | head -n 1 | awk '{print $1-0.1}'`
lon_max=`saclst stlo f $sacfolder/$sac_file|  awk '{print $2}' | sort -n | tail -n 1 | awk '{print $1+0.1}'`
lat_min=`saclst stla f $sacfolder/$sac_file|  awk '{print $2}' | sort -n | head -n 1 | awk '{print $1-0.1}'`
lat_max=`saclst stla f $sacfolder/$sac_file|  awk '{print $2}' | sort -n | tail -n 1 | awk '{print $1+0.1}'`

# Sichuan
#lon_min=97
#lon_max=108
#lat_min=26.0
#lat_max=34.1

# CheungChau earthquake
#lon_min=113.2
#lon_max=114.9
#lat_min=22
#lat_max=23

R=`echo ${lon_min}/${lon_max}/${lat_min}/${lat_max}`
echo $R

#112.5/115.8/21.5/24 Pearl delta
#100.1/99.8/34.1/25.0 Sichun


###### BEGIN PLOT ########
gmt begin firstmotion png

### Create basemap
gmt basemap -JM9i -R$R -Ba0.5 -BWeSn --MAP_FRAME_TYPE=plain

### Plot map with land filled, political boundaries and rivers
gmt coast -Ggrey -A1 -Df -Ia/0.1p,azure1 -Na/0.25p -W4/0.5p

### Plot station and station codes
saclst stlo stla KSTNM f $sacfolder/$sac_file | awk '{print $2,$3,$4}'| gmt text -F+f7p+jRB #-Dj-0.05/0.05
#awk '{print substr($0,53,7),substr($0,43,6),substr($0,1,4)}' HK.stations  | gmt text -Dj0.05/0.05 -F+f4p+jML
# plot symbol 
saclst stlo stla f $sacfolder/$sac_file |  gmt plot -Sc0.1c -i1,2
#awk '{print substr($0,53,7),substr($0,43,6)}' HK.stations|  gmt plot -St0.1

### Plot waveform 
# -T(startfrom p arrival) 
# -C(starttime from -0.2s to 1s of the reference time) 
# -G(fill with color) 
# -M(pen) -S(scale) -F(remove mean)

#gmt sac $sacfolder/$sac_file -M0.3i -S0.3c -C-0.2/0.5 -Fr -T+t-2
# Fill waveform with color at positive and negative 
gmt sac $sacfolder/*HZ* -M0.3i -Si1.5i -C-0.1/0.7 -T+t-2 -Gp+gblack -Gn+gred -Fr

### Plot Earthquake epicenter
saclst evlo evla f $sacfolder/$sac_file | head -n 1 |awk '{print $2,$3,0}'| gmt plot -Sa0.2 -Cred

### Plot FOCAL MECHANISM (RESULTS FROM HASH)
echo 113.846 22.072 12.0  253.618698       83.4670563      -138.425079 3.47 113.75 21.95 2020/01/04T22:55:38.387 | gmt meca -Sa0.2i+jML+f5 -C
#echo 113.846 22.072 12.0  260.889893 76.4130859 -127.247932 3.47 113.75 21.95 2020/01/04T22:55:38.387 | gmt meca -Sa0.2i+jML+f5 -C
echo 113.8478 22.086 10.0  278.461639 64.3977356 -105.76642  2.53 113.7 22.0 2020/01/04T23:01:05.641| gmt meca -Sa0.2i+jML+f5 -C

### Caption #> [X, Y] [font] [angle] [justify] linespace parwidth parjust(left,centered,right,justified)
lat_cap=`echo $lat_min| awk '{print $1 - 0.15}'`

#gmt text -F+f+a+j -B1 -M -N << EOF
#> $lon_min $lat_cap 8p,black 0 LT 10p 9i j
#@%5%Figure 1.@%% This figure shows the first motion of the arrivals.
#Star is the epiccenter of the recent earthquake on 2020/01/05 morning. 
#EOF

gmt end show
