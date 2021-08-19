"""
    Usage of basic_utils: 
    year,month,day = month_day(year,julday)
    days = julday(year,month,day)
    status, st = get_data(net,sta,chn,starttime,endtime,f_folder,time_zone_effect=0)
    degree = spherical_dist(lon_1,lat_1,lon_2,lat_2)
    idx,diff = find_nearest(array,value)
    hypoDD_ref_days(dd_file,ref_time,shift_hours), no returns
    hypoDD_mag_mapper(reloc_file,out_sum)
    eve_list,df = load_hypoDD(reloc_file='hypoDD.reloc',shift_hours)
    filt_dict = _region_subset(eve_dict,filt_condition="-999/-999/-999/-999/-999/-999")
    filt_dict = _time_subset(eve_dict,starttime,endtime)

    Usage of phase_convert:
    REAL2phs(input_file,phase_filt=8)
    SC2phs(file_list=[],region_condition="-9/-9/-9/-9",mag_condition=-9,res_condition=-9,starttime=0,endtime=0))
"""
