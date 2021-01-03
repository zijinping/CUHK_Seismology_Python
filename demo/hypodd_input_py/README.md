## Introduction of the script
```python
    """
    Create input event.dat, phase pair (dt.ct) and waveform cross-correlation
    (dt.cc) files for hypodd 
        Processes:
            0 Subset events
            1 ph2dt
            2.1 Generate event streams (Optional for prebuilt already)
            2.2 Cross-correlation
        Configure the script after __main__
        Input: 
            Catalog file (readable via obspy.read_events)
                - Formats: .pha, .xml, etc...
            Continuous Waveform
                - Store in daily base (year/doy/*.[sac/mseed])
        Output: (In the same directory of the script)
            Event.dat   - Subsetted event with assigned id
            dt.ct       - phase-pair
            dt.cc       - cross-correlation pair
            [ For intermediate step ]
            event_stream.p      - temporal storage of event wf
            event_stream_p.p    - temporal storage of P arrival wf
            event_stream_s.p    - temporal storage of S arrival wf
            event_id_mapper.p   - temporal storage of hypoDD id to catalog id mapper
            hypodd_subset.p     - temporal storage of subseted obspy.catalog object
    Created on Fri Nov  6 11:57:19 HKT 2020
    @author: jw
    """
```
## Future update of the script
- [ ] Local magnitude estimation
- [ ] Support catalog plotting (Map, time series, B value)
- [ ] **Parallel Processing**

Possible future development: 
- [ ] Convert to event sacfile
- [ ] Use ASDF to speed up file I/O processes 


## Configure eqcorr_hypdd.py at line 252-296
```python
    if __name__ == "__main__":
        ######## Configurations ########
        hypodd_pha = read_events('/NAS2/Sichuan_data/phasenet/event_waveform/HYPODD2/weiyuan_subset.pha')


        ### 0 Subset events with latlon and time ranges
        lat = [29.48, 29.73]
        lon = [104.65, 104.95]

        starttime = UTC(2018,9,1)
        endtime = UTC(2019,3,1)

        ### 1 Generate hypodd output (ph2dt)
        max_sep = 20            # Maximum seperations between event paris in km
        min_link = 8            # minimum links for an event to be paired,
                                # e.g. min. no. of picks with same sta& chan
                                # shared between two paired events


        ### 2 Cross-correlation Configurations
        load_prebuilt = False           # Select true if event waveforms are stored
                                        #   event_stream.p, ...

        ## Continuous wf directory
        # Waveform store in waveform_dir/year/doy
        waveform_dir = '/NAS2/Sichuan_data/event_waveform_PN'

        ## Waveform CC configurations
        before_arrival = 0.5
        after_arrival = 1.5

        filter_kwargs = {'type':'bandpass',
                         'freqmin': 2,
                         'freqmax': 15}

        print('Filter settings')

        ## Cross-correlation pair configurations
        max_sep = 20
        min_link_cc = 4
        min_cc = 0.5
        max_shift=1.0



        ######## End of Configurations ########
```
