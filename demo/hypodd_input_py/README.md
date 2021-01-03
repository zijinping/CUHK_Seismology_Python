# Configure eqcorr_hypdd.py at line 252-296

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
