README for Simon's W-band processing

VOCALS 2008 :: 2009-09-23 :: Simon de Szoeke

A number of m-files were written and run in a particular order, ensuring the data
flows properly from one stage to the next:

MAIN program line (in execution order)
proc_wband_1min_stat.m
    nc_varget_lapxm.m
    read_kongsberg.m
proc_wband_cloudtop.m
proc_wband_cloudtop_10min.m
compile_1min_stat.m

plot_1min_stat.m
cloud_drizzle_retrieval.m
    ave_small2large.m
    ave_large2small.m
    nc_varget_lapxm.m
    read_kongsberg.m

TEST programs
compare_top_sonde_wband.m

PRELIMINARY programs (in reverse chronological order)
proc_wband_Z_w.m
WbandPlot.m