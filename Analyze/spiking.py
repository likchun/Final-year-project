import spiking_lib


input_folder = [
    # 'OUT_DATA_7500_01',
    # 'OUT_DATA_7500_005',
    # 'OUT_DATA_7500_05',
    # 'OUT_DATA_7500_0025',
    # 'OUT_DATA_7500_25',
    'OUT_DATA_7500_125',
    # 'OUT_DATA_10000_01',
    # 'OUT_DATA_10000_05',
    # 'OUT_DATA_10000_125',
    # 'OUT_DATA_15000_01',
    # 'OUT_DATA_15000_05',
    # 'OUT_DATA_15000_125',
    # 'OUT_DATA_DIV66',
    # 'OUT_DATA_DIV66_enhance_exc_k1',
    # 'OUT_DATA_DIV66_enhance_exc_k075',
    # 'OUT_DATA_DIV66_enhance_exc_k05',
    # 'OUT_DATA_DIV66_enhance_exc_k025',
    # 'OUT_DATA_DIV66_suppress_inh_k1',
    # 'OUT_DATA_DIV66_suppress_inh_k075',
    # 'OUT_DATA_DIV66_suppress_inh_k05',
    # 'OUT_DATA_DIV66_suppress_inh_k025',
    # 'OUT_DATA_RANDOM',
    # 'OUT_DATA_RANDOM_suppress_inh_k1',
    # 'OUT_DATA_RANDOM_suppress_inh_k075',
    # 'OUT_DATA_RANDOM_suppress_inh_k05',
    # 'OUT_DATA_RANDOM_suppress_inh_k025',
    # 'OUT_DATA_RANDOM_enhance_exc_k1',
    # 'OUT_DATA_RANDOM_enhance_exc_k075',
    # 'OUT_DATA_RANDOM_enhance_exc_k05',
    # 'OUT_DATA_RANDOM_enhance_exc_k025',
]

for each_path in input_folder:
    s = Spiking(each_path+'\\OUT_SPIK.txt', each_path+'\\INI_CNFG', output_folder=each_path)
    s_calc = s.calculate
    s_plot = s.plot

    s_calc.AverageFiringRate()
    # s_calc.FiringRate_Stat()
    # s_calc.InterSpikeInterval_Stat()

    s_calc.identifyBurstingNode(output=True)
    s_plot.reformatSpikeData()
    s_plot.SpikeRaster()

    # s_plot.FiringRateDistribution()
    # s_plot.FiringRateDistribution(bins=[0, 65, 660], xrange=[0,20], yrange=[0, None], file_label='fit') # DIV66 10000
    # s_plot.FiringRateDistribution(bins=[-2,4,60], xrange=[-2,4], show_norm=True) # for RANDOM

    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,150], xrange=[0.0005,1], file_label='fit') # DIV66
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,110], xrange=[0.0005,1.5], file_label='fit') # RANDOM



    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_0.svg', bins=100, xrange=(0,5))
    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_1.svg', bins=80, xrange=(0,5))
    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_fit.svg', bins=78, xrange=(0,10), yrange=(0, 1.3)) #for DIV 7500
    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_fit.svg', bins=150, xrange=(0,10), yrange=(0, 1.3)) #for DIV 15000
    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_fit_0.svg', bins=150, xrange=(0,16), yrange=(0, 1.3)) #for DIV enhance exc
    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_fit.svg', bins=150, xrange=(0,100), show_norm=True) # for RANDOM network D enhance exc
    # s_plot.FiringRateDistribution(output_file='Firing_Rate_Density_lowfit.svg', bins=30, xrange=(0,100), show_norm=True) # for RANDOM network D enhance exc

    # s_plot.InterSpikeIntervalDensity(output_file='Interspike_Interval_Density_0.svg', bins=400)
