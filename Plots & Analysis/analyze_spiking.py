from tools import Spiking


input_folder = [
    # 'OUT_DATA_DIV66_7500_125',
    # 'OUT_DATA_DIV66_7500_125_seed1',
    # 'OUT_DATA_DIV66_7500_05',
    # 'OUT_DATA_DIV66_7500_05_seed0',
    # 'OUT_DATA_DIV66_7500_01',
    # 'OUT_DATA_DIV66_10000_125',
    # 'OUT_DATA_DIV66_10000_05',
    # 'OUT_DATA_DIV66_10000_05_seed0',
    'OUT_DATA_DIV66_10000_025',
    # 'OUT_DATA_DIV66_20000_125',
    # 'OUT_DATA_DIV66_20000_05',
    # 'OUT_DATA_DIV66_40000_125',
    # 'OUT_DATA_DIV66_40000_05',
    # 'OUT_DATA_DIV66_60000_125',
    # 'OUT_DATA_DIV66_60000_05',
    # 'OUT_DATA_DIV66_60000_01'
]

for each_path in input_folder:
    s = Spiking(each_path+'\\OUT_SPIK.txt', each_path+'\\INI_CNFG', output_folder=each_path)
    s_calc = s.calculate
    s_plot = s.plot

    print('\n> Displaying details for \'{}\' <\n'.format(each_path))

    # s_calc.FiringRate_Stat()
    # s_calc.InterSpikeInterval_Stat()
    s_calc.AverageFiringRate()
    s_calc.TotalNumberOfSpikes()
    s_calc.SpikeCountBounds()

    # s_calc.identifyBurstingNode(output=True)
    # s_calc.reformatSpikeData()
    # s_plot.SpikeRaster(plot_horizontal_area_stretch=1) # 7500: 1, 20000: 1.5, 40000: 3, 60000: 4.5

    # s_plot.FiringRateDistribution(bins=[0, 65, 500], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 7500 01
    # s_plot.FiringRateDistribution(bins=[0, 65, 500], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 7500 05
    # s_plot.FiringRateDistribution(bins=[0, 65, 660], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 10000 05
    s_plot.FiringRateDistribution(bins=[0, 65, 660], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 10000 025
    # s_plot.FiringRateDistribution(bins=[0, 65, 660], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 20000 05
    # s_plot.FiringRateDistribution(bins=[0, 65, 660], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 40000 05
    # s_plot.FiringRateDistribution(bins=[0, 65, 655], xrange=[0,10], yrange=[0, None], file_label='fit') # DIV66 60000 05
    # s_plot.FiringRateDistribution(bins=[-2,4,60], xrange=[-2,4], show_norm=True, file_label='fit')      # RANDOM 10000

    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,40,138], xrange=[0.0005,10], file_label='fit') # DIV66 7500 01
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,40,180], xrange=[0.0005,10], file_label='fit') # DIV66 7500 05
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,175], xrange=[0.0005,10], file_label='fit') # DIV66 10000 05
    s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,175], xrange=[0.0005,10], file_label='fit') # DIV66 10000 025
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,30,180], xrange=[0.0005,10], file_label='fit') # DIV66 20000 05
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,40,180], xrange=[0.0005,10], file_label='fit') # DIV66 40000 05
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,50,149], xrange=[0.0005,40], file_label='fit') # DIV66 60000 05
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,50,172], xrange=[0.0005,40], file_label='fit') # DIV66 60000 125
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,110], xrange=[0.0005,35], file_label='fit') # RANDOM 10000

    print('\n------------------------------next------------------------------')
