from tools import Spiking


input_folder = [
    '0', '1', '2', '3', '4'
]

for each_path in input_folder:
    s = Spiking(each_path+'\\OUT_SPIK.txt', each_path+'\\INI_CNFG', output_folder=each_path)
    s_calc = s.calculate
    s_plot = s.plot

    # s_calc.AverageFiringRate()
    # s_calc.FiringRate_Stat()
    # s_calc.InterSpikeInterval_Stat()

    # s_calc.identifyBurstingNode(output=True)
    s_plot.reformatSpikeData()
    # s_plot.SpikeRaster()

    # s_plot.FiringRateDistribution()
    s_plot.FiringRateDistribution(bins=[0, 55, 420], xrange=[0,20], yrange=[0, None], file_label='fit') # DIV66 7500
    # s_plot.FiringRateDistribution(bins=[0, 65, 660], xrange=[0,20], yrange=[0, None], file_label='fit') # DIV66 10000
    # s_plot.FiringRateDistribution(bins=[-2,4,60], xrange=[-2,4], show_norm=True) # for RANDOM

    s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,85], xrange=[0.0005,1], file_label='fit') # DIV66 7500
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,150], xrange=[0.0005,1], file_label='fit') # DIV66 10000
    # s_plot.InterSpikeIntervalDistribution(bins=[0.0005,10,110], xrange=[0.0005,1.5], file_label='fit') # RANDOM