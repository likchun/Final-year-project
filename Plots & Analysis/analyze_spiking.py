from tools import Spiking


input_data = [
    {'name':'DIV66',                'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    {'name':'DIV66_INH_SUP_k1',     'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    {'name':'DIV66_INH_SUP_k075',   'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    {'name':'DIV66_INH_SUP_k05',    'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    {'name':'DIV66_INH_SUP_k025',   'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},

    # {'name':'DIV66_7500_05',    'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,180], 'raster_stretch':1},
    # {'name':'DIV66_7500_01',    'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,180], 'raster_stretch':1},
    # {'name':'DIV66_10000_125',  'firing_rate_bins':[0, 65, 655], 'ISI_bins':[0.0005,10,135], 'raster_stretch':1},
    # {'name':'DIV66_10000_05',   'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,10,175], 'raster_stretch':1},
    # {'name':'DIV66_10000_025',  'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,10,175], 'raster_stretch':1},
    # {'name':'DIV66_10000_01',   'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,10,175], 'raster_stretch':1},
    # {'name':'DIV66_10000_005',  'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,10,175], 'raster_stretch':1},
    # {'name':'DIV66_20000_05',   'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,30,180], 'raster_stretch':1.5},
    # {'name':'DIV66_20000_125',  'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,30,150], 'raster_stretch':1.5},
    # {'name':'DIV66_40000_05',   'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,40,180], 'raster_stretch':3},
    # {'name':'DIV66_60000_125',  'firing_rate_bins':[0, 65, 655], 'ISI_bins':[0.0005,50,135], 'raster_stretch':4.5},
    # {'name':'DIV66_60000_05',   'firing_rate_bins':[0, 65, 655], 'ISI_bins':[0.0005,50,149], 'raster_stretch':4.5},
    # {'name':'DIV66_60000_01',   'firing_rate_bins':[0, 65, 655], 'ISI_bins':[0.0005,50,180], 'raster_stretch':4.5},

    # {'name':'RANDOM_10000_05',  'firing_rate_bins':[-2, 4,  60], 'ISI_bins':[0.0005,10,110], 'show_norm':True},
]

for datum in input_data:
    s = Spiking(datum['name']+'\\OUT_SPIK.txt', datum['name']+'\\INI_CNFG', output_folder=datum['name'])
    s_calc = s.calculate
    s_plot = s.plot

    print('\n> Generating results for \'{}\' <\n'.format(datum['name']))

    # s_calc.FiringRate_Stat()
    # s_calc.InterSpikeInterval_Stat()
    # s_calc.AverageFiringRate()
    # s_calc.TotalNumberOfSpikes()
    # s_calc.SpikeCountBounds()

    # s_calc.identifyBurstingNode(output=True)
    s_calc.reformatSpikeData_noindex()
    # s_plot.SpikeRaster(plot_horizontal_stretch=datum['raster_stretch'])

    s_plot.FiringRateDistribution(bins=datum['firing_rate_bins'], xrange=[0,10], yrange=[0, None])
    # s_plot.FiringRateDistribution(bins=datum['firing_rate_bins'], xrange=[0,10], yrange=[0, None], show_norm=datum['show_norm'])

    s_plot.InterSpikeIntervalDistribution(bins=datum['ISI_bins'], xrange=[0.0005,10], file_label=datum['tag'])

    print('\n------------------------------next------------------------------')
