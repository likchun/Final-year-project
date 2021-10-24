from tools import Dynamic


input_data = [
    # {'name':'DIV66',                'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    # {'name':'DIV66_INH_SUP_k1',     'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    # {'name':'DIV66_INH_SUP_k075',   'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    # {'name':'DIV66_INH_SUP_k05',    'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},
    # {'name':'DIV66_INH_SUP_k025',   'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,118], 'raster_stretch':1, 'tag':''},

    # {'name':'RANDOM_EXC_ENH_k025'},
    # {'name':'RANDOM_EXC_ENH_k075'},

    {'name':'DIV66_7500_05',    'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,180], 'raster_stretch':1},
    # {'name':'DIV66_7500_01',    'firing_rate_bins':[0, 65, 500], 'ISI_bins':[0.0005,40,180], 'raster_stretch':1},
    # {'name':'DIV66_10000_125',  'firing_rate_bins':[0, 65, 655], 'ISI_bins':[0.0005,10,135], 'raster_stretch':1},
    {'name':'DIV66_10000_05',   'firing_rate_bins':[0, 65, 660], 'ISI_bins':[0.0005,10,175], 'raster_stretch':1},
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
    dynamic = Dynamic(datum['name']+'\\OUT_SPIK.txt', datum['name']+'\\INI_CNFG', output_folder=datum['name'])

    print('\n> Generating results for \'{}\' <\n'.format(datum['name']))

    # dynamic.calculate.FiringRate_Stat()
    # dynamic.calculate.InterSpikeInterval_Stat()
    # dynamic.calculate.AverageFiringRate()
    # dynamic.calculate.TotalNumberOfSpikes()
    dynamic.calculate.SpikeCountBounds()

    # dynamic.calculate.identifyBurstingNode(output=True)
    # dynamic.calculate.reformatSpikeData_noindex()

    dynamic.plot.SpikeRaster(plot_horizontal_stretch=datum['raster_stretch'])

    dynamic.plot.FiringRateDistribution(bins=datum['firing_rate_bins'], xrange=[0,10], yrange=[0, None])
    # dynamic.plot.FiringRateDistribution(bins=datum['firing_rate_bins'], xrange=[0,10], yrange=[0, None], show_norm=datum['show_norm'])

    dynamic.plot.InterSpikeIntervalDistribution(bins=datum['ISI_bins'], xrange=[0.0005,10], file_label=datum['tag'])

    print('\n------------------------------next------------------------------')
