import spiking_compare_lib


### DIV66 Network C ###

input_folder = [
    'OUT_DATA_DIV66',
    'OUT_DATA_DIV66_suppress_inh_k025',
    'OUT_DATA_DIV66_suppress_inh_k05',
    'OUT_DATA_DIV66_suppress_inh_k075',
    'OUT_DATA_DIV66_suppress_inh_k1',
    # 'OUT_DATA_DIV66_enhance_exc_k025',
    # 'OUT_DATA_DIV66_enhance_exc_k05',
    # 'OUT_DATA_DIV66_enhance_exc_k075',
    # 'OUT_DATA_DIV66_enhance_exc_k1',
]

spiking_data = []
for each_path in input_folder:
    spiking_data.append(each_path+'\\OUT_SPIK.txt')

data_info = [
    # 'k = 0.25', 'k = 0.5', 'k = 1'
    'k = 0.25', 'k = 0.5', 'k = 0.75', 'k = 1'
]

file_info = [
    0.25, 0.5, 0.75, 1.0
]

suppressed_values = [
    0.25, 0.5, 0.75, 1.0
]; suppressed_values = [value * 0.00374561833693429 for value in suppressed_values]

enhanced_values = [
    0.25, 0.5, 0.75, 1.0
]; enhanced_values = [value * 0.006411179005874374 for value in enhanced_values]

sc = SpikingCompare(spiking_data, data_info, 'OUT_DATA_DIV66\\INI_CNFG', coupling='DIV66.txt', coupling_enhance_factor=2, coupling_delimiter='\t')
sc_calc = sc.calculate
sc_plot = sc.plot

# for file_idx in range(1, len(suppressed_value)+1):
#     sc_calc.findIncreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Increased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findDecreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Decreased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findUnchangedInFiringRate(file_idx, output=True, output_file='Nodes_of_Unchanged_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')

sc_plot.ChangeInFiringRateDistribution(bins=[-2,14,160], xrange=[-2,12], yaxis_logscale=True, file_label='DIV66_suppressed_logscale')
# sc_plot.ChangeInFiringRateDistribution(bins=[-2,14,150], xrange=[-2,12], yaxis_logscale=False, file_label='DIV66')
# sc_plot.ChangeInFiringRateDistribution(bins=[-2,14,150], xrange=[-2,2], yaxis_logscale=False, file_label='DIV66_nearzero')
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66_enhance.svg', xrange=(-2, 100), bins=100, yaxis_logscale=True)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66_enhance_nearzero.svg', xrange=(-1.5, 1.5), bins=32, yaxis_logscale=False)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66.svg', xrange=(-2, 10), bins=80, yaxis_logscale=False)

# sc_plot.RatioOfSuppression_vs_RatioOfIncreaseInFiringRate(suppressed_values, output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate_DIV66.svg', label='DIV66 (Network C)')
# sc_plot.RatioOfExcitation_vs_RatioOfIncreaseInFiringRate(excited_values, output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate_DIV66.svg', label='DIV66 (Network C)')


### RANDOM Network D ###

input_folder = [
    'OUT_DATA_RANDOM',
    'OUT_DATA_RANDOM_suppress_inh_k025',
    'OUT_DATA_RANDOM_suppress_inh_k05',
    'OUT_DATA_RANDOM_suppress_inh_k075',
    'OUT_DATA_RANDOM_suppress_inh_k1',
    # 'OUT_DATA_RANDOM_enhance_exc_k025',
    # 'OUT_DATA_RANDOM_enhance_exc_k05',
    # 'OUT_DATA_RANDOM_enhance_exc_k075',
    # 'OUT_DATA_RANDOM_enhance_exc_k1',
]

spiking_data = []
for each_path in input_folder:
    spiking_data.append(each_path+'\\OUT_SPIK.txt')

data_info = [
    # 'k = 0.25', 'k = 0.5', 'k = 1'
    'k = 0.25', 'k = 0.5', 'k = 0.75', 'k = 1'
]

file_info = [
    0.25, 0.5, 0.75#, 1.0
]

suppressed_values = [
    0.25, 0.5, 0.75, 1.0
]; suppressed_values = [value * 0.005611535353324736 for value in suppressed_values]

enhanced_values = [
    0.25, 0.5, 0.75#, 1.0
]; enhanced_values = [value * 0.005625560960642306 for value in enhanced_values]

sc = SpikingCompare(spiking_data, data_info, 'OUT_DATA_RANDOM\\INI_CNFG', coupling='Random.txt', coupling_enhance_factor=2, coupling_delimiter=' ')
sc_calc = sc.calculate
sc_plot = sc.plot

# for file_idx in range(1, len(suppressed_value)+1):
#     sc_calc.findIncreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Increased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findDecreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Decreased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findUnchangedInFiringRate(file_idx, output=True, output_file='Nodes_of_Unchanged_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')

sc_plot.ChangeInFiringRateDistribution(bins=[-1,1,24], xrange=[-1, 1], yaxis_logscale=False, file_label='RANDOM_suppressed')
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_RANDOM_suppress_0.svg', xrange=(-1, 1), bins=40, yaxis_logscale=False)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_RANDOM_enhance.svg', xrange=(-5, 100), bins=150, yaxis_logscale=True)

# sc_plot.RatioOfSuppression_vs_RatioOfIncreaseInFiringRate(suppressed_values, output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate_RANDOM.svg')
# sc_plot.RatioOfExcitation_vs_RatioOfIncreaseInFiringRate(enhanced_values, output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate_RANDOM.svg', label='RANDOM (Network D)')


### DIV66 & RANDOM ###

input_folder_1 = [
    'OUT_DATA_DIV66',
    'OUT_DATA_DIV66_suppress_inh_k025',
    'OUT_DATA_DIV66_suppress_inh_k05',
    'OUT_DATA_DIV66_suppress_inh_k075',
    'OUT_DATA_DIV66_suppress_inh_k1',
    # 'OUT_DATA_DIV66_enhance_exc_k025',
    # 'OUT_DATA_DIV66_enhance_exc_k05',
    # 'OUT_DATA_DIV66_enhance_exc_k075',
    # 'OUT_DATA_DIV66_enhance_exc_k1',
]

spiking_data_1 = []
for each_path in input_folder_1:
    spiking_data_1.append(each_path+'\\OUT_SPIK.txt')

input_folder_2 = [
    'OUT_DATA_RANDOM',
    'OUT_DATA_RANDOM_suppress_inh_k025',
    'OUT_DATA_RANDOM_suppress_inh_k05',
    'OUT_DATA_RANDOM_suppress_inh_k075',
    'OUT_DATA_RANDOM_suppress_inh_k1',
    # 'OUT_DATA_RANDOM_enhance_exc_k025',
    # 'OUT_DATA_RANDOM_enhance_exc_k05',
    # 'OUT_DATA_RANDOM_enhance_exc_k075',
    # 'OUT_DATA_RANDOM_enhance_exc_k1',
]

spiking_data_2 = []
for each_path in input_folder_2:
    spiking_data_2.append(each_path+'\\OUT_SPIK.txt')

suppressed_values_1 = [
    0.25, 0.5, 0.75, 1.0
]; suppressed_values_1 = [value * 0.00374561833693429 for value in suppressed_values_1]

suppressed_values_2 = [
    0.25, 0.5, 0.75, 1.0
]; suppressed_values_2 = [value * 0.005611535353324736 for value in suppressed_values_2]

# enhanced_values_1 = [
#     0.25*0.0128224, 0.5*0.0128224, 0.75*0.0128224, 1*0.0128224
# ]

# enhanced_values_2 = [
#     0.25*0.0112511, 0.5*0.0112511, 0.75*0.0112511#, 1*0.0112511
# ]

# SuppressionRatio_vs_FiringRateIncreaseRatio_Combined(spiking_data_1, spiking_data_2, 'OUT_DATA_DIV66\\INI_CNFG', 'DIV66.txt', 'Random.txt', suppressed_values_1, suppressed_values_2, coupling_enhance_factor_1=1, coupling_enhance_factor_2=1, coupling_delimiter_1='\t', coupling_delimiter_2=' ', label_1='DIV66 (Network C)', label_2='RANDOM (Network D)')
# RatioOfExcitation_vs_RatioOfIncreaseInFiringRate_Combined(spiking_data_1, spiking_data_2, 'OUT_DATA_DIV66\\INI_CNFG', 'DIV66.txt', 'Random.txt', enhanced_values_1, enhanced_values_2, coupling_enhance_factor_1=2, coupling_enhance_factor_2=2, coupling_delimiter_1='\t', coupling_delimiter_2=' ', label_1='DIV66 (Network C)', label_2='RANDOM (Network D)')


##########################
### Studying dt and Tn ###
##########################

input_folder = [
    'OUT_DATA_7500_0025',
    'OUT_DATA_7500_005',
    'OUT_DATA_7500_01',
    'OUT_DATA_7500_05',
    'OUT_DATA_7500_125',
    'OUT_DATA_7500_25',
    # 'OUT_DATA_10000_01',
    # 'OUT_DATA_10000_05',
    # 'OUT_DATA_10000_125',
    # 'OUT_DATA_15000_01',
    # 'OUT_DATA_15000_05',
    # 'OUT_DATA_15000_125'
]

spiking_data = []
for each_path in input_folder:
    spiking_data.append(each_path+'\\OUT_SPIK.txt')

data_info = [
    'dt = 0.0025', 'dt = 0.005', 'dt = 0.01', 'dt = 0.05', 'dt = 0.125', 'dt = 0.25', 
    # 'dt = 0.1', 'dt = 0.05', 'dt = 0.125'
    # 'Tn = 7500', 'Tn = 10000', 'Tn = 15000'
]

# sc = SpikingCompare(spiking_data, data_info, 'OUT_DATA_7500_01\\INI_CNFG')
# sc = SpikingCompare(spiking_data, data_info, 'OUT_DATA_15000_01\\INI_CNFG', simulataion_duration_override=[7500, 10000, 15000])
# sc_calc = sc.calculate
# sc_plot = sc.plot

# sc_plot.FiringRateDensity(file_label='T7500', bins=78, xrange=(0,10), yrange=(0, 1.3))
# sc_plot.FiringRateDensity(file_label='T10000', bins=105, xrange=(0,10), yrange=(0, 1.3))
# sc_plot.FiringRateDensity(file_label='T15000', bins=150, xrange=(0,10), yrange=(0, 1.3))

# sc_plot.InterSpikeIntervalDensity(file_label='T15000', xrange=(-3,1), bins=150)

# sc_plot.FiringRateDensity(file_label='dt125', bins=[78, 105, 150], xrange=(0,10), yrange=(0, 1.3))
# sc_plot.InterSpikeIntervalDensity(file_label='dt01', xrange=(-3,1), bins=150)
