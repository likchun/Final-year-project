from tools import Compare, SuppressionRatio_vs_FiringRateIncreaseRatio_Combined, EnhancementRatio_vs_RatioOfIncreaseInFiringRate_Combined

import os, csv, math
import numpy as np
import matplotlib.pyplot as plt

from tools import Spiking

current_directory = os.path.dirname(os.path.abspath(__file__))

### varying T ###

input_data = [
    ['OUT_DATA_DIV66_60000_05', [0,65,655], [0.0005,50,149], 'cD', 'T = 60000'],
    ['OUT_DATA_DIV66_40000_05', [0,65,660], [0.0005,40,180], 'mX', 'T = 40000'],
    ['OUT_DATA_DIV66_20000_05', [0,65,660], [0.0005,30,180], 'gs', 'T = 20000'],
    ['OUT_DATA_DIV66_10000_05', [0,65,660], [0.0005,10,175], 'ro', 'T = 10000'],
    ['OUT_DATA_DIV66_7500_05' , [0,65,500], [0.0005,40,180], 'b^', 'T = 7500' ]
]

fig_FiR, ax_FiR = plt.subplots(figsize=(9, 6), dpi=50)
fig_IsI, ax_IsI = plt.subplots(figsize=(9, 6), dpi=50)

for datum in input_data:
    s = Spiking(datum[0]+'\\OUT_SPIK.txt', datum[0]+'\\INI_CNFG')
    s_plot = s.plot

    ax_FiR = s_plot.FiringRateDistribution(bins=datum[1], plot_axes=ax_FiR, info_list=datum[3:])
    ax_IsI = s_plot.InterSpikeIntervalDistribution(bins=datum[2], plot_axes=ax_IsI, info_list=datum[3:])

    ax_FiR.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
    ax_FiR.set_xlim(0, 10)
    ax_FiR.set_ylim(0, 1.5)
    ax_FiR.grid(True)
    ax_FiR.legend()

    ax_IsI.set(xlabel='ISI (s)', ylabel='Probability density')
    ax_IsI.set_xlim(0.0005, 10)
    ax_IsI.set_ylim(0, 2.5)
    ax_IsI.grid(True)
    ax_IsI.legend()

fig_FiR.savefig(os.path.join(current_directory, 'Firing_Rate_Distribution_varT.svg'))
fig_IsI.savefig(os.path.join(current_directory, 'Interspike_Interval_Distribution_varT.svg'))


### varying dt ###

# input_data = [
#     ['OUT_DATA_DIV66_10000_125', [0,65,655], 'cD', 'dt = 0.125'],
#     ['OUT_DATA_DIV66_10000_05',  [0,65,660], 'mX', 'dt = 0.05'],
#     ['OUT_DATA_DIV66_10000_025', [0,65,660], 'gs', 'dt = 0.025'],
#     ['OUT_DATA_DIV66_10000_01',  [0,65,660], 'ro', 'dt = 0.01'],
#     ['OUT_DATA_DIV66_10000_005', [0,65,500], 'b^', 'dt = 0.005']
# ]


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

# sc = Compare(spiking_data, data_info, 'OUT_DATA_DIV66\\INI_CNFG', coupling='DIV66.txt', coupling_enhance_factor=2, coupling_delimiter='\t')
# sc_calc = sc.calculate
# sc_plot = sc.plot

# sc_plot.ChangeInFiringRateDistribution(bins=[-2,14,160], xrange=[-2,12], yaxis_logscale=True, file_label='DIV66_suppressed_logscale')
# sc_plot.ChangeInFiringRateDistribution(bins=[-2,14,150], xrange=[-2,12], yaxis_logscale=False, file_label='DIV66')
# sc_plot.ChangeInFiringRateDistribution(bins=[-2,14,150], xrange=[-2,2], yaxis_logscale=False, file_label='DIV66_nearzero')

# sc_plot.SuppressionRatio_vs_FiringRateIncreaseRatio(suppressed_values, output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate_DIV66.svg', label='DIV66 (Network C)')
# sc_plot.EnhancementRatio_vs_FiringRateIncreaseRatio(excited_values, output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate_DIV66.svg', label='DIV66 (Network C)')

# for file_idx in range(1, len(suppressed_value)+1):
#     sc_calc.findIncreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Increased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findDecreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Decreased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findUnchangedInFiringRate(file_idx, output=True, output_file='Nodes_of_Unchanged_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')


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
    'k = 0.25', 'k = 0.5', 'k = 0.75', 'k = 1'
]

file_info = [
    0.25, 0.5, 0.75, 1.0
]

suppressed_values = [
    0.25, 0.5, 0.75, 1.0
]; suppressed_values = [value * 0.005611535353324736 for value in suppressed_values]

enhanced_values = [
    0.25, 0.5, 0.75, 1.0
]; enhanced_values = [value * 0.005625560960642306 for value in enhanced_values]

# sc = Compare(spiking_data, data_info, 'OUT_DATA_RANDOM\\INI_CNFG', coupling='Random.txt', coupling_enhance_factor=2, coupling_delimiter=' ')
# sc_calc = sc.calculate
# sc_plot = sc.plot

# sc_plot.ChangeInFiringRateDistribution(bins=[-1,1,24], xrange=[-1, 1], yaxis_logscale=False, file_label='RANDOM_suppressed')

# sc_plot.SuppressionRatio_vs_FiringRateIncreaseRatio(suppressed_values, plot_label='RANDOM (Network D)', file_label='RANDOM')
# sc_plot.EnhancementRatio_vs_FiringRateIncreaseRatio(enhanced_values, plot_label='RANDOM (Network D)', file_label='RANDOM')

# for file_idx in range(1, len(suppressed_value)+1):
#     sc_calc.findIncreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Increased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findDecreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Decreased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findUnchangedInFiringRate(file_idx, output=True, output_file='Nodes_of_Unchanged_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')


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

# ]

# enhanced_values_2 = [

# ]

# SuppressionRatio_vs_FiringRateIncreaseRatio_Combined(spiking_data_1, spiking_data_2, 'OUT_DATA_DIV66\\INI_CNFG', 'DIV66.txt', 'Random.txt', suppressed_values_1, suppressed_values_2, coupling_enhance_factor_1=1, coupling_enhance_factor_2=1, coupling_delimiter_1='\t', coupling_delimiter_2=' ', label_1='DIV66 (Network C)', label_2='RANDOM (Network D)')
# EnhancementRatio_vs_FiringRateIncreaseRatio(spiking_data_1, spiking_data_2, 'OUT_DATA_DIV66\\INI_CNFG', 'DIV66.txt', 'Random.txt', enhanced_values_1, enhanced_values_2, coupling_enhance_factor_1=2, coupling_enhance_factor_2=2, coupling_delimiter_1='\t', coupling_delimiter_2=' ', label_1='DIV66 (Network C)', label_2='RANDOM (Network D)')
