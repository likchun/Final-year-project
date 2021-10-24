import os, csv, math
import matplotlib.pyplot as plt

from tools import Dynamic, Compare

current_directory = os.path.dirname(os.path.abspath(__file__))

def plot_Change_in_Firing_Rate_Distribution(bins=[-3,12,125], xrange=[-3,12], yrange=[0.001,None], xaxis_logscale=False, yaxis_logscale=True, file_label=''):
    original_spiking_data_path = oringinal_spiking + '\\OUT_SPIK.txt'
    altered_spiking_data_path = [path['name']+'\\OUT_SPIK.txt' for path in altered_spiking]

    original_config_path = oringinal_spiking + '\\INI_CNFG'
    altered_config_path = [path['name']+'\\INI_CNFG' for path in altered_spiking]

    label_list = [path['label'] for path in altered_spiking]
    style_list = [path['style'] for path in altered_spiking]

    compare = Compare(original_spiking_data_path, original_config_path, altered_spiking_data_path, altered_config_path)

    compare.plot.ChangeInFiringRateDistribution(bins=bins, label_list=label_list, style_list=style_list, xrange=xrange, yrange=yrange, xaxis_logscale=xaxis_logscale, yaxis_logscale=yaxis_logscale, file_label=file_label)


def plot_Firing_Rate_and_ISI_Distribution_Stacked_Graph(input_data, FiR_range=[(0,10), (0,1.5)], IsI_range=[(0.0005,10), (0,2.5)], file_name=['Firing_Rate_Distribution', 'Interspike_Interval_Distribution'], file_label=['', '']):
    fig_FiR, ax_FiR = plt.subplots(figsize=(9, 6), dpi=50)
    fig_IsI, ax_IsI = plt.subplots(figsize=(9, 6), dpi=50)

    for datum in input_data:
        dynamic = Dynamic(datum['name']+'\\OUT_SPIK.txt', datum['name']+'\\INI_CNFG')

        ax_FiR = dynamic.plot.FiringRateDistribution(bins=datum['firing_rate_bins'], plot_axes=ax_FiR, info_list=[datum['style'], datum['label']])
        ax_IsI = dynamic.plot.InterSpikeIntervalDistribution(bins=datum['ISI_bins'], plot_axes=ax_IsI, info_list=[datum['style'], datum['label']])

    ax_FiR.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
    ax_FiR.set_xlim(FiR_range[0][0], FiR_range[0][1])
    ax_FiR.set_ylim(FiR_range[1][0], FiR_range[1][1])
    ax_FiR.grid(True)
    ax_FiR.legend()

    ax_IsI.set(xlabel='ISI (s)', ylabel='Probability density')
    ax_IsI.set_xlim(IsI_range[0][0], IsI_range[0][1])
    ax_IsI.set_ylim(IsI_range[1][0], IsI_range[1][1])
    ax_IsI.grid(True)
    ax_IsI.legend()

    output_file = [file_name[count]+file_label[count]+'.svg' for count in range(len(file_name))]
    fig_FiR.savefig(os.path.join(current_directory, output_file[0]))
    fig_IsI.savefig(os.path.join(current_directory, output_file[1]))


### Change in Firing Rate Distribution ###
# input folders of data #

oringinal_spiking = 'RANDOM'
altered_spiking = [
    {'name':'DIV66_INH_SUP_k025',   'label':'suppression level k = 0.25'},
    {'name':'DIV66_INH_SUP_k05',    'label':'suppression level k = 0.5' },
    {'name':'DIV66_INH_SUP_k075',   'label':'suppression level k = 0.75'},
    {'name':'DIV66_INH_SUP_k1',     'label':'suppression level k = 1'   },
    
    # {'name':'RANDOM_INH_SUP_k025',     'label':'suppression level k = 0.25', 'style':'cD'},
    # {'name':'RANDOM_INH_SUP_k075',     'label':'suppression level k = 0.75', 'style':'b^'},
    # {'name':'RANDOM_INH_ENH_k025',     'label':'enhancement level k = 0.25', 'style':'ms'},
    # {'name':'RANDOM_INH_ENH_k075',     'label':'enhancement level k = 0.75', 'style':'ro'},

    # {'name':'RANDOM_EXC_SUP_k025',     'label':'suppression level k = 0.25', 'style':'cD'},
    # {'name':'RANDOM_EXC_SUP_k075',     'label':'suppression level k = 0.75', 'style':'b^'},
    # {'name':'RANDOM_EXC_ENH_k025',     'label':'enhancement level k = 0.25', 'style':'ms'},
    # {'name':'RANDOM_EXC_ENH_k075',     'label':'enhancement level k = 0.75', 'style':'ro'},
]


### varying T ###
# input folders of data #

input_data_T = [
    {'name':'DIV66_60000_05',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,149], 'style':'cD', 'label':'T = 60000'},
    {'name':'DIV66_40000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,40,180], 'style':'mX', 'label':'T = 40000'},
    {'name':'DIV66_20000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,30,180], 'style':'gs', 'label':'T = 20000'},
    {'name':'DIV66_10000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'ro', 'label':'T = 10000'},
    {'name':'DIV66_7500_05',    'firing_rate_bins':[0,65,500], 'ISI_bins':[0.0005,30,180], 'style':'b^', 'label':'T = 7500' },
    # FiR_range=[(0,10), (0,1.6)], IsI_range=[(0.0005,10), (0,2.5)]

    # {'name':'DIV66_60000_01',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,180], 'style':'gs', 'label':'T = 60000'},
    # {'name':'DIV66_10000_01',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'ro', 'label':'T = 10000'},
    # {'name':'DIV66_7500_01',    'firing_rate_bins':[0,65,500], 'ISI_bins':[0.0005,40,138], 'style':'b^', 'label':'T = 7500' },
]


### varying dt ###
# input folders of data #

input_data_dt = [
    # [input folder name, bins for firing rate plot, bins for IsI plot, marker & line style, legend]
    # # bins for firing rate plot = [lower bound of bins, upper bound, total number of bins]
    # # bins for IsI plot = [lower bound of bins, upper bound, total number of bins]
    {'name':'DIV66_10000_125',  'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,10,135], 'style':'cD', 'label':'dt = 0.125'},
    {'name':'DIV66_10000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'mX', 'label':'dt = 0.05' },
    {'name':'DIV66_10000_025',  'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'gs', 'label':'dt = 0.025'},
    {'name':'DIV66_10000_01',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'ro', 'label':'dt = 0.01' },
    {'name':'DIV66_10000_005',  'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'b^', 'label':'dt = 0.005'},
    # FiR_range=[(0,10), (0,1.3)], IsI_range=[(0.0005,10), (0,2.7)]
    
#     {'name':'DIV66_60000_125',  'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,135], 'style':'gs', 'label':'dt = 0.125'},
#     {'name':'DIV66_60000_05',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,149], 'style':'ro', 'label':'dt = 0.05' },
#     {'name':'DIV66_60000_01',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,180], 'style':'b^', 'label':'dt = 0.01' },
]

if __name__ == '__main__':

    plot_Change_in_Firing_Rate_Distribution(bins=[-3,12,125], yaxis_logscale=True, file_label='DIV66') # DIV66
    # plot_Change_in_Firing_Rate_Distribution(bins=[0,70,125],  xrange=[None,None], yrange=[0,None], xaxis_logscale=False, yaxis_logscale=False, file_label='RANDOM_EXC_ENH') # RANDOM EXC ENH
    # plot_Change_in_Firing_Rate_Distribution(bins=[-3,1.5,38],  xrange=[-1.5,1.5], yrange=[0,None], yaxis_logscale=False, file_label='RANDOM_EXC_SUP') # RANDOM EXC SUP
    # plot_Change_in_Firing_Rate_Distribution(bins=[-10,10,160], xrange=[-1,1], yrange=[0,None], yaxis_logscale=False, file_label='RANDOM_INH') # RANDOM INH

    plot_Firing_Rate_and_ISI_Distribution_Stacked_Graph(input_data_T, FiR_range=[(0,10), (0,1.7)], IsI_range=[(0.0005,10), (0,2.7)], file_label=['_dt01', '_dt01'])

    plot_Firing_Rate_and_ISI_Distribution_Stacked_Graph(input_data_dt, FiR_range=[(0,10), (0,1.6)], IsI_range=[(0.0005,10), (0,1.6)], file_label=['_T60000', '_T60000'])
