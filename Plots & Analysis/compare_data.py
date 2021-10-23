import os, csv, math
import matplotlib.pyplot as plt

from tools import Spiking, Compare

current_directory = os.path.dirname(os.path.abspath(__file__))

def plot_Change_in_Firing_Rate_Distribution(bins=[-3,12,125]):
    original_spiking_data_path = oringinal_spiking + '\\OUT_SPIK.txt'
    altered_spiking_data_path = [path['name']+'\\OUT_SPIK.txt' for path in altered_spiking]

    original_config_path = oringinal_spiking + '\\INI_CNFG'
    altered_config_path = [path['name']+'\\INI_CNFG' for path in altered_spiking]

    label_list = [path['label'] for path in altered_spiking]

    c = Compare(original_spiking_data_path, original_config_path, altered_spiking_data_path, altered_config_path)
    c_plot = c.plot

    c_plot.ChangeInFiringRateDistribution(bins=bins, label_list=label_list, xrange=[-3,12], yrange=[0.001,None])


def plot_Firing_Rate_and_ISI_Distribution_Stacked_Graph(FiR_range=None, IsI_range=None, file_name=['Firing_Rate_Distribution', 'Interspike_Interval_Distribution'], file_label=['', '']):
    fig_FiR, ax_FiR = plt.subplots(figsize=(9, 6), dpi=50)
    fig_IsI, ax_IsI = plt.subplots(figsize=(9, 6), dpi=50)

    for datum in input_data:
        s = Spiking(datum['name']+'\\OUT_SPIK.txt', datum['name']+'\\INI_CNFG')
        s_plot = s.plot

        ax_FiR = s_plot.FiringRateDistribution(bins=datum['firing_rate_bins'], plot_axes=ax_FiR, info_list=[datum['style'], datum['label']])
        ax_IsI = s_plot.InterSpikeIntervalDistribution(bins=datum['ISI_bins'], plot_axes=ax_IsI, info_list=[datum['style'], datum['label']])

        ax_FiR.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
        ax_FiR.set_xlim(xrange[0][0], 10)
        ax_FiR.set_ylim(0, 1.6)
        ax_FiR.grid(True)
        ax_FiR.legend()

        ax_IsI.set(xlabel='ISI (s)', ylabel='Probability density')
        ax_IsI.set_xlim(0.0005, 10)
        ax_IsI.set_ylim(0, 3)
        ax_IsI.grid(True)
        ax_IsI.legend()

    output_file = [file_name[count]+file_label[count]+'.svg' for count in range(file_name)]
    fig_FiR.savefig(os.path.join(current_directory, output_file[0]))
    fig_IsI.savefig(os.path.join(current_directory, output_file[1]))


### Change in Firing Rate Distribution ###
# input folders of data #

oringinal_spiking = 'DIV66'
altered_spiking = [
    {'name':'DIV66_INH_SUP_k025',   'label':'suppression level k = 0.25'},
    {'name':'DIV66_INH_SUP_k05',    'label':'suppression level k = 0.5' },
    {'name':'DIV66_INH_SUP_k075',   'label':'suppression level k = 0.75'},
    {'name':'DIV66_INH_SUP_k1',     'label':'suppression level k = 1'   },
]


### varying T ###
# input folders of data #

input_data = [
    {'name':'DIV66_60000_05',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,149], 'style':'cD', 'label':'T = 60000'},
    {'name':'DIV66_40000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,40,180], 'style':'mX', 'label':'T = 40000'},
    {'name':'DIV66_20000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,30,180], 'style':'gs', 'label':'T = 20000'},
    {'name':'DIV66_10000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'ro', 'label':'T = 10000'},
    {'name':'DIV66_7500_05',    'firing_rate_bins':[0,65,500], 'ISI_bins':[0.0005,30,180], 'style':'b^', 'label':'T = 7500' },

    # {'name':'DIV66_60000_01',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,180], 'style':'gs', 'label':'T = 60000'},
    # {'name':'DIV66_10000_01',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'ro', 'label':'T = 10000'},
    # {'name':'DIV66_7500_01',    'firing_rate_bins':[0,65,500], 'ISI_bins':[0.0005,40,138], 'style':'b^', 'label':'T = 7500' },
]


### varying dt ###
# input folders of data #

input_data = [
    # [input folder name, bins for firing rate plot, bins for IsI plot, marker & line style, legend]
    # # bins for firing rate plot = [lower bound of bins, upper bound, total number of bins]
    # # bins for IsI plot = [lower bound of bins, upper bound, total number of bins]
    {'name':'DIV66_10000_125',  'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,10,135], 'style':'cD', 'label':'dt = 0.125'},
    {'name':'DIV66_10000_05',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'mX', 'label':'dt = 0.05' },
    {'name':'DIV66_10000_025',  'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'gs', 'label':'dt = 0.025'},
    {'name':'DIV66_10000_01',   'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'ro', 'label':'dt = 0.01' },
    {'name':'DIV66_10000_005',  'firing_rate_bins':[0,65,660], 'ISI_bins':[0.0005,10,175], 'style':'b^', 'label':'dt = 0.005'},
    
    # {'name':'DIV66_60000_125',  'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,135], 'style':'gs', 'label':'dt = 0.05' },
    # {'name':'DIV66_60000_05',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,149], 'style':'ro', 'label':'dt = 0.025'},
    # {'name':'DIV66_60000_01',   'firing_rate_bins':[0,65,655], 'ISI_bins':[0.0005,50,180], 'style':'b^', 'label':'dt = 0.005'},
]

if __name__ == '__main__':

    plot_Change_in_Firing_Rate_Distribution(bins=[-3,12,125])

    # plot_Firing_Rate_and_ISI_Distribution_Stacked_Graph(FiR_range=[(0,10), (0,1.6)], IsI_range=[(0.0005,10), (0,2.6)], file_label=['_dt05', '_dt05'])

    # plot_Firing_Rate_and_ISI_Distribution_Stacked_Graph(FiR_range=[(0,10), (0,1.5)], IsI_range=[(0.0005,10), (0,3)], file_label=['_T10000', '_T10000'])
