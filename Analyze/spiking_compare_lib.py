import os, csv, math
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.style as mplstyle
import seaborn as sns

from coupling_lib import Coupling
from spiking_lib import Spiking


class Initialization:

    def __init__(self, spiking_data, simulation_config, delimiter, input_folder, output_folder):
        self.initSpikeData(spiking_data, simulation_config, delimiter, input_folder, output_folder)

    def initSpikeData(self, spiking_data, simulation_config, delimiter, input_folder, output_folder):
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if input_folder == None: input_folder = ''
        else: input_folder = input_folder + '\\'
        input_path = this_dir + '\\' + input_folder
        if output_folder == None: output_folder = ''
        else: output_folder = output_folder + '\\'
        self.output_path = this_dir + '\\' + output_folder
        try: os.mkdir(self.output_path)
        except FileExistsError: pass
        with open(input_path+simulation_config, 'r') as file_config:
            reader = csv.reader(file_config, delimiter=delimiter)
            self.Config = np.array(list(reader), dtype=object)
        Nsize = int(self.Config[0][0])
        self.SpikeTimes = [None] * len(spiking_data)
        self.SpikeCount = [None] * len(spiking_data)
        file_spike = [None] * len(spiking_data)
        count_file = 0
        for file_path in spiking_data:
            self.SpikeTimes[count_file] = np.empty((Nsize), dtype=object)
            self.SpikeCount[count_file] = np.zeros(Nsize)
            with open(input_path+file_path, 'r') as file_spike:
                reader = csv.reader(file_spike, delimiter=delimiter)
                count_node = 0
                for row in reader:
                    self.SpikeTimes[count_file][count_node] = np.delete(np.array(list(row)).astype('float'), [0,1], 0)
                    self.SpikeCount[count_file][count_node] = int(row[1])
                    count_node += 1
            count_file += 1
        self.SpikeCount = np.array(self.SpikeCount)
        self.SpikeTimes = np.array(self.SpikeTimes)

class SpikingCompare(Initialization):

    def __init__(self, spiking_data: list, data_info: list, simulation_config: str, delimiter='\t', input_folder=None, output_folder=None, coupling=None, coupling_enhance_factor=1, coupling_delimiter='\t', simulataion_duration_override=[]):
        super().__init__(spiking_data, simulation_config, delimiter, input_folder, output_folder)
        if len(simulataion_duration_override) == 0: FiringRate = self.SpikeCount / float(self.Config[0][2]) * 1000
        else:
            # simulation duration override: used when comparing data with diffrent simulation duration Tn
            FiringRate = self.SpikeCount
            for count in range(len(simulataion_duration_override)):
                FiringRate[count] = FiringRate[count] / float(simulataion_duration_override[count]) * 1000
        ChangeInFiringRate = np.zeros((len(FiringRate)-1, int(self.Config[0][0])))
        for count in range(1, len(FiringRate)):
            ChangeInFiringRate[count-1] = FiringRate[count] - FiringRate[0]
        self.calculate = self.Calculate(ChangeInFiringRate, self.Config, data_info, self.output_path, spiking_data, simulation_config)
        self.plot = self.Plot(self.SpikeTimes, FiringRate, ChangeInFiringRate, self.Config, data_info, self.output_path, spiking_data, simulation_config, coupling, coupling_enhance_factor, coupling_delimiter)

    class Calculate:

        def __init__(self, FiringRate, Config, data_info, output_path, spiking_data, simulation_config):
            self.FiringRate = FiringRate
            self.Config = Config
            self.data_info = data_info
            self.Nsize = int(self.Config[0][0])
            self.output_path = output_path
            self.spiking_data = spiking_data
            self.simulation_config = simulation_config

        def findIncreasedInFiringRate(self, file_index: int, output=False, output_file='Nodes_of_Increased_in_Firing_Rate.txt'):
            # file_index starts from 1
            index_of_increased = np.argwhere(self.ChangeInFiringRate[file_index-1] > 0).flatten()
            activity_of_increased = self.ChangeInFiringRate[file_index-1][index_of_increased].flatten()
            if output == True:
                with open(self.output_path+output_file, 'w') as file_output:
                    for count in range(len(index_of_increased)):
                        file_output.write('{:d}\t{}\n'.format(index_of_increased[count], activity_of_increased[count]))
            return index_of_increased, activity_of_increased

        def findDecreasedInFiringRate(self, file_index: int, output=False, output_file='Nodes_of_Decreased_in_Firing_Rate.txt'):
            # file_index starts from 1
            index_of_decreased = np.argwhere(self.ChangeInFiringRate[file_index-1] < 0).flatten()
            activity_of_decreased = self.ChangeInFiringRate[file_index-1][index_of_decreased].flatten()
            if output == True:
                with open(self.output_path+output_file, 'w') as file_output:
                    for count in range(len(index_of_decreased)):
                        file_output.write('{:d}\t{}\n'.format(index_of_decreased[count], activity_of_decreased[count]))
            return index_of_decreased, activity_of_decreased

        def findUnchangedInFiringRate(self, file_index: int, output=False, output_file='Nodes_of_Unchanged_in_Firing_Rate.txt'):
            # file_index starts from 1
            index_of_unchanged = np.argwhere(self.ChangeInFiringRate[file_index-1] == 0).flatten()
            activity_of_unchanged = self.ChangeInFiringRate[file_index-1][index_of_unchanged].flatten()
            if output == True:
                with open(self.output_path+output_file, 'w') as file_output:
                    for count in range(len(index_of_unchanged)):
                        file_output.write('{:d}\t{}\n'.format(index_of_unchanged[count], activity_of_unchanged[count]))
            return index_of_unchanged, activity_of_unchanged

        def RatioOfIncreaseInFiringRate(self):
            avg_firing_rate_original_network = Spiking(self.spiking_data[0], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            count = 0
            avg_firing_rate_altered_networks = np.zeros(len(self.spiking_data)-1)
            for each_datum in self.spiking_data[1:]:
                avg_firing_rate_altered_networks[count] = Spiking(each_datum, self.simulation_config).calculate.AverageFiringRate(console_print=False)
                count += 1
            return avg_firing_rate_altered_networks / avg_firing_rate_original_network-1

    class Plot:

        def __init__(self, SpikeTimes, FiringRate, ChangeInFiringRate, Config, data_info, output_path, spiking_data, simulation_config, coupling, coupling_enhance_factor, coupling_delimiter):
            self.SpikeTimes = SpikeTimes
            self.FiringRate = FiringRate
            self.ChangeInFiringRate = ChangeInFiringRate
            self.Config = Config
            self.data_info = data_info
            self.Nsize = int(self.Config[0][0])
            self.output_path = output_path
            self.spiking_data = spiking_data
            self.simulation_config = simulation_config
            self.coupling = coupling
            self.coupling_enhance_factor = coupling_enhance_factor
            self.coupling_delimiter = coupling_delimiter

        def FiringRateDensity(self, output_file='Firing_Rate_Density', bins=60, xrange=(0,10), yrange=(0,None), file_label=''):
            density = []; x_value = []
            if type(bins) == list: 
                for count in range(0, len(self.FiringRate)):
                    density_each, bin_edges = np.histogram(self.FiringRate[count], bins=np.linspace(xrange[0], xrange[1], bins[count]), density=True)
                    density.append(density_each)
                    x_value.append((bin_edges[1:] + bin_edges[:-1]) / 2)
                
                fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
                color_list = ['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', '']
                for count in range(len(density)):
                    ax.plot(x_value[count], density[count], str(color_list[count])+'-', lw=2, label='Firing rate density '+str(self.data_info[count]))
                ax.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
                ax.set_xlim(xrange[0], xrange[1])
                ax.set_ylim(yrange[0], yrange[1])
                ax.grid(True)
                ax.legend()
                if file_label == '': output_file += '.svg'
                else: output_file += '_' + file_label + '.svg'
                fig.savefig(os.path.join(self.output_path, output_file))
                plt.clf()
            else:
                for count in range(0, len(self.FiringRate)):
                    density_each, bin_edges = np.histogram(self.FiringRate[count], bins=np.linspace(xrange[0], xrange[1], bins), density=True)
                    density.append(density_each)
                x_value = (bin_edges[1:] + bin_edges[:-1]) / 2
                
                fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
                color_list = ['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', '']
                count = 0
                for each in density:
                    ax.plot(x_value, each, str(color_list[count])+'-', lw=2, label='Firing rate density '+str(self.data_info[count]))
                    count += 1
                ax.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
                ax.set_xlim(xrange[0], xrange[1])
                ax.set_ylim(yrange[0], yrange[1])
                ax.grid(True)
                ax.legend()
                if file_label == '': output_file += '.svg'
                else: output_file += '_' + file_label + '.svg'
                fig.savefig(os.path.join(self.output_path, output_file))
                plt.clf()
        
        def InterSpikeIntervalDensity(self, output_file='Interspike_Interval_Density', bins=150, xrange=(-3,1), file_label=''):
            density = []
            for count_1 in range(0, len(self.SpikeTimes)):
                IsI = np.empty(self.SpikeTimes[count_1].shape[0], dtype=object)
                count_2 = 0
                for row in self.SpikeTimes[count_1]:
                    try: IsI[count_2] = np.diff(row)
                    except ValueError: IsI[count_2] = np.diff(np.array([0]))
                    count_2 += 1
                IsI = np.concatenate([item for item in IsI.flatten()], 0) / 1000
                density_each, bin_edges = np.histogram(IsI, bins=np.logspace(xrange[0], xrange[1], bins))
                density_each = np.array(density_each, dtype=float)
                density_each /= np.dot(density_each, np.diff(bin_edges)) # normalization
                density.append(density_each)
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            color_list = ['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', '']
            count = 0
            for each in density:
                ax.semilogx(x_value, each, str(color_list[count])+'-', lw=2, label='Log ISI density '+str(self.data_info[count]))
                count += 1
            ax.set(xlabel='ISI (s)', ylabel='Probability density')
            # ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.grid(True)
            ax.legend()
            if file_label == '': output_file += '.svg'
            else: output_file += '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()

        def ChangeInFiringRateDistribution(self, bins=[-2,12,120], xrange=[-2,12], yrange=[None,None], yaxis_logscale=False, output_file='Change_in_Firing_Rate_Distribution', file_label=''):
            if np.amax(self.ChangeInFiringRate) > bins[1]: print('Warning! Maximum of Change in Firing Rate exceeds upper bound of bins range. Max Change: {}; Max bins range: {}'.format(np.amax(self.ChangeInFiringRate), bins[1]))
            if np.amin(self.ChangeInFiringRate) < bins[0]: print('Warning! Minimum of Change in Firing Rate subceeds lower bound of bins range. Min Change: {}; Min bins range: {}'.format(np.amin(self.ChangeInFiringRate), bins[0]))
            density = np.zeros(len(self.ChangeInFiringRate), dtype=object)
            for count in range(1, len(self.ChangeInFiringRate)+1):
                density_each, bin_edges = np.histogram(self.ChangeInFiringRate[count-1], bins=np.linspace(bins[0], bins[1], bins[2]), density=True)
                density[count-1] = density_each
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            color_list = ['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', '']
            if yaxis_logscale == False:
                for count in range(len(density)):
                    ax.plot(x_value, density[count], str(color_list[count])+'-', lw=2, label=''+str(self.data_info[count]))
            elif yaxis_logscale == True:
                for count in range(len(density)):
                    plot_data_corrected = np.vstack((density[count], x_value))
                    plot_data_corrected = np.delete(plot_data_corrected, np.argwhere(plot_data_corrected[0] == 0), 1)                          
                    ax.semilogy(plot_data_corrected[1], plot_data_corrected[0], str(color_list[count])+'-', lw=2, label=''+str(self.data_info[count]))
            ax.set(xlabel='Change in firing rate (Hz)', ylabel='Probability density')
            ax.set_xlim(xrange[0], xrange[1])
            ax.set_ylim(yrange[0], yrange[1])
            ax.grid(True)
            ax.legend()

            if file_label == '': output_file_plot = output_file + '.svg'
            else: output_file_plot = output_file + '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()
            if file_label == '': output_file_info = output_file + '.txt'
            else: output_file_info = output_file + '_' + file_label + '.txt'
            with open(self.output_path+output_file_info, 'w') as fp_info:
                fp_info.write('### Plot information ###\n\n')
                fp_info.write('Max change in firing rate: {}\n'.format(np.amax(self.ChangeInFiringRate)))
                fp_info.write('Min change in firing rate: {}\n'.format(np.amin(self.ChangeInFiringRate)))
                fp_info.write('\n### Plot settings ###\n\n')
                fp_info.write('Bin size: {}\n'.format((bins[1]-bins[0])/bins[2]))
                fp_info.write('Bin bounds: lower: {}, upper: {}\n'.format(bins[0], bins[1]))
                fp_info.write('Number of bins: {}\n'.format(bins[2]))

        def RatioOfSuppression_vs_RatioOfIncreaseInFiringRate(self, suppressed_values: list, output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate.svg', label=''):
            if self.coupling == None: print('Error. No coupling matrix input.'); return -1
            ratio_of_suppression = Coupling(self.coupling, coupling_enhance_factor=self.coupling_enhance_factor, delimiter=self.coupling_delimiter).calculate.RatioOfSuppression(suppressed_values)
            avg_firing_rate_original_network = Spiking(self.spiking_data[0], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            avg_firing_rate_altered_networks = np.zeros(len(suppressed_values))
            for count in range(0, len(suppressed_values)):
                avg_firing_rate_altered_networks[count] = Spiking(self.spiking_data[count+1], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            ratio_of_increase_in_firing_rate = avg_firing_rate_altered_networks / avg_firing_rate_original_network - 1
            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(ratio_of_suppression, ratio_of_increase_in_firing_rate, 'b^-', lw=2, label=label)
            ax.set(xlabel='Ratio of suppression in inhibitory synaptic weights', ylabel='Ratio of increase in avergae firing rate') 
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))
        
        def RatioOfExcitation_vs_RatioOfIncreaseInFiringRate(self, excited_values: list, output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate.svg', label=''):
            if self.coupling == None: print('Error. No coupling matrix input.'); return -1
            ratio_of_excitation = Coupling(self.coupling, coupling_enhance_factor=self.coupling_enhance_factor, delimiter=self.coupling_delimiter).calculate.RatioOfExcitation(excited_values)
            avg_firing_rate_original_network = Spiking(self.spiking_data[0], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            avg_firing_rate_altered_networks = np.zeros(len(excited_values))
            for count in range(0, len(excited_values)):
                avg_firing_rate_altered_networks[count] = Spiking(self.spiking_data[count+1], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            ratio_of_increase_in_firing_rate = avg_firing_rate_altered_networks / avg_firing_rate_original_network - 1
            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(ratio_of_excitation, ratio_of_increase_in_firing_rate, 'b^-', lw=2, label=label)
            ax.set(xlabel='Ratio of excitation in excitatory synaptic weights', ylabel='Ratio of increase in avergae firing rate') 
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))

def SuppressionRatio_vs_FiringRateIncreaseRatio_Combined(spiking_data_1, spiking_data_2, simulation_config, coupling_1, coupling_2, suppressed_values_1: list, suppressed_values_2: list, coupling_enhance_factor_1, coupling_enhance_factor_2, coupling_delimiter_1, coupling_delimiter_2, label_1='', label_2='', output_file='Suppression_Ratio_vs_Firing_Rate_Increase_Ratio', file_label=''):
    ratio_of_suppression_1 = Coupling(coupling_1, coupling_enhance_factor=coupling_enhance_factor_1, delimiter=coupling_delimiter_1).calculate.RatioOfSuppression(suppressed_values_1)
    ratio_of_suppression_2 = Coupling(coupling_2, coupling_enhance_factor=coupling_enhance_factor_2, delimiter=coupling_delimiter_2).calculate.RatioOfSuppression(suppressed_values_2)
    avg_firing_rate_original_network_1 = Spiking(spiking_data_1[0], simulation_config).calculate.AverageFiringRate(console_print=False)
    avg_firing_rate_original_network_2 = Spiking(spiking_data_2[0], simulation_config).calculate.AverageFiringRate(console_print=False)
    avg_firing_rate_altered_networks_1 = np.zeros(len(suppressed_values_1))
    avg_firing_rate_altered_networks_2 = np.zeros(len(suppressed_values_2))
    for count in range(0, len(suppressed_values_1)):
        avg_firing_rate_altered_networks_1[count] = Spiking(spiking_data_1[count+1], simulation_config).calculate.AverageFiringRate(console_print=False)
    for count in range(0, len(suppressed_values_2)):
        avg_firing_rate_altered_networks_2[count] = Spiking(spiking_data_2[count+1], simulation_config).calculate.AverageFiringRate(console_print=False)
    ratio_of_increase_in_firing_rate_1 = avg_firing_rate_altered_networks_1 / avg_firing_rate_original_network_1 - 1
    ratio_of_increase_in_firing_rate_2 = avg_firing_rate_altered_networks_2 / avg_firing_rate_original_network_2 - 1
    fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
    ax.plot(ratio_of_suppression_1, ratio_of_increase_in_firing_rate_1, 'b^-', lw=2, label=label_1)
    ax.plot(ratio_of_suppression_2, ratio_of_increase_in_firing_rate_2, 'ro-', lw=2, label=label_2)
    ax.set(xlabel='Ratio of suppression in inhibitory synaptic weights', ylabel='Ratio of increase in avergae firing rate')
    ax.set_ylim(-0.15, 1.05)
    ax.grid(True)
    ax.legend()
    
    if file_label == '': output_file_plot = output_file + '.svg'
    else: output_file_plot = output_file + '_' + file_label + '.svg'
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file_plot)); plt.clf()

def EnhancementRatio_vs_RatioOfIncreaseInFiringRate_Combined(spiking_data_1, spiking_data_2, simulation_config, coupling_1, coupling_2, enhanced_values_1: list, enhanced_values_2: list, coupling_enhance_factor_1, coupling_enhance_factor_2, coupling_delimiter_1, coupling_delimiter_2, label_1='', label_2='', output_file='Enhancement_Ratio_vs_Firing_Rate_Increase_Ratio', file_label=''):
    ratio_of_enhancement_1 = Coupling(coupling_1, coupling_enhance_factor=coupling_enhance_factor_1, delimiter=coupling_delimiter_1).calculate.RatioOfEnhancement(enhanced_values_1)
    ratio_of_enhancement_2 = Coupling(coupling_2, coupling_enhance_factor=coupling_enhance_factor_2, delimiter=coupling_delimiter_2).calculate.RatioOfEnhancement(enhanced_values_2)
    avg_firing_rate_original_network_1 = Spiking(spiking_data_1[0], simulation_config).calculate.AverageFiringRate(console_print=False)
    avg_firing_rate_original_network_2 = Spiking(spiking_data_2[0], simulation_config).calculate.AverageFiringRate(console_print=False)
    avg_firing_rate_altered_networks_1 = np.zeros(len(enhanced_values_1))
    avg_firing_rate_altered_networks_2 = np.zeros(len(enhanced_values_2))
    for count in range(0, len(enhanced_values_1)):
        avg_firing_rate_altered_networks_1[count] = Spiking(spiking_data_1[count+1], simulation_config).calculate.AverageFiringRate(console_print=False)
    for count in range(0, len(enhanced_values_2)):
        avg_firing_rate_altered_networks_2[count] = Spiking(spiking_data_2[count+1], simulation_config).calculate.AverageFiringRate(console_print=False)
    ratio_of_increase_in_firing_rate_1 = avg_firing_rate_altered_networks_1 / avg_firing_rate_original_network_1 - 1
    ratio_of_increase_in_firing_rate_2 = avg_firing_rate_altered_networks_2 / avg_firing_rate_original_network_2 - 1
    fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
    ax.plot(ratio_of_enhancement_1, ratio_of_increase_in_firing_rate_1, 'b^-', lw=2, label=label_1)
    ax.plot(ratio_of_enhancement_2, ratio_of_increase_in_firing_rate_2, 'ro-', lw=2, label=label_2)
    ax.set(xlabel='Ratio of enhancement in excitatory synaptic weights', ylabel='Ratio of increase in avergae firing rate')
    ax.grid(True)
    ax.legend()
    
    if file_label == '': output_file_plot = output_file + '.svg'
    else: output_file_plot = output_file + '_' + file_label + '.svg'
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file_plot)); plt.clf()