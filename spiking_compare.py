import os, csv, math
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.style as mplstyle
import seaborn as sns

from coupling import Coupling
from spiking import Spiking


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

        def ChangeInFiringRateDensity(self, output_file='Change_in_Firing_Rate_Density.svg', bins=100, xrange=(-2,12), show_norm=False, yaxis_logscale=False):
            density = np.zeros(len(self.ChangeInFiringRate), dtype=object)
            for count in range(1, len(self.ChangeInFiringRate)+1):
                density_each, bin_edges = np.histogram(self.ChangeInFiringRate[count-1], bins=np.linspace(xrange[0], xrange[1], bins), density=True)
                density[count-1] = density_each
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            color_list = ['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', '']
            count = 0
            if yaxis_logscale == False:
                for each in density:
                    ax.plot(x_value, each, str(color_list[count])+'-', lw=2, label='Change in Firing rate density '+str(self.data_info[count]))
                    count += 1
            elif yaxis_logscale == True:
                for each in density:
                    ax.semilogy(x_value, each, str(color_list[count])+'-', lw=2, label='Change in Firing rate density '+str(self.data_info[count]))
                    count += 1
            # if show_norm == True:
            #     mu = np.mean(FiringRate.flatten()); sigma = np.std(FiringRate.flatten())
            #     norm = stats.norm(loc=mu, scale=sigma)
            #     ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')
            ax.set(xlabel='Change in Firing rate (Hz)', ylabel='Probabili ty density')
            # ax.set_xlim(xrange[0], xrange[1])
            # ax.set_ylim(0.001, None)
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))

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

def RatioOfSuppression_vs_RatioOfIncreaseInFiringRate_Combined(spiking_data_1, spiking_data_2, simulation_config, coupling_1, coupling_2, suppressed_values_1: list, suppressed_values_2: list, coupling_enhance_factor_1, coupling_enhance_factor_2, coupling_delimiter_1, coupling_delimiter_2, label_1='', label_2='', output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate.svg'):
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
    ax.grid(True)
    ax.legend()
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file))

def RatioOfExcitation_vs_RatioOfIncreaseInFiringRate_Combined(spiking_data_1, spiking_data_2, simulation_config, coupling_1, coupling_2, enhanced_values_1: list, enhanced_values_2: list, coupling_enhance_factor_1, coupling_enhance_factor_2, coupling_delimiter_1, coupling_delimiter_2, label_1='', label_2='', output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate.svg'):
    ratio_of_excitation_1 = Coupling(coupling_1, coupling_enhance_factor=coupling_enhance_factor_1, delimiter=coupling_delimiter_1).calculate.RatioOfSuppression(enhanced_values_1)
    ratio_of_excitation_2 = Coupling(coupling_2, coupling_enhance_factor=coupling_enhance_factor_2, delimiter=coupling_delimiter_2).calculate.RatioOfSuppression(enhanced_values_2)
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
    ax.plot(ratio_of_excitation_1, ratio_of_increase_in_firing_rate_1, 'b^-', lw=2, label=label_1)
    ax.plot(ratio_of_excitation_2, ratio_of_increase_in_firing_rate_2, 'ro-', lw=2, label=label_2)
    ax.set(xlabel='Ratio of suppression in inhibitory synaptic weights', ylabel='Ratio of increase in avergae firing rate')
    ax.grid(True)
    ax.legend()
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file))


#######################
### Spiking Compare ###
#######################

input_folder = [
    'OUT_DATA_DIV66',
    # 'OUT_DATA_DIV66_suppress_inh_k025',
    # 'OUT_DATA_DIV66_suppress_inh_k05',
    # 'OUT_DATA_DIV66_suppress_inh_k075',
    # 'OUT_DATA_DIV66_suppress_inh_k1',
    'OUT_DATA_DIV66_enhance_exc_k025',
    'OUT_DATA_DIV66_enhance_exc_k05',
    'OUT_DATA_DIV66_enhance_exc_k075',
    'OUT_DATA_DIV66_enhance_exc_k1',
]

spiking_data = []
for each_path in input_folder:
    spiking_data.append(each_path+'\\OUT_SPIK.txt')

data_info = [
    # 'k = 0.25', 'k = 0.5', 'k = 1'
    'k = 0.25', 'k = 0.5', 'k = 0.75', 'k = 1'
]

file_info = [
    0.25, 0.5, 0.75, 1
]

suppressed_values = [
    0.25*0.0074912, 0.5*0.0074912, 0.75*0.0074912, 1*0.0074912
]

excited_values = [
    0.25*0.0128224, 0.5*0.0128224, 0.75*0.0128224, 1*0.0128224
]

# sc = SpikingCompare(spiking_data, data_info, 'OUT_DATA_DIV66\\INI_CNFG', coupling='DIV66.txt', coupling_enhance_factor=2, coupling_delimiter='\t')
# sc_calc = sc.calculate
# sc_plot = sc.plot

# for file_idx in range(1, len(suppressed_value)+1):
#     sc_calc.findIncreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Increased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findDecreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Decreased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findUnchangedInFiringRate(file_idx, output=True, output_file='Nodes_of_Unchanged_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')

# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66_suppress.svg', xrange=(-2, 45), bins=100, yaxis_logscale=True)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66_suppress_nearzero.svg', xrange=(-1.5, 1.5), bins=32, yaxis_logscale=False)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66_enhance.svg', xrange=(-2, 100), bins=100, yaxis_logscale=True)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66_enhance_nearzero.svg', xrange=(-1.5, 1.5), bins=32, yaxis_logscale=False)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_DIV66.svg', xrange=(-2, 10), bins=80, yaxis_logscale=False)

# sc_plot.RatioOfSuppression_vs_RatioOfIncreaseInFiringRate(suppressed_values, output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate_DIV66.svg', label='DIV66 (Network C)')
# sc_plot.RatioOfExcitation_vs_RatioOfIncreaseInFiringRate(excited_values, output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate_DIV66.svg', label='DIV66 (Network C)')


input_folder = [
    'OUT_DATA_RANDOM',
    # 'OUT_DATA_RANDOM_suppress_inh_k025',
    # 'OUT_DATA_RANDOM_suppress_inh_k05',
    # 'OUT_DATA_RANDOM_suppress_inh_k075',
    # 'OUT_DATA_RANDOM_suppress_inh_k1',
    'OUT_DATA_RANDOM_enhance_exc_k025',
    'OUT_DATA_RANDOM_enhance_exc_k05',
    'OUT_DATA_RANDOM_enhance_exc_k075',
    # 'OUT_DATA_RANDOM_enhance_exc_k1',
]

spiking_data = []
for each_path in input_folder:
    spiking_data.append(each_path+'\\OUT_SPIK.txt')

data_info = [
    # 'k = 0.25', 'k = 0.5', 'k = 1'
    'k = 0.25', 'k = 0.5', 'k = 0.75'#, 'k = 1'
]

file_info = [
    0.25, 0.5, 0.75#, 1
]

suppressed_values = [
    0.25*0.011223, 0.5*0.011223, 0.75*0.011223, 1*0.011223
]

enhanced_values = [
    0.25*0.0112511, 0.5*0.0112511, 0.75*0.0112511#, 1*0.0112511
]

sc = SpikingCompare(spiking_data, data_info, 'OUT_DATA_RANDOM\\INI_CNFG', coupling='Random.txt', coupling_enhance_factor=2, coupling_delimiter=' ')
sc_calc = sc.calculate
sc_plot = sc.plot

# for file_idx in range(1, len(suppressed_value)+1):
#     sc_calc.findIncreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Increased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findDecreasedInFiringRate(file_idx, output=True, output_file='Nodes_of_Decreased_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')
#     sc_calc.findUnchangedInFiringRate(file_idx, output=True, output_file='Nodes_of_Unchanged_in_Firing_Rate_k_'+str(file_info[file_idx-1])+'.txt')

# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_RANDOM_suppress_fit.svg', xrange=(-1, 1), bins=22, yaxis_logscale=False)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_RANDOM_suppress_0.svg', xrange=(-1, 1), bins=40, yaxis_logscale=False)
# sc_plot.ChangeInFiringRateDensity(output_file='Change_in_Firing_Rate_Density_RANDOM_enhance.svg', xrange=(-5, 100), bins=150, yaxis_logscale=True)

# sc_plot.RatioOfSuppression_vs_RatioOfIncreaseInFiringRate(suppressed_values, output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate_RANDOM.svg')
sc_plot.RatioOfExcitation_vs_RatioOfIncreaseInFiringRate(enhanced_values, output_file='Ratio_of_Excitation_vs_Ratio_of_Increase_in_Firing_Rate_RANDOM.svg', label='RANDOM (Network D)')




input_folder_1 = [
    'OUT_DATA_DIV66',
    # 'OUT_DATA_DIV66_suppress_inh_k025',
    # 'OUT_DATA_DIV66_suppress_inh_k05',
    # 'OUT_DATA_DIV66_suppress_inh_k075',
    # 'OUT_DATA_DIV66_suppress_inh_k1',
    'OUT_DATA_DIV66_enhance_exc_k025',
    'OUT_DATA_DIV66_enhance_exc_k05',
    'OUT_DATA_DIV66_enhance_exc_k075',
    'OUT_DATA_DIV66_enhance_exc_k1',
]

spiking_data_1 = []
for each_path in input_folder_1:
    spiking_data_1.append(each_path+'\\OUT_SPIK.txt')

input_folder_2 = [
    'OUT_DATA_RANDOM',
    # 'OUT_DATA_RANDOM_suppress_inh_k025',
    # 'OUT_DATA_RANDOM_suppress_inh_k05',
    # 'OUT_DATA_RANDOM_suppress_inh_k075',
    # 'OUT_DATA_RANDOM_suppress_inh_k1',
    'OUT_DATA_RANDOM_enhance_exc_k025',
    'OUT_DATA_RANDOM_enhance_exc_k05',
    'OUT_DATA_RANDOM_enhance_exc_k075',
    # 'OUT_DATA_RANDOM_enhance_exc_k1',
]

spiking_data_2 = []
for each_path in input_folder_2:
    spiking_data_2.append(each_path+'\\OUT_SPIK.txt')

# suppressed_values_1 = [
#     0.25*0.0074912, 0.5*0.0074912, 0.75*0.0074912, 1*0.0074912
# ]

# suppressed_values_2 = [
#     0.25*0.011223, 0.5*0.011223, 0.75*0.011223, 1*0.011223
# ]

enhanced_values_1 = [
    0.25*0.0128224, 0.5*0.0128224, 0.75*0.0128224, 1*0.0128224
]

enhanced_values_2 = [
    0.25*0.0112511, 0.5*0.0112511, 0.75*0.0112511#, 1*0.0112511
]

# RatioOfSuppression_vs_RatioOfIncreaseInFiringRate_Combined(spiking_data_1, spiking_data_2, 'OUT_DATA_DIV66\\INI_CNFG', 'DIV66.txt', 'Random.txt', suppressed_values_1, suppressed_values_2, coupling_enhance_factor_1=2, coupling_enhance_factor_2=2, coupling_delimiter_1='\t', coupling_delimiter_2=' ', label_1='DIV66 (Network C)', label_2='RANDOM (Network D)')
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