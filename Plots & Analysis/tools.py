import os, csv, math
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.style as mplstyle
import seaborn as sns


class Coupling:

    def __init__(self, coupling_matrix, coupling_enhance_factor=1, delimiter='\t', input_folder=None, output_folder=None):
        self.initCouplingMatrix(coupling_matrix, delimiter, input_folder, output_folder)
        self.Coupling = coupling_enhance_factor * self.Coupling
        self.calculate = self.Calculate(self.Coupling)
        self.plot = self.Plot(self.Coupling, self.output_path)

    def initCouplingMatrix(self, coupling_matrix, delimiter, input_folder, output_folder):
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if input_folder == None: input_folder = ''
        else: input_folder = input_folder + '\\'
        input_path = this_dir + '\\' + input_folder
        if output_folder == None: output_folder = ''
        else: output_folder = output_folder + '\\'
        self.output_path = this_dir + '\\' + output_folder
        try: os.mkdir(self.output_path)
        except FileExistsError: pass

        with open(input_path+coupling_matrix, 'r', newline='') as file_input:
            reader = csv.reader(file_input, delimiter=delimiter)
            self.Coupling = np.array(list(reader)).astype('float')

    class Calculate:

        def __init__(self, Coupling):
            self.Coupling = Coupling
            self.Nsize = np.shape(self.Coupling)[0]
            self.ExcitatoryWeight = np.array([np.array([pos_weight for pos_weight in self.Coupling[idx][self.Coupling[idx] > 0]]) for idx in range(self.Nsize)], dtype=object)
            self.InhibitoryWeight = np.array([np.array([neg_weight for neg_weight in self.Coupling[idx][self.Coupling[idx] < 0]]) for idx in range(self.Nsize)], dtype=object)
            self.ExcitatoryWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() > 0])
            self.InhibitoryWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() < 0])

        def CouplingWeightBound(self, print_console=True):
            min_coupling_weight = np.amin(self.Coupling)
            max_coupling_weight = np.amax(self.Coupling)
            if print_console == True:
                print('Minimum coupling weight: {}, Maximum: {}'.format(min_coupling_weight, max_coupling_weight))
            return min_coupling_weight, max_coupling_weight

        def ConnectionProbability(self, print_console=True):
            number_of_links = len(np.argwhere(self.Coupling.flatten()))
            connection_probability = number_of_links / ( self.Nsize * (self.Nsize-1) )
            if print_console == True:
                print('Connection probability: {}'.format(connection_probability))
            return connection_probability

        def SynapticWeight_Stat(self, print_console=True):
            synaptic_weight_mean = np.mean(self.Coupling.flatten()[np.nonzero(self.Coupling.flatten())])
            synaptic_weight_sd = np.std(self.Coupling.flatten()[np.nonzero(self.Coupling.flatten())])
            if print_console == True:
                print('Mean of coupling weight: {}'.format(synaptic_weight_mean))
                print('SD of coupling weight: {}'.format(synaptic_weight_sd))
            return synaptic_weight_mean, synaptic_weight_sd

        def ExcitatorySynapticWeight_Stat(self, print_console=True):
            excitatory_weight_mean = np.mean(self.ExcitatoryWeight_flatten)
            excitatory_weight_sd = np.std(self.ExcitatoryWeight_flatten)
            if print_console == True:
                print('Mean of excitatory node weight: {}'.format(excitatory_weight_mean))
                print('SD of excitatory node weight: {}'.format(excitatory_weight_sd))
            return excitatory_weight_mean, excitatory_weight_sd

        def InhibitorySynapticWeight_Stat(self, print_console=True):
            inhibitory_weight_mean = np.mean(self.InhibitoryWeight_flatten)
            inhibitory_weight_sd = np.std(self.InhibitoryWeight_flatten)
            if print_console == True:
                print('Mean of inhibitory node weight: {}'.format(inhibitory_weight_mean))
                print('SD of inhibitory node weight: {}'.format(inhibitory_weight_sd))
            return inhibitory_weight_mean, inhibitory_weight_sd

        def AverageSynapticWeight_ExcitatoryIncomingLinks(self, idx, print_console=True):
            # index of node starts from 1, end with 4095
            # s+_in(i)
            avg_syn_weight_exc_in = np.sum(self.ExcitatoryWeight[idx-1]) / len(self.ExcitatoryWeight[idx-1])
            if print_console == True:
                print('Average synaptic weight of excitatory incoming links s+_in({}): {}'.format(idx, avg_syn_weight_exc_in))
            return avg_syn_weight_exc_in

        def AverageSynapticWeight_InhibitoryIncomingLinks(self, idx, print_console=True):
            # s-_in(i)
            avg_syn_weight_inh_in = np.sum(self.InhibitoryWeight[idx-1]) / len(self.InhibitoryWeight[idx-1])
            if print_console == True:
                print('Average synaptic weight of inhibitory incoming links |s-_in({})|: {}'.format(idx, avg_syn_weight_inh_in))
            return avg_syn_weight_inh_in
        
        def AverageSynapticWeight_OutgoingLinks(self, idx, print_console=True):
            # s_out(i)
            avg_syn_weight_out = np.sum(self.Coupling.transpose()[idx-1]) / len(np.argwhere(self.Coupling.transpose()[idx-1]))
            if print_console == True:
                print('Average synaptic weight of outgoing links |s_out({})|: {}'.format(idx, avg_syn_weight_out))
            return avg_syn_weight_out

        def RatioOfSuppression(self, suppressed_values):
            NewWeight = np.zeros((len(suppressed_values), len(self.InhibitoryWeight_flatten)))
            avg_change_inh_weight = np.zeros(len(suppressed_values))
            for count in range(0, len(suppressed_values)):
                NewWeight[count] += self.InhibitoryWeight_flatten + suppressed_values[count]
                NewWeight[count][NewWeight[count] > 0] = 0
                avg_change_inh_weight[count] = np.mean(np.array(NewWeight[count]) - self.InhibitoryWeight_flatten)
            avg_magnitude_inh_weight = np.mean(self.InhibitoryWeight_flatten)
            ratio_of_suppression = np.abs(avg_change_inh_weight) / np.abs(avg_magnitude_inh_weight)
            return ratio_of_suppression

        def RatioOfEnhancement(self, excited_values):
            avg_change_exc_weight = np.full(len(excited_values), excited_values)
            avg_magnitude_exc_weight = np.mean(self.ExcitatoryWeight_flatten)
            ratio_of_enhancement = np.abs(avg_change_exc_weight) / np.abs(avg_magnitude_exc_weight)
            return ratio_of_enhancement

    class Plot:

        def __init__(self, Coupling, output_path):
            self.Coupling = Coupling
            self.Nsize = np.shape(self.Coupling)[0]
            self.ExcitatoryWeight = np.array([np.array([pos_weight for pos_weight in self.Coupling[idx][self.Coupling[idx] > 0]]) for idx in range(self.Nsize)], dtype=object)
            self.InhibitoryWeight = np.array([np.array([neg_weight for neg_weight in self.Coupling[idx][self.Coupling[idx] < 0]]) for idx in range(self.Nsize)], dtype=object)
            # self.ExcitatoryWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() > 0])
            # self.InhibitoryWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() < 0])
            self.output_path = output_path

        def AverageSynapticWeight_ExcitatoryIncomingLinks(self, output_file='Average_Synaptic_Weight_Excitatory_Incoming_Links.svg', bins=100, show_norm=False):
            AvgSynWeight_ExcIn = np.array([np.sum(self.ExcitatoryWeight[idx]) / len(self.ExcitatoryWeight[idx]) for idx in range(self.Nsize)], dtype=float)
            density, bin_edges = np.histogram(AvgSynWeight_ExcIn, bins=bins, density=True)
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(x_value, density, 'b^-', lw=2, label='|s+_in| density')
            if show_norm == True:
                mu = np.mean(AvgSynWeight_ExcIn.flatten()); sigma = np.std(AvgSynWeight_ExcIn.flatten())
                norm = stats.norm(loc=mu, scale=sigma)
                ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')
            ax.set(xlabel='s+_in', ylabel='Probability density')
            ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()

        def AverageSynapticWeight_InhibitoryIncomingLinks(self, output_file='Average_Synaptic_Weight_Inhibitory_Incoming_Links.svg', bins=100, show_norm=False):
            AvgSynWeight_InhIn = -1 * np.array([np.sum(self.InhibitoryWeight[idx]) / len(self.InhibitoryWeight[idx]) for idx in range(self.Nsize)], dtype=float)
            density, bin_edges = np.histogram(AvgSynWeight_InhIn, bins=bins, density=True)
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(x_value, density, 'b^-', lw=2, label='|s-_in| density')
            if show_norm == True:
                mu = np.mean(AvgSynWeight_InhIn.flatten()); sigma = np.std(AvgSynWeight_InhIn.flatten())
                norm = stats.norm(loc=mu, scale=sigma)
                ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')
            ax.set(xlabel='|s-_in|', ylabel='Probability density')
            ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()

        def AverageSynapticWeight_OutgoingLinks(self, output_file='Average_Synaptic_Weight_Outgoing_Links.svg', bins=100, show_norm=False):
            AvgSynWeight_Out = np.array([np.sum(self.Coupling.transpose()[idx]) / len(np.argwhere(self.Coupling.transpose()[idx])) for idx in range(self.Nsize)], dtype=float)
            density, bin_edges = np.histogram(np.abs(AvgSynWeight_Out.flatten()), bins=bins, density=True)
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(x_value, density, 'b^-', lw=2, label='|s_out| density')
            if show_norm == True:
                mu = np.mean(np.abs(AvgSynWeight_Out.flatten())); sigma = np.std(np.abs(AvgSynWeight_Out.flatten()))
                norm = stats.norm(loc=mu, scale=sigma)
                ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')
            ax.set(xlabel='|s_out|', ylabel='Probability density')
            ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()


class Spiking:

    def __init__(self, spiking_data, simulation_config, delimiter='\t', input_folder=None, output_folder=None):
        self.initSpikeData(spiking_data, simulation_config, delimiter, input_folder, output_folder)
        self.calculate = self.Calculate(self.SpikeCount, self.SpikeTimes, self.Config, self.output_path)
        self.plot = self.Plot(self.SpikeCount, self.SpikeTimes, self.Config, self.output_path)

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
            reader = csv.reader(file_config, delimiter='\t')
            self.Config = np.array(list(reader), dtype=object)
        
        Nsize = int(self.Config[0][0])
        self.SpikeTimes = np.empty((Nsize), dtype=object)
        self.SpikeCount = np.zeros(Nsize)

        with open(input_path+spiking_data, 'r') as file_spike:
            reader = csv.reader(file_spike, delimiter=delimiter)
            counter = 0
            for row in reader:
                try: self.SpikeTimes[counter] = np.delete(np.array(list(row)).astype('float'), [0,1], 0)
                except ValueError: pass
                self.SpikeCount[counter] = int(row[1])
                counter += 1


    class Calculate:

        def __init__(self, SpikeCount, SpikeTimes, Config, output_path):
            self.SpikeCount = SpikeCount
            self.SpikeTimes = SpikeTimes
            self.Config = Config
            # self.Config[0][0] = N
            # self.Config[0][1] = dt
            # self.Config[0][2] = tn
            self.output_path = output_path

        def AverageFiringRate(self, console_print=True):
            spike_count_average = np.mean(self.SpikeCount)
            firing_rate_average = spike_count_average / float(self.Config[0][2]) * 1000
            if console_print == True:
                print('Average firing rate: {}'.format(firing_rate_average))
            return firing_rate_average

        def FiringRate_Stat(self, console_print=True):
            FiringRate = self.SpikeCount / float(self.Config[0][2]) * 1000
            log_mean = np.mean(FiringRate.flatten())
            log_sd = np.std(FiringRate.flatten())
            if console_print == True:
                print('Mean of Firing Rate: {}'.format(log_mean))
                print('SD of Firing Rate: {}'.format(log_sd))
            return log_mean, log_sd

        def InterSpikeInterval_Stat(self, console_print=True):
            IsI = np.empty(self.SpikeTimes.shape[0], dtype=object)
            count = 0
            for row in self.SpikeTimes:
                IsI[count] = np.diff(row)
                count += 1
            IsI = np.log10(np.concatenate([item for item in IsI.flatten()], 0) / 1000)
            IsI_log_mean = np.mean(IsI)
            IsI_log_sd = np.std(IsI)
            if console_print == True:
                print('Mean of LOG Inter-spike Interval: {}'.format(IsI_log_mean))
                print('SD of LOG Inter-spike Interval: {}'.format(IsI_log_sd))
            return IsI_log_mean, IsI_log_sd

        def identifyBurstingNode(self, output=False, output_file='Bursting_Nodes.txt'):
            BurstingNode = []
            IsI = np.empty(self.SpikeTimes.shape[0], dtype=object)
            count = 0
            for row in self.SpikeTimes:
                IsI[count] = np.diff(row)
                count += 1
            for idx in range(int(self.Config[0][0])):
                if len(IsI[idx]) != 0 and len(IsI[idx][IsI[idx] < 10]) / len(IsI[idx]) > 0.5:
                    BurstingNode.append(idx)
            if output == True:
                with open(self.output_path+output_file, 'w') as file_output:
                    for idx in BurstingNode:
                        file_output.write('{}\n'.format(idx))
            return BurstingNode


    class Plot:

        def __init__(self, SpikeCount, SpikeTimes, Config, output_path):
            self.SpikeCount = SpikeCount
            self.SpikeTimes = SpikeTimes
            self.Config = Config
            self.Nsize = int(self.Config[0][0])
            self.output_path = output_path

        def reformatSpikeData(self, output_file='OUT_SPIK_reformatted.txt'):
            with open(self.output_path+output_file, 'w') as file_output:
                for idx in range(self.Nsize):
                    file_output.write('{:.0f}'.format(self.SpikeCount[idx]))
                    for timestamp in self.SpikeTimes[idx]:
                        file_output.write('\t{:.0f}'.format(timestamp/float(self.Config[0][1])))
                    file_output.write('\n')

        def SpikeRaster(self, output_file='Spiking_Raster_Plot.png'):
            start_node, end_node = 0, self.Nsize

            fig, ax = plt.subplots(figsize=(9, 6), dpi=250)
            count = 0
            for each_node in (self.SpikeTimes / 1000):
                count += 1
                ax.scatter(each_node, np.full(np.size(each_node), count), s=0.5, c='black')
            ax.set(xlabel='Time (s)', ylabel='Node index')
            ax.set_xlim(0, float(self.Config[0][2])/1000)
            ax.set_ylim(start_node-2, end_node+1)
            ax.grid(True)
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()
        
        def FiringRateDistribution(self, bins=[0, 80, 1000], xrange=[0,20], yrange=[0,None], show_norm=False, output_file='Firing_Rate_Distribution', file_label=''):
            FiringRate = self.SpikeCount / float(self.Config[0][2]) * 1000
            if np.amax(FiringRate) > bins[1]: print('Warning! Maximum of Firing Rate exceeds upper bound of bins range. Max Firing Rate: {}; Max bins range: {}'.format(np.amax(FiringRate), bins[1]))
            if np.amin(FiringRate) < bins[0]: print('Warning! Minimum of Firing Rate subceeds lower bound of bins range. Min Firing Rate: {}; Min bins range: {}'.format(np.amin(FiringRate), bins[0]))
            hist_density, bin_edges = np.histogram(FiringRate, bins=np.linspace(bins[0], bins[1], bins[2]), density=True)
            # bin1 = np.linspace(0, 1.5, 8); bin2 = np.linspace(1.5, 70, 150)[1:]
            # hist_density, bin_edges = np.histogram(FiringRate, bins=np.concatenate([bin1, bin2]), density=True)
            total_density = np.dot(hist_density, np.diff(bin_edges))
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(x_value, hist_density, 'b^-', lw=2, label='Firing rate distribution')
            if show_norm == True:
                mu = np.mean(FiringRate.flatten()); sigma = np.std(FiringRate.flatten())
                norm = stats.norm(loc=mu, scale=sigma)
                ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')
            ax.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
            ax.set_xlim(xrange[0], xrange[1]); ax.set_ylim(yrange[0], yrange[1])
            ax.grid(True)
            ax.legend()

            if file_label == '': output_file_plot = output_file + '.svg'
            else: output_file_plot = output_file + '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()
            if file_label == '': output_file_info = output_file + '.txt'
            else: output_file_info = output_file + '_' + file_label + '.txt'
            with open(self.output_path+output_file_info, 'w') as fp_info:
                fp_info.write('### Plot information ###\n')
                fp_info.write('Max firing rate: {}\n'.format(np.amax(FiringRate)))
                fp_info.write('Min firing rate: {}\n'.format(np.amin(FiringRate)))
                fp_info.write('Total density by areal summation: {}\n'.format(total_density))
                fp_info.write('\n### Plot settings ###\n')
                fp_info.write('Bin size: {}\n'.format((bins[1]-bins[0])/bins[2]))
                fp_info.write('Bin bounds: lower: {}, upper: {}\n'.format(bins[0], bins[1]))
                fp_info.write('Number of bins: {}\n'.format(bins[2]))
                fp_info.write('Show Gaussian distribution: {}'.format(str(show_norm)))

        def InterSpikeIntervalDistribution(self, bins=[0.001, 10, 150], xrange=[0.001,1], output_file='Interspike_Interval_Distribution', file_label=''):
            IsI = np.empty(self.SpikeTimes.shape[0], dtype=object)
            for count in range(len(self.SpikeTimes)):
                try: IsI[count] = np.diff(self.SpikeTimes[count])
                except ValueError: IsI[count] = np.diff(np.array([0]))
            IsI = np.concatenate([item for item in IsI.flatten()], 0) / 1000
            if np.amax(IsI) > bins[1]: print('Warning! Maximum of ISI exceeds upper bound of bins range. Max ISI: {}; Max bins range: {}'.format(np.amax(IsI), bins[1]))
            if np.amin(IsI) < bins[0]: print('Warning! Minimum of ISI subceeds lower bound of bins range. Min ISI: {}; Min bins range: {}'.format(np.amin(IsI), bins[0]))
            bins[0] = math.log10(bins[0]); bins[1] = math.log10(bins[1]); xrange[0] = math.log10(xrange[0]); xrange[1] = math.log10(xrange[1])
            hist_density, bin_edges = np.histogram(IsI, bins=np.logspace(bins[0], bins[1], bins[2]))
            hist_density = np.array(hist_density, dtype=float)
            hist_density /= np.dot(hist_density, np.diff(bin_edges)) # normalization
            total_density = np.dot(hist_density, np.diff(bin_edges))
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2
            
            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.semilogx(x_value, hist_density, 'b^-', lw=2, label='Log ISI distribution')
            ax.set(xlabel='ISI (s)', ylabel='Probability density')
            # ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.grid(True)
            # ax.legend()
            
            if file_label == '': output_file_plot = output_file + '.svg'
            else: output_file_plot = output_file + '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()
            if file_label == '': output_file_info = output_file + '.txt'
            else: output_file_info = output_file + '_' + file_label + '.txt'
            with open(self.output_path+output_file_info, 'w') as fp_info:
                fp_info.write('### Plot information ###\n')
                fp_info.write('Max ISI: {}\n'.format(np.amax(IsI)))
                fp_info.write('Min ISI: {}\n'.format(np.amin(IsI)))
                fp_info.write('Total density by areal summation: {}\n'.format(total_density))
                fp_info.write('\n### Plot settings ###\n')
                fp_info.write('Bin size in log scale: {}\n'.format((bins[1]-bins[0])/bins[2]))
                fp_info.write('Bin bounds: lower: {}, upper: {}\n'.format(math.pow(10,bins[0]), math.pow(10,bins[1])))
                fp_info.write('Number of bins: {}'.format(bins[2]))


class Compare:

    def __init__(self, spiking_data: list, data_info: list, simulation_config: str, delimiter='\t', input_folder=None, output_folder=None, coupling=None, coupling_enhance_factor=1, coupling_delimiter='\t', simulataion_duration_override=[]):
        self.initSpikeData(spiking_data, simulation_config, delimiter, input_folder, output_folder)
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

        def SuppressionRatio_vs_FiringRateIncreaseRatio(self, suppressed_values: list, plot_label='', output_file='Ratio_of_Suppression_vs_Ratio_of_Increase_in_Firing_Rate', file_label=''):
            if self.coupling == None: print('Error. No coupling matrix input.'); return -1
            ratio_of_suppression = Coupling(self.coupling, coupling_enhance_factor=self.coupling_enhance_factor, delimiter=self.coupling_delimiter).calculate.RatioOfSuppression(suppressed_values)
            avg_firing_rate_original_network = Spiking(self.spiking_data[0], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            avg_firing_rate_altered_networks = np.zeros(len(suppressed_values))
            for count in range(0, len(suppressed_values)):
                avg_firing_rate_altered_networks[count] = Spiking(self.spiking_data[count+1], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            ratio_of_increase_in_firing_rate = avg_firing_rate_altered_networks / avg_firing_rate_original_network - 1
            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(ratio_of_suppression, ratio_of_increase_in_firing_rate, 'b^-', lw=2, label=plot_label)
            ax.set(xlabel='Ratio of suppression in inhibitory synaptic weights', ylabel='Ratio of increase in avergae firing rate') 
            ax.grid(True)
            ax.legend()
            if file_label == '': output_file_plot = output_file + '.svg'
            else: output_file_plot = output_file + '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()
        
        def EnhancementRatio_vs_FiringRateIncreaseRatio(self, enhanced_values: list, plot_label = '', output_file='Ratio_of_Enhancement_vs_Ratio_of_Increase_in_Firing_Rate', file_label=''):
            if self.coupling == None: print('Error. No coupling matrix input.'); return -1
            ratio_of_enhancement = Coupling(self.coupling, coupling_enhance_factor=self.coupling_enhance_factor, delimiter=self.coupling_delimiter).calculate.RatioOfEnhancement(enhanced_values)
            avg_firing_rate_original_network = Spiking(self.spiking_data[0], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            avg_firing_rate_altered_networks = np.zeros(len(enhanced_values))
            for count in range(0, len(enhanced_values)):
                avg_firing_rate_altered_networks[count] = Spiking(self.spiking_data[count+1], self.simulation_config).calculate.AverageFiringRate(console_print=False)
            ratio_of_increase_in_firing_rate = avg_firing_rate_altered_networks / avg_firing_rate_original_network - 1
            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(ratio_of_enhancement, ratio_of_increase_in_firing_rate, 'b^-', lw=2, label=plot_label)
            ax.set(xlabel='Ratio of enhancement in excitatory synaptic weights', ylabel='Ratio of increase in avergae firing rate') 
            ax.grid(True)
            ax.legend()
            if file_label == '': output_file_plot = output_file + '.svg'
            else: output_file_plot = output_file + '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()


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

def EnhancementRatio_vs_FiringRateIncreaseRatio_Combined(spiking_data_1, spiking_data_2, simulation_config, coupling_1, coupling_2, enhanced_values_1: list, enhanced_values_2: list, coupling_enhance_factor_1, coupling_enhance_factor_2, coupling_delimiter_1, coupling_delimiter_2, label_1='', label_2='', output_file='Enhancement_Ratio_vs_Firing_Rate_Increase_Ratio', file_label=''):
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