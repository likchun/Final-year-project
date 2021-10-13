import os, csv, math
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.style as mplstyle
import seaborn as sns


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

class Spiking(Initialization):

    def __init__(self, spiking_data, simulation_config, delimiter='\t', input_folder=None, output_folder=None):
        super().__init__(spiking_data, simulation_config, delimiter, input_folder, output_folder)
        self.calculate = self.Calculate(self.SpikeCount, self.SpikeTimes, self.Config, self.output_path)
        self.plot = self.Plot(self.SpikeCount, self.SpikeTimes, self.Config, self.output_path)

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