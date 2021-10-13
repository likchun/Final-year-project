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
        
        def FiringRateDensity(self, output_file='Firing_Rate_Density.svg', bins=60, xrange=(0,10), yrange=(0,None), show_norm=False):
            FiringRate = self.SpikeCount / float(self.Config[0][2]) * 1000
            density, bin_edges = np.histogram(FiringRate, bins=np.linspace(xrange[0], xrange[1], bins), density=True)
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.plot(x_value, density, 'b^-', lw=2, label='Firing rate density')
            if show_norm == True:
                mu = np.mean(FiringRate.flatten()); sigma = np.std(FiringRate.flatten())
                norm = stats.norm(loc=mu, scale=sigma)
                ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')
            ax.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
            ax.set_xlim(xrange[0], xrange[1])
            ax.set_ylim(yrange[0], yrange[1])
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()

        def InterSpikeIntervalDensity(self, output_file='Interspike_Interval_Density.svg', bins=150, xrange=(-3,1)):
            IsI = np.empty(self.SpikeTimes.shape[0], dtype=object)
            count = 0
            for row in self.SpikeTimes:
                try: IsI[count] = np.diff(row)
                except ValueError: IsI[count] = np.diff(np.array([0]))
                count += 1
            IsI = np.concatenate([item for item in IsI.flatten()], 0) / 1000
            density, bin_edges = np.histogram(IsI, bins=np.logspace(xrange[0], xrange[1], bins))
            density = np.array(density, dtype=float)
            density /= np.dot(density, np.diff(bin_edges)) # normalization
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2
            
            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            ax.semilogx(x_value, density, 'b^-', lw=2, label='Log ISI density')
            ax.set(xlabel='ISI (s)', ylabel='Probability density')
            # ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.grid(True)
            ax.legend()
            fig.savefig(os.path.join(self.output_path, output_file))
            plt.clf()


###############
### Spiking ###
###############

input_folder = [
    # 'OUT_DATA_7500_01',
    # 'OUT_DATA_7500_005',
    # 'OUT_DATA_7500_05',
    # 'OUT_DATA_7500_0025',
    # 'OUT_DATA_7500_25',
    # 'OUT_DATA_7500_125',
    # 'OUT_DATA_10000_01',
    # 'OUT_DATA_10000_05',
    # 'OUT_DATA_10000_125',
    # 'OUT_DATA_15000_01',
    # 'OUT_DATA_15000_05',
    # 'OUT_DATA_15000_125',
    # 'OUT_DATA_DIV66',
    # 'OUT_DATA_DIV66_enhance_exc_k1',
    # 'OUT_DATA_DIV66_enhance_exc_k075',
    # 'OUT_DATA_DIV66_enhance_exc_k05',
    # 'OUT_DATA_DIV66_enhance_exc_k025',
    # 'OUT_DATA_DIV66_suppress_inh_k1',
    'OUT_DATA_DIV66_suppress_inh_k075',
    # 'OUT_DATA_DIV66_suppress_inh_k05',
    # 'OUT_DATA_DIV66_suppress_inh_k025',
    # 'OUT_DATA_RANDOM',
    # 'OUT_DATA_RANDOM_suppress_inh_k1',
    # 'OUT_DATA_RANDOM_suppress_inh_k075',
    # 'OUT_DATA_RANDOM_suppress_inh_k05',
    # 'OUT_DATA_RANDOM_suppress_inh_k025',
    # 'OUT_DATA_RANDOM_enhance_exc_k1',
    # 'OUT_DATA_RANDOM_enhance_exc_k075',
    # 'OUT_DATA_RANDOM_enhance_exc_k05',
    # 'OUT_DATA_RANDOM_enhance_exc_k025',
]

for each_path in input_folder:
    s = Spiking(each_path+'\\OUT_SPIK.txt', each_path+'\\INI_CNFG', output_folder=each_path)
    s_calc = s.calculate
    s_plot = s.plot

    # s_calc.AverageFiringRate()
    # s_calc.FiringRate_Stat()
    # s_calc.InterSpikeInterval_Stat()

    # s_calc.identifyBurstingNode(output=True)

    # s_plot.reformatSpikeData()
    # s_plot.SpikeRaster()
    # s_plot.FiringRateDensity()
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_fit.svg', bins=78, xrange=(0,10), yrange=(0, 1.3)) #for DIV 7500
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_fit.svg', bins=105, xrange=(0,10), yrange=(0, 1.3)) #for DIV 10000
    s_plot.FiringRateDensity(output_file='Firing_Rate_Density_fit_0.svg', bins=150, xrange=(0,16), yrange=(0, 1.3)) #for DIV enhance
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_fit.svg', bins=150, xrange=(0,10), yrange=(0, 1.3)) #for DIV 15000
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_0.svg', bins=100, xrange=(0,5))
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_1.svg', bins=80, xrange=(0,5))
    # s_plot.FiringRateDensity(bins=42, xrange=(0,4), show_norm=True) # for random network D
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_fit.svg', bins=150, xrange=(0,100), show_norm=True) # for random network D enhance exc
    # s_plot.FiringRateDensity(output_file='Firing_Rate_Density_lowfit.svg', bins=30, xrange=(0,100), show_norm=True) # for random network D enhance exc
    # s_plot.InterSpikeIntervalDensity(output_file='Interspike_Interval_Density_0.svg', bins=400)
    # s_plot.InterSpikeIntervalDensity(output_file='Interspike_Interval_Density_fit.svg', xrange=(-4,1), bins=150)
