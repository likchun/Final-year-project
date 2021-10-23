import os, csv, math
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.style as mplstyle
import seaborn as sns


class Coupling:
    """
    Analysis and graphings for coupling matrix.

    Parameters
    ----------
    coupling_matrix_path : str
        Path of the input coupling matrix file, e.g., DIV66.txt, exclude input folder if `input_folder` is used.
    coupling_enhance_factor : int or float, default: 1
        Multiply all synaptic weights by a factor.
    delimiter : str, default: '\\t'
        Delimiter to the input coupling matrix file.
    input_folder : str, default: None
        Path of input folder.
    output_folder : str, default: None
        Path of output folder.

    Notes
    -----
    Invoke object `calculate` to find basic features and properties of the coupling matrix.\n
    Invoke object `plot` to plot graphs of basic features and properties of the coupling matrix.

    Examples
    --------
    To start with, create a object for class `Coupling` by specifying the file path and delimiter (optional).
    Invoke object `calculate` / `plot` to find / plot basic features and properties of the coupling matrix.
    >>> c = Coupling('DIV66.txt', delimiter='\\t')
    >>> c_calc = c.calculate
    >>> c_plot = c.plot

    Find the *connection probability* of the coupling matrix:
    >>> c = Coupling('DIV66.txt', delimiter='\\t')
    >>> c.calculate.ConnectionProbability()
    'terminal> Connection probability: 0.015344'
    >>> result = c.calculate.ConnectionProbability(print_console=False)
    >>> print(result)
    'terminal> 0.015344'

    Plot the *average synaptic weight excitatory incoming links distribution* of the coupling matrix:
    >>> c = Coupling('DIV66.txt', delimiter='\\t')
    >>> c.plot.AverageSynapticWeight_ExcitatoryIncomingLinks()
    """

    def __init__(self, coupling_matrix_path: str, coupling_enhance_factor=1, delimiter='\t', input_folder=None, output_folder=None):
        self.initCouplingMatrix(coupling_matrix_path, delimiter, input_folder, output_folder)
        self.Coupling = coupling_enhance_factor * self.Coupling
        self.calculate = self.Calculate(self.Coupling)
        self.plot = self.Plot(self.Coupling, self.output_path)

    def initCouplingMatrix(self, coupling_matrix_path, delimiter, input_folder, output_folder):
        # Initialize coupling matrix
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if input_folder == None: input_folder = ''
        else: input_folder = input_folder + '\\'
        input_path = this_dir + '\\' + input_folder
        if output_folder == None: output_folder = ''
        else: output_folder = output_folder + '\\'
        self.output_path = this_dir + '\\' + output_folder
        try: os.mkdir(self.output_path)
        except FileExistsError: pass

        with open(input_path+coupling_matrix_path, 'r', newline='') as file_input:
            reader = csv.reader(file_input, delimiter=delimiter)
            self.Coupling = np.array(list(reader)).astype('float')

    class Calculate:
        """
        Find basic features and properties of the coupling matrix.
        """

        def __init__(self, Coupling):
            # Initialize instance attributes
            self.Coupling = Coupling
            self.Nsize = np.shape(self.Coupling)[0]
            self.PositiveWeight = np.array([np.array([pos_weight for pos_weight in self.Coupling[idx][self.Coupling[idx] > 0]]) for idx in range(self.Nsize)], dtype=object)
            self.NegativeWeight = np.array([np.array([neg_weight for neg_weight in self.Coupling[idx][self.Coupling[idx] < 0]]) for idx in range(self.Nsize)], dtype=object)
            self.PositiveWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() > 0])
            self.NegativeWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() < 0])

        def ConnectionProbability(self, print_console=True):
            """
            Find the connection probability.

            Parameters
            ----------
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            float
                Connection probability of the network.
            """
            number_of_links = len(np.argwhere(self.Coupling.flatten()))
            connection_probability = number_of_links / ( self.Nsize * (self.Nsize-1) )
            if print_console == True:
                print('Connection probability:\t{}'.format(connection_probability))
            return connection_probability
            
        def SynapticWeightBounds(self, print_console=True):
            """
            Find the minimum and maximum synaptic weights.

            Parameters
            ----------
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            min_synaptic_weight : float
                Minimum synaptic weight of the network.
            max_synaptic_weight : float
                Maximum synaptic weight of the network.
            """
            min_synaptic_weight = np.amin(self.Coupling)
            max_synaptic_weight = np.amax(self.Coupling)
            if print_console == True:
                print('Minimum synaptic weight:\t{}\nMaximum synaptic weight:\t{}'.format(min_synaptic_weight, max_synaptic_weight))
            return min_synaptic_weight, max_synaptic_weight

        def SynapticWeight_Stat(self, print_console=True):
            """
            Find the mean and standard deviation of synaptic weights.

            Parameters
            ----------
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            synaptic_weight_mean : float
                Mean of synaptic weight of all links in the network.
            synaptic_weight_sd : float
                Population standard deviation of synaptic weight of all links in the network.
            """
            synaptic_weight_mean = np.mean(self.Coupling.flatten()[np.nonzero(self.Coupling.flatten())])
            synaptic_weight_sd = np.std(self.Coupling.flatten()[np.nonzero(self.Coupling.flatten())])
            if print_console == True:
                print('Mean of all synaptic weights:\t{}'.format(synaptic_weight_mean))
                print('S.D. of all synaptic weights:\t{}'.format(synaptic_weight_sd))
            return synaptic_weight_mean, synaptic_weight_sd

        def PositiveSynapticWeight_Stat(self, print_console=True):
            """
            Find the mean and population standard deviation of positive synaptic weights.

            Parameters
            ----------
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            positive_weight_mean : float
                Mean of synaptic weight of all positive links of the network.
            positive_weight_sd : float
                Standard deviation of synaptic weight of all positive links in the network.
            """
            positive_weight_mean = np.mean(self.PositiveWeight_flatten)
            positive_weight_sd = np.std(self.PositiveWeight_flatten)
            if print_console == True:
                print('Mean of synaptic weight of positive links:\t{}'.format(positive_weight_mean))
                print('S.D. of synaptic weight of positive links:\t{}'.format(positive_weight_sd))
            return positive_weight_mean, positive_weight_sd

        def NegativeSynapticWeight_Stat(self, print_console=True):
            """
            Find the mean and population standard deviation of negative synaptic weights.

            Parameters
            ----------
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            negative_weight_mean : float
                Mean of synaptic weight of all negative links of the network.
            negative_weight_sd : float
                Standard deviation of synaptic weight of all negative links in the network.
            """
            negative_weight_mean = np.mean(self.NegativeWeight_flatten)
            negative_weight_sd = np.std(self.NegativeWeight_flatten)
            if print_console == True:
                print('Mean of synaptic weight of negative links:\t{}'.format(negative_weight_mean))
                print('S.D. of synaptic weight of negative links:\t{}'.format(negative_weight_sd))
            return negative_weight_mean, negative_weight_sd

        def AverageSynapticWeight_ExcitatoryIncomingLinks(self, idx, print_console=True):
            """
            Find the average synaptic weight of excitatory incoming links of a specific node, s+_in(i).

            Parameters
            ----------
            idx : int
                The index of node to be calculated, starts from 1.
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            float
                Average synaptic weight of excitatory incoming links of node `idx`.
            """
            # index of node starts from 1, end with 4095
            # s+_in(i)
            avg_syn_weight_exc_in = np.sum(self.PositiveWeight[idx-1]) / len(self.PositiveWeight[idx-1])
            if print_console == True:
                print('Average synaptic weight of excitatory incoming links s+_in({}): {}'.format(idx, avg_syn_weight_exc_in))
            return avg_syn_weight_exc_in

        def AverageSynapticWeight_InhibitoryIncomingLinks(self, idx, print_console=True):
            """
            Find the average synaptic weight of inhibitory incoming links of a specific node, s-_in(i).

            Parameters
            ----------
            idx : int
                The index of node to be calculated, starts from 1.
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            float
                Average synaptic weight of inhibitory incoming links of node `idx`.
            """
            # s-_in(i)
            avg_syn_weight_inh_in = np.sum(self.NegativeWeight[idx-1]) / len(self.NegativeWeight[idx-1])
            if print_console == True:
                print('Average synaptic weight of inhibitory incoming links |s-_in({})|: {}'.format(idx, avg_syn_weight_inh_in))
            return avg_syn_weight_inh_in
        
        def AverageSynapticWeight_OutgoingLinks(self, idx, print_console=True):
            """
            Find the average synaptic weight of outgoing links of a specific node, s_out(i).

            Parameters
            ----------
            idx : int
                The index of node to be calculated, starts from 1.
            print_console : bool
                Print result in terminal if `True`.

            Returns
            -------
            float
                Average synaptic weight of outgoing links of node `idx`.
            """
            # s_out(i)
            avg_syn_weight_out = np.sum(self.Coupling.transpose()[idx-1]) / len(np.argwhere(self.Coupling.transpose()[idx-1]))
            if print_console == True:
                print('Average synaptic weight of outgoing links |s_out({})|: {}'.format(idx, avg_syn_weight_out))
            return avg_syn_weight_out

        def RatioOfSuppression(self, suppressed_node_type: str, suppressed_values: list):
            """
            Calculate the ratio of suppression for each suppressed value.

            Parameters
            ----------
            suppressed_node_type : {'inhibitory', 'excitatory'}
                The node type, inhibitory or excitatory, to be suppressed.
            suppressed_values : list of float
                A list of suppressed value.

            Returns
            -------
            list of float
                A list of suppression ratios.
            
            Notes
            -----
            Suppressed value = (suppression level k) * (S.D. of the suppressed node type).

            Suppressing *inhibitory*:
                suppressed value added to the synaptic weight of every *negative* links, which cannot be greater than zero.

            Suppressing *excitatory*:
                suppressed value subtracted from the synaptic weight of every *positive* links, which cannot be smaller than zero.
            """
            avg_change_weight = np.zeros(len(suppressed_values))

            if suppressed_node_type == 'inhibitory':
                avg_magnitude_weight = np.mean(self.NegativeWeight_flatten)
                NewWeight = np.zeros((len(suppressed_values), len(self.NegativeWeight_flatten)))
                for count in range(0, len(suppressed_values)):
                    NewWeight[count] += self.NegativeWeight_flatten + suppressed_values[count]
                    NewWeight[count][NewWeight[count] > 0] = 0
                    avg_change_weight[count] = np.mean(np.array(NewWeight[count]) - self.NegativeWeight_flatten)

            elif suppressed_node_type == 'excitatory':
                avg_magnitude_weight = np.mean(self.PositiveWeight_flatten)
                NewWeight = np.zeros((len(suppressed_values), len(self.PositiveWeight_flatten)))
                for count in range(0, len(suppressed_values)):
                    NewWeight[count] += self.PositiveWeight_flatten - suppressed_values[count]
                    NewWeight[count][NewWeight[count] < 0] = 0
                    avg_change_weight[count] = np.mean(np.array(NewWeight[count]) - self.PositiveWeight_flatten)
            
            else:
                print('Error. the valid input for argument \'suppressed_node_type\' is either \'inhibitory\' or \'excitatory\'.')
            
            ratio_of_suppression = np.abs(avg_change_weight) / np.abs(avg_magnitude_weight)
            return ratio_of_suppression

        def RatioOfEnhancement(self, enhanced_node_type: str, enhanced_values: list):
            """
            Calculate the ratio of enhancement for each enhanced value.

            Parameters
            ----------
            enhanced_node_type : {'inhibitory', 'excitatory'}
                The node type, inhibitory or excitatory, to be enhanced.
            enhanced_values : list of float
                A list of enhanced value.

            Returns
            -------
            list of float
                A list of enhancement ratios.
            
            Notes
            -----
            Enhanced value = (enhancement level k) * (S.D. of the enhanced node type).

            Enhancing *inhibitory*:
                enhanced value subtracted from the synaptic weight of every *negative* links.

            Enhancing *excitatory*:
                enhanced value added to the synaptic weight of every *positive* links.
            """
            avg_change_weight = np.full(len(enhanced_values), enhanced_values)
            if enhanced_node_type == 'inhibitory': avg_magnitude_weight = np.mean(self.NegativeWeight_flatten)
            if enhanced_node_type == 'excitatory': avg_magnitude_weight = np.mean(self.PositiveWeight_flatten)
            ratio_of_enhancement = np.abs(avg_change_weight) / np.abs(avg_magnitude_weight)
            return ratio_of_enhancement

    class Plot:
        """
        Plot graphs of basic features and properties of the coupling matrix.
        """

        def __init__(self, Coupling, output_path):
            self.Coupling = Coupling
            self.Nsize = np.shape(self.Coupling)[0]
            self.PositiveWeight = np.array([np.array([pos_weight for pos_weight in self.Coupling[idx][self.Coupling[idx] > 0]]) for idx in range(self.Nsize)], dtype=object)
            self.NegativeWeight = np.array([np.array([neg_weight for neg_weight in self.Coupling[idx][self.Coupling[idx] < 0]]) for idx in range(self.Nsize)], dtype=object)
            # self.PositiveWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() > 0])
            # self.NegativeWeight_flatten = np.array(self.Coupling.flatten()[self.Coupling.flatten() < 0])
            self.output_path = output_path

        def AverageSynapticWeight_ExcitatoryIncomingLinks(self, output_file='Average_Synaptic_Weight_Excitatory_Incoming_Links.svg', bins=100, show_norm=False):
            AvgSynWeight_ExcIn = np.array([np.sum(self.PositiveWeight[idx]) / len(self.PositiveWeight[idx]) for idx in range(self.Nsize)], dtype=float)
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
            AvgSynWeight_InhIn = -1 * np.array([np.sum(self.NegativeWeight[idx]) / len(self.NegativeWeight[idx]) for idx in range(self.Nsize)], dtype=float)
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
    """
    Analysis and graphings for spiking data.

    Parameters
    ----------
    spiking_data_path : str
        Path of the input spiking data file, OUT_SPIK.txt, exclude input folder if `input_folder` is used.
    config_path : str
        Path of the input spiking data config file, INI_CNFG, exclude input folder if `input_folder` is used.
    delimiter : str, default: '\\t'
        Delimiter to the input spiking data file.
    input_folder : str, default: None
        Path of input folder.
    output_folder : str, default: None
        Path of output folder.

    Notes
    -----
    Invoke object `calculate` to find basic features and properties of the spiking data.\n
    Invoke object `plot` to plot graphs of basic features and properties of the spiking data.

    Examples
    --------
    To start with, create a object for class `Spiking` by specifying the spiking data file path and the configuration file path.
    Invoke object `calculate` / `plot` to find / plot basic features and properties of the spiking data.
    >>> s = Spiking('INPUT_FOLDER\\\OUT_SPIK.txt, INPUT_FOLDER\\\INI_CNFG')
    >>> s_calc = s.calculate
    >>> s_plot = s.plot

    Find the *average firing rate* of the neuron dynamics from spiking data:
    >>> s = Spiking('INPUT_FOLDER\\\OUT_SPIK.txt, INPUT_FOLDER\\\INI_CNFG')
    >>> s.calculate.AverageFiringRate()
    'terminal> Average firing rate: 2.1246'
    >>> result = s.calculate.AverageFiringRate(print_console=False)
    >>> print(result)
    'terminal> 2.1246'

    Reformat spiking data, reformat spike times from *ms* to *timestep*:
    >>> s = Spiking('INPUT_FOLDER\\\OUT_SPIK.txt, INPUT_FOLDER\\\INI_CNFG')
    >>> s.calculate.reformatSpikeData_timesteps()

    Plot the *inter-spike interval distribution* of the neuron dynamics from spiking data:
    >>> s = Spiking('INPUT_FOLDER\\\OUT_SPIK.txt, INPUT_FOLDER\\\INI_CNFG')
    >>> s.plot.InterSpikeIntervalDistribution(bins=[0.0005,50,180], xrange=[0.0005,10])
    """

    def __init__(self, spiking_data_path, config_path, delimiter='\t', input_folder=None, output_folder=None):
        self.initSpikeData(spiking_data_path, config_path, delimiter, input_folder, output_folder)
        self.calculate = self.Calculate(self.SpikeCount, self.SpikeTimes, self.Config, self.output_path)
        self.plot = self.Plot(self.SpikeCount, self.SpikeTimes, self.Config, self.output_path)

    def initSpikeData(self, spiking_data_path, config_path, delimiter, input_folder, output_folder):
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if input_folder == None: input_folder = ''
        else: input_folder = input_folder + '\\'
        input_path = this_dir + '\\' + input_folder
        if output_folder == None: output_folder = ''
        else: output_folder = output_folder + '\\'
        self.output_path = this_dir + '\\' + output_folder
        try: os.mkdir(self.output_path)
        except FileExistsError: pass

        with open(input_path+config_path, 'r') as file_config:
            reader = csv.reader(file_config, delimiter='\t')
            self.Config = np.array(list(reader), dtype=object)
        
        Nsize = int(self.Config[0][0])
        self.SpikeTimes = np.empty((Nsize), dtype=object)
        self.SpikeCount = np.zeros(Nsize)

        with open(input_path+spiking_data_path, 'r') as file_spike:
            reader = csv.reader(file_spike, delimiter=delimiter)
            counter = 0
            for row in reader:
                try: self.SpikeTimes[counter] = np.delete(np.array(list(row)).astype('float'), [0,1], 0)
                except ValueError: pass
                self.SpikeCount[counter] = int(row[1])
                counter += 1


    class Calculate:
        """
        Find basic features and properties of the spiking data.
        """

        def __init__(self, SpikeCount, SpikeTimes, Config, output_path):
            self.SpikeCount = SpikeCount
            self.SpikeTimes = SpikeTimes
            self.Nsize = int(Config[0][0])
            self.dt = float(Config[0][1])
            self.Tn = float(Config[0][2])
            self.output_path = output_path

        def SpikeCountBounds(self, console_print=True):
            max_spike_count = np.amax(self.SpikeCount)
            min_spike_count = np.amin(self.SpikeCount)
            if console_print == True:
                print('Maximum number of spike for nodes: {}'.format(max_spike_count))
                print('Minimum number of spike for nodes: {}'.format(min_spike_count))
            return min_spike_count, max_spike_count

        def TotalNumberOfSpikes(self, console_print=True):
            sum_spike_count = np.sum(self.SpikeCount)
            if console_print == True:
                print('Total number of spikes: {}'.format(sum_spike_count))
            return sum_spike_count

        def AverageFiringRate(self, console_print=True):
            average_spike_count = np.mean(self.SpikeCount)
            average_firing_rate = average_spike_count / float(self.Tn) * 1000
            if console_print == True:
                print('Average firing rate: {}'.format(average_firing_rate))
            return average_firing_rate

        def FiringRate_Stat(self, console_print=True):
            FiringRate = self.SpikeCount / float(self.Tn) * 1000
            log_mean = np.mean(FiringRate.flatten())
            log_sd = np.std(FiringRate.flatten())
            if console_print == True:
                print('Mean of Firing Rate: {}'.format(log_mean))
                print('SD of Firing Rate: {}'.format(log_sd))
            return log_mean, log_sd

        def InterSpikeInterval_Stat(self, console_print=True):
            IsI = np.empty(self.Nsize, dtype=object)
            for idx in range(len(self.SpikeTimes)):
                IsI[idx] = np.diff(self.SpikeTimes[idx])
            IsI = np.log10(np.concatenate([item for item in IsI.flatten()], 0) / 1000)
            IsI_log_mean = np.mean(IsI)
            IsI_log_sd = np.std(IsI)
            if console_print == True:
                print('Mean of LOG Inter-spike Interval: {}'.format(IsI_log_mean))
                print('SD of LOG Inter-spike Interval: {}'.format(IsI_log_sd))
            return IsI_log_mean, IsI_log_sd

        def identifyBurstingNode(self, output=False, file_name='Bursting_Nodes.txt'):
            BurstingNode = []
            IsI = np.empty(self.Nsize, dtype=object)
            for idx in range(len(self.SpikeTimes)):
                IsI[idx] = np.diff(self.SpikeTimes[idx])
            for idx in range(self.Nsize):
                if len(IsI[idx]) != 0 and len(IsI[idx][IsI[idx] < 10]) / len(IsI[idx]) > 0.5:
                    BurstingNode.append(idx)
            if output == True:
                with open(self.output_path+file_name, 'w') as file_output:
                    for idx in BurstingNode:
                        file_output.write('{}\n'.format(idx))
            return BurstingNode

        def reformatSpikeData_timesteps(self, file_name='OUT_SPIK_timesteps.txt'):
            """
            Reformat the data from spiking timestamp file, OUT_SPIK.txt.

            Remove index column and record timestamps *t*1, *t*2, ..., *tn* in simulation timesteps.

            Reformat as\:
            - N rows correspond to N spiking data for N nodes
            - Column 1: the number of spikes *n(i)* of the corresponding node *i*
            - Column 2 and onwards: *n(i)* timestamps *t*1, *t*2, ..., *tn* **(in timestep)** at which the spikes occur for that node

            Parameters
            ----------
            file_name : str, default: 'OUT_SPIK_reformatted.txt'
                Defines the output file name.

            Notes
            -----
            Raw data format\:
            - N rows correspond to N spiking data for N nodes
            - Column 1: the index of node *i*. (*i* runs from 1 to N)
            - Column 2: the number of spikes *n(i)* of the corresponding node *i*
            - Column 3 and onwards: *n(i)* timestamps *t*1, *t*2, ..., *tn* **(in ms)** at which the spikes occur for that node
            """
            with open(self.output_path+file_name, 'w') as file_output:
                for idx in range(self.Nsize):
                    file_output.write('{:.0f}'.format(self.SpikeCount[idx]))
                    for timestamp in self.SpikeTimes[idx]:
                        file_output.write('\t{:.0f}'.format(timestamp/float(self.dt)))
                    file_output.write('\n')
        
        def reformatSpikeData_noindex(self, file_name='OUT_SPIK_noindex.txt'):
            """
            Reformat the data from spiking timestamp file, OUT_SPIK.txt.

            Remove index column.

            Reformat as\:
            - N rows correspond to N spiking data for N nodes
            - Column 1: the number of spikes *n(i)* of the corresponding node *i*
            - Column 2 and onwards: *n(i)* timestamps *t*1, *t*2, ..., *tn* **(in ms)** at which the spikes occur for that node

            Parameters
            ----------
            file_name : str, default: 'OUT_SPIK_reformatted.txt'
                Defines the output file name.

            Notes
            -----
            Raw data format\:
            - N rows correspond to N spiking data for N nodes
            - Column 1: the index of node *i*. (*i* runs from 1 to N)
            - Column 2: the number of spikes *n(i)* of the corresponding node *i*
            - Column 3 and onwards: *n(i)* timestamps *t*1, *t*2, ..., *tn* **(in ms)** at which the spikes occur for that node
            """
            with open(self.output_path+file_name, 'w') as file_output:
                for idx in range(self.Nsize):
                    file_output.write('{:.0f}'.format(self.SpikeCount[idx]))
                    for timestamp in self.SpikeTimes[idx]:
                        file_output.write('\t{:.2f}'.format(timestamp))
                    file_output.write('\n')


    class Plot:
        """
        Plot graphs of basic features and properties of the spiking data.
        """

        def __init__(self, SpikeCount, SpikeTimes, Config, output_path):
            self.SpikeCount = SpikeCount
            self.SpikeTimes = SpikeTimes
            self.Nsize = int(Config[0][0])
            self.dt = float(Config[0][1])
            self.Tn = float(Config[0][2])
            self.output_path = output_path

        def SpikeRaster(self, plot_horizontal_stretch=1, file_name='Spiking_Raster_Plot', file_label=''):
            """
            Plot a raster graph of neuron spiking activity.

            Parameters
            ----------
            plot_horizontal_stretch : float or int
                Stretches the plot horizontally by the `plot_horizontal_stretch` factor. Useful for plotting raster graphs for simulations with large simulation duration T.
            file_name : str
                Defines the output file name. Default: Spiking_Raster_Plot.svg.
            file_label : str
                Appends a label / tag at the end of the file name.

            Notes
            -----
            Export a PNG file by default.
            """
            fig, ax = plt.subplots(figsize=(9*plot_horizontal_stretch, 6), dpi=250)
            count = 0
            for each_node in (self.SpikeTimes / 1000):
                count += 1
                ax.scatter(each_node, np.full(np.size(each_node), count), s=0.5, c='black')
            ax.set(xlabel='Time (s)', ylabel='Node index')
            ax.set_xlim(0, float(self.Tn)/1000)
            start_node, end_node = 0, self.Nsize
            ax.set_ylim(start_node-2, end_node+1)
            ax.grid(True)

            if file_label == '': output_file_plot = file_name + '.png'
            else: output_file_plot = file_name + '_' + file_label + '.png'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()
        
        def FiringRateDistribution(self, bins=[0, 80, 800], xrange=[0,10], yrange=[0,None], show_norm=False, file_name='Firing_Rate_Distribution', file_label='', plot_axes=None, info_list=[]):
            """
            Plot a distribution graph of neuron firing rate.

            Parameters
            ----------
            bins : list
                Defines the lower, upper bounds of bins and the number of bins for histogram, in the following format: `[lower bound of bins, upper bound of bins, number of bins]`.
            xrange : list
                Defines the lower, upper limits for x-axis. Set corresponding element to `None` to remove limit.
            yrange : list
                Defines the lower, upper limits for y-axis. Set corresponding element to `None` to remove limit.
            show_norm : bool
                Displays a Gaussian distribution curve with mean and S.D. extracted from input data.
            file_name : str
                Defines the output file name. Default: Firing_Rate_Distribution.svg.
            file_label : str
                Appends a label / tag at the end of the file name.
            plot_axes
                Takes the MATPLOTLIB 'axes' as input. Suppressed all file output functions when used. Plots will be return directly by the method, instead of exporting a file.
            info_list : list
                Specifies the style, legend of the current plot, in the following format: `[line style, legend]`.
            
            Returns
            -------
            `matplotlib.axes.Axes`
                A matplotlib axes returned.

            Notes
            -----
            Export a SVG file by default.\n
            When exporting a SVG file, an TXT file with the same name containing details of firing rate and bin size will also be exported.\n
            `show_norm` and file output are disabled when `plot_axes` takes an matplotlib axes as an argument.
            """
            FiringRate = self.SpikeCount / float(self.Tn) * 1000

            if np.amax(FiringRate) > bins[1]: print('Warning! Maximum of Firing Rate exceeds upper bound of bins range. Max Firing Rate: {}; Max bins range: {}'.format(np.amax(FiringRate), bins[1]))
            if np.amin(FiringRate) < bins[0]: print('Warning! Minimum of Firing Rate subceeds lower bound of bins range. Min Firing Rate: {}; Min bins range: {}'.format(np.amin(FiringRate), bins[0]))
            hist_density, bin_edges = np.histogram(FiringRate, bins=np.linspace(bins[0], bins[1], bins[2]), density=True)
            total_density = np.dot(hist_density, np.diff(bin_edges))
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            # Export plot as SVG file
            if plot_axes == None:
                fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
                ax.plot(x_value, hist_density, 'b^-', lw=2, label='Firing rate distribution\n------------------------------\nT: {:.7}, dt: {:.5}\nmin firing rate: {:.5}, max: {:.5}\n------------------------------\nbin size: {:.5}\nbins range: ({:.5}, {:.5})\nnumber of bins: {:d}'.format( float(self.Tn), float(self.dt), np.amin(FiringRate), np.amax(FiringRate), float((bins[1]-bins[0])/bins[2]), float(bins[0]), float(bins[1]), int(bins[2]) ))

                if show_norm == True:
                    mu = np.mean(FiringRate.flatten()); sigma = np.std(FiringRate.flatten())
                    norm = stats.norm(loc=mu, scale=sigma)
                    ax.plot(x_value, norm.pdf(x_value), 'r--', lw=2, label='Normal distribution with same mean and sd')

                ax.set(xlabel='Firing rate (Hz)', ylabel='Probability density')
                ax.set_xlim(xrange[0], xrange[1]); ax.set_ylim(yrange[0], yrange[1])
                ax.grid(True)
                ax.legend()

                if file_label == '': output_file_plot = file_name + '.svg'
                else: output_file_plot = file_name + '_' + file_label + '.svg'
                fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()

                if file_label == '': output_file_info = file_name + '.txt'
                else: output_file_info = file_name + '_' + file_label + '.txt'
                with open(self.output_path+output_file_info, 'w') as fp_info:
                    fp_info.write('### Plot information ###\n')
                    fp_info.write('Max firing rate: {0}\nMin firing rate: {1}\n'.format(np.amax(FiringRate), np.amin(FiringRate)))
                    fp_info.write('Total density by areal summation: {0}\n'.format(total_density))
                    fp_info.write('\n### Plot settings ###\n')
                    fp_info.write('Bin size: {0}\nBin bounds: lower: {1}, upper: {2}\nNumber of bins: {3}\n'.format( (bins[1]-bins[0])/bins[2], bins[0], bins[1], bins[2] ))
                    fp_info.write('Show Gaussian distribution: {0}'.format(str(show_norm)))

            # Return plot as matplotlib.axes
            else:
                ax = plot_axes
                ax.plot(x_value, hist_density, info_list[0]+'-', lw=2, label=info_list[1])
                return ax

        def InterSpikeIntervalDistribution(self, bins=[0.0005,50,180], xrange=[0.0005,10], yrange=[0,None], file_name='Interspike_Interval_Distribution', file_label='', plot_axes=None, info_list=[]):
            """
            Plot a distribution graph of inter-spike interval.

            Parameters
            ----------
            bins : list
                Defines the lower, upper bounds of bins and the number of bins for histogram, in the following format: `[lower bound of bins, upper bound of bins, number of bins]`.
            xrange : list
                Defines the lower, upper limits for x-axis. Set corresponding element to `None` to remove limit.
            yrange : list
                Defines the lower, upper limits for y-axis. Set corresponding element to `None` to remove limit. If \'yrange\' contains zero, the condition will be ingored for log scale.
            file_name : str
                Defines the output file name. Default: Interspike_Interval_Distribution.svg.
            file_label : str
                Appends a label / tag at the end of the file name.
            plot_axes
                Takes the MATPLOTLIB 'axes' as input. Suppressed all file output functions when used. Plots will be return directly by the method, instead of exporting a file.
            info_list : list
                Specifies the style, legend of the current plot, in the following format: `[line style, legend]`.
            
            Returns
            -------
            `matplotlib.axes.Axes`
                A matplotlib axes returned.

            Notes
            -----
            Export a SVG file by default.\n
            When exporting a SVG file, an TXT file with the same name containing details of firing rate and bin size will also be exported.\n
            `show_norm` and file output are disabled when `plot_axes` takes an matplotlib axes as an argument.
            """
            IsI = np.empty(self.Nsize, dtype=object)
            for count in range(self.Nsize):
                try: IsI[count] = np.array(np.diff(self.SpikeTimes[count]), dtype=float)
                except ValueError: IsI[count] = np.diff(np.array([0])); print('Value error.')
            IsI = np.concatenate([item for item in IsI.flatten()], 0) / 1000

            if np.amax(IsI) > bins[1]: print('Warning! Maximum of ISI exceeds upper bound of bins range. Max ISI: {}; Max bins range: {}'.format(np.amax(IsI), bins[1]))
            if np.amin(IsI) < bins[0]: print('Warning! Minimum of ISI subceeds lower bound of bins range. Min ISI: {}; Min bins range: {}'.format(np.amin(IsI), bins[0]))
            bins[0] = math.log10(bins[0]); bins[1] = math.log10(bins[1])
            hist_density, bin_edges = np.histogram(IsI, bins=np.logspace(bins[0], bins[1], bins[2]))
            hist_density = np.array(hist_density, dtype=float)
            hist_density /= np.dot(hist_density, np.diff(bin_edges)) # normalization
            total_density = np.dot(hist_density, np.diff(bin_edges))
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2
            
            # Export plot as SVG file
            if plot_axes == None:
                fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
                ax.semilogx(x_value, hist_density, 'b^-', lw=2, label='Log ISI distribution\n------------------------------\nT: {:.7}, dt: {:.5}\nmin ISI: {:.5}, max: {:.5}\n------------------------------\nbin size (log scale): {:.5}\nbins range: ({:.5}, {:.5})\nnumber of bins: {:d}'.format( float(self.Tn), float(self.dt), np.amin(IsI), np.amax(IsI), float((bins[1]-bins[0])/bins[2]), float(math.pow(10,bins[0])), float(math.pow(10,bins[1])), int(bins[2]) ))

                ax.set(xlabel='ISI (s)', ylabel='Probability density')
                ax.set_xlim(xrange[0], xrange[1])
                if any(yrange) < 0: print('Warning! \'yrange\' contains zero, the condition will be ingored for log scale.')
                ax.set_ylim(yrange[0], yrange[1])
                ax.grid(True)
                ax.legend()
            
                if file_label == '': output_file_plot = file_name + '.svg'
                else: output_file_plot = file_name + '_' + file_label + '.svg'
                fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()
                
                if file_label == '': output_file_info = file_name + '.txt'
                else: output_file_info = file_name + '_' + file_label + '.txt'
                with open(self.output_path+output_file_info, 'w') as fp_info:
                    fp_info.write('# Plot information #\n')
                    fp_info.write('Max ISI: {0}\nMin ISI: {1}\n'.format(np.amax(IsI), np.amin(IsI)))
                    fp_info.write('Total density by areal summation: {0}\n'.format(total_density))
                    fp_info.write('\n# Plot settings #\n')
                    fp_info.write('Bin size in log scale: {0}\nBin bounds: lower: {1}, upper: {2}\nNumber of bins: {3}'.format( (bins[1]-bins[0])/bins[2], math.pow(10,bins[0]), math.pow(10,bins[1]), bins[2] ))
        
            # Return plot as matplotlib.axes
            else:
                ax = plot_axes
                ax.semilogx(x_value, hist_density, info_list[0]+'-', lw=2, label=info_list[1])
                return ax


class Compare:
    """
    Compare spiking data from different simulations.

    Especially, methods like Plot.ChangeInFiringRateDistribution() can plot the changes in firing rate for neuron dynamics from altered networks.

    Parameters
    ----------
    original_spiking_data_path : str
        Path of the original spiking data file, OUT_SPIK.txt, exclude input folder if `input_folder` is used.
    original_config_path : str
        Path of the original spiking data config file, OUT_SPIK.txt, exclude input folder if `input_folder` is used.
    altered_spiking_data_path : list
        Path of the altered spiking data file, OUT_SPIK.txt, exclude input folder if `input_folder` is used.
    altered_config_path : list
        Path of the altered spiking data config file, OUT_SPIK.txt, exclude input folder if `input_folder` is used.
    delimiter : str, default: '\\t'
        Delimiter to the input spiking data file.
    input_folder : str, default: None
        Path of input folder.
    output_folder : str, default: None
        Path of output folder.

    Notes
    -----
    Invoke object `plot` to plot graphs to compare spiking data.

    Examples
    --------
    To start with, create a object for class `Compare` by specifying the spiking data file paths and the configuration file paths for both the original and altered networks.
    Invoke object `plot` to plot the comparison.
    >>> original_spiking_data_path = 'DIV66\\\OUT_SPIK.txt'
    >>> original_config_path = 'DIV66\\\INI_CNFG'
    >>> # a list of spiking data from networks with inhibitory weights suppressed
    >>> altered_spiking_data_path = ['DIV66_INH_SUP_k05\\\OUT_SPIK.txt', 'DIV66_INH_SUP_k1\\\OUT_SPIK.txt']
    >>> altered_config_path = ['DIV66_INH_SUP_k05\\\INI_CNFG', 'DIV66_INH_SUP_k1\\\INI_CNFG']
    >>> sc = Compare(original_spiking_data_path, original_config_path, altered_spiking_data_path, altered_config_path)
    >>> sc_plot = sc.plot

    Plot the *change in firing rate distribution* of the neuron dynamics from spiking data:
    >>> sc.plot.ChangeInFiringRateDistribution(bins=[-3,12,125], label_list=['k = 0.5','k = 1'])
    """

    def __init__(self, original_spiking_data_path: str, original_config_path: str, altered_spiking_data_path: list, altered_config_path: list, delimiter='\t', input_folder=None, output_folder=None):
        spiking_data_path = [original_spiking_data_path] + altered_spiking_data_path
        config_path = [original_config_path] + altered_config_path
        SpikeCount, SpikeTimes, Config, number_of_files, output_path = self.initSpikeData(spiking_data_path, config_path, delimiter, input_folder, output_folder)
        ChangeInFiringRate = self.initChangeInFiringRate(SpikeCount, Config, number_of_files)
        self.plot = self.Plot(ChangeInFiringRate, Config, number_of_files, output_path)

    def initSpikeData(self, spiking_data_path, config_path, delimiter, input_folder, output_folder):
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if input_folder == None: input_folder = ''
        else: input_folder = input_folder + '\\'
        input_path = this_dir + '\\' + input_folder
        if output_folder == None: output_folder = ''
        else: output_folder = output_folder + '\\'
        output_path = this_dir + '\\' + output_folder
        try: os.mkdir(output_path)
        except FileExistsError: pass

        number_of_files = len(spiking_data_path)
        Config = []

        for each_config in config_path:
            with open(input_path+each_config, 'r') as file_config:
                reader = csv.reader(file_config, delimiter='\t')
                this_Config = np.array(list(reader), dtype=object)
                Config.append(np.array(this_Config[0], dtype=float))
                # Read N, dt, T from each data file
                # # [0]: N, [1]: dt, [2]: T
                # # e.g., self.Config[2][1] --> dt of the simulation in the 2nd data file
        
        Nsize = int(Config[0][0])

        SpikeTimes = [None] * number_of_files
        SpikeCount = [None] * number_of_files

        for count_file in range(number_of_files):
            SpikeTimes[count_file] = np.empty((Nsize), dtype=object)
            SpikeCount[count_file] = np.zeros(Nsize)
            with open(input_path+spiking_data_path[count_file], 'r') as file_spike:
                reader = csv.reader(file_spike, delimiter=delimiter)
                count_node = 0
                for row in reader:
                    SpikeTimes[count_file][count_node] = np.delete(np.array(list(row)).astype('float'), [0,1], 0)
                    SpikeCount[count_file][count_node] = int(row[1])
                    count_node += 1

        SpikeCount = np.array(SpikeCount)
        SpikeTimes = np.array(SpikeTimes)

        return SpikeCount, SpikeTimes, Config, number_of_files, output_path

    def initChangeInFiringRate(self, SpikeCount, Config, number_of_files):
        FiringRate = SpikeCount
        for count in range(number_of_files):
            FiringRate[count] = FiringRate[count] / Config[count][2] * 1000     # self.Config[0][2] = Tn
        ChangeInFiringRate = np.zeros((number_of_files-1, int(Config[0][0])))   # self.Config[0][0] = Nsize
        for count in range(1, number_of_files):
            ChangeInFiringRate[count-1] = FiringRate[count] - FiringRate[0]
        return ChangeInFiringRate

    class Plot:

        def __init__(self, ChangeInFiringRate, Config, number_of_files, output_path):
            self.ChangeInFiringRate = ChangeInFiringRate
            self.number_of_plots = number_of_files - 1
            self.output_path = output_path
        
        def ChangeInFiringRateDistribution(self, bins=[-2,12,120], label_list=[], style_list = ['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', ''], xrange=[-2,12], yrange=[None,None], yaxis_logscale=True, file_name='Change_in_Firing_Rate_Distribution', file_label=''):
            """
            Plot a distribution graph of changes in neuron firing rate.

            Parameters
            ----------
            bins : list
                Defines the lower, upper bounds of bins and the number of bins for histogram, in the following format: `[lower bound of bins, upper bound of bins, number of bins]`.
            label_list : list of str
                Specifies the legend of the plots.
            style_list : list of str
                Specifies the style of the plots. Default: `['gs', 'ro', 'b^', 'mX', 'cD', 'yP', '', '', '', '']`.
            xrange : list
                Defines the lower, upper limits for x-axis. Set corresponding element to `None` to remove limit.
            yrange : list
                Defines the lower, upper limits for y-axis. Set corresponding element to `None` to remove limit. If \'yrange\' contains zero, the condition will be ingored for log scale.
            file_name : str
                Defines the output file name. Default: Interspike_Interval_Distribution.svg.
            file_label : str
                Appends a label / tag at the end of the file name.
            plot_axes
                Takes the MATPLOTLIB 'axes' as input. Suppressed all file output functions when used. Plots will be return directly by the method, instead of exporting a file.
            """
            if np.amax(self.ChangeInFiringRate) > bins[1]: print('Warning! Maximum of Change in Firing Rate exceeds upper bound of bins range. Max Change: {:.5}; Max bins range: {:.5}'.format(np.amax(self.ChangeInFiringRate), float(bins[1])))
            if np.amin(self.ChangeInFiringRate) < bins[0]: print('Warning! Minimum of Change in Firing Rate subceeds lower bound of bins range. Min Change: {:.5}; Min bins range: {:.5}'.format(np.amin(self.ChangeInFiringRate), float(bins[0])))
            hist_density = np.zeros(self.number_of_plots, dtype=object)
            for count in range(self.number_of_plots):
                each_density, bin_edges = np.histogram(self.ChangeInFiringRate[count], bins=np.linspace(bins[0], bins[1], bins[2]), density=True)
                hist_density[count] = each_density
            x_value = (bin_edges[1:] + bin_edges[:-1]) / 2

            fig, ax = plt.subplots(figsize=(9, 6), dpi=50)
            if yaxis_logscale == False:
                for count in range(self.number_of_plots):
                    ax.plot(x_value, hist_density[count], str(style_list[count])+'-', lw=2, label=''+str(label_list[count]))
            elif yaxis_logscale == True:
                for count in range(self.number_of_plots):
                    # # The following lines removes the points with zero density in the plot
                    plot_data_corrected = np.vstack((hist_density[count], x_value))
                    plot_data_corrected = np.delete(plot_data_corrected, np.argwhere(plot_data_corrected[0] == 0), 1)
                    ax.semilogy(plot_data_corrected[1], plot_data_corrected[0], str(style_list[count])+'-', lw=2, label=''+str(label_list[count]))
                    # # Or to use the following line without removing zero points
                    # ax.semilogy(x_value, hist_density[count], str(style_list[count])+'-', lw=2, label=''+str(label_list[count]))
            ax.set(xlabel='Change in firing rate (Hz)', ylabel='Probability density')
            ax.set_xlim(xrange[0], xrange[1])
            ax.set_ylim(yrange[0], yrange[1])
            ax.grid(True)
            ax.legend()

            if file_label == '': output_file_plot = file_name + '.svg'
            else: output_file_plot = file_name + '_' + file_label + '.svg'
            fig.savefig(os.path.join(self.output_path, output_file_plot)); plt.clf()

            if file_label == '': output_file_info = file_name + '.txt'
            else: output_file_info = file_name + '_' + file_label + '.txt'
            with open(self.output_path+output_file_info, 'w') as fp_info:
                fp_info.write('### Plot information ###\n\n')
                fp_info.write('Max change in firing rate: {}\nMin change in firing rate: {}\n'.format(np.amax(self.ChangeInFiringRate), np.amin(self.ChangeInFiringRate)))
                fp_info.write('\n### Plot settings ###\n\n')
                fp_info.write('Bin size: {}\nBin bounds: lower: {}, upper: {}\nNumber of bins: {}\n'.format( (bins[1]-bins[0])/bins[2], bins[0], bins[1], bins[2] ))
