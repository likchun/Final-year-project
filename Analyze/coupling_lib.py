import os, csv, math
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.style as mplstyle
import seaborn as sns


class Initialization:

    def __init__(self, coupling_matrix, delimiter, input_folder, output_folder):
        self.initCouplingMatrix(coupling_matrix, delimiter, input_folder, output_folder)

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

class Coupling(Initialization):

    def __init__(self, coupling_matrix, coupling_enhance_factor=1, delimiter='\t', input_folder=None, output_folder=None):
        super().__init__(coupling_matrix, delimiter, input_folder, output_folder)
        self.Coupling = coupling_enhance_factor * self.Coupling
        self.calculate = self.Calculate(self.Coupling)
        self.plot = self.Plot(self.Coupling, self.output_path)

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