from tools import Coupling


networks = [('DIV66','\t'), ('RANDOM',' ')]

for each in networks:
    c = Coupling(each[0]+'.txt', delimiter=each[1], coupling_enhance_factor=1)
    c_calc = c.calculate
    c_plot = c.plot

    c_plot.AverageSynapticWeight_ExcitatoryIncomingLinks(show_norm=True)
    c_plot.AverageSynapticWeight_InhibitoryIncomingLinks(show_norm=True)
    c_plot.AverageSynapticWeight_OutgoingLinks(show_norm=True)

    coupling_bound = c_calc.CouplingWeightBound()
    connection_probability = c_calc.ConnectionProbability()
    coupling_stat = c_calc.SynapticWeight_Stat()
    exc_node_stat = c_calc.ExcitatorySynapticWeight_Stat()
    inh_node_stat = c_calc.InhibitorySynapticWeight_Stat()

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), each[0]+'.txt'), 'w') as fp:
        fp.write('### Statistics of {} Network ###\n\n'.format(each[0]))
        fp.write('Connection probability: {}\n'.format(connection_probability))
        fp.write('Minimum: {}\nMaximum: {}\n'.format(coupling_bound[0], coupling_bound[1]))
        fp.write('Mean: {}\n'.format(coupling_stat[0]))
        fp.write('SD: {}\n'.format(coupling_stat[1]))
        fp.write('Mean (excitatory): {}\n'.format(coupling_stat[0]))
        fp.write('SD (excitatory): {}\n'.format(exc_node_stat[1]))
        fp.write('Mean (inhibitory): {}\n'.format(inh_node_stat[0]))
        fp.write('SD (inhibitory): {}'.format(inh_node_stat[1]))