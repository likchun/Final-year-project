import os
from tools import Network


networks = [
    {'name':'DIV66',    'delimiter':'\t'},
    {'name':'Random',   'delimiter':' ' },
]

for network in networks:
    net = Network(network['name']+'.txt', delimiter=network['delimiter'], coupling_enhance_factor=1)

    print('\n> Displaying details for network \'{}\' <\n'.format(network['name']))

    # net.plot.AverageSynapticWeight_ExcitatoryIncomingLinks(show_norm=True)
    # net.plot.AverageSynapticWeight_InhibitoryIncomingLinks(show_norm=True)
    # net.plot.AverageSynapticWeight_OutgoingLinks(show_norm=True)

    connection_probability = net.calculate.ConnectionProbability()
    synaptic_weight_bound = net.calculate.SynapticWeightBounds()
    synaptic_weight_stat = net.calculate.SynapticWeight_Stat()
    pos_link_stat = net.calculate.PositiveSynapticWeight_Stat()
    neg_link_stat = net.calculate.NegativeSynapticWeight_Stat()

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), network['name']+'_stat.txt'), 'w') as fp:
        fp.write('### Statistics of Network \'{}\' ###\n\n'.format(network['name']))
        fp.write('Connection probability: {}\n'.format(connection_probability))
        fp.write('\nStatistics of synaptic weights:\n')
        fp.write('Minimum weight:\t{}\nMaximum weight:\t{}\n'.format(synaptic_weight_bound[0], synaptic_weight_bound[1]))
        fp.write('Mean (all links):\t{}\n'.format(synaptic_weight_stat[0]))
        fp.write('S.D. (all links):\t{}\n'.format(synaptic_weight_stat[1]))
        fp.write('Mean (negative links):\t{}\n'.format(neg_link_stat[0]))
        fp.write('S.D. (negative links):\t{}\n'.format(neg_link_stat[1]))
        fp.write('Mean (positive links):\t{}\n'.format(pos_link_stat[0]))
        fp.write('S.D. (positive links):\t{}'.format(pos_link_stat[1]))

    print('\n------------------------------next------------------------------')
