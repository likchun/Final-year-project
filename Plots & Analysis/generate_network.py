import os
from tools import Network


coupling_matrix = 'DIV66.txt'

net = Network(coupling_matrix)

net.generate.ShuffleColumns(output_file=True)
net.generate.ShuffleRows(output_file=True)