###### A Final Year Project in CUHK, Autumn 2021

## **Network Dynaimcs Simulation**

### **Files**

_param.h_
- edit all the variables & settings here

_simulate.c_
- the main program to run the network dynaimcs

### **How to use**
1. edit variables in param.h
2. place param.h and simulate.c in the same folder
3. compile simulate.c
4. wait for results

### **Output**
export up to 4 files

_OUT_SPIK_
- stores all the spiking data
  - column 1: index of nodes, starting from 1
  - column 2: number of spikes of the corresponding node
  - remining columns: time-stamps of each spikes

_OUT_POTV_
- stores the time series of membrane potential v(t) for the network dynamics

_OUT_INFO_
- stores all the variables and settings as well as execution time for a simulation, for later reference

_INI_CNFG_
- same as _OUT_INFO_, designated for easy computer program importation

### **Notes**
1. results will be output in the same folder as the codes, i.e., next to them

## **Analysing Network and Their Dynamics**

### **Files**

_coupling.py_
- calculate
  - connection probability
  - statistics of synaptic weight
  - average synaptic weight
  - ratio of suppression & enhancement
- plot
  - average synaptic weight distribution

_spiking.py_
- calculate
  - average firing rate and its statistics
  - statistics of ISI (inter-spike interval)
  - identifying bursting nodes (work in progess)
  - statistics of synaptic weight of a network
  - average synaptic weight
  - ratio of suppression & enhancement
- plot
  - reformat spiking data
  - spike raster plot
  - firing rate distribution
  - ISI distribution

_spiking_compare.py_
- calculate
  - changes in firing rate
  - ratio of change in firing rate
- plot
  - firing rate distribution (compared)
  - ISI distribution (compared)
  - change in firing rate distribution (&combined)
  - Ratio of suppression / enhancement vs ratio of increase in firing rate (&combined)
