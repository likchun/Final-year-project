###### A Final Year Project in CUHK, Autumn 2021

## Share Portal

https://www.dropbox.com/s/bxs5s3tfq9zsrhk/results-data-plots-on-22-Oct.zip?dl=0

</br>

# Network Dynaimcs Simulation

</br>

### **Files**

_param.h_
- edit all the variables & settings here

_simulate.c_
- the main program to run the network dynaimcs

</br>

### **How to use**
1. edit variables in param.h
   - set input matrix file in _input_data_
   - set delimiter in _Delimiter_
   - set time step in _dt_ and simulation duration in _T_
3. place param.h and simulate.c in the same folder
4. compile simulate.c in Windows CMD by ```gcc -O3 simulate.c -o simulate```
5. wait for results

</br>

### **Advance Settings & Debug**
_output_potential_enabled_
- set to **true** if you need the time series for membrane potential
- Caution: this function is work in progress, don't use this for very small time step _dt_ or very large simulation duration _T_.

_InputBuffer_
- if the program terminates unexpectedly right after it's executed, it's possible that _InputBuffer_ is too small, so that the matrix cannot be read.
- for **DIV66** and **RANDOM** matrix, 15000 will be enough. Increase _InputBuffer_ if it's neccessary.

_SpikeBuffer_
- if the program terminates unexpectedly after it's executed for some time, especially for lengthy simulation duration or strongly coupled networks, it's possible that _SpikeBuffer_ is too small.
- _SpikeBuffer_ must be greater than the total number of spikes of any node.

</br>

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

</br>

### **Notes**
1. results will be output in the same folder as the codes, i.e., next to them

</br>

### **Optimization**

#### Choice of compiler
After compiling the source code with several C compilerson Windows system, MinGW TDM-GCC 64 seems to be a good choice. Its execution time is shorter.
You can find MinGW TDM-GCC 64 here: https://jmeubank.github.io/tdm-gcc/

#### Compiling flag
I recommend using the -O3 flag when compiling, e.g., >gcc -O3 simulate.c -o simulate
It turns on all the -O3 optimization flags, which reduce the running time significantly.
Visit here for more details: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html

#### Notes
This program creates multiple 1-/2-dimensional arrays when running. It accesses the array elements in the tightest loops. Fast memory is essential as the program freqently reads from / writes into RAM.
Also, if you enable output for time series, try to write the file on a fast drive, such as SSD, it will be substantially faster. You can change the output path for time series data file in 'param.h'.

___

</br>

# Analysing Network and Their Dynamics

</br>

### **Functions in _tool.py_**

_Coupling_ ( )
- calculate
  - connection probability
  - statistics of synaptic weight
  - average synaptic weight
  - ratio of suppression & enhancement
- plot
  - average synaptic weight distribution

_Spiking_ ( )
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

_Compare_ ( )
- calculate
  - changes in firing rate
  - ratio of change in firing rate
- plot
  - firing rate distribution (compared)
  - ISI distribution (compared)
  - change in firing rate distribution (&combined)
  - Ratio of suppression / enhancement vs ratio of increase in firing rate (&combined)
