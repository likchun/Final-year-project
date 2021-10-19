/* header.h */

// Set all the variables & settings here.
// This file defines all the variables, such as total number of nodes N, time step dt, etc. as well as settings like MaxBuffer.
// For more details, find README.

// Use this command line in CMD (in Windows) to compile: gcc -O3 simulate.c -o simulate


/* includes & defines */

#include <stdbool.h>


/* import */

// Paths for file import
char const* input_data = "DIV66.txt";
// Delimiter for matrix file input; '\t' for DIV66.txt, ' ' for Random.txt
#define Delimiter "\t"


/* simulation parameters */

// Number of nodes to be simulated
// when this value is smaller than the size of coupling matrix in the input file, top-left sub-matrix will be used
#define N 4095
// Time step size
float dt = 0.05;
// Total simulation time
float tn = 10000;
// Strength of noise
float sigma = 3;


/* suppress, enhance */

bool suppress_inhibitory = false;
bool enhance_inhibitory = false;
bool suppress_excitatory = false;
bool enhance_excitatory = false;
float alter_inh_k = 0;
float alter_inh_sd = 0.00374561833693429;
float alter_exc_k = 0;
float alter_exc_sd = 0.006411179005874374;


/* spiking neuron model parameters */

// Izhikevich's Spiking Neuron Model Parameters
float a_exc = 0.02, b_exc = 0.2, c_exc = -65, d_exc = 8;
float a_inh = 0.1, b_inh = 0.2, c_inh = -65, d_inh = 2;
float thresh_v_exc = 0.0, thresh_v_inh = -80.0, tau_exc = 5, tau_inh = 6, beta = 2;

// Initial values
float v0 = -65, u0 = 8;
// Method to generate randomize initial values; 0: fixed values, 1: gaussian, 2: draw from an interval uniformly
int init_val_rnd_method = 2;
float v0_rnd_sd = 20, u0_rnd_sd = 2.5;
int v0_rnd_max = -20, v0_rnd_min = -75, u0_rnd_max = 12, u0_rnd_min = 0;

// All units are in ms or mV


/* simulation settings */

// Seed for random number generation
float seed_for_random = 0;
// Speeding up by truncation of spiking history, recommend ON
// when enabled, execution time grows approximately linearly with total timesteps
bool history_truncation = true;
float trunc_time_exc = 120;
float trunc_time_inh = 120;


/* controls for I/O and execution */

// Enabling this will output a single file storing the number of spikes and their corresponding time-stamps for each node
bool output_spikeinfo_enabled = true;
// Enabling this will output a single file storing the time series for membrane potential, useful for plotting graphs, file size is typically of several GB.
bool output_potential_enabled = false;
// Export paths
char const* output_data = "OUT_POTV.txt";
char const* output_info = "OUT_INFO.txt";
char const* output_spik = "OUT_SPIK.txt";
char const* output_ini_cnfg = "INI_CNFG";


/* buffers */

// Max buffer for total number of spikes of each node
#define SpikeBuffer 1500
// Max buffer for matrix file input, increase it if the program terminates unexpectedly or 'Function conversion error' is encountered
#define InputBuffer 15000
