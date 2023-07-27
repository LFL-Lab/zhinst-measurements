'''
dd_sequences.py

description:


'''

import numpy as np
from scipy.signal.windows import gaussian
from zhinst.toolkit.waveform import Wave, OutputType
from zhinst.toolkit import CommandTable,Sequence,Waveforms
from utilities import roundToBase
from sequence_setup import make_wave
import matplotlib.pyplot as plt

#%% something descriptive here... 
def seq_time(dd_seq, delay_time,num_pulses, sample_rate,sym,qb_pars={}):
    ''' 
    Calculates time (in seconds) for one dd sequence
    '''
    if dd_seq == 'pure-x': 
        pulse_time = 2 * ( qb_pars['pi_len'] * sample_rate ) 
        if sym is False: # x - d - x - d
            total_delay_time = 2 * delay_time 
        else: # d - x - 2d - x - d 
            total_delay_time = 4 * delay_time
    elif dd_seq == 'xy4':     
        pulse_time = 4 * ( qb_pars['pi_len'] * sample_rate ) 
        if sym is False: # y - d - x - d - y - d - x - d
            total_delay_time = 4 * delay_time
        else: # d - y - 2d - x - 2d- y - 2d -x -d 
            total_delay_time = 8 * delay_time
    #elif dd_seq == 'udd':
    #    pulse_time = num_pulses *  ( qb_pars['pi_len'] * sample_rate ) 
    sequence_time = total_delay_time + pulse_time
    return sequence_time

def max_sequences(experiment_time, sequence_time):
    ''' 
    Calculates the maximum number of dd sequences given an experiment time
    '''
    return int(experiment_time / sequence_time)


#%% main dd sequence stuff
def run_dd(dd_seq,num_pulses, experiment_time, delay_time, sample_rate, qb_pars={}):

    # compute time to run one dd sequence
    time_seq = seq_time(dd_seq, delay_time,num_pulses, sample_rate,sym = True,qb_pars=qb_pars) 
    # compute maximum number of dd sequences given experiment time
    n_max = max_sequences(experiment_time, time_seq)

    for ii in range (0, n_max):
        for jj in range (0,ii):
            # run specific dd sequence here...
            return
        

        # tomography
    return

#%% main dd sequence generator
def gen_dd_seq_code(dd_seq,sym):
    '''
    generates dd sequence code
    '''
    if dd_seq ==  'pure-x':
        code = purex_sequence(sym)
    elif dd_seq == 'xy4': 
        code = xy4_sequence(sym)
    elif dd_seq == 'xy8':
        code = xy8_sequence(sym)
    elif dd_seq == 'udd':
        code = udd_sequence(sym)
    code = finalize_sequence(code)

    return code

#%% dd awg sequence functions

def purex_sequence(sym):
    '''
    generates pure x dd sequence
    '''
    if sym  is False:
        # plays pure-x sequence: x - d - x - d
        awg_program = '''
            repeat(n_seq){
            // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);
            // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);
            }
        '''
    else:
        # plays symmetric pure-x sequence: d - x - 2d - x - d 
        awg_program = '''
            repeat(n_seq){
            // Wait time delay 
            playZero(delay_samples);
           // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(2 * delay_samples);
            // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);
            }
        '''
    

    return awg_program

def xy4_sequence(sym):
    if sym is False:
        # y - d - x - d - y - d - x - d
        awg_program = '''
            repeat(n_seq){
            // Add a y pulse
            executeTableEntry(1);
            // Wait Time Delay
            playZero(delay_samples);
             // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);
            // Add a y pulse
            executeTableEntry(1);
            // Wait Time Delay
            playZero(delay_samples);
             // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);
            }
        '''
    else:
        # d - y - 2d - x - 2d- y - 2d -x -d 
        awg_program = '''
            repeat(n_seq){
            // Wait time delay 
            playZero(delay_samples);
           // Execute a pi y-pulse
            executeTableEntry(1);
             // Wait time delay 
            playZero(2* delay_samples);
           // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(2 * delay_samples);
            // Execute a pi y-pulse
            executeTableEntry(1);
            // Wait time delay 
            playZero(2 * delay_samples);
            // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);
            }
        '''
    return awg_program



def udd_sequence(sym):
    if sym is False:
        awg_program = ''' 
        '''
    return awg_program

#%% Make pulses

def setup_dd_waveforms(sequence,qb_pars={}):
    ''' Creates the waveforms necessary for dd'''
    
    sequence.waveforms = Waveforms()
    N = qb_pars['pi_len']
    # x - pulse
    sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'pi',
                        wave_name = 'x_pulse',
                        amplitude = 1,
                        output_order = '12')),
        Wave(*make_wave(pulse_type = 'zero,',
                        wave_name = 'x_zero',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
                        )
    
    sequence.waveforms[1] = (
        Wave(*make_wave(pulse_type = 'zero,',
                        wave_name = 'y_zero',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12')),
        Wave(*make_wave(pulse_type = 'pi',
                       wave_name = 'y_pulse',
                       amplitude = 1,
                       pulse_length = N,
                       output_order = '12'))
                        )
#%% something something

def setup_dd_pars(sequence,dd_seq, n_seq,delay_time=10e-6,sample_rate = 2.4e9, T = 0):
    if dd_seq == 'udd':
        sequence.constants['n_seq'] = n_seq
    else: 
        delay_samples = delay_time * sample_rate
        sequence.constants['delay_samples'] = delay_samples
        sequence.constants['n_seq'] = n_seq


#%% some kinda function to make command table

def make_command_table():
   return


def construct_time_array(experiment_time,num_pulses):
    """ 
    Constructs an array of time intervals for pi
    pulse sequence to be applied
    """
    if (np.mod(num_pulses,2)!= 0):
        j = np.linspace(1, num_pulses+1, num_pulses+1, dtype = int)
    else: 
        j = np.linspace(1, num_pulses, num_pulses, dtype = int)
    tj = experiment_time * np.sin(j * np.pi / (2 * num_pulses + 2))**2

    plt.plot(j,tj)
    return tj

#%% UDD Time Array Functions (TO BE EDITED)
def tj_k(k, m, tau_j, tj_m_1):
    """
    Returns tj_[k] for QDD_n_[m] where inner pulses occur over
    time [tau_j] and last outer pulse is at [tj_m_1].
    """
    frac = (k * np.pi) / (2 * m + 2)
    return tau_j * (np.sin(frac))**2 + tj_m_1

def make_tj_k_list(m, tau_j, tj_m_1):
    """
    Returns list of tj_k times for QDD_n_[m] inner pulses
    over time [tau_j] where last outer pulse is at [tj_m_1].
    """
    if m % 2 == 0:
        k_list = [k for k in range(1, m + 1)]
    else:
        k_list = [k for k in range(1, m + 2)]
    return np.array([tj_k(k, m, tau_j, tj_m_1) for k in k_list])

def make_udd_sched(n, T, pulse_width=0):
    """
    Make a UDD_[n] time schedule over total DD time [T] where non-ideal
    pulses have finite [pulse_width]. Makes idealized tj times and then
    substracts off [pulse_width] then checks no pulse overlaps.

    Returns physical_times when pulses should START.
    """
    # get idealized tj times
    tj_list = make_tj_list(n, T)
    # turn physical times into integers
    tj_list = np.array(list(map(int, np.floor(tj_list))))
    # subtract off pulse_width to get correct time to begin gates
    phys_tj_list = list(map(int, tj_list - pulse_width))

    # ensure that no gates must be applied at "negative" times
    if phys_tj_list[0] < 0:
        e = (f"Either n too large or T too small to accomodate pulses with\
        width {pulse_width} since first gate must be applied at t1\
        = {phys_tj_list[0]}.\n")
        raise ValueError(e)
    # ensure no diff in t between gates is smaller than pulse_width
    diffs = []
    for idx in range(1, len(phys_tj_list)):
        diffs.append(phys_tj_list[idx] - phys_tj_list[idx - 1])
    min_diff = min(diffs)
    if min_diff < pulse_width:
        e = (f"Minimum pulse spacing required for n={n} and T={T} is\
            {min_diff}, but pulse_width is {pulse_width}.")

    return phys_tj_list

def find_min_udd_time(n, x_width):
    """
    Finds smallest acceptable UDD time which
    corresponds to smallest pulse delay.
    """
    # first, get order of magnitude guess for T
    T = x_width
    success = False
    while success is False:
        try:
            make_udd_sched(n, T, x_width)
            success = True
        except:
            prev_T = T
            T *= 10

    # now try guesses incrementally until minimum T found
    for Tg in range(prev_T, T):
        try:
            make_udd_sched(n, Tg, x_width)
            break
        except:
            continue

    return Tg

def make_mid_udd_sched(n, T, pulse_width=0):
    """
    Make a UDD_[n] time schedule over total DD time [T] where non-ideal
    pulses have finite [pulse_width]. Makes idealized tj times and then
    substracts off [pulse_width]/2 then checks no pulse overlaps.

    Returns physical_times when pulses should START.
    """
    # get idealized tj times
    tj_list = make_tj_list(n, T)
    # subtract off pulse_width to get correct time to begin gates
    phys_tj_list = tj_list - (pulse_width / 2)
    # schedule can only be specified to nearest integer
    phys_tj_list = list(map(int, np.ceil(phys_tj_list)))

    # ensure that no gates must be applied at "negative" times
    if phys_tj_list[0] < 0:
        e = (f"Either n too large or T too small to accomodate pulses with\
        width {pulse_width} since first gate must be applied at t1\
        = {phys_tj_list[0]}.\n")
        raise ValueError(e)
    # ensure no diff in t between gates is smaller than pulse_width
    diffs = []
    for idx in range(1, len(phys_tj_list)):
        diffs.append(phys_tj_list[idx] - phys_tj_list[idx - 1])
    min_diff = min(diffs)
    if min_diff < pulse_width:
        e = (f"Minimum pulse spacing required for n={n} and T={T} is\
            {min_diff}, but pulse_width is {pulse_width}.")

    return phys_tj_list

def make_end_udd_sched(n, T, pulse_width=0):
    """
    Make a UDD_[n] time schedule over total DD time [T] where non-ideal
    pulses have finite [pulse_width]. Makes idealized tj times and uses
    these as start times of pulses.

    Returns physical_times when pulses should START.
    """
    # get idealized tj times
    tj_list = make_tj_list(n, T)
    # subtract off pulse_width to get correct time to begin gates
    phys_tj_list = tj_list - 0
    # schedule can only be specified to nearest integer
    phys_tj_list = list(map(int, np.ceil(phys_tj_list)))

    # ensure that no gates must be applied at "negative" times
    if phys_tj_list[0] < 0:
        e = (f"Either n too large or T too small to accomodate pulses with\
        width {pulse_width} since first gate must be applied at t1\
        = {phys_tj_list[0]}.\n")
        raise ValueError(e)
    # ensure no diff in t between gates is smaller than pulse_width
    diffs = []
    for idx in range(1, len(phys_tj_list)):
        diffs.append(phys_tj_list[idx] - phys_tj_list[idx - 1])
    min_diff = min(diffs)
    if min_diff < pulse_width:
        e = (f"Minimum pulse spacing required for n={n} and T={T} is\
            {min_diff}, but pulse_width is {pulse_width}.")

    return phys_tj_list