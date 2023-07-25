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

#%% something descriptive here... 
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
        awg_program = '''
            repeat(n_pulses){
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
        awg_program = '''
            repeat(n_pulses){
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
        awg_program = '''
            repeat(n_pulses){
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
        awg_program = '''
            repeat(n_pulses){
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
                        output_order = '21')),
        Wave(*make_wave(pulse_type = 'pi',
                       wave_name = 'y_pulse',
                       amplitude = 1,
                       pulse_length = N,
                       output_order = '21'))
                        )
#%% something something

def setup_dd_pars(sequence,dd_seq, n_pulses, T = 0)