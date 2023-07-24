from sequence_setup import make_wave

def gen_dd_seq_code(dd_seq):
    '''
    generates dd sequence code
    '''
    if dd_seq ==  'pure-x':
        code = purex_sequence()
    elif dd_seq == 'xy4': 
        code = xy4_sequence()
    elif dd_seq == 'xy8':
        code = xy8_sequence()
    elif dd_seq == 'udd':
        code = udd_sequence()
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
            // Add a y pulse
            executeTableEntry(1);
            // Wait Time Delay
            playZero(delay_samples);
             // Execute a pi x-pulse
            executeTableEntry(0);
            // Wait time delay 
            playZero(delay_samples);

        '''
    else:
        awg_program = '''
        '''
    return awg_program