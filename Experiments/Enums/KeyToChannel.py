from Utilities.Utils import Utils

# --------------------------------------------------------------------------------------------------------
# KeyToChannel Map
#
# This map indicates for each key, what channels in the QuadRF device are affected by a value change.
# This way, only channels which we actually changed are updated on the QuadRF.
# (e.g. send MOT_duration to both channels 1 & 4)
# --------------------------------------------------------------------------------------------------------

KeyToChannel = {
    'Operation_Mode': [],
    'Triggering_Phase': [1, 4],
    'MOT_freq': [1],
    'PGC_final_freq': [1],
    'Fountain_initial_freq': [1],
    'Fountain_final_freq': [1],
    'Flash_freq': [1],
    'AOM_Repump_freq': [4],
    'MOT_AOM_freq': [2],
    'MOT_duration': [1, 4],        # [msec]
    'MOT_rep': [1, 4],
    'PGC_duration': [1, 4],
    'Post_MOT_delay': [1, 4],
    'PGC_prep_duration': [1, 4],      # [msec]
    'PGC_final_amp': [1],
    'PGC_final_amp_0': [1],       # Relative AOM amplitude between 0 to 1
    'PGC_final_amp_plus': [1],       # Relative AOM amplitude between 0 to 1
    'Repump_PGC_freq': [4],     # The final frequency of the Repump in the PGC sequence
    'Fountain_duration': [1, 4],    # [msec]
    'Fountain_prep_duration': [1, 4],          # [msec]

    # magnetic fountain:
    'PGC_beams_0_off_duration': [],  # [msec]
    'AOM_0_att_3rd_stage':[],

    'Fountain_initial_amp_0': [],  # (also minus) Relative AOM amplitude between 0 to 1
    'Fountain_initial_amp_plus': [],
    'Fountain_initial_amp_minus': [],
    'Fountain_final_amp_0': [],      # (also minus) Relative AOM amplitude between 0 to 1
    'Fountain_final_amp_plus': [],   # Relative AOM amplitude between 0 to 1
    'PGC_final_Delta_freq': [],        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    'PGC_initial_Delta_freq': [],        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    'Fountain_Delta_freq': [],        # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].
    'Fountain_final_Delta_freq': [],
    # The frequency in the frame of reference where the atoms have just enough velocity to reach the toroids [Hz].

    # Free-Fall parameters:
    # The time take the atoms to reach the Toroids and fall back = 2 * sqrt(2 * (Toroid height[m] - MOT height[m]) / g[m/sec^2]) * 1e3[msec]
    'FreeFall_duration': [1, 4],  # We use 1.4 because there is a limit for the pulse duration. = 63ms
    'OD_FS_sleep': [1, 4],
    # Imaging parameters:
    'AOM_Off_Detuning': [4],      # [Hz]; Frequency for essentially turnning light coming from AONs off
    'Snapshot_Intervals': [1, 4],   # [msec]
    'PrePulse_duration': [1, 4],            # [msec]
    'PrePulse_CH1_freq': [1],  # [Hz]
    'Trigger_delay': [1, 4],  # [msec]
    'Pulse_1_duration': [1, 4],       # [msec]
    'InterPulses_duration': [1, 4], # [msec]; Time between pulse 1 & 2.
    'Pulse_2_duration': [1, 4],      # [msec]; If 0 -> no snd pulse
    'Zeeman_delay': [1, 4],          # [msec]
    'N_Snaps': [1, 4],
    # 'Buffer_Cycles': [1,2,3,4]
    'Pulse_1_decay_duration': [1, 4],
    'Pulse_1_amp_i': [1, 4],
    'Pulse_1_amp_f': [1, 4],
    'Pulse_1_CH1_Freq_i': [1],
    'Pulse_1_CH1_Freq_f': [1],
    'Pulse_2_CH1_Freq': [1],
    'Pulse_2_CH_2_3_Freq': [2, 3],
    'Pulse_2_CH4_Freq': [4],
    'PostPulses_freq': [1],
    # Continuous params
    'CH1_freq': [1],
    'CH2_freq': [2],
    'CH3_freq': [3],
    'CH4_freq': [4],
    'CH1_amp': [1],
    'CH2_amp': [2],
    'CH3_amp': [3],
    'CH4_amp': [4]
}

# Run a sanity check on keys in all dictionaries
Utils.verify_keys_for_case_sensitivity([KeyToChannel])
