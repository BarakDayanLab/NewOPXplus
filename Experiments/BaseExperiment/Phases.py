Phases_Names = ['MOT', 'Fountain', 'PGC', 'Free_Fall', 'Pulse_1', 'Inter_Pulses', 'Pulse_2', 'Post_Pulse', 'Push Beam']


# Phases used in experiments.
# Note - we cannot use the Python Enum class as we use the Phase values in JSONs which we serialize and save.
#        JSON save does not allow Enum values.
class Phases:
    UNDEFINED = -1
    MOT = Phases_Names.index('MOT')
    FOUNTAIN = Phases_Names.index('Fountain')
    PGC = Phases_Names.index('PGC')
    FREE_FALL = Phases_Names.index('Free_Fall')
    PULSE_1 = Phases_Names.index('Pulse_1')
    INTER_PULSES = Phases_Names.index('Inter_Pulses')
    PULSE_2 = Phases_Names.index('Pulse_2')
    POST_PULSE = Phases_Names.index('Post_Pulse')
    PUSH_BEAM = Phases_Names.index('Push Beam')
