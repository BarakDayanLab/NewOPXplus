from Experiments.SNSPDCounts import SNSPDCounts_Config_Experiment as Config
from qm.qua import *


def opx_control(obj, qm):

    with program() as opx_control_prog:

        counts1 = declare(int)
        counts2 = declare(int)
        counts3 = declare(int)
        counts4 = declare(int)
        counts5 = declare(int)
        counts6 = declare(int)
        counts7 = declare(int)
        counts8 = declare(int)
        counts9 = declare(int)
        counts10 = declare(int)

        counts_st1 = declare_stream()
        counts_st2 = declare_stream()
        counts_st3 = declare_stream()
        counts_st4 = declare_stream()
        counts_st5 = declare_stream()
        counts_st6 = declare_stream()
        counts_st7 = declare_stream()
        counts_st8 = declare_stream()
        counts_st9 = declare_stream()
        counts_st10 = declare_stream()

        n = declare(int)

        m_window = Config.MOT_pulse_len  # [nsec]
        # diff = declare(int)
        # g2 = declare(fixed, size=m_window)
        # g2_idx = declare(int)
        # g2_st = declare_stream()
        measuring_time = 100 * 1e6  # [nsec]
        rep = int(measuring_time / m_window)

        with infinite_loop_():
            with for_(n, 0, n < rep, n + 1):
                play("Const_open_triggered", "PULSER_N")
                # play("Const_open", "PULSER_N")
                # play("Const_open", "PULSER_S")
                play("Const_open_triggered", "PULSER_S")

                # playing early and late AOM's
                play("Const_open", "PULSER_L")
                # play("Const_high_open", "PULSER_E")

                # play("Square_Pulse", "PULSER_LO")
                # play("Const_open"*amp(0.7), "PULSER_LO")
                play("AntiHelmholtz_MOT", "AntiHelmholtz_Coils")
                # play("Spectrum_pulse", "AOM_Spectrum")
                # play("CRUS_pulse", "Pulser_CRUS")

                # measure("OD_measure", "digital_detectors_S", None,
                measure("OD_measure", "digital_detectors_N", None,
                        counting.digital(counts1, m_window, element_outputs="out1"),
                        counting.digital(counts2, m_window, element_outputs="out2"),
                        counting.digital(counts3, m_window, element_outputs="out3"),
                        counting.digital(counts4, m_window, element_outputs="out4"),
                        counting.digital(counts5, m_window, element_outputs="out5"),
                        counting.digital(counts6, m_window, element_outputs="out6"),
                        counting.digital(counts7, m_window, element_outputs="out7"),
                        counting.digital(counts8, m_window, element_outputs="out8"),
                        counting.digital(counts9, m_window, element_outputs="out9"),
                        counting.digital(counts10, m_window, element_outputs="out10"),
                        )

                ## Save Data: ##
                save(counts1, counts_st1)
                save(counts2, counts_st2)
                save(counts3, counts_st3)
                save(counts4, counts_st4)
                save(counts5, counts_st5)
                save(counts6, counts_st6)
                save(counts7, counts_st7)
                save(counts8, counts_st8)
                save(counts9, counts_st9)
                save(counts10, counts_st10)

        with stream_processing():
            counts_st1.buffer(rep).save("Detector_1_Avg_Counts")
            counts_st2.buffer(rep).save("Detector_2_Avg_Counts")
            counts_st3.buffer(rep).save("Detector_3_Avg_Counts")
            counts_st4.buffer(rep).save("Detector_4_Avg_Counts")
            counts_st5.buffer(rep).save("Detector_5_Avg_Counts")
            counts_st6.buffer(rep).save("Detector_6_Avg_Counts")
            counts_st7.buffer(rep).save("Detector_7_Avg_Counts")
            counts_st8.buffer(rep).save("Detector_8_Avg_Counts")
            counts_st9.buffer(rep).save("Detector_9_Avg_Counts")
            counts_st10.buffer(rep).save("Detector_10_Avg_Counts")

    job = qm.execute(opx_control_prog)
    #job = qm.execute(opx_control_prog, flags=['auto-element-thread'])

    return job
