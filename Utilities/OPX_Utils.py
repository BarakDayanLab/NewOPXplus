from qm.qua import *

# -------------------------------------------------------------------
# OPX Stuff - can be used from different OPX_Code in all experiments
# -------------------------------------------------------------------


class OPX_Utils:

    def __init__(self):
        raise Exception('No need to initialize OPX_Utils class. It used for statis methods only')

    @staticmethod
    def parameters_update(param):
        return

    # TODO: Q:
    # TODO: Need to check what this does. It seems to do nothing as it multiplies by ZERO.
    # TODO: Dor says it's a hack by Quantum Machine that is suppose to tie variables to elements
    # TODO: Need to check it.
    #
    @staticmethod
    def assign_variables_to_element(element, *variables):
        """
        Forces the given variables to be used by the given element thread. Useful as a workaround for when the compiler
        wrongly assigns variables which can cause gaps.
        To be used at the beginning of a program, will add a 16ns wait to the given element. Use an `align()` if needed.

        Example::

            >>> with program() as program_name:
            >>>     a = declare(int)
            >>>     b = declare(fixed)
            >>>     assign_variables_to_element('resonator', a, b)
            >>>     align()
            >>>     ...

        """
        _exp = variables[0]
        for variable in variables[1:]:
            _exp += variable
        wait(4 + 0 * _exp, element)
