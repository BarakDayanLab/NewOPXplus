import time

# TODO:
#  1) implement logging into file
#  3) Debug - filtering

class BDLogger:
    HEADER = '\033[95m'
    OK_BLUE_COLOR = '\033[94m'
    OK_CYAN_COLOR = '\033[96m'
    INFO_COLOR = '\033[92m'  # Green
    WARNING_COLOR = '\033[93m'  # Read
    FAIL_COLOR = '\033[91m'
    DEBUG_COLOR = '\033[96m'  # Cyan
    END_COLOR = '\033[0m'
    BOLD_FORMAT = '\033[1m'
    UNDERLINE_FORMAT = '\033[4m'

    def __init__(self):
        pass

    def blue(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        str_r = self._resolve_format(f'{BDLogger.OK_BLUE_COLOR}{time_str} - {str}{BDLogger.END_COLOR}')
        print(str_r)

    def error(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        str_r = self._resolve_format(f'{BDLogger.FAIL_COLOR}{time_str} - Error: {str}{BDLogger.END_COLOR}')
        print(str_r)

    def warn(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        str_r = self._resolve_format(f'{BDLogger.WARNING_COLOR}{time_str} - Warning: {str}{BDLogger.END_COLOR}')
        print(str_r)

    def info(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        str_r = self._resolve_format(f'{BDLogger.INFO_COLOR}{time_str} - Info: {str}{BDLogger.END_COLOR}')
        print(str_r)

    def debug(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        str_r = self._resolve_format(f'{BDLogger.DEBUG_COLOR}{time_str} - Debug: {str}{BDLogger.END_COLOR}')
        print(str_r)

    def _resolve_format(self, str):
        """ TODO: fix bug where we do not return the color format"""
        str = str.replace('**', BDLogger.BOLD_FORMAT, 1)
        str = str.replace('**', BDLogger.END_COLOR, 2)
        str = str.replace('__', BDLogger.UNDERLINE_FORMAT, 1)
        str = str.replace('__', BDLogger.END_COLOR, 2)
        return str