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
        print(f'{BDLogger.OK_BLUE_COLOR}{time_str} - {str}{BDLogger.END_COLOR}')

    def error(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        print(f'{BDLogger.FAIL_COLOR}{time_str} - Error: {str}{BDLogger.END_COLOR}')

    def warn(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        print(f'{BDLogger.WARNING_COLOR}{time_str} - Warning: {str}{BDLogger.END_COLOR}')

    def info(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        print(f'{BDLogger.INFO_COLOR}{time_str} - Info: {str}{BDLogger.END_COLOR}')

    def debug(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        print(f'{BDLogger.DEBUG_COLOR}{time_str} - Debug: {str}{BDLogger.END_COLOR}')
