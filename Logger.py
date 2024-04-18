import logging
from colorama import Fore, Back, Style, init

# Initialize colorama
init()

class ColorizingStreamHandler(logging.StreamHandler):
    color_map = {
        logging.DEBUG: Fore.BLUE,
        logging.INFO: Fore.GREEN,
        logging.WARNING: Fore.YELLOW,
        logging.ERROR: Fore.RED,
        logging.CRITICAL: Fore.RED + Back.WHITE + Style.BRIGHT,
    }

    def emit(self, record):
        try:
            message = self.format(record)
            self.stream.write(self.color_map.get(record.levelno) + message + Style.RESET_ALL + "\n")
            self.flush()
        except Exception:
            self.handleError(record)

def setup_logger():
    # Create a custom logger
    logger = logging.getLogger(__name__)

    if not logger.handlers:
        # Set the level of this logger. 
        logger.setLevel(logging.DEBUG)

        # Create handlers
        console_handler = ColorizingStreamHandler()
        file_handler = logging.FileHandler('logfile.log')

        # Set level for handlers
        console_handler.setLevel(logging.DEBUG)

        # Create formatters and add it to handlers
        console_format = logging.Formatter('[%(levelname)s (%(asctime)s)] ' + Style.RESET_ALL + '- %(message)s', datefmt='%I:%M:%S')
        file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_format)
        file_handler.setFormatter(file_format)

        # Add handlers to the logger
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

    return logger