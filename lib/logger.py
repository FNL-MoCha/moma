# -*- coding: utf-8 -*-
# Custom logger class for package
import sys
import datetime
from termcolor import colored
from textwrap import indent, fill


class Logger(object):
    def __init__(self, loglevel='debug', colored_output=False, dest=None, 
            quiet=True):
        """
        Basic Logger class for the MOMA package.
        This logger will take in the loglevel and destination for logging, and
        set up a system to write messages to a logfile (or STDERR if no
        destination indicated.

        Args:
            loglevel (str): Maximum level to which log messages should be
                            written. Valid log levels are, from lowest to 
                            highest: `debug`, `warning`, `error`, `info`. We
                            will output anything on a lower tier than the log
                            level entered.
            colored_output (bool): Output log messages in color or in black and 
                            white.
            dest (str):    Log file to which the data should be written. If not
                           set, the log messages will be written to STDOUT.
            quiet (bool):  If set to true, do not output anything to stdout;
                           only output to file.
        """

        info_tag = colored(" INFO ", 'green', 
                attrs=['bold']) if colored_output else "INFO"
        warn_tag = colored(" WARN ", 'yellow', 
                attrs=['bold']) if colored_output else "WARN"
        error_tag = colored(" ERROR ", 'red', 
                attrs=['bold']) if colored_output else "ERROR"
        debug_tag = colored(" DEBUG ", 'cyan', 
                attrs=['bold']) if colored_output else "DEBUG"
        note_tag = colored(" NOTE ", 'magenta', 
                attrs=['bold']) if colored_output else "NOTE"

        self.colored_output = colored_output
        self.quiet = quiet

        self.log_levels = {
            'info'   : (0, info_tag), 
            'note'   : (0, note_tag),
            'warn'   : (1, warn_tag), 
            'error'  : (2, error_tag), 
            'debug'  : (3, debug_tag),
        }
        self.log_key = loglevel
        try:
            self.log_level, self.log_tag = self.log_levels[self.log_key]
        except KeyError:
            sys.stderr.write("ERROR: '%s' is not a valid log level. Valid "
                    "choices are: 'info', 'warn', 'error', or 'debug'. Can not "
                    "create logger.\n\n" % loglevel)
            raise

        # Create a filehandle
        if dest is not None:
            self.outfh = open(dest, "w")
        else:
            self.outfh = sys.stderr
        
    def __enumerate_loglevel(self, loglevel):
        try:
            level_tag = self.log_levels[loglevel]
        except KeyError:
            sys.stderr.write("ERROR: '%s' is not a valid loglevel. Can not "
                    "write this entry to log file.\n" % loglevel)
            return (None, None)
        return level_tag

    def __get_time(self):
        return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def write_log(self, logtype, message):
        """
        Main log writer method. Will write the log message if the logtype is
        below the logger threshold.

        Args:
            logtype (str): Type of log message to be considered for printing.
                           Can be any one of the legitimate levels, such as
                           'info', 'warn', 'error', or 'debug'. If a blank
                           string is pass, then print a '\\t' char followed by
                           the message; intended as a continuation of the
                           previous.
            message (str): Message to be printed.

        Returns:
            Output to class filehandle, which is either STDERR or a logfile.

        Examples:

            >>> Logger.write_log('debug', "This is a test message")
            2019-06-11 13:53:52  [ DEBUG ]: This is a test message

        """
        if logtype:
            if logtype == 'header' or logtype == 'unformatted':
                self.outfh.write(f'{message}\n')
                self.outfh.flush()
                if self.quiet is False:
                    sys.stderr.write(f'{message}\n')
                    sys.stderr.flush()
            else: 
                level, tag = self.__enumerate_loglevel(logtype)
                
                if level is not None and level <= self.log_level:
                    if self.colored_output:
                        label = (colored('[', 'white', attrs=['bold']) + tag
                            + colored(']', 'white', attrs=['bold']) + ':')
                        outstr = '{:33} {:<49}  {:<40}\n'.format(
                            colored(self.__get_time(), 'white', attrs=['bold']),
                            label,
                            message
                         )
                    else:
                        outstr = '{:20} {}{:^8}{}:  {}\n'.format(
                            self.__get_time(), '[', tag, ']', message
                        )
                    self.outfh.write(outstr)
                    self.outfh.flush()
                    if self.quiet is False:
                        sys.stderr.write(outstr)
                        sys.stderr.flush()
        else:
            outstr = fill(message, width=80)
            self.outfh.write(indent(f'{outstr}\n', '    '))
            self.outfh.flush()
            if self.quiet is False:
                sys.stderr.write(indent(f'{outstr}\n', '    '))
                sys.stderr.flush()

