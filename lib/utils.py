# -*- coding: utf-8 -*-
# Misc global functions for the package
import os
import sys
import inspect
import datetime

from termcolor import cprint
from pprint import pprint


def __exit__(line=None, msg=None, color=None):
    """
    Exit, outputting a colored, formatted, line indicating where we stopped in
    the script with an optional message telling us why we stopped.  Useful for
    debugging and writing new code.
    """
    filename, lineno, function = inspect.stack()[1][1:4]
    outcolor = 'on_green' if color is None else 'on_%s' % color

    line = lineno if line is None else line
    output = ('Script "{}" stopped in `{}()` at line: {} with message: '
            '"{}".'.format(os.path.basename(filename), function, line, msg))
    sys.stderr.write('\n')
    cprint(output, "white", outcolor, attrs=['bold'], file=sys.stderr)
    sys.exit()

def today():
        return datetime.datetime.today().strftime('%Y%m%d')

def pp(data):
    pprint(data, stream=sys.stderr)
