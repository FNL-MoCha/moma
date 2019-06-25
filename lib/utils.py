# -*- coding: utf-8 -*-
# Misc global functions for the package
import os
import sys
import inspect

from termcolor import cprint


def __exit__(line=None, msg=None):
    """
    Exit, outputting a colored, formatted, line indicating where we stopped in
    the script with an optional message telling us why we stopped.  Useful for
    debugging and writing new code.
    """
    if line is None:
        filename, line, function = inspect.stack()[1][1:4]
    output = ('Script "{}" stopped in `{}()` at line: {} with message: '
            '"{}".'.format(os.path.basename(filename), function, line, msg))
    sys.stderr.write('\n')
    cprint(output, "white", 'on_green', attrs=['bold'], file=sys.stderr)
    sys.exit()
