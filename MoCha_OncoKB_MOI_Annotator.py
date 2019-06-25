#!/usr/bin/env python3
import sys
import os
import subprocess

from ion.plugin import *

class MoCha_OncoKB_MOI_Annotator(IonPlugin):
    """
    Plugin to generate an MOI report from Oncomine data using Annovar and the
    OncoKB databases.
    """
    with open(os.path.join(os.environ['DIRNAME'], '_version.py')) as fh:
        version = fh.readline().rstrip("'\n").split("'")[1]
    # version = '0.10.20180921'
    major_block = False
    runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]
    runlevels = [RunLevel.DEFAULT]
    depends = ["variantCaller"]

    def launch(self, data=None):
        cmd = [
            os.path.join(os.environ['DIRNAME'], 'moma_plugin.py'),
            '-V', self.version, 
            'startplugin.json',
            'barcodes.json'
        ]
        plugin = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
        plugin.communicate()
        sys.exit(plugin.poll())

if __name__ == '__main__':
    PluginCLI()
