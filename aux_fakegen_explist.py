"""
Script to perform auxiliary tasks to the Fake Generation code by Stephanie
Hamilton
"""

import os
import sys
import shutil
import socket
import datetime
import time
import gc
import logging
import argparse
import itertools
import errno
import shlex
import subprocess
import collections
import uuid
import numpy as np
import pandas as pd

class Aid():
    def merge_tab(sef, path, outnm, nite=None, exp="newdata_", DEPTH=0):
        """ Method that merges all tables inside the given folder, having the
        input expression in the filename into one file. Only goes to the given
        DEPTH
        """
        cnt = 0
        for root, dirs, files in os.walk(path):
            if root.count(os.sep) >= DEPTH:
                del dirs[:]
            for index,item in enumerate(files):
                fnm = os.path.join(root, item)
                if (exp in item) and (os.access(fnm, os.R_OK)):
                    if (c == 0):
                        df0 = pd.read_csv(fnm, engine="python")
                        c += 1
                    else:
                        df_i = pd.read_csv(fnm, engine="python")
                        df0 = pd.concat([df0, df_i])
                        c += 1
        # Reset index
        df0.reset_index(drop=True, inplace=True)
        print df0.info()
        print df0[1]
        df0.write_csv(outnm, header=True, index=False)
        return df0

if __name__ == "__main__":
    # For simplicity, use sys.argv

    args = sys.argv

    A = Aid()
    A.merge_tab(args[1], args[2], nite=20171119)
