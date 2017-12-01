""" Script to check md5 of the transferred files

The filename of the desired log file will be input by the supra script,
as well as the source from where to extract the md5sum coming from the DB

Each time this code runs will generate a new LOG file, then run it only
once per set of files belonging to one night
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
import fitsio

class Care():
    """ Class devoted to check for integrity of the transferred file
    """

    def __init__(self, logfile=None):
        if (logfile is None):
            hhmmss = datetime.datetime.today().strftime("%H:%M:%S")
            aux = "failCheck_{0}_{1}.txt".format(hhmmss, str(uuid.uuid4()))
            self.logfile = aux

    def write_info(self, line_info):
        """ Method to write line by line to the logfile.
        First it checks if logfile exists, if not exists, creates it
        """
        try:
            with open(self.logfile, "a+") as objfile:
                # Open to append
                objfile.write(line_info)
        except IOError as ioe:
            with open(self.logfile, "a+") as objfile:
                # Create to write
                objfile.info(line_info)
        return True

    def integrity_checksum(self, filename):
        """ Method to check CHECKSUM/DATASUM.
        If want to use CHECKSUM/DATASUM: use fitsio.verify_checksum() for
        each extension. If verification fails, then raise ValueError being
        'data checksum failed' or 'hdu checksum failed'
        Note this option is HIGHLY dependent of fitsio error message
        """
        hdu = fitsio.FITS(filename)
        gonext = True
        ext = 0
        error_check = False
        while gonext:
            try:
                verify = hdu[ext].verify_checksum()
            except Exception as e: #ValueError as ve:
                if ("checksum failed" in str(e)):
                    print 123
                    # save occurence to logfile
                    error_check = True
                gonext = False
            ext += 1
        hdu.close()
        # If the file checksum failed, then save only 1 occurrence on the logs
        if error_check:
            # Format: reason, filename
            line = "CHECKSUM,{0}\n".format(filename)
            self.write_info(line)


    def integrity_md5(self, filename):
        """ Method to compare MD5 from DB against local md5 output
        If want to use md5sum: use MD5SUM from DESFILE and compare with
        value obtained from md5 locally
        """
        pass

class Perm():
    """ Class devoted to change permissions to a set of files belonging to the
    same night, to a folder, or to a file
    """

    def dir_up(self, path, level_up=1):
        """ Method to look for parent directory, to a desired level
        Inputs:
        - path: full path over which operate
        - level_up: number of levels to go up before run chmod
        """
        # If last element is '/', get rid of it
        if (path[-1] == "/"):
            path = path[:-1]
        # If want to go up
        if (level_up > 0):
            # If file, iterate an additional time
            if os.path.isfile(path):
                level_up += 1
            # Get the parent up to a desired level
            for l in range(level_up):
                path = os.path.dirname(path)
        # Check if exists
        if os.path.exists(path):
            return path
        else:
            raise ValueError("Path don't exists: {0}".format(path))
            return False

    def mode_change(self, path, code=775, recursive=True):
        """ Method to change permissions on file or directories
        Inputs:
        - path: full path over which operate
        - code: chmod code to be used
        - recursive: True for using '-R' option in chmod
        """
        # Call of chmod
        if recursive:
            cmd = "chmod -R {0} {1}".format(code, path)
        else:
            cmd = "chmod {0} {1}".format(code, path)
        cA = subprocess.call(shlex.split(cmd))
        cA.wait()
        return True


if __name__ == "__main__":

    # Iterate over the files listed in a text file?
    # Or call once for parent folder?
    # In this scenario, will use the text file containing a list of paths to
    # call integrity_checksum, and then call mode_change() for the set of
    # parent folders
    args = sys.argv

    C = Care()
    P = Perm()

    # Load the table, assuming it contains only full paths
    tab = np.loadtxt(args[0])
    path = []
    for fnm in tab:
        print fnm
        C.integrity_checksum(fnm)
        path.append( P.dir_up(fnm) )
    # Avoid duplicates
    path = np.unique( np.array(path) )
    for upath in path:
        print "chmod working on {0}".format(upath)
        P.mode_change(upath)

    print "Ended!"
