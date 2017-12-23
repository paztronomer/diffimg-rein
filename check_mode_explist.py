"""
Script to change permissions and checksum of the transferred files
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
imprt errno
import shlex
import subprocess
import collections
import uuid
import numpy as np
import pandas as pd
import fitsio

# =====================================
# Pending:
# - get md5sum from DB and compare with the results from bash md5
# =====================================

class Keep_Logs():
    def __init__(self, nite, dir_log=None):
        """ Method to setup the log output and start with some information
        """
        # Check/assign LOG directory
        if (dir_log is None):
            self.dir_log = os.path.join(os.getcwd(), "logs/")
        else:
            self.dir_log = dir_log
        # Check/create directory
        try:
            self.dir_log = os.path.join(self.dir_log, nite)
            os.makedirs(self.dir_log)
        except OSError as exception:
            if (exception.errno != errno.EEXIST):
                raise
                logging.error("ERROR when creating {0}".format(self.dir_log))
        # Setup write out
        slef.hhmmss = datetime.datetime.today().strftime("%Hh%Mm%Ss")
        lognm = "explist_and_copy_{0}_{1}.log".format(nite, self.hhmmss)
        logpath = os.path.join(self.dir_log, lognm)
        logging.basicConfig(filename=logpath, level=logging.DEBUG,
                            format="%(asctime)s - %(levelname)s - %(message)s")
        # First information
        logging.info("\nRunning on: {0}\n".format(socket.gethostname()))
        logging.info("Script: {0}\n".format(os.path.basename(__file__)))


class Care():
    def __init__(self, log2file=None):
        if (log2file is None):
            hhmmss = datetime.datetime.today().strftime("%Hh%Mm%Ss")
            aux = "failCheck_{0}_{1}.log".format(hhmmss, str(uuid.uuid4()))
            self.log2file = aux
        else:
            self.log2file = log2file

    def write_info(self, line_info):
        """ Method to write line by line to the log2file.
        First it checks if log2file exists, if not exists, creates it
        Inputs:
        - line_info: line to be appended to file, including newline character
        """
        try:
            with open(self.log2file, "a+") as objfile:
                # Open to append
                objfile.write(line_info)
        except IOError as ioe:
            with open(self.log2file, "w+") as objfile:
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
        if os.path.exists(filename):
            hdu = fitsio.FITS(filename)
            gonext = True
            ext = 0
            error_check = False
            while gonext:
                try:
                    verify = hdu[ext].verify_checksum()
                except Exception as e: #ValueError as ve:
                    if ("checksum failed" in str(e)):
                        # save occurence to log2file
                        error_check = True
                    gonext = False
                ext += 1
            hdu.close()
            # If the file checksum failed, then save only 1 occurrence on log
            if error_check:
                # Format: reason, filename
                line = "CHECKSUM,{0}\n".format(filename)
                self.write_info(line)
        else:
            # This error means file doesn't exists or user has not read 
            # permission
            line = "IOERROR,{0}\n".format(filename)
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

    def mode_change(self, path, code=None, recursive=True):
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
        # The subprocess.call() waits for the process to finish and returns
        # returncode attribute
        cA = subprocess.call(shlex.split(cmd))
        if (cA != 0):
            raise ValueError("Error on chmod call")
            exit(1)
        else:
            return True


if __name__ == "__main__":
    txt = "Script to checksum and datasum from FITS header, using fitsio. It"
    txt += " also change permissions on the files"
    arg = argparse.ArgumentParser(description=txt)
    # 
    txt1 = "Text file containing the list of full paths of the files over"
    txt1 += " which to checksum/datasum and change permissions."
    arg.add_argument("file_list", help=txt1, metavar="")
    #
    txt1b = "Night on which we are working"
    arg.add_argument("night", help=txt1b, metavar="", type=int)
    # 
    aux_mod = 764
    txt2 = "Permission to be recursively set on the file."
    txt2 += "  Default: {0}".format(aux_mod)
    arg.add_argument("--perm", help=txt2, metavar="", type=int)
    #
    txt3 = "Directory where to store the LOGs. One folder per night."
    txt3 = " Default: <current_directory>/logs"
    abc.add_argument("--d_log", help=txt3, metavar="")
    
    # Setup the logs
    Keep_Logs(arg.night, )

    # Load the table, assuming it contains only full paths
    tab = pd.read_table(arg.file_list, header=None, names=["path"])
    tab = tab["path"].values
    
    # Initiate classes
    P = Perm()
    C = Care(log2file=arg.perm)
    
    # Run checksum
    x_path = []
    for fnm in tab:
        C.integrity_checksum(fnm)
        x_path.append( P.dir_up(fnm) )
    
    # Change permissions, avoiding duplicates
    x_path = np.unique( np.array(x_path) )
    for upath in x_path:
        logging.info("chmod working on {0}".format(upath))
        P.mode_change(upath)
    
