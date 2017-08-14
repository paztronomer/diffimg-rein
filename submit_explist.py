""" Script derived from loopsubmit.sh, adding some extra features
"""
import os
import sys
import errno
import time
import datetime
import socket
import shlex
import subprocess
import logging
import argparse


class Toolbox():
    def split_path(self, path):
        """Method to return the relative root folder (one level upper),
        given a path.
        Inputs
        - path: complete path
        Returns
        - 2 strings, one being the parent folder, and the filename
        """
        #relat_root = os.path.abspath(os.path.join(path,os.pardir))
        relroot, filename = os.path.split(path)
        return (relroot, filename)


class Work():
    def __init__(self, im_list=None, im_table=None, dir_log=None, nite=None):
        """ Method to initialize Work()
        Inputs
        - im_list: list. List of pull paths to be used as the images to be
        processed
        - im_table: str. Filename of the plain text (1 column) file harboring
        the full paths to the images to be processed
        - dir_log: str, folder to harbor the logs for the processing
        - nite: int. Night of observation to be processed
        """
        if (im_list is None) and (im_table is not None):
            # Read the im_table and save in a list
            with open(im_table) as f:
                aux_lines = f.read().splitlines()
            self.im_path = aux_lines
        elif (im_list is not None) and (im_table is None):
            self.im_path = im_list
        else:
            logging.error("Input image list OR image table, not both")
            exit(1)
        self.dir_log = dir_log
        self.nite = nite
        # String containing the exec time of this script, used for naming logs
        self.isotime = datetime.datetime.today().strftime("%Y%m%dT%H:%M:%S") 

    def runproc(self):
        """ Method to iterate over the files and submit the work
        """
        T = Toolbox()
        for idx, im in enumerate(self.im_path):
            print im
            rootpath, fname = T.split_path(im)
            expnum = int(fname[1 : fname.find("_")])
            #
            msg1 = "Submitting exposure {0}".format(fname)
            msg1 += " through {0}.".format(os.path.basename(__file__))
            msg1 += " On {0}\n".format(time.ctime())
            msg1 += " N={0} of {1}".format(idx + 1, len(self.im_path))
            logging.info(msg1)
            #
            cmdA = "./DAGMaker.sh {0:08}".format(expnum)
            # cmdA = "./DAGMaker.sh {0}".format(fname)
            cmdA += " > logs/DAGMaker_{0:08}.log 2>&1".format(expnum)
            logging.info("Running: {0}".format(cmdA))
            cmdA = shlex.split(cmdA)
            prA = subprocess.Popen(cmdA,
                                   shell=True,
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
            prA.wait()
            #
            cmdB = "jobsub_submit_dag -G des --role=DESGW"
            cmdB += " file://desgw_pipeline_{0:08}.dag".format(expnum)
            cmdB += " > logs/jobsub_submit_dag_{0:08}.log 2>&1".format(expnum)
            logging.info("Running: {0}".format(cmdB))
            cmdB = shlex.split(cmdB)
            prB = subprocess.Popen(cmdB,
                                   shell=True,
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
            prB.wait()
        logging.info("Submission of {0} has ended".format(len(self.im_path)))
        return True

    def setup_log(self):
        """ Method to setup the log output and start with some information
        """
        # Check/create directory
        if self.dir_log is None:
            self.dir_log = os.path.join(os.getcwd(), "logs/")
        try:
            self.dir_log = os.path.join(self.dir_log, self.nite)
            os.makedirs(self.dir_log)
        except OSError as exception:
            if (exception.errno != errno.EEXIST):
                raise
                logging.error("ERROR when creating {0}".format(self.dir_log))
        # Setup write out
        lognm = "proc_nite_{0}_run_{1}.log".format(self.nite, self.isotime)
        logpath = os.path.join(self.dir_log, lognm)
        logging.basicConfig(filename=logpath, level=logging.DEBUG, 
                            format="%(asctime)s - %(levelname)s - %(message)s")
        # First information
        logging.info("\nRunning on: {0}\n".format(socket.gethostname()))
        logging.info("Script: {0}\n".format(os.path.basename(__file__)))
        return True

if __name__ == "__main__":
    # Parse of arguments
    h = "Script to iterate over a list of exposures (full path) and submit"
    h += " the diffimg processing. The LOG files will be named as"
    h += " proc_nite_<nite>_run_<actual time at which proc started>.log"
    abc = argparse.ArgumentParser(description=h)
    # Positional
    txt0 = "Filename of the list containing the path to the exposures"
    abc.add_argument("tab_exp", help=txt0)
    # Optional
    txt1 = "Directory where to store the LOGs. One folder per night."
    txt1 += " Default: <current_directory>/logs"
    abc.add_argument("--dir_log", help=txt1, metavar="")
    #
    txt2 = "Night of observing to which the current images belongs."
    txt2 += " Default: last night"
    abc.add_argument("--nite", help=txt2, metavar="")
    # Recover args
    val = abc.parse_args()
    kw = dict()
    kw["im_table"] = val.tab_exp
    kw["dir_log"] = val.dir_log
    kw["nite"] = val.nite
    #
    W = Work(**kw)
    # Setup the logs
    W.setup_log()
    # Run processing
    W.runproc()
