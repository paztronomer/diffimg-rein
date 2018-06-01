"""
Script to perform auxiliary tasks to the Fake Generation code by Stephanie
Hamilton
The input must be a table or dataframe object containing as a minimum:
telra, teldec, band, mjd_obs
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
#
import shamilton_faketno_gen2obs_y5diffim as fakeTNO

class Aid():
    def __init__(self, mcol=["telra", "teldec", "band", "mjd_obs"]):
        """
        Inputs
        - mcol: minimum columns to be contained on the table used as input
        to fake TNO generation
        """
        self.mcol =  map(lambda x: x.upper(), mcol)

    def is_subset(self, child_1D, parent_1D):
        """ Check if all the elements in the child list are contained in the
        parent list. Returns a  boolean.
        Inputs
        - child_1D: list to check if is contained in the parent one
        - parent_1D: list to be used as the parent list, on which to search for
        the smaller list
        Returns
        - boolean
        """
        child_1D = list(map(str.upper, child_1D))
        parent_1D = list(map(str.upper, parent_1D))
        return set(child_1D).issubset(set(parent_1D))

    def merge_tab(self, path, outnm=None, nite=None, exp=None, DEPTH=None):
        """ exp="newdata", DEPTH=0
        Method to construct an unique table from the files generated
        as output of mk_explist.py. There are two types of generated files:
        those accumulating the results 'explist_{nite}_{hhmmss}.csv', and
        those having only new added data named 'newdata_{nite}_{hhmmss}.csv'
        Inputs
        - path: path to the folder containing the output files from mk_explist
        - outnm: full path for the table used as input of SHamilton code. It
        is not necessary to write out, as the dataframe is managed internally
        - nite: night to be processed, used as a double check when looking at
        the files
        - exp: string to match in the filenames. Two options has been coded,
        explist and newdata
        - DEPTH: depth at which to walk through files. Default is zero
        Return
        - df0: output merged dataframe containing the necessary information
        to run fake TNO generation. This dataframe is also saved to disk,
        using filename in arg outname
        """
        logging.info("Merging tables \'{0}\' from mk_explist.py".format(exp))
        c = 0
        # Walk through the files on the input path
        for root, dirs, files in os.walk(path):
            # This is the step checking for the depth of root path
            if root.count(os.sep) >= DEPTH:
                del dirs[:]
            for index, item in enumerate(files):
                fnm = os.path.join(root, item)
                if ((exp in item) and (str(nite) in item)
                     and (os.access(fnm, os.R_OK))):
                    if (c == 0):
                        df0 = pd.read_csv(fnm, engine="python")
                        c += 1
                        logging.info("Table {0} was loaded".format(item))
                    else:
                        df_i = pd.read_csv(fnm, engine="python")
                        df0 = pd.concat([df0, df_i])
                        c += 1
                        logging.info("Table {0} was loaded".format(item))
        # Check the minimum amount of columns is in place
        checkCol = self.is_subset(self.mcol, df0.columns.values)
        if checkCol:
            pass
        else:
            err0 = "Minimum set of columns was not satisfied. Exiting"
            logging.error(err0)
            exit(1)
        # Here use criteria to remove duplicates, using EXPNUM. We expect
        # exposure numbers from the mk_explist.py code
        df0.drop_duplicates(keep="first", inplace=True)
        # Reset index
        df0.reset_index(drop=True, inplace=True)
        if (outnm is not None):
            # Write out
            df0.to_csv(outnm, header=True, index=False)
            msg0 = "Dataframe saved to disk ({0})".format(outnm)
            msg0 += " Duplicates were removed"
            logging.info(msg0)
        return df0

    def file_tab(self, file_fnm):
        """ Load the table given as input, checking if the minimal set of cols
        are in place. To give a file to be loaded is an option to use the
        outputs of mk_explist. Note that no drop of duplicates will be made.
        Inputs
        - file_fnm: complete path to the file
        Returns
        - df1: loaded data frame
        """
        logging.info("Loading file {0}".format(file_fnm))
        try:
            # In pandas, python engine should be able to detect the
            # column separator
            df1 = pd.read_table(file_fnm, sep=None, engine="python")
        except:
                # sys.exc_info() returns a tuple with info of the exception
                # handled now. Values are (type, value, traceback)
            e = sys.exc_info()[0]
            logging.error("Error loading {0}".format(file_fnm))
            logging.error(e)
        # Check for minimal set of columns
        checkCol = self.is_subset(self.mcol, df1.columns.values)
        if checkCol:
            pass
        else:
            err1 = "Minimum set of columns was not satisfied. Exiting"
            logging.error(err1)
            exit(1)
        # Note no duplicate drop will be performed
        logging.info("No duplicate drop was performed on {0}".format(file_fnm))
        return df1

    def remove_old(self, path_rm=None):
        """ Remove old file generated by SHamilton code before create the
        new table of fakes TNO. This is because the filename does not
        change from run to run (convenient for DAGmaker bash script). This
        code ONLY removes files, for directories see shutil.rmtree(), and
        os.removedirs() os.rmdir() for empty directories
        Inputs
        - path_rm: path to the TNO fake file to be removed. If None, then
        the default will be used
        Returns
        - boolean: indicates wether the file was removed
        """
        if (path_rm is None):
            path_rm = "/cvmfs/des.osgstorage.org/stash/wsdiff/fakecats/"
            path_rm += "TNOFAKEOBS_SVTOY4.csv"
        # Here need to keep the distinction between error and warning
        try:
            if (not os.path.exists(path_rm)):
                w2 = "File {0} does not exists. Not removed".format(path_rm)
                logging.warning(w2)
                return False
            os.remove(path_rm)
            i2 = "File {0} was successfully removed".format(path_rm)
            logging.info(i2)
        except:
            # sys.exc_info() returns a tuple with info of the exception
            # handled now. Values are (type, value, traceback)
            e = sys.exc_info()[0]
            logging.error("Error removing {0}".format(path_rm))
            logging.error(e)
            exit(1)
        return True

    def wait_update_storage(self, MM=40):
        """ Method to wait for MM minutes
        """
        time.sleep(MM * float(60))

    def main_tno(self,
                 file_tab=None,
                 exp_path=None, outnm=None, nite=None, exp="newdata", depth=0,
                 KBO_file=None, fakeout_file=None):
        """ Wrapper to: (1) create/read input table; (2) remove old file;
        (3) call SHamilton code, (4) copy generated file to pnfs, (5) wait
        until file is copied to cvmfs, (6) success?/fail?
        """
        # merge_tab(self, path, outnm=None, nite=None, exp="newdata", DEPTH=0)
        # file_tab(self, file_fnm)
        if (file_tab is not None) and (exp_path is None):
            # Load the input table
            df = self.file_tab(file_tab)
        elif (file_tab is None) and (exp_path is not None):
            # Merge tables
            dict1 = {
                "outnm" : outnm,
                "nite" : nite,
                "exp" : exp,
                "DEPTH" : depth,
            }
            df = self.merge_tab(exp_path, **dict1)
        else:
            errN = "Must select either a path to search for tables or give an"
            errN += " table as input. Exiting"
            logging.error(errN)
            exit(1)
        # Using this dataframe, run SHamilton code
        # KBO_file is the table needed by the code to constraint orbital
        # parameters
        if (KBO_file is None):
            KBO_file = "TNOfakes_KBO_distant_P9_gen.csv"
        if (fakeout_file is None):
            fakeout_file = "TNOFAKEOBS_SVTOY4.csv"
        kw_fakeTNO = {
            "exp_df" : df,
            "fakegen_file" : KBO_file,
            "fakeobs_outfile" : fakeout_file,
        }
        #
        # SHamilton code
        fakeTNO.main_caller(**kw_fakeTNO)
        #
        # Check file was generated
        if os.path.exists(fakeout_file):
            pass
        else:
            err = "Fake TNO output ({0}) was not saved.".format(fakeout_file)
            err += " Exiting"
            logging.error(err)
            exit(1)
        # Move file to /pnfs/des/persistent/stash/wsdiff/fakecats/

        # Wait 30min and keep waiting until the file is in place


if __name__ == "__main__":
    # For simplicity, use sys.argv

    ## 1) Add file input
    ## 2) Setup logging!!!
    ## 3) Manage copy/remove files

    ## NEED TO TEST THIS CODE

    # Pending
    # - add option to give a file as input
    # - remove previous file before copy
    # - add a waiting time of 30 minutes

    aux0 = "explist/20171207/newdata_20171207_22h11m11s.csv"

    A = Aid()
    A.file_tab(aux0)
    exit()

    args = sys.argv

    A = Aid()
    A.merge_tab(args[1], args[2], nite=20171119)
