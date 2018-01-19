""" Script to modify the header of the immasked images, to meet the
DESWS criteria for DiffImg processing
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
import fitsio

# PENDING
# - add argparse
# - run in parallel

class Change_Keys():
    """ Class to change the keys of the header
    """
    def __init__(self, logname=None):
        # Setup the file on whic to save the probable failed outputs
        # IMPORTANT! Check how to append and not overwrite the LOG
        if not (logname is None):
            logging.basicConfig(
                filename=logname,
                level=logging.DEBUG,
                format="%(asctime)s - %(levelname)s - %(message)s"
            )

    def modify(self, fits, extname="SCI"):
        """ Method to change/add header keys
        """
        if os.path.exists(fits):
            try:
                # Copy from archive to local home, to have both files is
                # a safe measure
                cmdA = "cp {0} {1}".format(fits, os.path.expanduser("~"))
                cA = subprocess.call(shlex.split(cmdA))
                logging.info("{0} copied to home".format(fits))
                #
                local_fits = os.path.join(os.path.expanduser("~"),
                                          os.path.basename(fits))
                hdu = fitsio.FITS(local_fits, "rw")
                sci_ext = hdu[extname]
                h = sci_ext.read_header()
                # Get RA, DEC, transform them and use
                ra_hh = [np.float(x) for x in h["RA"].split(":")]
                ra_hh = ra_hh[0] + ra_hh[1] / 60. + ra_hh[2] / 3600.
                ra_d = ra_hh * 15.
                dec_d = [np.float(y) for y in h["DEC"].split(":")]
                sign = np.sign(dec_d[0])
                dec_d = np.abs(dec_d[0]) + dec_d[1] / 60. + dec_d[2] / 3600.
                # Construct the strings to go into the header
                aux_ra_d = int(np.floor(ra_d * 10.))
                aux_dec_d = int(np.floor(dec_d * 10))
                if (sign == -1):
                    s = "-"
                elif (sign == 1):
                    s = "+"
                else:
                    s = ""
                XY = "{0}{1}{2}".format(aux_ra_d, s, aux_dec_d)
                k_object = "DESWS hex WS{0} tiling 1".format(XY)
                k_field = "WS{0}".format(XY)
                k_tiling = 1
                # Construct the list of dictionaries to be used as input for
                # write new header values
                d0 = {
                    "name": "OBJECT",
                    "value": k_object,
                    "comment": "Modified for DiffImg on {0}".format(time.ctime()),
                }
                d1 = {
                    "name": "FIELD",
                    "value": k_field,
                    "comment": "Modified on {0}".format(time.ctime()),
                }
                d2 = {
                    "name": "TILING",
                    "value": k_tiling,
                    "comment": "Modified on {0}".format(time.ctime()),
                }
                new_keys = [d0, d1, d2]
                sci_ext.write_keys(new_keys)
                # Close HDU
                hdu.close()
                # Move back to original location, erasing from home
                cmdB = "mv -v --force {0} {1}".format(local_fits, fits)
                cB = subprocess.Popen(shlex.split(cmdB),
                                      stdout=subprocess.PIPE)
                out_cB = cB.communicate()
                cB.wait()
                logging.info(out_cB)
                logging.info("{0} moved back to {1}".format(fits, local_fits))
            except Exception as e:
                logging.error(e)
                logging.error("Failed HEADER modification on {0}".format(fits))
        else:
            logging.error("Inexistent of not readable file: {0}".format(fits))

if __name__ == "__main__":
    t0 = "Script to modify the header (and thus change the image itself)"
    t0 += " to match DiffImg requirements of keywords"
    ax = argparse.ArgumentParser(description=t0)
    h1 = "File containing full paths to the FITS, one per line"
    ax.add_argument("fits_list", help=h1, metavar="")
    h2 = "Full path of the log to be created. Suggested: logs/{night}/{N.log}"
    ax.add_argument("logfile", help=h2, metavar="", type=str)
    # Parse arguments
    arg = ax.parse_args()
    
    # Initialize
    CK = Change_Keys(logname=arg.logfile)

    logging.info("Starting on {0}".format(time.ctime()))

    # Run over a list of FITS files
    # Load the table, assuming it contains only full paths
    tab = pd.read_table(arg.fits_list, header=None, names=["path"])
    for idx, row in tab.iterrows():
        CK.modify(row["path"])

    txt_end = "Ended on {0}, {1} FITS files".format(time.ctime(), tab.size)
    logging.info(txt_end)
