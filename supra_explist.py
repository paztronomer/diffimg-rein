""" Script to manager the execution of mk_explist.py and submit_explist.py

Important
=========
Before run this script, please do source of a file containing the below:
-----------------------------------------------------------------
#!/bin/bash

source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
setup -v firstcut Y5Ndev+3
setup -v easyaccess
setup -v pandas 0.15.2+5
setup -v astropy
setup -v numpy 1.9.1+11
setup -v python 2.7.9+1

kx509
voms-proxy-init -rfc -noregen -voms des:/des/Role=DESGW
-----------------------------------------------------------------
"""

import os
import logging
import argparse
import subprocess
# my modules
import mk_explist
import submit_explist


if __name__=="__main__":
    #
    # Using the same argparse as in mk_explist, but changing the gral help
    #
    # Parse of arguments
    hlp = "Script to detect the last night (or other) already processed"
    hlp += " exposures, and then create bash executable files for remote"
    hlp += " copy. Then, the remote copy is done and exposures are used"
    hlp += " for diffimg"
    abc = argparse.ArgumentParser(description=hlp)
    # Optional
    txt1 = "Night to be queried"
    abc.add_argument("--nite", help=txt1, metavar="")
    #
    immask_aux = "/pnfs/des/persistent/wsdiff/exp"
    txt2 = "Parent folder to harbor immask files per CCD."
    txt2 += " Default: {0}".format(immask_aux)
    abc.add_argument("--d_msk", help=txt2, metavar="", default=immask_aux)
    #
    txt3 = "Parent folder to harbor the exposure lists (various files)."
    txt3 += " Each night has a different folder."
    txt3 += " Default: <current_folder>/explist/"
    abc.add_argument("--d_exp", help=txt3, metavar="")
    #
    txt4 = "Parent folder to harbor the bash files for remote copy."
    txt4 += " Each night in its own folder"
    txt4 += " Default: <current_folder>/bash_scp/"
    abc.add_argument("--d_bash", help=txt4, metavar="")
    #
    txt5 = "Number of exposures to be included in each bash file to be copied."
    txt5 += " Default: 25"
    abc.add_argument("--N", help=txt5, metavar="", default=25, type=int)
    #
    user_aux = os.getlogin()
    txt6 = "Username to be employed for connect to DESDM."
    txt6 += " Default: actual username, {0}".format(user_aux)
    abc.add_argument("--user", help=txt6, metavar="", default=user_aux)
    #
    txt7 = "Prefix to be used on the written exposure lists. The default is"
    txt7 += " explist, so the final names are explist_<nite>_<hh:mm:ss>.csv"
    txt7 += " where hh:mm:ss is the time at which the query was saved"
    abc.add_argument("--pref", help=txt7, metavar="")
    #
    time_aux = 30.
    txt8 = "Minimum exposure time for the images."
    txt8 += " Default: {0}".format(time_aux)
    abc.add_argument("--exptime", help=txt8, metavar="", default=time_aux,
                     type=float)
    #
    teff_g_aux = 0.2
    txt9 = "Minimum T_EFF for g-band. Default: {0}".format(teff_g_aux)
    abc.add_argument("--teff_g", help=txt9, metavar="", default=teff_g_aux,
                     type=float)
    #
    teff_riz_aux = 0.3
    txt10 = "Minimum T_EFF for r, i, and z-bands."
    txt10 += " Default: {0}".format(teff_riz_aux)
    abc.add_argument("--teff_riz", help=txt10, metavar="", 
                     default=teff_riz_aux, type=float)
    #
    txt11 = "Is this a test run? This flag allows to query only 5 exposures"
    abc.add_argument("--test", help=txt11, action="store_true")
    #
    txt12 = "Directory where to store the LOGs. One folder per night."
    txt12 = " Default: <current_directory>/logs"
    abc.add_argument("--dir_log", help=txt12, metavar="")
    # Recover args
    val = abc.parse_args()
    kw = dict()
    kw["username"] = val.user
    kw["nite"] = val.nite
    kw["exptime"] = val.exptime
    kw["Nexpnum"] = val.N
    kw["dir_bash"] = val.d_bash
    kw["dir_exp"] = val.d_exp
    kw["dir_immask"] = val.d_msk
    kw["dir_log"] = val.dir_log
    kw["prefix"] = val.pref
    kw["teff_g"] = val.teff_g
    kw["teff_riz"] = val.teff_riz
    kw["testing"] = val.test
    #
    # Calling the methods from mk_explist
    #
    DB = mk_explist.DBInfo(**kw)
    if False:
        print ', '.join("%s: %s" % item for item in vars(DB).items())
    #Setup log
    DB.setup_log()
    # Make queries, save tables and files
    ended_ok = DB.exp_mask()
    # Remote copy
    if ended_ok:
        DB.run_scp()
    # Here the remote copy finished

    #
    # Calling methods from submit_explist
    #
    # Save a plain text list of the copied files, to be used in case the 
    # sumbission fails and we need to only run diffimaging again. 
    backup_list = "immaskFiles_{0}_{1}.txt".format(DB.nite1, DB.hhmmss)
    with open(backup_list, "w+") as b:
        for im in DB.immask_files:
            b.write("{0}\n".format(im))
    # 
    kwin["dir_log"] = val.dir_log
    kwin["nite"] = DB.nite1
    kwin["im_list"] = DB.immask_files
    Proc = submit_explist.Work(**kwin)
    # Setup the logs
    Proc.setup_log()
    # Run processing
    Proc.runproc()
    logging.info("Processing ended ok. At{0}".format(time.ctime()))

    #
    # Calling deletion of tmp list sved in case of failure
    #
    # If all went ok, delete the backup file of the immask full paths
    subprocess.call(shlex.split("rm -v {0}".format(backup_list)), shell=True)
    logging.info("File {0} was deleted".format(backup_list))

    # Final message
    logging.info("Exiting. Work done.")

