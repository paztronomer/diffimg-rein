""" Script to make a list of the new exposures, from the last night (or 
another).
The condition is to be already processed by DESDM. As the processing goes 
through the day, this code needs to run managed by a CRON.
Mind to be respectful with the files transfer, to not get into the processing 
path.
Francisco Paz-Chinchon
"""

import os
import sys
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
import numpy as np
import pandas as pd
from astropy.time import Time
try:
    import easyaccess as ea
    ea_import = True
except:
    try:
        import despydb.desdbi as desdbi
        ea_import = False
    except:
        errmsg = "No easyaccess neither despydb.desdbi available. Exiting"
        logging.error(errmsg)
        exit(1)

# =============================
# PENDINGS
# - keep a log
# - write a supra code to manage it
# =============================

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

    def to_path(self, parent=None, nite=None, expnum=None, reqnum=None,
                attnum=None, fnm=None, modify_fnm=False, str_run=None):
        """ Method to check for the existence of the destination folder, and
        to modify the filename used for save files
        Inputs
        - parent: root folder, /pnfs/des/persistent/wsdiff/exp/NITE/EXPNUM/
        - nite: night
        - expnum
        - reqnum
        - attnum
        - fnm: actual filename of immask files
        - modify_fnm: whether to change the ReqnumAttnum string on filename
        Returns
        - string with the destination path
        """
        # Modify filename
        if modify_fnm:
            aux = "r{0}p{1:02}".format(reqnum, attnum)
            fnm = fnm.replace(aux, str_run)
        # Check if directory exists, if not, then create it
        folder = os.path.join(parent, "{0}/{1:08}".format(nite, expnum))
        try:
            os.makedirs(folder)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
                logging.error("ERROR when creating {0}".format(folder))
        # Output the string with the final filename
        return os.path.join(folder, fnm)

    def chunk_N(self, y, size, fill_val=np.nan):
        """ Method to divide the list or 1D array in chunks of size N
        Inputs
        - data: list, tuple, array 1D to be chunked
        - size_n: number of elements of each chunk
        Returns
        - list of tuples containing the elements. When the last tuple has
        remaining spaces, fill with fill_val
        """
        args = [iter(y)] * size
        return list(itertools.izip_longest(*args, fillvalue=fill_val))

    def isot2mjd(self, isot, time_format="isot", time_scale="utc"):
        """ Method to transform from something like 2015-09-18T07:19:27.935065
        to MJD float. If a list is inputed, a list is returned
        Inputs
        - isot: string, in ISOT format
        - time_format: format for the time, as listed in the documentation
        - time_scale: scale of time (barycentric, utc, etc)
        Returns
        - float containing the mjd (or a list if a the input is multiple)
        """
        aux_t = Time(isot, format=time_format, scale=time_scale).mjd
        return aux_t

    def db_query(self, to_query, outdtype=None):
        """ Method to query the DB
        Inputs
        - to_query: str, contains the query. Do it changes the final SEMICOLON
        from easyaccess to desdbi?
        Returns
        - structured array
        """
        # What happens in the case of no rows? Test it
        if ea_import:
            # Needs: to_query
            connect = ea.connect("desoper")
            cursor = connect.cursor()
            try:
                df_obj = connect.query_to_pandas(to_query)
            except:
                logging.error("Error in querying\n\n\t{0}\n\n".format(to_query))
                exit(1)
            connect.close()
            return df_obj
            # Test if dtype works fine, if not, use zip and construct
            # the structured array scratch
            # return df_obj.to_records(index=False)
        else:
            # Needed variables: to_query, outdtype
            desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
            section = "db-desoper"
            dbi = desdbi.DesDbi(desfile, section)
            cursor = dbi.cursor()
            cursor.execute(to_query)
            cols = [line[0].lower() for line in cursor.description]
            rows = cursor.fetchall()
            #
            # change to pandas! is easier for collabs
            #
            exit(0)
            outtab = np.rec.array(rows, dtype=zip(cols, outdtype))
            return outtab


class DBInfo():
    def __init__(self, username=None, nite=None, exptime=None, Nexpnum=None,
                 dir_bash=None, dir_exp=None, dir_immask=None, dir_log=None,
                 prefix=None, teff_g=None, teff_riz=None, testing=None):
        """ Method to feed relevant info
        Inputs
        - username: str, user to connect to DESDM
        - nite: str or int, consider it as last night, not today
        - exptime: float, min time in seconds for select exposures
        - Nexpnum: int, number of exposures to be included in each separate
        bash file to copy
        - dir_bash: str, folder to save the bash files
        - dir_exp: str, folder to save the explist tables
        - dir_immask: str, folder where immask.fits.fz files will be stored
        - dir_log: str, folder where LOGs will be stored
        - prefix: str, prefix to the filename of the explists
        - teff_g, teff_riz: float, minimum values for T_EFF, g-band and
        r, i, z-bands
        - testing: boolean, if True, the only use 50 exposures
        """
        if nite is None:
            d1 = datetime.date.today() - datetime.timedelta(days=1)
            d2 = datetime.date.today()
            self.nite1 = d1.strftime("%Y%m%d")
            self.nite2 = d2.strftime("%Y%m%d")
        else:
            if isinstance(nite, str):
                pass
            else:
                nite = str(nite)
            d1 = datetime.datetime.strptime(nite, "%Y%m%d")
            d2 = d1 + datetime.timedelta(days=1)
            self.nite1 = d1.strftime("%Y%m%d")
            self.nite2 = d2.strftime("%Y%m%d")
        self.hhmmss = datetime.datetime.today().strftime("%H:%M:%S")
        self.username = username
        self.exptime = exptime
        self.Nexpnum = Nexpnum
        self.dir_bash = dir_bash
        self.dir_exp = dir_exp
        self.dir_immask = dir_immask
        self.dir_log = dir_log
        self.prefix = prefix
        self.teff_g = teff_g
        self.teff_riz = teff_riz
        self.bash_files = []
        self.testing = testing

    def setup_log(self):
        """ Method to setup the log output and start with some information
        """
        # Check/crete directory
        if self.dir_log is None:
            self.dir_log = os.path.join(os.getcwd(), "logs/")
        try:
            self.dir_log = os.path.join(self.dir_log, self.nite1)
            os.makedirs(self.dir_log)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
                logging.error("ERROR when creating {0}".format(self.dir_log))
        # Setup write out
        lognm = "explist_and_copy_{0}_{1}.log".format(self.nite1, self.hhmmss)
        logpath = os.path.join(self.dir_log, lognm)
        logging.basicConfig(filename=logpath, level=logging.DEBUG, 
                            format="%(asctime)s - %(levelname)s - %(message)s")
        # First information
        logging.info("\nRunning on: {0}\n".format(socket.gethostname()))
        logging.info("Script: {0}\n".format(os.path.basename(__file__)))

    def exp_info(self, minEXPTIME=None, minTEFF_g=None, minTEFF_riz=None,
                 parent_explist=None, outnm=None):
        """ Method to get information from the exposure, related to assessments
        from the firstcut processing and from the initial values coming from
        the telescope.
        The table results will be saves as *.csv. But remember! there are gonna
        be multiple tables per day, so make filename meaningful
        Inputs
        - minEXPTIME: minimum exposure time
        - minTEFF_g: minimum value of T_EFF for g-band
        - minTEFF_riz: minimum value of T_EFF for i, r, z bands
        - outnm: string used as filename for the output exposure list
        Returns
        - the constructed dataframe
        """
        minEXPTIME = self.exptime
        minTEFF_g = self.teff_g
        minTEFF_riz = self.teff_riz
        parent_explist = self.dir_exp
        outnm = self.prefix
        # Define the query which assumes last try to process as the valid
        qi = "with z as ("
        qi += "  select fcut.expnum, max(fcut.lastchanged_time) as evaltime"
        qi += "  from firstcut_eval fcut"
        qi += "  where fcut.analyst!='SNQUALITY'"
        qi += "  group by fcut.expnum"
        qi += "  )"
        qi += "  select e.expnum, e.nite, e.airmass, e.obstype, e.date_obs,"
        qi += "  e.mjd_obs, e.telra, e.teldec, e.radeg, e.decdeg, e.band,"
        qi += "  e.exptime, val.pfw_attempt_id,"
        qi += "  att.reqnum, att.attnum,"
        qi += "  fcut.fwhm_asec, fcut.t_eff, fcut.skybrightness"
        qi += "  from z, exposure e, firstcut_eval fcut, pfw_attempt_val val,"
        qi += "  pfw_attempt att"
        qi += "  where e.obstype='object'"
        qi += "  and e.exptime>={0}".format(minEXPTIME)
        qi += "  and e.nite between"
        qi += "  {0} and {1}".format(int(self.nite1), int(self.nite2))
        qi += "  and fcut.expnum=z.expnum"
        qi += "  and fcut.expnum=e.expnum"
        qi += "  and fcut.lastchanged_time=z.evaltime"
        qi += "  and fcut.program='survey'"
        qi += "  and fcut.accepted='True'"
        qi += "  and fcut.processed='True'"
        qi += "  and val.key='expnum'"
        qi += "  and to_number(val.val,'999999')=e.expnum"
        qi += "  and val.pfw_attempt_id=fcut.pfw_attempt_id"
        qi += "  and att.id=val.pfw_attempt_id"
        if self.testing:
            qi += "  and rownum<6"
        qi += "  order by e.nite"
        #
        T = Toolbox()
        df0 = T.db_query(qi)
        if (len(df0.index) == 0):
            noexp = "\tNo exposures were found for nite={0}".format(self.nite1)
            noexp += "\n\tExiting\n{0}".format("="*80)
            logging.warning(noexp)
            exit(0)
        # Add a T_EFF condition for g>=0.2 and r,i,z>=0.3
        c1 = (df0["BAND"] == "g") & (df0["T_EFF"] < minTEFF_g)
        if np.any(c1.values):
            df0.drop(c1, inplace=True)
        for b in ["r", "i", "z"]:
            cx = (df0["BAND"] == b) & (df0["T_EFF"] < minTEFF_riz)
            if np.any(cx.values):
                df0.drop(cx, inplace=True)
        # Drop rows where T_EFF is NaN
        df0.dropna(axis=0, subset=["T_EFF"], inplace=True)
        # Add a new column with a more precise MJD value
        mjd_aux = map(T.isot2mjd, df0["DATE_OBS"])
        # Pandas assign appeared on version 0.16
        if (float(str(pd.__version__)[:4]) >= 0.16):
            df0 = df0.assign(MJD_OBS_ADD=mjd_aux)
        else:
            df0.loc[:, "MJD_OBS_ADD"] = pd.Series(mjd_aux, index=df0.index)
        # Pandas sort_values appeared on v 0.17
        if (float(str(pd.__version__)[:4]) >= 0.17):
            # Re-sort for nite and band
            df0.sort_values(["NITE", "BAND", "EXPNUM"], ascending=True,
                            inplace=True)
        else:
            df0.sort(["NITE", "BAND", "EXPNUM"], ascending=True,
                     inplace=True)
        # Re-index
        df0 = df0.reset_index(drop=True)
        # Folder to save the explist csv files
        if parent_explist is None:
            parent_explist = os.path.join(os.getcwd(), "explist/")
        try:
            parent_explist = os.path.join(parent_explist, self.nite1)
            os.makedirs(parent_explist)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
                logging.error("ERROR when creating {0}".format(parent_explist))
        # Write out the table, one slightly different filename for each
        # time we query the DB
        if outnm is None:
            outnm = "explist_{0}".format(self.nite1)
            outnm += "_{0}.csv".format(self.hhmmss)
        else:
            outnm += "_{0}_{1}.csv".format(self.nite1, self.hhmmss)
        outnm = os.path.join(parent_explist, outnm)
        df0.to_csv(outnm, index=False, header=True)
        return df0

    def exp_mask(self, size_copy=None, root_path="/archive_data/desarchive",
                 parent_immask=None, parent_scp=None, desdm_user=None):
        """ Method to get the filenames and paths for the red_immask
        files associated to each CCD, for the previously selected exposures.
        After that, bash files are written.
        Inputs
        - size_copy: integer, number of EXPNUM to be used for transfer on each
        bash files
        - root_path: string, parent root path to the DESDM files
        - parent_immask: string, parent root onf Fermi where immask files will
        be stored
        - parent_scp: string, parent folder where bash files for SCP will be
        saved
        - desdm_user: user on the DESDM side, which will access the files.  If
        user is None, get the one in the session
        Outputs
        - boolean when end
        """
        size_copy = self.Nexpnum
        parent_immask = self.dir_immask
        parent_scp =self.dir_bash
        desdm_user = self.username
        #
        TT = Toolbox()
        # Get the explist
        dfexp = self.exp_info()
        # Need to query every pair pfw_attempt_id-expnum
        dfpath = pd.DataFrame()
        for index, row in dfexp.iterrows():
            gc.collect()
            qp = "select im.expnum, im.pfw_attempt_id, fai.path,"
            qp += "  fai.filename, fai.compression"
            qp += "  from image im, file_archive_info fai"
            qp += "  where im.pfw_attempt_id={0}".format(row["PFW_ATTEMPT_ID"])
            qp += "  and im.filetype='red_immask'"
            qp += "  and im.expnum={0}".format(row["EXPNUM"])
            qp += "  and fai.filename=im.filename"
            qp += "  order by fai.filename"
            dfaux = TT.db_query(qp)
            dfpath = dfpath.append(dfaux)
        # Save for testing
        # dfpath.to_csv("path.csv", index=False, header=True)
        # dfpath = pd.read_csv("path.csv")
        #
        # Folder to save the bash SCP files
        if parent_scp is None:
            parent_scp = os.path.join(os.getcwd(), "bash_scp/")
        try:
            parent_scp = os.path.join(parent_scp, self.nite1)
            os.makedirs(parent_scp)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
                logging.error("ERROR when creating {0}".format(parent_scp))
        # Write bash SCP in packs of N exposures each. chunk_N() returns a
        # list of tuples
        expnum = dfexp["EXPNUM"].values
        Nexp = map(np.array, TT.chunk_N(expnum, size_copy))
        lineout = ["#!/bin/bash \n"]
        str0 = "scp {0}@deslogin.cosmology.illinois.edu:".format(desdm_user)
        for write_exp in Nexp:
            # Remove the filling NaN
            write_exp = write_exp[np.logical_not(np.isnan(write_exp))]
            # Account for possible float
            write_exp = np.array(map(int, write_exp))
            for idx, row in dfpath.iterrows():
                if row["EXPNUM"] in write_exp:
                    cond = ((dfexp["EXPNUM"] == row["EXPNUM"]) &
                            (dfexp["PFW_ATTEMPT_ID"] == row["PFW_ATTEMPT_ID"]))
                    # Pending: check for unique req, att
                    req = dfexp.loc[cond, "REQNUM"].values[0]
                    att = dfexp.loc[cond, "ATTNUM"].values[0]
                    aux_fnm = row["FILENAME"] + row["COMPRESSION"]
                    destin = TT.to_path(parent=parent_immask,
                                        nite=self.nite1,
                                        expnum=row["EXPNUM"],
                                        fnm=aux_fnm,
                                        reqnum=req,
                                        attnum=att,
                                        modify_fnm=True,
                                        str_run="r4p4")
                    argu = [root_path, row["PATH"], row["FILENAME"]]
                    tmp = str0
                    tmp += os.path.join(*argu)
                    tmp += row["COMPRESSION"] + " "
                    tmp += destin
                    tmp += "\n"
                    lineout.append(tmp)
            # Write out the chunk files to be copied
            outfnm = "copy_{0}_{1}t{2}.sh".format(self.nite1, write_exp[0],
                                                    write_exp[-1])
            outfnm = os.path.join(parent_scp, outfnm)
            with open(outfnm, "w+") as f:
                f.writelines(lineout)
            logging.info("\twritten bash file {0}".format(outfnm))
            self.bash_files.append(outfnm)
            lineout = ["#!/bin/bash \n"]
        return True

    def run_scp(self):
        """ Method to run the created bash files for remote copy
        """
        Tbox = Toolbox()
        if (len(self.bash_files) > 0):
            for sh in self.bash_files:
                # Make the file executable
                root_path, fname= Tbox.split_path(sh)
                cmd = "chmod +x {0}".format(sh)
                pA = subprocess.call(shlex.split(cmd))
                # Run the bash file
                cmdscp = shlex.split("bash {0}".format(sh))
                pB = subprocess.Popen(cmdscp,
                                      shell=False,
                                      stdin=subprocess.PIPE,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True)
                # outM, errM = pB.communicate()
                pB.wait()
                logging.info("Ended run of {0}".format(fname))
        return True


if __name__ == "__main__":
    # Parse of arguments
    intro = "Script to detect the last night (or other) already processed"
    intro += " exposures, and then create bash executable files for remote"
    intro += " copy."
    abc = argparse.ArgumentParser(description=intro)
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
    # Calling
    DB = DBInfo(**kw)
    #Setup log
    DB.setup_log()
    # Make queries, save tables and files
    ended_ok = DB.exp_mask()
    # Remote copy
    if ended_ok:
        DB.run_scp()
