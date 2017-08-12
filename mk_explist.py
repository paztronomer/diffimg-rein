""" Script to make a list of the new exposures, from the last night (or another). 
The condition is to be already processed by DESDM. As the processing goes through
the day, this code needs to run managed by a CRON.
Mind to be respectful with the files transfer, to not get into the processing path. 
Francisco Paz-Chinchon
"""

import os
import sys
import socket
import time
import datetime
import gc
import logging
import argparse
import itertools
import errno
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
# - make scp script for transfer
# - divide in chunks
# - keep a log
# - write a supra code to manage it
# =============================

class Toolbox():
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
            df_obj = connect.query_to_pandas(to_query)
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
    def __init__(self, nite=None):
        """ Method to feed relevant info
        Inputs
        - nite: str or int, consider it as last night, not today
        """
        self.hhmmss = datetime.datetime.today().strftime("%H:%M:%S")
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

    def exp_info(self, minEXPTIME=30, minTEFF_g=0.2, minTEFF_riz=0.3,
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
        qi += "  order by e.nite"
        #
        T = Toolbox()
        df0 = T.db_query(qi)
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
        df0 = df0.assign(MJD_OBS_ADD=mjd_aux)
        # Re-sort for nite and band
        df0.sort_values(["NITE", "BAND", "EXPNUM"], ascending=True, 
                        inplace=True)
        # Re-index
        df0 = df0.reset_index(drop=True)
        # Folder to save the explist csv files
        if parent_explist is None:
            parent_explist = os.path.join(os.getcwd(), "explist/")
        try:
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
            outnm = os.path.join(parent_explist, outnm)
        df0.to_csv(outnm, index=False, header=True)
        return df0

    def exp_mask(self, size_copy=25, root_path="/archive_data/desarchive",
                 parent_immask="/pnfs/des/persistent/wsdiff/exp",
                 parent_scp=None):
        """ Method to get the filenames and paths for the red_immask
        files associated to each CCD, for the previously selected exposures.
        After that, bash files are written.

        
        """
        # Folder to save the bash SCP files
        if parent_scp is None:
            parent_scp = os.path.join(os.getcwd(), "bash_scp/")
        try:
            os.makedirs(parent_scp)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
                logging.error("ERROR when creating {0}".format(parent_scp))
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
        #
        #
        dfpath = pd.read_csv("path.csv")
        #
        # Write bash SCP in packs of N exposures each. chunk_N() returns a
        # list of tuples
        expnum = dfexp["EXPNUM"].values
        Nexp = map(np.array, TT.chunk_N(expnum, size_copy))
        lineout = ["#!/bin/bash \n"]
        str0 = "scp fpazch@deslogin.cosmology.illinois.edu:" 
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
            lineout = ["#!/bin/bash \n"]
        return True
            

if __name__ == "__main__":
    print socket.gethostname()

    DB = DBInfo(nite=20170201)
    # DB.exp_mask(parent_immask="/Users/fco/Code/diffimg_des/des-diffimg-small")
    DB.exp_mask(parent_immask="/home/s1/fpazchin")
