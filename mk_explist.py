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
                 outnm=None):
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
        qi += "  fcut.fwhm_asec, fcut.t_eff, fcut.skybrightness"
        qi += "  from z, exposure e, firstcut_eval fcut, pfw_attempt_val val"
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
        qi += "  order by e.nite"
        T = Toolbox()
        df0 = T.db_query(qi)
        # Add a T_EFF condition for g>=0.2 and r,i,z>=0.3
        c1 = (df0["BAND"] == "g") & (df0["T_EFF"] < minTEFF_g)
        print np.any(c1.values)
        if np.any(c1.values):
            df0.drop(c1, inplace=True)
        for b in ["r", "i", "z"]:
            cx = (df0["BAND"] == b) & (df0["T_EFF"] < minTEFF_riz)
            print np.any(cx.values)
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
        # Write out the table, one slightly different filename for each 
        # time we query the DB
        if outnm is None:
            hhmmss = datetime.datetime.today().strftime("%H:%M:%S")
            outnm = "explist_{0}".format(self.nite1)
            outnm += "_{0}.csv".format(hhmmss)
        df0.to_csv(outnm, index=False, header=True)
        return df0


if __name__ == "__main__":
    print socket.gethostname()

    DB = DBInfo(nite=20170201)
    DB.exp_info()
