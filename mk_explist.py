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
import numpy as np
import pandas as pd
from astropy import Time
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
# - Do the query, incorporating FCUT_EVAL. 
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
            # Test if dtype works fine, if not, use zip and construct 
            # the structured array scratch
            return df_obj.to_records(index=False)
        else:
            # Needed variables: to_query, outdtype
            desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
            section = "db-desoper"
            dbi = desdbi.DesDbi(desfile, section)
            cursor = dbi.cursor()
            cursor.execute(to_query)
            cols = [line[0].lower() for line in cursor.description]
            rows = cursor.fetchall()
            outtab = np.rec.array(rows, dtype=zip(cols, outdtype))
            return outtab


class DBInfo():
    def __init__(self, nite=None):
        """ Method to feed relevant info
        Inputs
        - nite: str or int, consider it as last night, not today
        """
        if nite is None:
            d1 = datetime.date.today() - datetime.timedelta(day=1)
            d2 = datetime.date.today()
            self.nite1 = d1.strftime("%Y%m%d")
            self.nite2 = d2.strftime("%Y%m%d")
        else:
            if isinstance(nite, str):
               aux_tday = nite
            else:
                aux_tday = str(nite)
            d1 = datetime.datetime.strptime(nite, "%Y%m%d") 
            d2 = d1 + datetime.timedelta(days=1)
            self.nite1 = d1.strftime("%Y%m%d")
            self.nite2 = d2.strftime("%Y%m%d")

    def query_img(self, ):
        """ 
        * Mind ea_import ddboolean flag
        * Save in *.csv
        """
        pass

if __name__ == "__main__":
    print socket.gethostname()
