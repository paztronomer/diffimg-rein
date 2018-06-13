''' Script to make a list of the new exposures, from the last night (or
another).
The condition is to be already processed by DESDM. As the processing goes
through the day, this code needs to run managed by a CRON.
Mind to be respectful with the files transfer, to not get into the processing
path.
Francisco Paz-Chinchon
'''

import os
import sys
import stat
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
        errmsg = 'No easyaccess neither despydb.desdbi available. Exiting'
        logging.error(errmsg)
        exit(1)

# =============================
# PENDINGS
# - keep a log of processed files
# - create a list of the used self.* variables, per class
# - capture all exceptions, print them besides
#   except: # catch *all* exceptions
#       e = sys.exc_info()[0]
# =============================

class Toolbox():
    def split_path(self, path):
        '''Method to return the relative root folder (one level upper),
        given a path.
        Inputs
        - path: complete path
        Returns
        - 2 strings, one being the parent folder, and the filename
        '''
        #relat_root = os.path.abspath(os.path.join(path,os.pardir))
        relroot, filename = os.path.split(path)
        return (relroot, filename)

    def to_path(self, parent=None, nite=None, expnum=None, reqnum=None,
                attnum=None, fnm=None, modify_fnm=False, str_run=None,
                rwx_mode=774):
        ''' Method to check for the existence of the destination folder, and
        to modify the filename used for save files
        Inputs
        - parent: root folder, /pnfs/des/persistent/wsdiff/exp/NITE/EXPNUM/
        - nite: night
        - expnum
        - reqnum
        - attnum
        - fnm: actual filename of immask files
        - modify_fnm: whether to change the ReqnumAttnum string on filename
        - str_run: string to replace the ReqnumAttnum if modify_fn=True. It is
        'r4p4' for Y5N
        Returns
        - string with the destination path
        '''
        # Modify filename
        if modify_fnm:
            # 1) Change reqnum, attnum
            #
            # Check if really attnum is zero-padded
            #
            aux = 'r{0}p{1:02}'.format(reqnum, attnum)
            # First, check if the string exists in the filename
            if (fnm.find(aux) >= 0):
                fnm = fnm.replace(aux, str_run)
            else:
                err_aux = 'Error in filename modification: {0}'.format(fnm)
                err_aux += ' String {0} was not found'.format(aux)
                logging.error(err_aux)
                exit(1)
                #
                # Important: here I'm exiting, but other solutions should be
                # better
                #
            # 2) Change _c{ccdnum}_ by _{ccdnum}_
            if (fnm.find('_c') >= 0):
                fnm = fnm.replace('_c', '_')
            else:
                err_aux = 'Error in filename modification: {0}'.format(fnm)
                err_aux += ' String {0} was not found'.format('_c{CCDNUM}')
                logging.error(err_aux)
                exit(1)
            # 3) Change immaked by immask
            if (fnm.find('immasked') >= 0):
                fnm = fnm.replace('immasked', 'immask')
            else:
                err_aux = 'Error in filename modification: {0}'.format(fnm)
                err_aux += ' String {0} was not found'.format('immasked')
                logging.error(err_aux)
                exit(1)
        # Check if directory exists, if not, then create it
        folder = os.path.join(parent, '{0}/{1}'.format(nite, expnum))
        # If directory has rwx permissions, do nothing. If doesn't have 
        # permissions, try change them
        if (os.path.exists(folder) and not os.access(folder, os.R_OK)
            and not os.access(folder, os.W_OK) 
            and not os.access(folder, os.X_OK)):
            # If path exists, force the folder to have the permissions we need
            try:
                # chmod 774
                cmd = 'chmod {0} {1}'.format(rwx_mode, folder)
                #subprocess.call(shlex.split(cmd))
                pA = subprocess.Popen(shlex,split(cmd), stderr=subprocess.PIPE)
                err = pA.communicate()
                if err:
                    logging.error(err)
                    logging.error(pA.returncode)
            except:
                e = sys.exc_info()[0]
                logging.error(e)
        else:
            try:
                os.makedirs(folder)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
                    logging.error('ERROR when creating {0}'.format(folder))
        # Output the string with the final filename
        return os.path.join(folder, fnm)

    def chunk_N(self, y, size, fill_val=np.nan):
        ''' Method to divide the list or 1D array in chunks of size N
        Inputs
        - data: list, tuple, array 1D to be chunked
        - size_n: number of elements of each chunk
        Returns
        - list of tuples containing the elements. When the last tuple has
        remaining spaces, fill with fill_val
        '''
        args = [iter(y)] * size
        return list(itertools.izip_longest(*args, fillvalue=fill_val))

    def isot2mjd(self, isot, time_format='isot', time_scale='utc'):
        ''' Method to transform from something like 2015-09-18T07:19:27.935065
        to MJD float. If a list is inputed, a list is returned
        Inputs
        - isot: string, in ISOT format
        - time_format: format for the time, as listed in the documentation
        - time_scale: scale of time (barycentric, utc, etc)
        Returns
        - float containing the mjd (or a list if a the input is multiple)
        '''
        aux_t = Time(isot, format=time_format, scale=time_scale).mjd
        return aux_t

    def db_query(self, to_query, outdtype=None):
        ''' Method to query the DB
        Inputs
        - to_query: str, contains the query. Do it changes the final SEMICOLON
        from easyaccess to desdbi?
        Returns
        - structured array
        '''
        # What happens in the case of no rows? Test it
        if ea_import:
            # Needs: to_query
            connect = ea.connect('desoper')
            cursor = connect.cursor()
            try:
                df_obj = connect.query_to_pandas(to_query)
            except:
                t_e = 'Error in querying\n\n\t{0}\n\n'.format(to_query)
                logging.error(t_e)
                e = sys.exc_info()[0]
                logging.error(e)
                exit(1)
            connect.close()
            return df_obj
            # Test if dtype works fine, if not, use zip and construct
            # the structured array scratch
            # return df_obj.to_records(index=False)
        else:
            logging.warning('No easyaccess, will exit')
            # Needed variables: to_query, outdtype
            desfile = os.path.join(os.getenv('HOME'), '.desservices.ini')
            section = 'db-desoper'
            dbi = desdbi.DesDbi(desfile, section)
            cursor = dbi.cursor()
            cursor.execute(to_query)
            cols = [line[0].lower() for line in cursor.description]
            rows = cursor.fetchall()
            #
            # change to pandas! is easier for collabs
            #
            t_w = 'DESDBI Not implemented!'
            logging.warning(t_w)
            #
            #
            #
            exit(0)
            outtab = np.rec.array(rows, dtype=zip(cols, outdtype))
            return outtab


class DBInfo():
    def __init__(self, username=None, nite=None, expnum_fnm=None,
                 exptime=None, Nexpnum=None,
                 dir_bash=None, dir_exp=None, dir_immask=None, dir_log=None,
                 prefix=None, teff_g=None, teff_riz=None, ra_range=None,
                 dec_range=None, testing=None, rNpN=None):
        ''' Method to feed relevant info
        Inputs
        - username: str, user to connect to DESDM
        - nite: str or int, consider it as last night, not today
        - expnum_fnm: filename of the list of exposure numbers
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
        - ra_range: list of 2 RA min and max, in degrees
        - dec_range: list of 2 DEC min and max, in degrees
        - testing: boolean, if True, the only use 50 exposures
        '''
        if (nite is None) and (expnum_fnm is None):
            d1 = datetime.date.today() - datetime.timedelta(days=1)
            d2 = datetime.date.today()
            self.nite1 = d1.strftime('%Y%m%d')
            self.nite2 = d2.strftime('%Y%m%d')
            self.expnum_df = None
        elif (nite is not None) and (expnum_fnm is None):
            # Check if nite was string, if not the case, convert it
            if isinstance(nite, str):
                pass
            else:
                nite = str(nite)
            d1 = datetime.datetime.strptime(nite, '%Y%m%d')
            d2 = d1 + datetime.timedelta(days=1)
            self.nite1 = d1.strftime('%Y%m%d')
            self.nite2 = d2.strftime('%Y%m%d')
            self.expnum_df = None
        elif (nite is None) and (expnum_fnm is not None):
            self.nite1 = None
            self.nite2 = None
            explist = pd.read_table(expnum_fnm, names=['EXPNUM'], comment='#')
            if explist.isnull().values.any():
                explist = explist.dropna(how='all')
                explist.reset_index(drop=True, inplace=True)
                t_w = 'Input explist file {0} contains NaN'.format(expnum_fnm)
                t_w += ' . NaN were removed'
                logging.warning(t_w)
            self.expnum_df = explist
            # Define auxiliary naming where exposure list is input
            xmin = self.expnum_df['EXPNUM'].min()
            xmax = self.expnum_df['EXPNUM'].max()
            xlen = len(self.expnum_df['EXPNUM'].index)
            self.fnm_expnum = '{0}t{1}n{2}'.format(xmin, xmax, xlen)
        self.hhmmss = datetime.datetime.today().strftime('%Hh%Mm%Ss')
        self.username = username
        self.exptime = exptime
        self.Nexpnum = Nexpnum
        self.dir_immask = dir_immask
        # Folder to save the bash SCP files
        if (dir_bash is None):
            self.dir_bash = os.path.join(os.getcwd(), 'bash_scp/')
        else:
            self.dir_bash = dir_bash
        # Folder to save the explist csv files
        if (dir_exp is None):
            self.dir_exp = os.path.join(os.getcwd(), 'explist/')
        else:
            self.dir_exp = dir_exp
        # Check/assign LOG directory
        if (dir_log is None):
            self.dir_log = os.path.join(os.getcwd(), 'logs/')
        else:
            self.dir_log = dir_log
        # Check/assign EXPLIST directory
        if (prefix is None):
            self.prefix = 'explist'
        else:
            self.prefix = prefix
        self.teff_g = teff_g
        self.teff_riz = teff_riz
        self.ra_range = ra_range
        self.dec_range = dec_range
        self.testing = testing
        self.rNpN = rNpN
        # Lists to keep track of the bash files and of the copied immask.fits
        self.bash_files = []
        self.immask_files = []
        self.aux_parent_explist = None

    def setup_log(self):
        ''' Method to setup the log output and start with some information
        '''
        # Need to take care of the 2 options (so far), to use a single night
        # each time, or to input a list of expnum
        if ((self.nite1 is not None) and (self.expnum_df is None)):
            # Check/create directory
            try:
                self.dir_log = os.path.join(self.dir_log, self.nite1)
                os.makedirs(self.dir_log)
            except OSError as exception:
                if (exception.errno != errno.EEXIST):
                    raise
                    t_e = 'ERROR when creating {0}'.format(self.dir_log)
                    logging.error(t_e)
            # Setup write out
            lognm = 'explist_and_copy_{0}_{1}.log'.format(self.nite1,
                                                          self.hhmmss)
            logpath = os.path.join(self.dir_log, lognm)
        elif ((self.nite1 is None) and (self.expnum_df is not None)):
            # Check/create directory, using minimum and maximum expnum
            try:
                self.dir_log = os.path.join(self.dir_log, self.fnm_expnum)
                os.makedirs(self.dir_log)
            except OSError as exception:
                if (exception.errno != errno.EEXIST):
                    raise
                    t_e = 'ERROR when creating {0}'.format(self.dir_log)
                    logging.error(t_e)
            # Setup write out
            lognm = 'explist_and_copy_{0}_{1}.log'.format(self.fnm_expnum,
                                                          self.hhmmss)
            logpath = os.path.join(self.dir_log, lognm)
        logging.basicConfig(
            filename=logpath, level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        # First information
        logging.info('\nRunning on: {0}\n'.format(socket.gethostname()))
        logging.info('Script: {0}\n'.format(os.path.basename(__file__)))

    def exp_info(self, minEXPTIME=None, minTEFF_g=None, minTEFF_riz=None,
                 parent_explist=None, outnm=None, query_len=100):
        ''' Method to get information from the exposure, related to assessments
        from the firstcut processing and from the initial values coming from
        the telescope.
        The table results will be saves as *.csv. But remember! there are gonna
        be multiple tables per day, so make filename meaningful
        Important: each time the query is performed, it saves cummulative
        results, so each time is bigger
        Inputs
        - minEXPTIME: minimum exposure time
        - minTEFF_g: minimum value of T_EFF for g-band
        - minTEFF_riz: minimum value of T_EFF for i, r, z bands
        - outnm: string used as filename for the output exposure list
        Returns
        - the constructed dataframe
        '''
        Tbox = Toolbox()
        minEXPTIME = self.exptime
        minTEFF_g = self.teff_g
        minTEFF_riz = self.teff_riz
        parent_explist = self.dir_exp
        # Auxiliary function to write query for given exposures
        def query_exp(exp, test):
            exp = ','.join(map(str, exp))
            qi = 'with z as ('
            qi += '  select fcut.expnum,'
            qi += '  max(fcut.lastchanged_time) as evaltime'
            qi += '  from firstcut_eval fcut'
            qi += '  where fcut.analyst!=\'SNQUALITY\''
            qi += '  group by fcut.expnum'
            qi += '  )'
            qi += '  select e.expnum, e.nite, e.airmass, e.obstype,'
            qi += '  e.date_obs,'
            qi += '  e.mjd_obs, e.telra, e.teldec, e.radeg, e.decdeg,'
            qi += '  e.band,'
            qi += '  e.exptime, val.pfw_attempt_id,'
            qi += '  att.reqnum, att.attnum,'
            qi += '  fcut.fwhm_asec, fcut.t_eff, fcut.skybrightness'
            qi += '  from z, exposure e, firstcut_eval fcut,'
            qi += '  pfw_attempt_val val, pfw_attempt att'
            qi += '  where e.obstype=\'object\''
            qi += '  and e.expnum in ({0})'.format(exp)
            qi += '  and fcut.expnum=z.expnum'
            qi += '  and fcut.expnum=e.expnum'
            qi += '  and fcut.lastchanged_time=z.evaltime'
            qi += '  and fcut.processed=\'True\''
            qi += '  and val.key=\'expnum\''
            qi += '  and to_number(val.val,\'999999\')=e.expnum'
            qi += '  and val.pfw_attempt_id=fcut.pfw_attempt_id'
            qi += '  and att.id=val.pfw_attempt_id'
            if test:
                qi += '  and rownum<6'
            qi += '  order by e.nite'
            return qi
        # Define the query which assumes last try to process as the valid
        # 2 cases: nite1, nite2 are input, and the other when expnum_df
        #
        #
        # HERE! do the query in a more elegant way, for the 2 cases. Loop for
        # when N > 1000
        #
        if ((self.nite1 is not None) and (self.nite2 is not None)
            and (self.expnum_df is None)):
            qi = 'with z as ('
            qi += '  select fcut.expnum,'
            qi += '  max(fcut.lastchanged_time) as evaltime'
            qi += '  from firstcut_eval fcut'
            qi += '  where fcut.analyst!=\'SNQUALITY\''
            qi += '  group by fcut.expnum'
            qi += '  )'
            qi += '  select e.expnum, e.nite, e.airmass, e.obstype,'
            qi += '  e.date_obs,'
            qi += '  e.mjd_obs, e.telra, e.teldec, e.radeg, e.decdeg, e.band,'
            qi += '  e.exptime, val.pfw_attempt_id,'
            qi += '  att.reqnum, att.attnum,'
            qi += '  fcut.fwhm_asec, fcut.t_eff, fcut.skybrightness'
            qi += '  from z, exposure e, firstcut_eval fcut,'
            qi += '  pfw_attempt_val val, pfw_attempt att'
            qi += '  where e.obstype=\'object\''
            qi += '  and e.exptime>={0}'.format(minEXPTIME)
            qi += '  and e.nite between'
            qi += ' {0} and {1}'.format(int(self.nite1), int(self.nite2))
            if not (self.ra_range is None):
                qi += ' and e.radeg between'
                qi += ' {0} and {1}'.format(*self.ra_range)
            if not (self.dec_range is None):
                qi += ' and e.decdeg between'
                qi += ' {0} and {1}'.format(*self.dec_range)
            qi += '  and fcut.expnum=z.expnum'
            qi += '  and fcut.expnum=e.expnum'
            qi += '  and fcut.lastchanged_time=z.evaltime'
            qi += '  and fcut.program=\'survey\''
            qi += '  and fcut.accepted=\'True\''
            qi += '  and fcut.processed=\'True\''
            qi += '  and val.key=\'expnum\''
            qi += '  and to_number(val.val,\'999999\')=e.expnum'
            qi += '  and val.pfw_attempt_id=fcut.pfw_attempt_id'
            qi += '  and att.id=val.pfw_attempt_id'
            if self.testing:
                qi += '  and rownum<6'
            qi += '  order by e.nite'
            df0 = T.db_query(qi)
        elif ((self.expnum_df is not None) and (self.nite1 is None)
              and (self.nite2 is None)):
            t_w = 'For the list of EXPNUM, RA and DEC constraints will not be'
            t_w += ' applied. Neither minimum EXPTIME, nor ACCEPTED,'
            t_w += ' nor PROGRAM'
            logging.warning(t_w)
            # If number of exposures is larger than 100, split in small pieces
            if (len(self.expnum_df.index) >= query_len):
                # Produce a list of arrays
                Nx = int(np.ceil(len(self.expnum_df.index) / float(query_len)))
                listN = np.array_split(self.expnum_df['EXPNUM'].values, Nx)
                # Query each set
                for idx, x in enumerate(listN):
                    qn = query_exp(x, self.testing)
                    if (idx == 0):
                        df0 = Tbox.db_query(qn)
                    df0 = pd.concat([df0, Tbox.db_query(qn)])
                df0.reset_index(drop=True, inplace=True)
            else:
                qj = query_exp(self.expnum_df['EXPNUM'].values, self.testing)
                df0 = Tbox.db_query(qj)
        #
        if (len(df0.index) == 0):
            noexp = '\tNo exposures were found for nite={0}'.format(self.nite1)
            noexp += '\n\tExiting\n{0}'.format('='*80)
            logging.warning(noexp)
            exit(0)
        # T_EFF condition for g>=min_teff_g and r,i,z>=min_teff_riz
        # g-band condition
        c1 = (df0['BAND'] == 'g') & (df0['T_EFF'] < minTEFF_g)
        if np.any(c1.values):
            df0.drop(df0[c1].index, inplace=True)
        for b in ['r', 'i', 'z']:
            cx = (df0['BAND'] == b) & (df0['T_EFF'] < minTEFF_riz)
            if np.any(cx.values):
                df0.drop(df0[cx].index, inplace=True)
        # if ((self.nite1 is not None) and (self.nite2 is not None)
        #     and (self.expnum_df is None)):
        #     # T_EFF condition for g>=min_teff_g and r,i,z>=min_teff_riz
        #     c1 = (df0['BAND'] == 'g') & (df0['T_EFF'] < minTEFF_g)
        #     if np.any(c1.values):
        #         df0.drop(c1, inplace=True)
        #     for b in ['r', 'i', 'z']:
        #         cx = (df0['BAND'] == b) & (df0['T_EFF'] < minTEFF_riz)
        #         if np.any(cx.values):
        #             df0.drop(cx, inplace=True)
        # elif ((self.expnum_df is not None) and (self.nite1 is None)
        #       and (self.nite2 is None)):
        #     t_w = 'For the list of EXPNUM, no cut on T_EFF will be performed'
        #     logging.warning(t_w)
        t_i = 'In case T_EFF is NaN, EXPNUM will be removed'
        logging.info(t_i)
        # Drop rows where T_EFF is NaN
        df0.dropna(axis=0, subset=['T_EFF'], inplace=True)
        # Add a new column with a more precise MJD value
        mjd_aux = map(Tbox.isot2mjd, df0['DATE_OBS'])
        # Pandas assign appeared on version 0.16
        if (float(str(pd.__version__)[:4]) >= 0.16):
            df0 = df0.assign(MJD_OBS_ADD=mjd_aux)
        else:
            df0.loc[:, 'MJD_OBS_ADD'] = pd.Series(mjd_aux, index=df0.index)
        # Pandas sort_values appeared on v 0.17
        if (float(str(pd.__version__)[:4]) >= 0.17):
            # Re-sort for nite and band
            df0.sort_values(['NITE', 'BAND', 'EXPNUM'], ascending=True,
                            inplace=True)
        else:
            df0.sort(['NITE', 'BAND', 'EXPNUM'], ascending=True,
                     inplace=True)
        # Re-index
        df0 = df0.reset_index(drop=True)
        # Check/create the EXPLIST directory
        if ((self.nite1 is not None) and (self.nite2 is not None)
            and (self.expnum_df is None)):
            try:
                parent_explist = os.path.join(parent_explist, self.nite1)
                os.makedirs(parent_explist)
            except OSError as exception:
                if (exception.errno != errno.EEXIST):
                    raise
                    t_e = 'ERROR when creating {0}'.format(parent_explist)
                    logging.error(t_e)
            # Store parent_explist for later
            self.aux_parent_explist = parent_explist
            # Write out the table, one slightly different filename for each
            # time we query the DB
            # The condition for prefix=None is now on the __init__
            aux_out= '{0}_{1}_{2}.csv'.format(self.prefix, self.nite1, 
                                              self.hhmmss)
        elif ((self.expnum_df is not None) and (self.nite1 is None)
              and (self.nite2 is None)):
            try:
                parent_explist = os.path.join(parent_explist, self.fnm_expnum)
                os.makedirs(parent_explist)
            except OSError as exception:
                if (exception.errno != errno.EEXIST):
                    raise
                    t_e = 'ERROR when creating {0}'.format(parent_explist)
                    logging.error(t_e)
            # Store parent_explist for later
            self.aux_parent_explist = parent_explist
            # Write out the table, one slightly different filename for each
            # time we query the DB
            # The condition for prefix=None is now on the __init__
            aux_out= '{0}_{1}_{2}.csv'.format(self.prefix, self.fnm_expnum, 
                                              self.hhmmss)
        #
        # Double check here past:  os.path.join(parent_explist, outnm)
        #
        outnm = os.path.join(parent_explist, aux_out)
        df0.to_csv(outnm, index=False, header=True)
        return df0

    def exp_mask(self, size_copy=None, root_path='/archive_data/desarchive',
                 parent_immask=None, parent_scp=None, desdm_user=None):
        ''' Method to get the filenames and paths for the red_immask
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
        - boolean when ends
        '''
        size_copy = self.Nexpnum
        parent_immask = self.dir_immask
        parent_scp = self.dir_bash
        desdm_user = self.username
        #
        TT = Toolbox()
        # Get the explist
        dfexp = self.exp_info()
        # Here check for which pairs of pfw_attempt_id-expnum has been already
        # saved in the BASH scp scripts from previous runs
        # Important: here a strong assumption is made, if an exposure has
        # been successfully processed, then all their CCDs has been successful
        # Check all previously saved tables for the night we're working, and
        # then continue with only the ones that are new
        #
        # Iterate over the folder harboring the tables from the past queries.
        # Construct a tmp dataframe to compare against
        #
        # NOTE: for the input of data we have implemented 2 options
        # 1) by the input night, in which case the naming schema uses the
        # nite string as part of the naming
        # 2) by the list of exposure, in which case the maximum and minimum
        # exposure numbers are used for naming
        if (self.nite1 is not None) and (self.expnum_df is None):
            aux_naming = self.nite1
        elif (self.nite1 is None) and (self.expnum_df is not None):
            aux_naming = self.fnm_expnum

        #if (self.nite1 is not None) and (self.expnum_df is None):
        dir_depth = 0
        dfcomp = pd.DataFrame()
        counter_files = 0
        for root, dirs, files in os.walk(self.aux_parent_explist):
            if root.count(os.sep) >= dir_depth:
                del dirs[:]
            for fnm in files:
                if ('{0}_{1}'.format(self.prefix, aux_naming) in fnm):
                    counter_files += 1
                    try:
                        aux_fnm = os.path.join(root, fnm)
                        tmp = pd.read_csv(aux_fnm, engine='python')
                        dfcomp = dfcomp.append(tmp)
                    except:
                        e = sys.exc_info()[0]
                        logging.error(e)
                        msg = 'Cannot load {0}'.format(fnm)
                        logging.error(msg)
        # If ONLY ONE file exists, then that file is copied to be
        # a newdata_ file
        #
        if (counter_files == 1):
            newaux = 'newdata_{0}'.format(aux_naming)
            newaux += '_{0}.csv'.format(self.hhmmss)
            newaux = os.path.join(root, newaux)
            orig = os.path.join(root, fnm)
            # Copy: explist to newdata
            # NOTE: copy2 preserves metadata
            shutil.copy2(orig, newaux)
            logging.info('Saving newdata: {0}'.format(newaux))
            #
        # Compare dfexp (list of retrieved expnum from the querys) with
        # dfcomp (all the info in the already saved files) and if there
        # are new entries, continue. If not, then stop and exit
        #
        # There is no need to concatenate [dfexp, dfcomp] because dfexp was
        # already written on disk, so it is already contained inside dfcomp
        df_tmp = dfcomp.reset_index(drop=True)
        t_i = 'Only {0} exposures passed the criteria'.format(len(df_tmp.index))
        logging.info(t_i)
        dfcomp = None
        #
        if df_tmp.empty:
            txt_emp = 'An error occurred when checking for new data. None was'
            txt_emp += ' encountered'
            logging.error(txt_emp)
            exit(1)
            '''
            Really exit here? can I use a less drastic way to end the script?
            '''
        else:
            # EUPS has only up to pandas 0.15, then I can not use
            # drop_duplicates with keep=False. The option keep=False
            # removes all duplicted items, keeping only the items with 1
            # occurence
            # Therefore, use a new method for version <= 0.17
            # Here I need to walk ONLY through the NON-duplicated expnum
            if (float(pd.__version__[:4]) >= 0.17):
                t_i = 'pandas version {0}'.format(pd.__version__)
                t_i += ' supports drop_duplicates with keep=False'
                logging.info(t_i)
                new_exp = df_tmp['EXPNUM'].drop_duplicates(keep=False)
                new_exp = list(df_tmp.values)
            else:
                t_i = 'pandas version {0}'.format(pd.__version__)
                t_i += ' does not supports drop_duplicates with'
                t_i += ' keep=False. Using alternative method instead'
                logging.info(t_i)
                counter = collections.Counter(df_tmp['EXPNUM'].values)
                uni_exp = np.unique(df_tmp['EXPNUM'].values)
                new_exp = []
                for e in uni_exp:
                    if (counter[e] == 1):
                        new_exp.append(e)
            if (len(new_exp) > 0):
                # Replace dataframe by the one containing only new entries
                dfexp = None
                dfexp = df_tmp[df_tmp['EXPNUM'].isin(new_exp)]
                dfexp = dfexp.reset_index(drop=True)
                # Save this because is containing ONLY new data
                # The first table to be saved was about 35 lines above
                newdt = 'newdata_{0}_{1}.csv'.format(aux_naming, self.hhmmss)
                newdt = os.path.join(self.aux_parent_explist, newdt)
                dfexp.to_csv(newdt, index=False, header=True)
                logging.info('Saving newdata: {0}'.format(newdt))
            else:
                dfexp = None
                msg_nonew = 'No new exposures at this time, compared with'
                msg_nonew += ' previous queries. Exiting'
                logging.warning(msg_nonew)
                exit(0)
                # Keep the above exit(0)
            # The above allows us to get the unique elements and thus use them
            # in the generation of bash files
            #
            #
            '''
            HERE!!!

            The below should be contained in another method, just devoted to
            create the bash scp files

            '''
            # Generation of BASH SCP scripts
            #
            #
            # Need to query every pair pfw_attempt_id-expnum
            dfpath = pd.DataFrame()
            for index, row in dfexp.iterrows():
                gc.collect()
                qp = 'select im.nite, im.expnum, im.pfw_attempt_id, fai.path,'
                qp += '  fai.filename, fai.compression'
                qp += '  from image im, file_archive_info fai'
                qp += '  where '
                qp += 'im.pfw_attempt_id={0}'.format(row['PFW_ATTEMPT_ID'])
                qp += '  and im.filetype=\'red_immask\''
                qp += '  and im.expnum={0}'.format(row['EXPNUM'])
                qp += '  and fai.filename=im.filename'
                qp += '  order by fai.filename'
                dfaux = TT.db_query(qp)
                if (len(dfaux.index) > 0):
                    dfpath = dfpath.append(dfaux)
                else:
                    t_w = 'Missing immask file,'
                    t_w += ' from EXPNUM={0}, with'.format(row['EXPNUM'])
                    t_w += ' PFW_ATTEMPT={0}'.format(row['PFW_ATTEMPT_ID'])
                    logging.warning(t_w)
            # Check if the result is empty or not
            if (dfpath.empty):
                t_e = 'No immask files were found. Exiting'
                logging.error(t_e)
                exit(1)
            #
            # Check/create BASH directory to save scp scripts
            try:
                parent_scp = os.path.join(parent_scp, aux_naming)
                os.makedirs(parent_scp)
            except OSError as exception:
                if (exception.errno != errno.EEXIST):
                    raise
                    logging.error('ERROR when creating {0}'.format(parent_scp))
            # Write bash SCP in packs of N exposures each. chunk_N() returns a
            # list of tuples
            expnum = dfexp['EXPNUM'].values
            Nexp = map(np.array, TT.chunk_N(expnum, size_copy))
            lineout = ['#!/bin/bash \n']
            str0 = 'scp'
            str0 += ' {0}@deslogin.cosmology.illinois.edu:'.format(desdm_user)
            for write_exp in Nexp:
                # Remove the filling NaN
                write_exp = write_exp[np.logical_not(np.isnan(write_exp))]
                # Account for possible float
                write_exp = np.array(map(int, write_exp))
                for idx, row in dfpath.iterrows():
                    if (row['EXPNUM'] in write_exp):
                        cond = (
                            (dfexp['EXPNUM'] == row['EXPNUM']) &
                            (dfexp['PFW_ATTEMPT_ID'] == row['PFW_ATTEMPT_ID'])
                        )
                        # Pending: check for unique req, att
                        req = dfexp.loc[cond, 'REQNUM'].values[0]
                        att = dfexp.loc[cond, 'ATTNUM'].values[0]
                        aux_fnm = row['FILENAME'] + row['COMPRESSION']
                        #
                        # Here the the change of immask filename is done
                        #
                        '''
                        HERE!! Need to change the way to obtain nite
                        '''
                        destin = TT.to_path(parent=parent_immask,
                                            nite=row['NITE'], #self.nite1,
                                            expnum=row['EXPNUM'],
                                            fnm=aux_fnm,
                                            reqnum=req,
                                            attnum=att,
                                            modify_fnm=True,
                                            str_run=self.rNpN)
                        argu = [root_path, row['PATH'], row['FILENAME']]
                        tmp = str0
                        tmp += os.path.join(*argu)
                        tmp += row['COMPRESSION'] + ' '
                        tmp += destin
                        tmp += '\n'
                        # Before write out the linr on the bash, check if
                        # file already exists. If true, then save the line,
                        # but commented
                        if os.path.exists(destin):
                            tmp = '# ' + tmp
                        lineout.append(tmp)
                        # Store the immask file entire path
                        self.immask_files.append(destin)
                # How many are already on disk?
                n_disk = [x[0] for x in lineout]
                n_disk = n_disk.count('#') - 1
                if (n_disk > 0):
                    comment_n = '# ' + '=' * 78 + '\n'
                    comment_n += '# ' + '=' * 78 + '\n'
                    comment_n += '# ' + '=' * 78 + '\n'
                    comment_n += '# Note: {0} files are'.format(n_disk)
                    comment_n += ' already of disk, from a total of'
                    comment_n += ' {0}'.format(len(lineout) - 1)
                    lineout.append(comment_n)
                # Write out the chunk files to be copied
                outfnm = 'copy_{0}_{1}t{2}.sh'.format(aux_naming, write_exp[0],
                                                      write_exp[-1])
                outfnm = os.path.join(parent_scp, outfnm)
                with open(outfnm, 'w+') as f:
                    f.writelines(lineout)
                logging.info('written bash file {0}'.format(outfnm))
                self.bash_files.append(outfnm)
                # Refresh lineout variable for the next iteration
                lineout = ['#!/bin/bash \n']
        return True

    def run_scp(self):
        ''' Method to run the created bash files for remote copy
        This method should not have an exit statement, and if fails should
        be able to record the failure on a file
        '''
        #
        # HERE
        # Issue, there is no checking if the copy was successfully completed.
        # Already occurred that the copy was not completed because has no
        # needed privileges, and no error was triggered.
        #
        Tbox = Toolbox()
        if (len(self.bash_files) > 0):
            for sh in self.bash_files:
                try:
                    # Make the file executable
                    root_path, fname= Tbox.split_path(sh)
                    cmd = 'chmod +x {0}'.format(sh)
                    pA = subprocess.call(shlex.split(cmd))
                except:
                    t_e = sys.exc_info()[0]
                    logging.error(t_e)
                try:
                    # Run the bash file
                    cmdscp = shlex.split('bash {0}'.format(sh))
                    # stdout=subprocess.PIPE,
                    # stdin=subprocess.PIPE,
                    pB = subprocess.Popen(cmdscp,
                                          shell=False,
                                          stderr=subprocess.PIPE,
                                          universal_newlines=True)
                    errM = pB.communicate()
                    pB.wait()
                    if errM[1]:
                        logging.error(errM[1])
                    logging.info('Ended run of {0}'.format(fname))
                except:
                    t_e = sys.exc_info()[0]
                    logging.error(t_e)
        return True


if __name__ == '__main__':
    # Parse of arguments
    intro = 'Script to detect the last night (or other) already processed'
    intro += ' exposures, and then create bash executable files for remote'
    intro += ' copy.'
    abc = argparse.ArgumentParser(description=intro)
    # Optional
    # Set the way the main set of science files will be selected: night based
    # or exposure number based. I defined a mutual-exclusive group for the
    # 2 options, 'breed'
    breed = abc.add_mutually_exclusive_group(required=True)
    txt1a = 'Night to be queried. RA and DEC constraints are applied. Default'
    txt1a += ' value is last night'
    breed.add_argument('--nite', help=txt1a, metavar='YYYYMMDD')
    txt1b = 'Set of exposure numbers to be used. Input a 1-columns text file'
    txt1b += ' having one EXPNUM per line.'
    txt1b += ' RA and DEC constraints are not applied.'
    breed.add_argument('--exp', help=txt1b, metavar='filename')
    #
    immask_aux = '/pnfs/des/persistent/wsdiff/exp'
    txt2 = 'Parent folder to harbor immask files per CCD.'
    txt2 += ' Default: {0}'.format(immask_aux)
    abc.add_argument('--d_msk', help=txt2, metavar='', default=immask_aux)
    #
    txt3 = 'Parent folder to harbor the exposure lists (various files).'
    txt3 += ' Each night has a different folder.'
    txt3 += ' Default: <current_folder>/explist/'
    abc.add_argument('--d_exp', help=txt3, metavar='')
    #
    txt4 = 'Parent folder to harbor the bash files for remote copy.'
    txt4 += ' Each night in its own folder'
    txt4 += ' Default: <current_folder>/bash_scp/'
    abc.add_argument('--d_bash', help=txt4, metavar='')
    #
    txt12 = 'Directory where to store the LOGs. One folder per night.'
    txt12 = ' Default: <current_directory>/logs'
    abc.add_argument('--d_log', help=txt12, metavar='')
    #
    Nmax = 20
    txt5 = 'Number of exposures to be included in each bash file to be copied.'
    txt5 += ' Default: {0}'.format(Nmax)
    abc.add_argument('--N', help=txt5, metavar='', default=Nmax, type=int)
    #
    # Reverse to os.getlogin()!
    #
    # user_aux = os.getlogin()
    user_aux = 'fpazch'
    txt6 = 'Username to be employed for connect to DESDM.'
    txt6 += ' Default: actual username, {0}'.format(user_aux)
    abc.add_argument('--user', help=txt6, metavar='', default=user_aux)
    #
    txt7 = 'Prefix to be used on the written exposure lists. The default is'
    txt7 += ' explist, so the final names are explist_<nite>_<HHhMMmSSs>.csv'
    txt7 += ' where HHh MMm SSs is the time at which the query was saved'
    abc.add_argument('--pref', help=txt7, metavar='')
    #
    time_aux = 30.
    txt8 = 'Minimum exposure time for the images.'
    txt8 += ' Default: {0}'.format(time_aux)
    abc.add_argument('--exptime', help=txt8, metavar='', default=time_aux,
                     type=float)
    #
    teff_g_aux = 0.2
    txt9 = 'Minimum T_EFF for g-band. Default: {0}'.format(teff_g_aux)
    abc.add_argument('--teff_g', help=txt9, metavar='', default=teff_g_aux,
                     type=float)
    #
    teff_riz_aux = 0.3
    txt10 = 'Minimum T_EFF for r, i, and z-bands.'
    txt10 += ' Default: {0}'.format(teff_riz_aux)
    abc.add_argument('--teff_riz', help=txt10, metavar='',
                     default=teff_riz_aux, type=float)
    #
    tmp_ra = [20., 40.]
    txt11 = 'Minimum and maximum RA (degrees), separated by space.'
    txt11 += ' Default: {0} {1}'.format(*tmp_ra)
    abc.add_argument('--ra', help=txt11, nargs=2, type=float, default=tmp_ra)
    #
    tmp_dec = [-90., -10.]
    txt12 = 'Minimum and maximum DEC (degrees), separated by space.'
    txt12 += ' Default: {0} {1}'.format(*tmp_dec)
    abc.add_argument('--dec', help=txt12, nargs=2, type=float, default=tmp_dec)
    #
    aux_rp = 'r4p4'
    txt13 = 'For renaming the \'filetype=red_immask\' files, which rNpN to'
    txt13 += ' use. Default: {0}'.format(aux_rp)
    abc.add_argument('--rp', help=txt13, default=aux_rp)
    #
    txtN = 'Is this a test run? This flag allows to query only 5 exposures'
    abc.add_argument('--test', help=txtN, action='store_true')
    # Recover args
    val = abc.parse_args()
    kw = dict()
    kw['username'] = val.user
    kw['nite'] = val.nite
    kw['expnum_fnm'] = val.exp
    kw['exptime'] = val.exptime
    kw['Nexpnum'] = val.N
    kw['dir_bash'] = val.d_bash
    kw['dir_exp'] = val.d_exp
    kw['dir_immask'] = val.d_msk
    kw['dir_log'] = val.d_log
    kw['prefix'] = val.pref
    kw['teff_g'] = val.teff_g
    kw['teff_riz'] = val.teff_riz
    kw['ra_range'] = val.ra
    kw['dec_range'] = val.dec
    kw['rNpN'] = val.rp
    kw['testing'] = val.test
    #
    # Calling
    DB = DBInfo(**kw)
    #Setup log
    DB.setup_log()
    t_i = 'Input parameters: {0}'.format(val)
    logging.info(t_i)
    #
    # Make queries, save tables and files
    ended_ok = DB.exp_mask()
    # Remote copy
    if ended_ok:
        logging.info('Run scp transfer')
        DB.run_scp()
    # Save a plain text list of the copied files, to be used in case the
    # sumbission fails and we need to only run query/scp again.
    if (DB.nite1 is None):
        aux_nm = DB.fnm_expnum
    else:
        aux_nm = DB.nite1
    backup_list = 'immaskFiles_{0}_{1}.txt'.format(aux_nm, DB.hhmmss)
    # Add a safe step here
    try:
        backup_list = os.path.join(DB.dir_log, backup_list)
        with open(backup_list, 'w+') as b:
            for im in DB.immask_files:
                b.write('{0}\n'.format(im))
    except Exception as e:
        logging.error(str(e))
        info1 = 'Trying to write {0}'.format(backup_list)
        info1 += ' in the current directory'
        logging.info(info1)
        with open(backup_list, 'w+') as b:
            for im in DB.immask_files:
                b.write('{0}\n'.format(im))
    logging.info('Backup file saved: {0}'.format(backup_list))
