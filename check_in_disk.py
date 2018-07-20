''' Simple code to check if the files for a expnum are in disk or not
It can also check if files in a ist are or not in disk
'''

import os
import sys
import glob
import argparse
import logging
import pandas as pd
import numpy as np

# Setup logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO,
)


def get_args():
    t_gral = 'Simple code to see if exposures from a table are in the'
    t_gral += ' directory tree structure. Option 1)  search will be performed'
    t_gral += ' under {root path}/{observing night}/{exposure number}/'
    t_gral += ' were NITE and EXPNUM must be given in an input table'
    t_gral += ' Option 2) Check the full path given as input in a list.'
    t_gral += ' These options are mutually exclusive'
    abc = argparse.ArgumentParser(description=t_gral)
    h0 = '<op1> Table containing NITE and EXPNUM. Header shoud contain'
    h0 += ' NITE and EXPNUM strings (capitalized and not commented)'
    abc.add_argument('-t', help=h0, metavar='filename')
    ppath = '/pnfs/des/persistent/wsdiff/exp'
    h1 = '<op1> Parent path to look for the files. Default: {0}'.format(ppath)
    abc.add_argument('-p', help=h1, metavar='/path', default=ppath)
    expr = 'immask'
    h2 = '<op1> String to look for into the filenames.'
    h2 += ' Default: {0}'.format(expr)
    abc.add_argument('-f', help=h2, metavar='filetype-like', default=expr)
    h3 = '<op2> Filename of the table containing the paths to search for'
    abc.add_argument('-l', help=h3, metavar='filename')
    h4 = 'Output filename. Default for op1/op2:'
    h4 += ' checkFiles_{summary/missing}_PID{nnnn}.csv'
    abc.add_argument('-o', help=h4, metavar='filename')
    abc = abc.parse_args()
    return abc

def aux_main(tab_fnm=None, root=None, expr=None, path_list=None, out_fnm=None):
    ''' This method works for 2 scenarios: either the input of NITE, EXPNUM
    or by a list of full paths
    '''
    if (path_list is None):
        kw = {
            'sep' : None,
            'comment' : '#',
            'engine' : 'python',
        }
        try:
            df = pd.read_table(tab_fnm, **kw)
        except:
            t_w = 'pandas{0} doesn\'t support guess sep'.format(pd.__version__)
            logging.warning(t_w)
            df = pd.read_csv(tab_fnm)
        # Loading csv table. Remenber it must have a header
        logging.info('Table {0} successfully loaded'.format(tab_fnm))
        f1 = lambda r, n, e: os.path.join(r, str(n), str(e))
        # aux = map(f1, [root] * len(df.index), df['NITE'], df['EXPNUM'])
        aux_path, aux_exist, aux_nite, aux_expnum, aux_N = [], [], [], [], []
        c1, c2 = 0, 0
        for idx, row in df.iterrows():
            p = f1(root, row['NITE'], row['EXPNUM'])
            # First check for existence
            aux_nite.append(row['NITE'])
            aux_expnum.append(row['EXPNUM'])
            aux_path.append(p)
            if os.path.exists(p):
                # The following gives a list of all the full path matches
                m = glob.glob(os.path.join(p, '*{0}*'.format(expr)))
                aux_exist.append('True')
                aux_N.append(len(m))
                c1 += 1
            else:
                logging.warning('Path doesn\'t exists: {0}'.format(p))
                aux_exist.append('False')
                aux_N.append(0)
                c2 += 1
        # Save results
        t_i = '{0} paths exists, {1} does not'.format(c1, c2)
        logging.info(t_i)
        res = pd.DataFrame({'path' : aux_path, 'nite' : aux_nite, 
                            'expnum' : aux_expnum, 'exists' : aux_exist,
                            'n_files' : aux_N})
        if (out_fnm is None):
            out_fnm = 'checkFiles_summary_PID{0}.csv'.format(os.getpid())
        res.to_csv(out_fnm, index=False, header=True)
        t_i = 'Saved file {0} with the results of the checking'.format(out_fnm)
        logging.info(t_i)
    elif (path_list is not None):
        kw = {
            'sep' : None,
            'comment' : '#',
            'engine' : 'python',
            'names' : ['PATH'],
        }
        try:
            df = pd.read_table(path_list, **kw)
        except:
            t_w = 'pandas{0} doesn\'t support guess sep'.format(pd.__version__)
            logging.warning(t_w)
            df = pd.read_csv(path_list, names=['PATH'])
        #
        df['PATH'].replace('', np.nan, inplace=True)
        df.dropna(subset=['PATH'], inplace=True)
        df.to_csv(path_list, index=False, header=False)
        #
        aux_miss = []
        for idx, row in df.iterrows():
            if os.path.exists(row['PATH']):
                pass
            else:
                t_w = '{0} does not exists'.format(row['PATH'])
                logging.warning(t_w)
                aux_miss.append(row['PATH'])
        dfz = pd.DataFrame({'path': aux_miss})
        if (out_fnm is None):
            out_fnm = 'checkFiles_missing_PID{0}.csv'.format(os.getpid())
        dfz.to_csv(out_fnm, index=False, header=False)
    return True

if __name__ == '__main__':
    args = get_args()
    aux_main(tab_fnm=args.t, root=args.p, expr=args.f, path_list=args.l,
             out_fnm=args.o)
