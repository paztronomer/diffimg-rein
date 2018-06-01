''' Simple code to check if the files for a expnum are in disk or not
'''

import os
import sys
import glob
import argparse
import logging
import pandas as pd

# Setup logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO,
)


def get_args():
    t_gral = 'Simple code to see if exposures from a table are in the'
    t_gral += ' directory tree structure. The search will be performed'
    t_gral += ' under {root path}/{observing night}/{exposure number}/'
    abc = argparse.ArgumentParser(description=t_gral)
    h0 = 'Table containing NITE and EXPNUM'
    abc.add_argument('tab', help=h0, metavar='filename')
    ppath = '/pnfs/des/persistent/wsdiff/exp'
    h1 = 'Parent path to look for the files. Default: {0}'.format(ppath)
    abc.add_argument('-p', help=h1, metavar='/path', default=ppath)
    expr = 'immask'
    h2 = 'String to look for into the filenames. Default: {0}'.format(expr)
    abc.add_argument('-f', help=h2, metavar='filetype-like', default=expr)
    abc = abc.parse_args()
    return abc

def aux_main(tab_fnm=None, root=None, expr=None):
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
            aux_exist.append('False')
            aux_N.append(0)
            c2 += 1
    # Save results
    t_i = '{0} paths exists, {1} does not'.format(c1, c2)
    logging.info(t_i)
    res = pd.DataFrame({'path' : aux_path, 'nite' : aux_nite, 
                        'expnum' : aux_expnum, 'exists' : aux_exist,
                        'n_files' : aux_N})
    outnm = 'Results_CheckInDisk_PID{0}.csv'.format(os.getpid())
    res.to_csv(outnm, index=False, header=True)
    t_i = 'Saved file {0} with the results of the checking'.format(outnm)
    logging.info(t_i)


if __name__ == '__main__':
    args = get_args()
    aux_main(tab_fnm=args.tab, root=args.p, expr=args.f)
