#!/usr/bin/env python
'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''


import logging
import multiprocessing as mp
import signal
import subprocess as sb
from functools import partial

def run_crispresso(crispresso_cmds,descriptor,idx):
    """
    Runs a specified crispresso command specified by idx
    Used for multiprocessing by run_crispresso_cmds
    input:
    crispresso_cmds: list of commands to run
    descriptor: label printed out describing a command e.g. "Could not process 'region' 5" or "Could not process 'batch' 5"
    idx: index of the command to run
    """
    crispresso_cmd=crispresso_cmds[idx]

    logging.info('Running CRISPResso on %s #%d/%d: %s' % (descriptor, idx+1, len(crispresso_cmds), crispresso_cmd))

    return_value = sb.call(crispresso_cmd,shell=True)


    if return_value != 0:
        logging.warn('CRISPResso command failed (return value %d) on %s #%d: "%s"'%(return_value,descriptor,idx,crispresso_cmd))
    else:
        logging.info('Finished CRISPResso %s #%d' %(descriptor,idx))
    return return_value

def run_crispresso_cmds(crispresso_cmds,n_processes=1,descriptor = 'region',continue_on_fail=False):
    """
    input: crispresso_cmds: list of crispresso commands to run
    info_cmd: info command from parent (logger)
    descriptor: label printed out describing a command e.g. "Could not process 'region' 5" or "Could not process 'batch' 5"
    """
    logging.info("Running CRISPResso with %d processes" % n_processes)
    pool = mp.Pool(processes = n_processes)
    idxs = range(len(crispresso_cmds))
    pFunc = partial(run_crispresso,crispresso_cmds,descriptor)

    #handle signals -- bug in python 2.7 (https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        res = pool.map_async(pFunc,idxs)
        ret_vals = res.get(60*60*10) # Without the timeout this blocking call ignores all signals.
        for idx, ret in enumerate(ret_vals):
            if ret != 0 and not continue_on_fail:
                raise Exception('CRISPResso %s #%d failed',descriptor,idx)
    except KeyboardInterrupt:
        pool.terminate()
        logging.warn('Caught SIGINT. Program Terminated')
        raise Exception('CRISPResso2 Terminated')
        exit (0)
    except Exception as e:
        print('CRISPResso2 failed')
        raise e
    else:
        plural = descriptor+"s"
        if descriptor.endswith("ch") or descriptor.endswith("sh"):
            plural = descriptor+"es"
        logging.info("Finished all " + plural)
        pool.close()
    pool.join()
