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
from inspect import getmodule, stack
import numpy as np
import pandas as pd
import traceback

def get_max_processes():
    return mp.cpu_count()

def run_crispresso(crispresso_cmds, descriptor, idx):
    """
    Runs a specified crispresso command specified by idx
    Used for multiprocessing by run_crispresso_cmds
    input:
    crispresso_cmds: list of commands to run
    descriptor: label printed out describing a command e.g. "Could not process 'region' 5" or "Could not process 'batch' 5"
    idx: index of the command to run
    """
    crispresso_cmd=crispresso_cmds[idx]

    logging.info('Running CRISPResso on %s #%d/%d: %s' % (descriptor, idx, len(crispresso_cmds), crispresso_cmd))

    return_value = sb.call(crispresso_cmd, shell=True)

    if return_value == 137:
        logging.warn('CRISPResso was killed by your system (return value %d) on %s #%d: "%s"\nPlease reduce the number of processes (-p) and run again.'%(return_value, descriptor, idx, crispresso_cmd))
    elif return_value != 0:
        logging.warn('CRISPResso command failed (return value %d) on %s #%d: "%s"'%(return_value, descriptor, idx, crispresso_cmd))
    else:
        logging.info('Finished CRISPResso %s #%d' %(descriptor, idx))
    return return_value


def wrapper(func, args):
    idx, args = args
    return (idx, func(args))


def run_crispresso_cmds(crispresso_cmds, n_processes="1", descriptor = 'region', continue_on_fail=False, start_end_percent=None, logger=None):
    """Run multiple CRISPResso commands in parallel.

    Parameters
    ----------
    crispresso_cmds: list
        The list of CRISPResso commands to run.
    n_processes: str
        The number of processes to use, if `max` it will use the maximum number
        of processors.
    descriptor: str
        The label printed out describing a command e.g. "Could not process
        'region' 5" or "Could not process 'batch' 5".
    continue_on_fail: bool
        If True, CRISPResso will continue to run even if one of the commands fails.
    start_end_percent: tuple
        The start and end percent of the commands to run. For example, if you
        have 10 commands to run and `start_end_percent` is (0, 100), then the
        percent_complete will be 10 when the first command is finished and 20
        when the second command is finished, etc.
    logger: logging.Logger | None
        The logger to use for logging. If None, the logger of the calling module
        is used.

    Returns
    -------
    None
    """
    if not crispresso_cmds:
        return

    if logger is None:
        # this line noise will get the name of the module from which this
        # function was called, and thereby the correct logger
        logger = logging.getLogger(getmodule(stack()[1][0]).__name__)

    int_n_processes = 1
    if n_processes == "max":
        int_n_processes = get_max_processes()
    else:
        int_n_processes = int(n_processes)

    logger.info("Running CRISPResso with %d processes" % int_n_processes)
    pool = mp.Pool(processes=int_n_processes)
    idxs = range(len(crispresso_cmds))
    ret_vals = [None] * len(crispresso_cmds)
    pFunc = partial(run_crispresso, crispresso_cmds, descriptor)
    p_wrapper = partial(wrapper, pFunc)
    if start_end_percent is not None:
        percent_complete_increment = start_end_percent[1] - start_end_percent[0]
        percent_complete_step = percent_complete_increment / len(crispresso_cmds)
        percent_complete = start_end_percent[0]
    else:
        percent_complete_step = 1 / len(crispresso_cmds)
        percent_complete = 0

    #handle signals -- bug in python 2.7 (https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        completed = 0
        for idx, res in pool.imap_unordered(p_wrapper, enumerate(idxs)):
            ret_vals[idx] = res
            completed += 1
            percent_complete += percent_complete_step
            logger.info(
                "Completed {0}/{1} runs".format(completed, len(crispresso_cmds)),
                {'percent_complete': percent_complete},
            )
        for idx, ret in enumerate(ret_vals):
            if ret == 137:
                raise Exception('CRISPResso %s #%d was killed by your system. Please decrease the number of processes (-p) and run again.'%(descriptor, idx))
            if ret != 0 and not continue_on_fail:
                raise Exception('CRISPResso %s #%d failed. For more information, try running the command: "%s"'%(descriptor, idx, crispresso_cmds[idx]))
    except KeyboardInterrupt:
        pool.terminate()
        logger.warn('Caught SIGINT. Program Terminated')
        raise Exception('CRISPResso2 Terminated')
        exit (0)
    except Exception as e:
        print('CRISPResso2 failed')
        raise e
    else:
        plural = descriptor+"s"
        if descriptor.endswith("ch") or descriptor.endswith("sh"):
            plural = descriptor+"es"
        logger.info("Finished all " + plural)
        pool.close()
    pool.join()

def run_pandas_apply_parallel(input_df, input_function_chunk, n_processes=1):
    """
    Runs a function on chunks of the input_df
    This is a little clunky, but seems to work better than serial runs
    The input_function should be a wrapper that takes in df chunks and applies a real function on the rows
    For example, the input_function should be:
    def input_function_chunk(df):
        return df.apply(input_function,axis=1)
    or this (apply runs the function twice on the first row -- this loop doesn't)
    def input_function_chunk(df):
        new_df = pd.DataFrame(columns=df.columns)
        for i in range(len(df)):
            new_df = new_df.append(input_function(df.iloc[i].copy()))
        return(new_df)

    The wrapped function should add to rows and return the row
    e.g.  def input_function(row):
            row['newval'] = row['oldval']*row['oldval2']
            return row
    This is useful for functions for which the operation on the row takes a significant amount of time and the overhead
        of dataframe manipulation is relatively small
    """
    #shuffle the dataset to avoid finishing all the ones on top while leaving the ones on the bottom unfinished
    n_splits = min(n_processes, len(input_df))
    df_split = np.array_split(input_df.sample(frac=1), n_splits)
    pool = mp.Pool(processes = n_splits)

    #handle signals -- bug in python 2.7 (https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        r = pool.map_async(input_function_chunk, df_split)
        results = r.get(60*60*60) # Without the timeout this blocking call ignores all signals.
        df_new = pd.concat(results)
    except KeyboardInterrupt:
        pool.terminate()
        logging.warn('Caught SIGINT. Program Terminated')
        raise Exception('CRISPResso2 Terminated')
        exit (0)
    except Exception as e:
        print('CRISPResso2 failed')
        raise e
    else:
        pool.close()
    pool.join()
    return df_new


def run_function_on_array_chunk_parallel(input_array, input_function, n_processes=1):
    """
    Runs a function on chunks of input_array
    input_array: array of values
    input_function: function to run on chunks of the array
        input_function should take in a smaller array of objects
    """
    pool = mp.Pool(processes = n_processes)

    #handle signals -- bug in python 2.7 (https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        n = int(max(10, len(input_array)/n_processes)) #don't parallelize unless at least 10 tasks
        input_chunks = [input_array[i * n:(i + 1) * n] for i in range((len(input_array) + n - 1) // n )]
        r = pool.map_async(input_function, input_chunks)
        results = r.get(60*60*60) # Without the timeout this blocking call ignores all signals.
    except KeyboardInterrupt:
        pool.terminate()
        logging.warn('Caught SIGINT. Program Terminated')
        raise Exception('CRISPResso2 Terminated')
        exit (0)
    except Exception as e:
        print('CRISPResso2 failed')
        raise e
    else:
        pool.close()
    pool.join()
    return [y for x in results for y in x]



def run_subprocess(cmd):
    return sb.call(cmd, shell=True)

def run_parallel_commands(commands_arr,n_processes=1,descriptor='CRISPResso2',continue_on_fail=False):
    """
    input: commands_arr: list of shell commands to run
    descriptor: string to print out to user describing run
    """
    pool = mp.Pool(processes = n_processes)

    #handle signals -- bug in python 2.7 (https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python)
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        res = pool.map_async(run_subprocess, commands_arr)
        ret_vals = res.get(60*60*60) # Without the timeout this blocking call ignores all signals.
        for idx, ret in enumerate(ret_vals):
            if ret == 137:
                raise Exception('%s #%d was killed by your system. Please decrease the number of processes (-p) and run again.'%(descriptor, idx))
            if ret != 0 and not continue_on_fail:
                raise Exception('%s #%d failed'%(descriptor, idx))
    except KeyboardInterrupt:
        pool.terminate()
        logging.warn('Caught SIGINT. Program Terminated')
        raise Exception('CRISPResso2 Terminated')
        exit (0)
    except Exception as e:
        print('CRISPResso2 failed')
        raise e
    else:
        pool.close()
    pool.join()


def run_plot(plot_func, plot_args, num_processes, process_futures, process_pool):
    """Run a plot in parallel if num_processes > 1, otherwise in serial.

    Parameters
    ----------
    plot_func: function
        The plotting function to call.
    plot_args: dict
        The arguments to pass to the plotting function.
    num_processes: int
        The number of processes to use in parallel.
    process_futures: List
        The list of futures that submitting the parallel job will return.
    process_pool: ProcessPoolExecutor or ThreadPoolExecutor
        The pool to submit the job to.

    Returns
    -------
    None
    """
    logger = logging.getLogger(getmodule(stack()[1][0]).__name__)
    try:
        if num_processes > 1:
            process_futures[process_pool.submit(plot_func, **plot_args)] = (plot_func, plot_args)
        else:
            plot_func(**plot_args)
    except Exception as e:
        logger.warn(f"Plot error {e}, skipping plot \n")
        logger.debug(traceback.format_exc())

