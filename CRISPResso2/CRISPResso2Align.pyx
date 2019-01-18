#with help from friends at https://github.com/brentp/align for cython implementation (no thanks for bug-ridden algorithm)
# https://github.com/dnase/affine-gap-sequence-alignment/blob/master/alignment.py for affine gap algorithm

from cython.view cimport array as cvarray
import numpy as np
cimport numpy as np



cimport cython
import sys
import os.path

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject *PyString_FromStringAndSize(char *, size_t)
    int _PyString_Resize(PyObject **, size_t)
    char * PyString_AS_STRING(PyObject *)

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.int8_t DTYPE_BOOL

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4
cdef size_t MARRAY = 1, IARRAY = 2, JARRAY = 3

def read_matrix(path):
    """
    Read a matrix in the NCBI format
    The score for a 'C' changing to an 'A' is stored in the matrix as:
        mat[ord('C'), ord('A')] = score
    """
    cdef np.ndarray[DTYPE_INT, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v, mat_size

    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [ord(x) for x in line.split(' ') if x]
    mat_size = max(headers) + 1

    a = np.zeros((mat_size, mat_size), dtype=np.int)

    line = fh.readline()
    while line:
        line_vals = [int(x) for x in line[:-1].split(' ')[1:] if x]
        for ohidx, val in zip(headers, line_vals):
            a[headers[ai], ohidx] = val
        ai += 1
        line = fh.readline()

    return a

@cython.boundscheck(False)
@cython.nonecheck(False)
def global_align(char* seqj, char* seqi, np.ndarray[DTYPE_INT, ndim=2] matrix,
          np.ndarray[DTYPE_INT,ndim=1] gap_incentive, int gap_open=-1,
          int gap_extend=-1):
    """
    Global sequence alignment (needleman-wunsch) on seq i and j.
    Reference is seq_i, sequenced read is seq_j
    Match and mismatch values are read from matrix object
    where matrix is of the format provided in the ncbi/data directory.
    gap_incentive is the incentive for having a gap at each position in seqi -
      this allows for the preferential location of gaps to be at the predicted
      cut site in genome editing experiments.

    """

    cdef size_t max_j = strlen(seqj)
    cdef size_t max_i = strlen(seqi)
    if len(gap_incentive) != max_i + 1:
        print('\nERROR: Mismatch in gap_incentive length (gap_incentive: ' + str(len(gap_incentive)) + ' ref: '+str(max_i+1) + '\n')
        return 0

    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, seqlen, align_counter = 0, p
    cdef int diag_score, up_score, left_score, tscore

    cdef char *align_j
    cdef char *align_i
    cdef char ci
    cdef char cj
    cdef PyObject *ai
    cdef PyObject *aj

    #create 3 arrays of scores and 3 arrays of pointers
    # M array - best alignment so far ending with a match
    # I array - best alignment so far ending with a gap in Read (J) (insertion in ref, deletion in read)
    # J array - best alignment so far ending with a gap in Ref (I) (deletion in ref, insertion in read)

    cdef int [:,:] mScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] iScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] jScore = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] mPointer = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] iPointer = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))
    cdef int [:,:] jPointer = np.empty((max_i + 1, max_j + 1), dtype=np.dtype("i"))


    cdef int min_score = gap_open * max_j * max_i

    #init match matrix
    mScore[0,1:] = min_score
    mScore[1:,0] = min_score
    mScore[0,0] = 0
    mPointer[0,1:] = IARRAY
    mPointer[1:,0] = JARRAY
    mPointer[0,0] = 0

# no penalty for gaps starting at beginning
#    score[0, 1:] = 0
#    score[1:, 0] = 0
# gap extension penalty for gaps starting at beginning
    #init i matrix
    for i in range(1,max_j+1):
      iScore[0,i] = gap_extend * i
#    iScore[0,1:] = [gap_extend * np.arange(1, max_j+1, dtype=np.int)]
    iScore[0:,0] = min_score
    iPointer[0,1:] = IARRAY

    #init j matrix
    for i in range(1,max_i+1):
      jScore[i,0] = gap_extend * i
    #jScore[1:,0] = np.vectorize(gap_extend * np.arange(1, max_i+1, dtype=np.int))
    jScore[0,0:] = min_score
    jPointer[1:,0] = JARRAY

#    print('gap penalty is'+str(gap_incentive))

    cdef int iFromMVal
    cdef int iExtendVal
    cdef int jFromMVal
    cdef int jExtendVal
    cdef int mVal, iVal, jVal

    #apply NW algorithm for inside squares (not last row or column)
    for i in range(1, max_i):
        ci = seqi[i - 1] #char in i

        for j in range(1, max_j):
            cj = seqj[j - 1] #char in j

            iFromMVal = gap_open + mScore[i, j - 1] + gap_incentive[i]
            iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
            if iFromMVal > iExtendVal:
                    iScore[i,j] = iFromMVal
                    iPointer[i,j] = MARRAY
            else:
                    iScore[i,j] = iExtendVal
                    iPointer[i,j] = IARRAY

            jFromMVal = gap_open + mScore[i - 1, j] + gap_incentive[i-1]
	    #no gap incentive here -- J already got the gap incentive when it transitioned from M, so don't add it again if we're extending.
            jExtendVal = gap_extend + jScore[i - 1, j]
            if jFromMVal > jExtendVal:
                    jScore[i,j] =  jFromMVal
                    jPointer[i,j] = MARRAY
            else:
                    jScore[i,j] = jExtendVal
                    jPointer[i,j] = JARRAY

            mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
            iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
            jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
            if mVal > jVal:
                if mVal > iVal:
                    mScore[i, j] = mVal
                    mPointer[i, j] = MARRAY
                else:
                    mScore[i, j]   = iVal
                    mPointer[i, j] = IARRAY
            else:
                if jVal > iVal:
                    mScore[i, j]  = jVal
                    mPointer[i, j] = JARRAY
                else:
                    mScore[i, j] = iVal
                    mPointer[i, j] = IARRAY

#            print('mScore['+str(i) + ',' + str(j) +']: ' + str(mScore[i,j]) + ': max(' + str(mScore[i - 1, j - 1])+ '+ (' + str(ci)+ ',' + str(cj) + ') ' + str(matrix[ci,cj]) + ', i:'+str(iVal) + ',j:' + str(jVal))

    #for last column and last row, ignore gap opening penalty
    #last column
    j = max_j
    cj = seqj[j-1]
    for i in range(1, max_i):
        ci = seqi[i-1]

        iFromMVal = gap_extend + mScore[i, j - 1] + gap_incentive[i]
        iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
        if iFromMVal > iExtendVal:
                iScore[i,j] =  iFromMVal
                iPointer[i,j] = MARRAY
        else:
                iScore[i,j] = iExtendVal
                iPointer[i,j] = IARRAY

        jFromMVal = gap_extend + mScore[i - 1, j] + gap_incentive[i-1]
        jExtendVal = gap_extend + jScore[i - 1, j]
        if jFromMVal > jExtendVal:
                jScore[i,j] =  jFromMVal
                jPointer[i,j] = MARRAY
        else:
                jScore[i,j] = jExtendVal
                jPointer[i,j] = JARRAY

        mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
        iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
        jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
        if mVal > jVal:
            if mVal > iVal:
                mScore[i, j] = mVal
                mPointer[i, j] = MARRAY
            else:
                mScore[i, j]   = iVal
                mPointer[i, j] = IARRAY
        else:
            if jVal > iVal:
                mScore[i, j]  = jVal
                mPointer[i, j] = JARRAY
            else:
                mScore[i, j] = iVal
                mPointer[i, j] = IARRAY
#        print('lastCol: mScore['+str(i) + ',' + str(j) +']: ' + str(mScore[i,j]) + ': max(' + str(mScore[i - 1, j - 1])+ '+ (' + str(ci)+ ',' + str(cj) + ') ' + str(matrix[ci,cj]) + ', i:'+str(iVal) + ',j:' + str(jVal))

    #last row
    i = max_i
    ci = seqi[i - 1]
    for j in range(1, max_j+1):
        cj = seqj[j - 1]

        iFromMVal = gap_extend + mScore[i, j - 1] + gap_incentive[i]
        iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
        if iFromMVal > iExtendVal:
                iScore[i,j] =  iFromMVal
                iPointer[i,j] = MARRAY
        else:
                iScore[i,j] = iExtendVal
                iPointer[i,j] = IARRAY

        jFromMVal = gap_extend + mScore[i - 1, j] + gap_incentive[i-1]
        jExtendVal = gap_extend + jScore[i - 1, j]
        if jFromMVal > jExtendVal:
                jScore[i,j] =  jFromMVal
                jPointer[i,j] = MARRAY
        else:
                jScore[i,j] = jExtendVal
                jPointer[i,j] = JARRAY


        mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
        iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
        jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
        if mVal > jVal:
            if mVal > iVal:
                mScore[i, j] = mVal
                mPointer[i, j] = MARRAY
            else:
                mScore[i, j]   = iVal
                mPointer[i, j] = IARRAY
        else:
            if jVal > iVal:
                mScore[i, j]  = jVal
                mPointer[i, j] = JARRAY
            else:
                mScore[i, j] = iVal
                mPointer[i, j] = IARRAY
#        print('lastRow: mScore['+str(i) + ',' + str(j) +']: ' + str(mScore[i,j]) + ': max(' + str(mScore[i - 1, j - 1])+ '+ (' + str(ci)+ ',' + str(cj) + ') ' + str(matrix[ci,cj]) + ', i:'+str(iVal) + ',j:' + str(jVal))


    seqlen = max_i + max_j
    ai = PyString_FromStringAndSize(NULL, seqlen)
    aj = PyString_FromStringAndSize(NULL, seqlen)

    # had to use this and PyObject instead of assigning directly...
    align_j = PyString_AS_STRING(aj)
    align_i = PyString_AS_STRING(ai)
#
#    print('mScore')
#    for ii in range(mScore.shape[0]):
#      for jj in range(mScore.shape[1]):
#        print(str(mScore[ii,jj]) + '.' + str(mPointer[ii,jj])+ ","),
#      print("\n"),
#    print('iScore')
#    for ii in range(iScore.shape[0]):
#      for jj in range(iScore.shape[1]):
#        print(str(iScore[ii,jj]) + '.' + str(iPointer[ii,jj])+ ","),
#      print("\n"),
#    print('jScore')
#    for ii in range(jScore.shape[0]):
#      for jj in range(jScore.shape[1]):
#        print(str(jScore[ii,jj]) + '.' + str(jPointer[ii,jj])+ ","),
#      print("\n"),

    cdef int matchCount = 0


    i = max_i
    j = max_j
    ci = seqi[i - 1]
    cj = seqj[j - 1]
    cdef int currMatrix
    currMatrix = MARRAY
    if mScore[i,j] > jScore[i,j]:
        if mScore[i,j] > iScore[i,j]:
            currMatrix = MARRAY
        else:
            currMatrix = IARRAY
    else:
        if jScore[i,j] > iScore[i,j]:
            currMatrix = JARRAY
        else:
            currMatrix = IARRAY
#    print('seqi' + str(seqi))
#    print('seqj' + str(seqj))
    while i > 0 or j > 0:
#        print("i: " + str(i) + " j: " + str(j) + " currMatrix: " + str(currMatrix) + " match score: " + str(mScore[i,j]) + " last match: " +  str(mScore[i-1,j-1]) + " matrix[" + str(ci) + "," + str(cj) + "]: " + str(matrix[ci,cj]) + " last j " + str(jScore[i,j]) + " last i: " + str(iScore[i,j]) + " mpointer: " + str(mPointer[i,j]) + " ipointer: " + str(iPointer[i,j]) + " jpointer: " + str(jPointer[i,j]))

        currVal = mScore[i,j]
        currPtr = mPointer[i,j]
        if currMatrix == 2:
          currVal = iScore[i,j]
          currPtr = iPointer[i,j]
        if currMatrix == 3:
          currVal = jScore[i,j]
          currPtr = jPointer[i,j]
#        print("i: " + str(i) + " j: " + str(j) + " " + str(currMatrix) +':' + str(currVal) + ' > ' + str(currPtr))
        if i > 0 or j > 0:
          if currMatrix == MARRAY:
            currMatrix = mPointer[i,j]
            i -= 1
            j -= 1
            align_j[align_counter] = cj
            align_i[align_counter] = ci
            if cj == ci:
                 matchCount += 1

            ci = seqi[i-1]
            cj = seqj[j-1]

#            print('in M set to ' + str(currMatrix))
          elif currMatrix == JARRAY:
            currMatrix = jPointer[i,j]
            i -= 1
            align_j[align_counter] = c"-"
            align_i[align_counter] = ci
            ci = seqi[i-1]
          elif currMatrix == IARRAY:
            currMatrix = iPointer[i,j]
            j -= 1
            align_j[align_counter] = cj
            align_i[align_counter] = c"-"
            cj = seqj[j-1]
          else:
            print('i: ' + str(i) + ' j: ' + str(j))
            print('currMatrix:' + str(currMatrix))
            print('seqj: ' + str(seqj) + ' seqi: ' + str(seqi))
            raise Exception('wtf4!:pointer: %i', i)
#          print('at end, currMatrix is ' + str(currMatrix))

        align_counter += 1

    _PyString_Resize(&aj, align_counter)
    _PyString_Resize(&ai, align_counter)

#    print(str(matchCount) + " aln: " + str(align_counter))
    final_score = 100*matchCount/float(align_counter)
    return (<object>aj)[::-1], (<object>ai)[::-1],final_score
