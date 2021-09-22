import cython
import numpy as np
cimport numpy as np
import re

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char* s)


cdef extern from "Python.h":
    ctypedef void PyObject
    int _PyBytes_Resize(PyObject **, size_t)
    char * PyBytes_AS_STRING(PyObject *)


re_find_indels = re.compile("(-*-)")

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def find_indels_substitutions(read_seq_al, ref_seq_al, _include_indx):

    cdef char* sub_seq=''

    cdef int st
    cdef int en

    cdef int idx_c
    cdef int idx

    #ref_positions holds the indices for which positions map back to the original reference
    # for example,
    #     1 2 3 4 5 6 7 8
    # ref A A T T G G C C
    #
    # and for a given alignment
    # ref A A T T - G G C C
    # aln A - T T T G G C C
    #     1 2 3 4-4 5 6 7 8 <ref positions. Note that the negative values/indices represent places that don't map back to the original reference
    ref_positions=[]
    all_substitution_positions=[]
    substitution_positions=[]
    all_substitution_values=[]
    substitution_values=[]

    nucSet = set(['A', 'T', 'C', 'G', 'N'])
    idx=0
    for idx_c, c in enumerate(ref_seq_al):
        if c in nucSet:
            ref_positions.append(idx)
            if ref_seq_al[idx_c]!=read_seq_al[idx_c] and read_seq_al[idx_c] != '-' and read_seq_al[idx_c] != 'N':
                all_substitution_positions.append(idx)
                all_substitution_values.append(read_seq_al[idx_c])
                if idx in _include_indx:
                    substitution_positions.append(idx)
                    substitution_values.append(read_seq_al[idx_c])

            idx+=1

        else:
            if idx==0:
                ref_positions.append(-1)
            else:
                ref_positions.append(-idx)

    substitution_n = len(substitution_positions)

    #the remainder of positions are with reference to the original reference sequence indexes we calculated above
    all_deletion_positions=[]
    deletion_positions=[]
    deletion_coordinates=[]
    deletion_sizes=[]

    all_insertion_positions=[]
    all_insertion_left_positions=[]
    insertion_positions=[]
    insertion_coordinates = []
    insertion_sizes=[]

    include_indx_set = set(_include_indx)
    for p in re_find_indels.finditer(read_seq_al):
        st,en=p.span()
        ref_st = 0
        if st-1 > 0:
          ref_st = ref_positions[st]
        ref_en = idx-1
        if en < len(ref_positions):
          ref_en = ref_positions[en]
        all_deletion_positions.extend(range(ref_st,ref_en))
        inc_del_pos = include_indx_set.intersection(range(ref_st,ref_en))
        if(len(inc_del_pos)>0):
          deletion_positions.extend(range(ref_st,ref_en))
          deletion_coordinates.append((ref_st,ref_en))
          deletion_sizes.append(en-st)

    deletion_n = np.sum(deletion_sizes)

    for p in re_find_indels.finditer(ref_seq_al):
        st,en=p.span()
        #sometimes insertions run off the end of the reference
        if st == 0: # if insertion happened before ref
          continue
        if en == len(ref_seq_al): # if insertion happened after ref
          continue
        ref_st = ref_positions[st-1]
        ref_en = ref_positions[en]

        all_insertion_left_positions.append(ref_st)
        all_insertion_positions.append(ref_st)
        all_insertion_positions.append(ref_en)
        if(ref_st in _include_indx and ref_en in _include_indx):
          insertion_coordinates.append((ref_st,ref_en))
          insertion_positions.append(ref_st)
          insertion_positions.append(ref_en)
          insertion_sizes.append(en-st)

    insertion_n = np.sum(insertion_sizes)


    retDict = {
	    'all_insertion_positions':all_insertion_positions,
	    'all_insertion_left_positions':all_insertion_left_positions,
	    'insertion_positions':insertion_positions,
	    'insertion_coordinates':insertion_coordinates,
	    'insertion_sizes':insertion_sizes,
	    'insertion_n':insertion_n,

	    'all_deletion_positions':all_deletion_positions,
	    'deletion_positions':deletion_positions,
	    'deletion_coordinates':deletion_coordinates,
	    'deletion_sizes':deletion_sizes,
	    'deletion_n':deletion_n,

	    'all_substitution_positions':all_substitution_positions,
	    'substitution_positions':substitution_positions,
	    'all_substitution_values':np.array(all_substitution_values),
	    'substitution_values':np.array(substitution_values),
	    'substitution_n':substitution_n,

	    'ref_positions':ref_positions,
    }
    return retDict

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def find_indels_substitutions_legacy(read_seq_al, ref_seq_al, _include_indx):

    cdef char* sub_seq=''

    cdef int st
    cdef int en

    cdef int idx_c
    cdef int idx

    #ref_positions holds the indices for which positions map back to the original reference
    # for example,
    #     1 2 3 4 5 6 7 8
    # ref A A T T G G C C
    #
    # and for a given alignment
    # ref A A T T - G G C C
    # aln A - T T T G G C C
    #     1 2 3 4-4 5 6 7 8 <ref positions. Note that the negative values/indices represent places that don't map back to the original reference
    ref_positions=[]
    all_substitution_positions=[]
    substitution_positions=[]
    all_substitution_values=[]
    substitution_values=[]

    nucSet = set(['A', 'T', 'C', 'G', 'N'])
    idx=0
    for idx_c, c in enumerate(ref_seq_al):
        if c in nucSet:
            ref_positions.append(idx)
            if ref_seq_al[idx_c]!=read_seq_al[idx_c] and read_seq_al[idx_c] != '-' and read_seq_al[idx_c] != 'N':
                all_substitution_positions.append(idx)
                all_substitution_values.append(read_seq_al[idx_c])
                if idx in _include_indx:
                    substitution_positions.append(idx)
                    substitution_values.append(read_seq_al[idx_c])

            idx+=1

        else:
            if idx==0:
                ref_positions.append(-1)
            else:
                ref_positions.append(-idx)

    substitution_n = len(substitution_positions)

    #the remainder of positions are with reference to the original reference sequence indexes we calculated above
    all_deletion_positions=[]
    deletion_positions=[]
    deletion_coordinates=[]
    deletion_sizes=[]

    all_insertion_positions=[]
    all_insertion_left_positions=[]
    insertion_positions=[]
    insertion_coordinates = []
    insertion_sizes=[]

    include_indx_set = set(_include_indx)
    for p in re_find_indels.finditer(read_seq_al):
        st,en=p.span()
        ref_st = 0
        if st-1 > 0:
          ref_st = ref_positions[st]
        ref_en = idx-1
        if en < len(ref_positions):
          ref_en = ref_positions[en]
        all_deletion_positions.extend(range(ref_st,ref_en))
        inc_del_pos = include_indx_set.intersection(range(ref_st,ref_en))
        if(len(inc_del_pos)>0):
          deletion_positions.extend(range(ref_st,ref_en))
          deletion_coordinates.append((ref_st,ref_en))
          deletion_sizes.append(en-st)

    deletion_n = np.sum(deletion_sizes)

    for p in re_find_indels.finditer(ref_seq_al):
        st,en=p.span()
        #sometimes insertions run off the end of the reference
        if st == 0: # if insertion happened before ref
          continue
        if en == len(ref_seq_al): # if insertion happened after ref
          continue
        ref_st = ref_positions[st-1]
        ref_en = ref_positions[en]

        all_insertion_left_positions.append(ref_st)
        all_insertion_positions.append(ref_st)
        all_insertion_positions.append(ref_en)
        if(ref_st in _include_indx or ref_en in _include_indx):
          insertion_coordinates.append((ref_st,ref_en))
          insertion_positions.append(ref_st)
          insertion_positions.append(ref_en)
          insertion_sizes.append(en-st)

    insertion_n = np.sum(insertion_sizes)


    retDict = {
	    'all_insertion_positions':all_insertion_positions,
	    'all_insertion_left_positions':all_insertion_left_positions,
	    'insertion_positions':insertion_positions,
	    'insertion_coordinates':insertion_coordinates,
	    'insertion_sizes':insertion_sizes,
	    'insertion_n':insertion_n,

	    'all_deletion_positions':all_deletion_positions,
	    'deletion_positions':deletion_positions,
	    'deletion_coordinates':deletion_coordinates,
	    'deletion_sizes':deletion_sizes,
	    'deletion_n':deletion_n,

	    'all_substitution_positions':all_substitution_positions,
	    'substitution_positions':substitution_positions,
	    'all_substitution_values':np.array(all_substitution_values),
	    'substitution_values':np.array(substitution_values),
	    'substitution_n':substitution_n,

	    'ref_positions':ref_positions,
    }
    return retDict


def calculate_homology(a, b):
    cdef char *al = a
    cdef char *bl = b
    cdef size_t l = strlen(al)
    cdef float score = 0.0

    for i in range(l):
        if al[i] == bl[i]:
            score+=1
    return score/l
