'''
CRISPResso2 - Kendell Clement and Luca Pinello 2018
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2018 The General Hospital Corporation. All Rights Reserved.
'''

import gzip
import argparse
import io
import os
import datetime
import numpy
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def main():
	parser = argparse.ArgumentParser(description='Filter fastqs for quality')
	parser.add_argument('--fastq_r1',required=True)
	parser.add_argument('--fastq_r2')
	parser.add_argument('--min_bp_qual_in_read',type=int)
	parser.add_argument('--min_av_read_qual',type=int)
	parser.add_argument('--min_bp_qual_or_N',type=int)
	parser.add_argument('--fastq_r1_out')
	parser.add_argument('--fastq_r2_out')

	args = parser.parse_args()

	filterFastqs(args.fastq_r1,args.fastq_r2,args.fastq_r1_out,args.fastq_r2_out,args.min_bp_qual_in_read,args.min_av_read_qual,args.min_bp_qual_or_N)

def filterFastqs(fastq_r1=None,fastq_r2=None,fastq_r1_out=None,fastq_r2_out=None,min_bp_qual_in_read=None,min_av_read_qual=None,min_bp_qual_or_N=None,debug=False):

	if (debug):
		print('--fastq_r1:'+str(fastq_r1))
		print('--fastq_r2:'+str(fastq_r2))
		print('--min_bp_qual_in_read:'+str(min_bp_qual_in_read))
		print('--min_av_read_qual:'+str(min_av_read_qual))
		print('--min_bp_qual_or_N:'+str(min_bp_qual_or_N))
		print('--fastq_r1_out:'+str(fastq_r1_out))
		print('--fastq_r2_out:'+str(fastq_r2_out))

	startTime = datetime.datetime.now()

	if not os.path.exists(fastq_r1):
		raise Exception("fastq_r1 file '"+fastq_r1+"' does not exit.")

	if fastq_r2 is not None and not os.path.exists(fastq_r2):
		raise Exception("fastq_r2 file '"+fastq_r2+"' does not exit.")

	##CREATION OF FILEHANDLES##
	if fastq_r1.endswith('.gz'):
		f1_in = io.BufferedReader(gzip.open(fastq_r1,'rb'))
		f1_out_filename=fastq_r1.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'
		if fastq_r1_out:
			f1_out_filename = fastq_r1_out

	else:
		f1_in = open(fastq_r1,'r')
		f1_out_filename=fastq_r1.replace('.fastq','')+'_filtered.fastq'
		if fastq_r1_out:
			f1_out_filename = fastq_r1_out

	if f1_out_filename.endswith('.gz'):
		f1_out = gzip.open(f1_out_filename, 'wb')
	else:
		f1_out = open(f1_out_filename,'w')



	if fastq_r2:
		if fastq_r2.endswith('.gz'):
			f2_in = io.BufferedReader(gzip.open(fastq_r2,'rb'))
			f2_out_filename=fastq_r2.replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'
			if fastq_r2_out:
				f2_out_filename = fastq_r2_out
		else:
			f2_in = open(fastq_r2,'r')
			f2_out_filename=fastq_r2.replace('.fastq','')+'_filtered.fastq'
			if fastq_r2_out:
				f2_out_filename = fastq_r2_out

		if f2_out_filename.endswith('.gz'):
			f2_out = gzip.open(f2_out_filename, 'wb')
		else:
			f2_out = open(f2_out_filename,'w')
	###END CREATION OF FILEHANDLES##

	if not fastq_r2:
		if min_bp_qual_in_read:
			if min_av_read_qual:
				if min_bp_qual_or_N:
					run_mBP_mRQ_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					run_mBP_mRQ(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
			else:
				if min_bp_qual_or_N:
					run_mBP_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					run_mBP(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
		else:
			if min_av_read_qual:
				if min_bp_qual_or_N:
					run_mRQ_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					run_mRQ(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
			else:
				if min_bp_qual_or_N:
					run_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					exit('Finished -- No modifications requested')
	else:#paired reads
		if min_bp_qual_in_read:
			if min_av_read_qual:
				if min_bp_qual_or_N:
					run_mBP_mRQ_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					run_mBP_mRQ_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
			else:
				if min_bp_qual_or_N:
					run_mBP_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					run_mBP_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
		else:
			if min_av_read_qual:
				if min_bp_qual_or_N:
					run_mRQ_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					run_mRQ_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
			else:
				if min_bp_qual_or_N:
					run_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N)
				else:
					exit('Finished -- No modifications requested')

	endTime = datetime.datetime.now()
	timeDiff = endTime - startTime
	print("Completed in %d seconds\n"%timeDiff.total_seconds())


def run_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		npSeqLine = numpy.fromstring(seqLine,'c')
		npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
		f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))

def run_mRQ(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean = numpy.mean(npQualLine)
		if mean >= min_av_read_qual:
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,seqLine,"+",qualLine))

def run_mBP(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		min = numpy.min(npQualLine)
		if min >= min_bp_qual_in_read:
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,seqLine,"+",qualLine))

def run_mBP_mRQ(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean = numpy.mean(npQualLine)
		if mean >= min_av_read_qual:
			min = numpy.min(npQualLine)
			if min >= min_bp_qual_in_read:
				f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,seqLine,"+",qualLine))

def run_mBP_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		min = numpy.min(npQualLine)
		if min >= min_bp_qual_in_read:
			npSeqLine = numpy.fromstring(seqLine,'c')
			npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))

def run_mRQ_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean = numpy.mean(npQualLine)
		if mean >= min_av_read_qual:
			npSeqLine = numpy.fromstring(seqLine,'c')
			npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))

def run_mBP_mRQ_mBPN(f1_in,f1_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	for (idLine,seqLine,qualLine) in iter1:
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		min = numpy.min(npQualLine)
		if min >= min_bp_qual_in_read:
			mean = numpy.mean(npQualLine)
			if mean >= min_av_read_qual:
				npSeqLine = numpy.fromstring(seqLine,'c')
				npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
				f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))


#PAIRED
def run_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		npSeqLine = numpy.fromstring(seqLine,'c')
		npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
		f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))
		npSeqLine2 = numpy.fromstring(seqLine2,'c')
		npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
		f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,npSeqLine2.tostring().decode('utf-8'),"+",qualLine2))

def run_mRQ_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean = numpy.mean(npQualLine)
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean2 = numpy.mean(npQualLine2)
		if mean >= min_av_read_qual and mean2 > min_av_read_qual:
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,seqLine,"+",qualLine))
			f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,seqLine2,"+",qualLine2))

def run_mBP_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		min = numpy.min(npQualLine)
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		min2 = numpy.min(npQualLine2)
		if min >= min_bp_qual_in_read and min2 > min_bp_qual_in_read:
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,seqLine,"+",qualLine))
			f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,seqLine2,"+",qualLine2))

def run_mBP_mRQ_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean = numpy.mean(npQualLine)
		mean2 = numpy.mean(npQualLine2)
		if mean >= min_av_read_qual and mean2 >= min_av_read_qual:
			min = numpy.min(npQualLine)
			min2 = numpy.min(npQualLine2)
			if min >= min_bp_qual_in_read and min2 >= min_bp_qual_in_read:
				f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,seqLine,"+",qualLine))
				f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,seqLine2,"+",qualLine2))

def run_mBP_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		min = numpy.min(npQualLine)
		min2 = numpy.min(npQualLine2)
		if min >= min_bp_qual_in_read and min2 >= min_bp_qual_in_read:
			npSeqLine = numpy.fromstring(seqLine,'c')
			npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))
			npSeqLine2 = numpy.fromstring(seqLine2,'c')
			npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
			f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,npSeqLine2.tostring().decode('utf-8'),"+",qualLine2))

def run_mRQ_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		mean = numpy.mean(npQualLine)
		mean2 = numpy.mean(npQualLine2)
		if mean >= min_av_read_qual and mean2 >= min_av_read_qual:
			npSeqLine = numpy.fromstring(seqLine,'c')
			npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
			f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))
			npSeqLine2 = numpy.fromstring(seqLine2,'c')
			npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
			f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,npSeqLine2.tostring().decode('utf-8'),"+",qualLine2))

def run_mBP_mRQ_mBPN_pair(f1_in,f1_out,f2_in,f2_out,min_bp_qual_in_read,min_av_read_qual,min_bp_qual_or_N):
	iter1 = FastqGeneralIterator(f1_in)
	iter2 = FastqGeneralIterator(f2_in)
	for (idLine,seqLine,qualLine) in iter1:
		(idLine2,seqLine2,qualLine2) = next(iter2)
		npQualLine = numpy.fromstring(qualLine,dtype=numpy.uint8)-33 #assume illumina 1.7
		npQualLine2 = numpy.fromstring(qualLine2,dtype=numpy.uint8)-33 #assume illumina 1.7
		min = numpy.min(npQualLine)
		min2 = numpy.min(npQualLine2)
		if min >= min_bp_qual_in_read and min2 >= min_bp_qual_in_read:
			mean = numpy.mean(npQualLine)
			mean2 = numpy.mean(npQualLine2)
			if mean >= min_av_read_qual and mean2 >= min_av_read_qual:
				npSeqLine = numpy.fromstring(seqLine,'c')
				npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
				f1_out.write("@%s\n%s\n%s\n%s\n"%(idLine,npSeqLine.tostring().decode('utf-8'),"+",qualLine))
				npSeqLine2 = numpy.fromstring(seqLine2,'c')
				npSeqLine2[npQualLine2 < min_bp_qual_or_N] = 'N'
				f2_out.write("@%s\n%s\n%s\n%s\n"%(idLine2,npSeqLine2.tostring().decode('utf-8'),"+",qualLine2))



if __name__ == "__main__":
	# execute only if run as a script
	main()
