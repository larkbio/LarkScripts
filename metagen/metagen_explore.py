from IPython.display import FileLinks, FileLink, Image, HTML
from Bio import SeqIO
import pysam
from collections import Counter
import pylab
import random
import math
import numpy as np
import pandas as pd
import MySQLdb as mdb

def aln_stat(bam_fname):
	left = 0
	right = 0
	unpaired = 0
	samfile = pysam.Samfile(bam_fname, "rb")
	for r in samfile.fetch():
		readpos = r.pos+r.alen/2
		readtype = 'n'
		if r.is_proper_pair:
			if r.is_read2:
				right = right + 1
				readtype = 'r2'
			else:
				left = left + 1
				readtype = 'r1'
		else:
			unpaired = unpaired + 1
	print "all reads: %d, paired: %d, unpaired: %d"%(left+right+unpaired, left, unpaired)

def write_coverage(bam_fname, out_fname):
	samfile = pysam.Samfile(bam_fname, "rb")
	poscnt = Counter()
	with open(out_fname, 'w') as out:
		for r in samfile.fetch():
			readpos = r.pos+r.alen/2
			height = poscnt[readpos]
			readtype = 'n'
			if r.is_proper_pair:
				readtype = 'r1'
				if r.is_read2:
					readtype = 'r2'
			up = [r.pos+i for i in range(r.alen)]
			poscnt.update(up)
			out.write("%s,%d,%d,%s\n"%(r.qname, readpos, height, readtype))

def get_coverage(bam_fname):
	samfile = pysam.Samfile(bam_fname, "rb")
	poscnt = Counter()
	for r in samfile.fetch(until_eof = True):
		readpos = r.pos+r.alen/2
		height = poscnt[readpos]
		readtype = 'n'
		if r.is_proper_pair:
			readtype = 'r1'
			if r.is_read2:
				readtype = 'r2'
		up = [r.pos+i for i in range(r.alen)]
		poscnt.update(up)
	return(poscnt)

def get_coverage_fast(bam_fname, binsize):
	samfile = pysam.Samfile(bam_fname, "rb")
	poscnt = Counter()
	for r in samfile.fetch(until_eof = True):
		readpos = r.pos+r.alen/2
		poscnt.update([readpos/binsize])
	return(poscnt)
	
def plot_coverage(bam_fname, out_fname, genome_size):
	pal = ["#999999", "#e69696", "#9696e6"]
	binsize = 10000
	poscnt = get_coverage_fast(bam_fname, binsize)
	x = []
	y = []
	for i in range(genome_size/binsize):
		summed = 0
		x.append(100.0*i*binsize/genome_size)
		y.append(1.0*poscnt[i]/binsize)
	pylab.figure()
	pylab.plot(x,y)
	pylab.ylim([0, np.mean(y)*3.0])
	pylab.savefig("results/"+out_fname)

def plot_insert_size_hist(bam_fname, out_fname):
	samfile = pysam.Samfile(bam_fname, "rb")
	isizes = []
	for r in samfile.fetch(until_eof = True):
		if r.is_proper_pair and r.tlen>0: 
			isizes.append(r.tlen) 
	pylab.figure()
	pylab.hist(isizes, 100)
	pylab.savefig("results/"+out_fname)

def get_tax_names(lst):
	con = mdb.connect('localhost', 'root', 'root', 'ncbi_tax')
	result = []
	with con:
		for gi in lst:
			cur = con.cursor()
			cur.execute("select name_txt from gi_taxid_nucl as gn, ncbi_names as nn where gi = \
							%s and nn.tax_id = gn.taxid and nn.name_class='scientific name'"%(gi))
			if cur.rowcount<=0:
				result.append(None)
			else:
				for i in range(cur.rowcount):
					row = cur.fetchone()
					result.append(row[0])
	return(result)

def annot_tax_counts(cnt):
	tax_tuples = cnt.most_common()
	gis = [rec[0].split("|")[1] for rec in cnt.most_common()]
	tax_names = get_tax_names(gis)
	df = pd.DataFrame({"gi": gis, "tax_names": tax_names, "count": [rec[1] for rec in tax_tuples]}, columns=["gi", "tax_names", "count"])
	return(df)

def get_tax_by_name(name):
	con = mdb.connect('localhost', 'root', 'root', 'ncbi_tax')
	result = []
	with con:
		for gi in lst:
			cur = con.cursor()
			cur.execute("select name_txt, taxid....")

def get_tax_by_gi(gi):
	con = mdb.connect('localhost', 'root', 'root', 'ncbi_tax')
	cur = con.cursor()
	cur.execute("select gt.taxid, nn.name_txt, nn.name_class from gi_taxid_nucl as gt, ncbi_names as nn \
		where gt.gi = '%d' and nn.tax_id = gt.taxid \
		and nn.name_class = 'scientific name'"%(gi))
	result = None
	assert cur.rowcount == 0 or cur.rowcount == 1
	if cur.rowcount:
		row = cur.fetchone()
		if row[2] == 'scientific name':
			result= (row[0], row[1])
	cur.close()
	return(result)
	
def get_parent_by_taxid(taxid):
	con = mdb.connect('localhost', 'root', 'root', 'ncbi_tax')
	cur = con.cursor()
	cur.execute("select parent_tax_id from ncbi_nodes where tax_id = %d"%(taxid))
	result = None
	assert cur.rowcount == 0 or cur.rowcount == 1
	if cur.rowcount:
		row = cur.fetchone()
		result = row[0]
	return(result)
	
def get_tax(taxid):
	con = mdb.connect('localhost', 'root', 'root', 'ncbi_tax')
	cur = con.cursor()
	cur.execute("select nn.tax_id, nn.name_txt from ncbi_names as nn \
		where nn.tax_id = '%d' \
		and nn.name_class = 'scientific name'"%(taxid))
	result = None
	assert cur.rowcount == 0 or cur.rowcount == 1
	if cur.rowcount:
		row = cur.fetchone()
		result = (row[0], row[1])
	cur.close()
	return(result)
	
def get_tree(lst):
	con = mdb.connect('localhost', 'root', 'root', 'ncbi_tax')
	result = []
	with con:
		for gi in lst:
			path = []
			tax = get_tax_by_gi(gi)
			path.append(tax)
			# print " > "
			parent = get_parent_by_taxid(tax[0]) 
			while parent:
				tax = get_tax(parent)
				# print tax
				path.append(tax)
				if parent==1:
					break;
				parent = get_parent_by_taxid(tax[0]) 
			result.append(path)
	con.close()
	return(result)
