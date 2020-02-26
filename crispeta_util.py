#sys and regular functions
import argparse, tempfile, os, re, math, sys, resource, time, datetime, _mysql, itertools, datetime
from subprocess import call
from os import listdir
from config import *

#Doench Score Rule_set2.0
sys.path.insert(0,'analysis')
from rs2_score_calculator import *

#Biopithon
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from operator import itemgetter

#graphics
import numpy as np
try:
	import plotly
	_plotly=1
except ImportError:
	_plotly=0
if _plotly==1:
	import plotly.plotly as py
	import plotly.graph_objs as go
	from plotly.graph_objs import Scatter, Layout


gibson_5 = 'atcttGTGGAAAGGACGAAACACCg'
constant_h1 = 'GTTTTAGAGCTAGAAGAGACGGAATTCCTAGGATCCCGTCTCTCTGTATGAGACCACTCTTTCCC'
gibson_3 = 'gttttagagctaGAAAtagcaagttaaaataaggc'


def get_IDs(infile):
	print '-'*40
	print "Getting input ID's"
	ids = []
	fi = open(infile, 'r')
	for line in fi:
		ids.append(line.strip())
	
	return ids


def get_transcript_ids(apprs, ids):
	print "Getting principal isoforms"

	transcripts = {}	
	fi = open(apprs, 'r')

	for line in fi:

		if not 'PRINCIPAL' in line:
			continue
		
		line = line.strip().split("\t")
		gene_name = line[0]
		gene_id = line[1]

		if gene_name in ids:
			ID = gene_name
		elif gene_id in ids:
			ID = gene_id
		else:
			continue	
		
		keys = transcripts.keys()
		enst = line[2]

		if ID in keys:
			if not enst in transcripts[ID]:
				transcripts[ID].append(enst)
		else:
			transcripts[ID] = [enst]

	keys = transcripts.keys()
	for id_ in ids:
		if id_ in keys:
			continue
		else:
			transcripts[id_] = ['empty']

	return transcripts


def get_regions(gtf, transcripts):
	print "Getting genomic regions for principal isoforms"
	exons = []
	regions = {}
	regions_non_coding = {}
	non_coding_ids = []
	fi = open(gtf, 'r')
	tn = 0
	t_ids = []

	for i in transcripts:
		regions[i] = []
		regions_non_coding[i] = {}
	# 	for h in transcripts[i]:
	# 		t_ids.append(h)

	for line in fi:

		if line.startswith("#"):
			continue
		
		line = line.strip().split("\t")
		tags = get_tags(line)

		try:
			tags['gene_name']
		except KeyError:
			tags['gene_name'] = ''

		if tags['gene_name'] in transcripts:
			gene_id_in_line = tags['gene_name']
		elif tags['gene_id'] in transcripts:
			gene_id_in_line = tags['gene_id']
		else:
			continue

		if not 'gene_biotype' in tags:
			tags['gene_biotype'] = tags['gene_type']

		if tags['gene_biotype'] == 'protein_coding':

			fetch_line_protein_coding(regions, transcripts, line, tags, gene_id_in_line)

#		else:
#			non_coding_ids.append(gene_id_in_line)
#			fetch_line_non_coding(regions_non_coding, transcripts, line, tags, gene_id_in_line)

#	append_non_coding_to_regions(regions, regions_non_coding, non_coding_ids)
	return regions


def fetch_line_protein_coding(regions, transcripts, line, tags, gene_id_in_line):

	if not line[2] == "CDS":
		return
		
	enst_id_in_line = tags['transcript_id']
	enst_check = 0
	for enst in transcripts[gene_id_in_line]:
		if enst == ['empty'] or enst in enst_id_in_line:
			enst_check = 1			
	if enst_check == 0:
		return

	en = tags['exon_number']
	inf_ = [tags['chrom'], int(tags['start'])-6, int(tags['end'])+6, gene_id_in_line, en, tags['strand']]

	v = 0
	for r in regions[gene_id_in_line]:
		if inf_[:3] == r[:3]:
			v = 1

	if v == 0: 	regions[gene_id_in_line].append(inf_)


def fetch_line_non_coding(regions_non_coding, transcripts, line, tags, gene_id_in_line):

	if not line[2] == "exon":
		return	

	enst = tags['transcript_id']
	en = tags['exon_number']
	inf_ = [tags['chrom'], int(tags['start'])-6, int(tags['end'])+6, gene_id_in_line, en, tags['strand']]	

	if enst in regions_non_coding[gene_id_in_line]:
		regions_non_coding[gene_id_in_line][enst].append(inf_)
	else:
		regions_non_coding[gene_id_in_line][enst] = [inf_]


def append_non_coding_to_regions(regions, regions_non_coding, non_coding_ids):

	for gene_id in non_coding_ids:
		longest = ['id', 0]
		for transcript in regions_non_coding[gene_id]:
			length = 0
			for exon in regions_non_coding[gene_id][transcript]:
				length += exon[2] - exon[1]
			if length > longest[1]:
				longest = [transcript, length]
		regions[gene_id] = regions_non_coding[gene_id][longest[0]]



def get_tags(line):

	n, tags = 0, {}
	tags['chrom'] = line[0]
	tags['type'] = line[2]
	tags['start'] = line[3]
	tags['end'] = line[4]
	tags['strand'] = line[6]

	info = line[-1].split(" ")
	for i in info:
		if n == 0:
			n, tag = 1, i
		else:
			if "." in i:
				n, tags[tag] = 0, i.replace(";","").replace("\"","")
			else:
				n, tags[tag] = 0, i.replace(";","").replace("\"","").split(".")[0]
	return tags


def get_sequences(dic, genome):
	print 'Reading Fasta'
	chrs = {}
	for seq in SeqIO.parse(genome, "fasta"):
		chrs[seq.id] = seq.seq
	print 'Getting sequences'
	new_dic = {}
	for ids in dic.keys():
		new_dic[ids] = []
		length = []
		for regions in dic[ids]:
			if not 'chr' in regions[0]:
				regions[0] = 'chr'+regions[0]
			regions.append(int(regions[2])-int(regions[1]))
			regions.append(chrs[regions[0]][regions[1]:regions[2]].upper())
			new_dic[ids].append(regions)
	return new_dic


def get_sgRNAs(regions, si, sp, t, model, s_method):
	print "Searching for sgRNAs"
	sgrnas = {}
	for ids in regions.keys():
		sgrnas[ids] = []
		length = 0
		total_length = 0

		for region in regions[ids]:
			total_length += region[-2]

		for region in regions[ids]:

			seq = region[-1]
			seq_length = region[-2]
			seq_strand = region[-3]

			for m in re.finditer(r'(?=(.{24}.GG...))',str(seq)):
				
				sgrna_info = get_sgrna_info(region, m, '+', length, seq_length, seq_strand)
				
				if filter_sgrna(m.group(1), t, si, sgrna_info, model, s_method, total_length):
					continue
				sgrnas[ids].append(sgrna_info)
				
			for m in re.finditer(r'(?=(...CC.{25}))',str(seq)):

				sgrna_info = get_sgrna_info(region, m, '-', length, seq_length, seq_strand)
				
				if filter_sgrna(m.group(1), t, si, sgrna_info, model, s_method, total_length):
					continue
				sgrnas[ids].append(sgrna_info)

			length += region[-2]

	return sgrnas


def filter_sgrna(sgrna, t, si, sgrna_info, model, s_method, total_length):

	if 'TTTTT' in sgrna[4:-6]:
		return 1

	if 'N' in sgrna:
		return 1

	aa = sgrna_info[4]
	per = (aa*100)/(total_length/3)
	if s_method == 'Doench':
		if calc_score_rule_set_2(sgrna, aa, per, si, sgrna_info, model) < float(si):
			return 1
	else:
		if calc_score_han_xu(sgrna[4:-6], matrix, si, sgrna_info) < float(si):
			return 1	

	if check_grna_off_targets(sgrna[4:-6], t) == 0:
		return 1


def check_grna_off_targets(sgRNA, t):
	'''
	Checks the number of offtargets for a sgRNA
	'''
	from config import config
	t=t.split(',')
	for i,j in enumerate(t):
		if j=='x':
			t[i]=float("inf")
		else:
			t[i]=int(j)
	table = config['table']
	sgrna = sgRNA
	conn = use_database()

	conn.query("SELECT * FROM "+table+" where grna = '"+sgrna+"';")
	result = conn.store_result()
	off = result.fetch_row()
	if len(off) == 0:
		return 0
	if int(off[0][1]) <= t[0]+1 and int(off[0][2]) <= t[1] and int(off[0][3]) <= t[2] and int(off[0][4]) <= t[3] and int(off[0][5]) <= t[4]:
		return 1
	else:
		return 0


def use_database(host='localhost', dbuser='crispeta', database='crispeta', dbpass='pwd'):

	from config import config
	dic = config['mysql_db']
	conn = _mysql.connect(user=dic['user'],passwd=dic['passwd'],host=dic['host'],db=dic['db'])
	return conn


def get_sgrna_info(region, m, strand, length, seq_length, seq_strand):

	chrom = region[0]
	start = get_start_position(m.start(), region[1], strand)
	end = start + 23
	cut_position = get_cut_position(m.start(), strand, length, seq_length, seq_strand)

	sgrna_info = [chrom, start, end, m.group(1), cut_position, strand]
	if strand == '-':
		sgrna_info[3] = reverse_complement(m.group(1))
	return sgrna_info

def get_start_position(grna_start, region_start, strand):

	if strand == '+':
		start = region_start + grna_start + 4
	else:
		start = region_start + grna_start + 3
	return start


def get_cut_position(start, strand, length, seq_length, seq_strand):

	if seq_strand == '+':
		if strand == '+':
			cut_position = int(length+start+21)
		else:
			cut_position = int(length+start+9)

	else:
		if strand == '+':
			cut_position = int(length+seq_length-start-21)
		else:
			cut_position = int(length+seq_length-start-9)
	
	return cut_position


def calc_score_rule_set_1(s):
	'''
	Doench old score rule set 1
	'''
	s_list = list(s)
	if len(s_list) != 30:
		return 0
	s_20mer = s[4:24]
	nuc_hash = {'A':0, 'T':1, 'C':2, 'G':3}
	score = 0.597636154
	gc = s_20mer.count('G')+s_20mer.count('C')
	gc_low = -0.202625894
	gc_high = -0.166587752
	if gc < 10:
		gc_val = abs(gc-10)
		score = score+(gc_val*gc_low)
	elif gc > 10:
		gc_val = gc-10
		score = score+(gc_val*gc_high)
	#rows[1-30]cols['ATCG']
	sing_nuc_hash = {'G2':-0.275377128,'A3':-0.323887456,'C3':0.172128871,'C4':-0.100666209,'C5':-0.20180294, \
					'G5':0.245956633,'A6':0.036440041,'C6':0.098376835,'C7':-0.741181291,\
					'G7':-0.393264397,'A12':-0.466099015,'A15':0.085376945,'C15':-0.013813972,\
					'A16':0.272620512,'C16':-0.119022648,'T16':-0.285944222,'A17':0.097454592,\
					'G17':-0.17554617,'C18':-0.345795451,'G18':-0.678096426,'A19':0.22508903,\
					'C19':-0.507794051,'G20':-0.417373597,'T20':-0.054306959,'G21':0.379899366,\
					'T21':-0.090712644,'C22':0.057823319,'T22':-0.530567296,'T23':-0.877007428,\
					'C24':-0.876235846,'G24':0.278916259,'T24':-0.403102218,'A25':-0.077300704,\
					'C25':0.287935617,'T25':-0.221637217,'G28':-0.689016682,'T28':0.117877577,\
					'C29':-0.160445304,'G30':0.386342585}
	dinuc_hash = {'GT2':-0.625778696,'GC5':0.300043317,'AA6':-0.834836245,'TA6':0.760627772,'GG7':-0.490816749,'GG12':-1.516907439,'TA12':0.7092612,'TC12':0.496298609,'TT12':-0.586873894,'GG13':-0.334563735,'GA14':0.76384993,'GC14':-0.53702517,'TG17':-0.798146133,'GG19':-0.66680873,'TC19':0.353183252,'CC20':0.748072092,'TG20':-0.367266772,'AC21':0.568209132,'CG21':0.329072074,'GA21':-0.836456755,'GG21':-0.782207584,'TC22':-1.029692957,'CG23':0.856197823,'CT23':-0.463207679,'AA24':-0.579492389,'AG24':0.649075537,'AG25':-0.077300704,'CG25':0.287935617,'TG25':-0.221637217,'GT27':0.117877577,'GG29':-0.697740024}
	for i,nuc in enumerate(s_list):
		key = nuc+str(i+1)
		if sing_nuc_hash.has_key(key):
			nuc_score = sing_nuc_hash[key]
		else:
			nuc_score = 0
		score = score+nuc_score
		if i<29:
			dinuc = nuc+s[i+1]+str(i+1)
			if dinuc in dinuc_hash.keys():
				score = score+dinuc_hash[dinuc]
	partial_score = math.e**-score
	final_score = 1/(1+partial_score)
	return final_score

def create_model():
	model_file_2 = 'saved_models/V3_model_full.pickle'
	model_file = model_file_2
	try:
		with open(model_file, 'rb') as f:
			model= pickle.load(f)    
	except:
		raise Exception("could not find model stored to file %s" % model_file)

	return model


def calc_score_rule_set_2(seq, aa_cut, per_peptide, si, sgrna_info, model):

	if sklearn.__version__ != '0.16.1':
		print 'Incorrect scikit-learn version. Use 0.16.1'
		sys.exit(1)
	seq = seq
	if len(seq)!=30: 
		print "Please enter a 30mer sequence."
		sys.exit(1)
	aa_cut = int(aa_cut/3)
	per_peptide = per_peptide

	if seq[25:27] == 'GG':
		score = model_comparison.predict(seq, aa_cut, per_peptide, model=model)
		if score > float(si):
			sgrna_info.append(str(score))
	else:
		score = 0

	return score


def get_pairs(sgrnas, c_dist, rank, g_dist):

	pairs = {}

	for ids in sgrnas.keys():

		if len(sgrnas[ids]) <= 1:
			continue
		pairs[ids] = []

		for i in range(0,len(sgrnas[ids])):
			for j in range(0,len(sgrnas[ids])):

				if j <= i:
					continue

				a = sgrnas[ids][i]
				b = sgrnas[ids][j]
				c = [str(float(a[6])+float(b[6]))]
				pair = a+b+c
				cut = max(pair[4],pair[11])-min(pair[4],pair[11])
				genome_dist = max(pair[1],pair[8])-min(pair[2],pair[9])

				if int(g_dist) == 0:
					pass
				elif genome_dist > int(g_dist) or genome_dist < int(c_dist):
					continue
					
				if cut % 3 == 0:	#remove pairs where cutting region is multiple of 3nt's or are closer than dist
					continue

				pair.append(int(pair[4])+int(pair[11]))						#order of appearance of the pair
				pair.append(max(pair[4],pair[11])-min(pair[4],pair[11]))	#distance between elements of the pair
				pairs[ids].append(pair)

		for el in pairs[ids]:
			el[14]=float(el[14])

		if rank == "score":
			pairs[ids].sort(key=itemgetter(14),reverse=True)	#paired score
		elif rank == "dist":
			pairs[ids].sort(key=itemgetter(16),reverse=True)	#distance between sgRNAs
		else:
			pairs[ids].sort(key=itemgetter(15))					#appearance order
	return pairs


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    bases = list(seq) 
    letters = [complement[base] for base in bases] 
    return ''.join(letters)


def reverse_complement(seq):
    return complement(seq[::-1])


def get_output(pairs, n , v, method):

	output_lines = []
	sgrna1_dict = dict()
	sgrna2_dict = dict()

	for i in pairs.keys():
		z=0
		
		for el in pairs[i]:

			if el[3] in sgrna1_dict:
				sgrna1_dict[el[3]] += 1.0
			else:
				sgrna1_dict[el[3]] = 1.0

			if el[10] in sgrna2_dict:
				sgrna2_dict[el[10]] += 1.0
			else:
				sgrna2_dict[el[10]] = 1.0

			if (sgrna1_dict[el[3]] / n) > v or (sgrna2_dict[el[10]] / n) > v:
				continue

			z+=1
			if z > n:
				break

			sgRNA1 = el[3][4:-6]
			sgRNA2 = el[10][4:-6]

			if method == "DECKO":
				if sgRNA1[0] == "G":
					oligo = gibson_5+sgRNA1+constant_h1+sgRNA2+gibson_3
				elif sgRNA2[0] == "G":
					oligo = gibson_5+sgRNA2+constant_h1+sgRNA1+gibson_3
				else:
					oligo = gibson_5+"G"+sgRNA1+constant_h1+sgRNA2+gibson_3
			else:
				oligo = '.'
			
			pair1 = [el[0], str(el[1]), str(el[2]), sgRNA1, str(el[6]), str(el[5]), str(el[4]/3)]
			pair2 = [el[7], str(el[8]), str(el[9]), sgRNA2, str(el[13]), str(el[12]), str(el[11]/3)]
			pscore = str(el[14])
			id_ = i+"("+str(z)+")"
			line = '\t'.join([id_]+pair1+pair2+[pscore]+[oligo]+["\n"])
			output_lines.append(line)

	return output_lines


def get_browser(pairs, n, v):

	browser_lines = []
	sgrna1_dict = dict()
	sgrna2_dict = dict()
	
	for i in pairs.keys():
		z=0
#		max_, min_ = 0, 2
#		for el in pairs[i]:
#			max_ = el[14] if el[14] > max_ else max_
#			min_ = el[14] if el[14] < min_ else min_ 
		max_, min_ = 1.3, 0.7

		for el in pairs[i]:

			OldValue = el[14]
			OldRange = (max_ - min_)  
			NewRange = (255.0 - 0.0)  
			if OldValue > max_:
				NewValue = 255
			elif OldValue < min_:
				NewValue = 0
			else:
				NewValue = (((OldValue - min_) * NewRange) / OldRange)
			el.append(NewValue)	
			
		for el in pairs[i]:

			if el[3] in sgrna1_dict:
				sgrna1_dict[el[3]] += 1.0
			else:
				sgrna1_dict[el[3]] = 1.0

			if el[10] in sgrna2_dict:
				sgrna2_dict[el[10]] += 1.0
			else:
				sgrna2_dict[el[10]] = 1.0

			if (sgrna1_dict[el[3]] / n) > v or (sgrna2_dict[el[10]] / n) > v:
				continue
			
			z+=1
			if z > n:
				break
			pair = el
			start = str(min(pair[1],pair[8]))
			end = str(max(pair[2],pair[9]))
			browser_lines.append(pair[0]+'\t'+start+'\t'+end+'\t'+i+"("+str(z)+")"+'\t'+str(pair[14])+'\t'+'.'+
					#'\t'+start+'\t'+end+'\t'+str(round((float(pair[14])*255.0)/2.0,0))+',0,0'+'\t'+'2'+'\t'+
					'\t'+start+'\t'+end+'\t'+str(pair[-1])+',0,0'+'\t'+'2'+'\t'+
					'23,23'+'\t'+'0,'+str(int(end)-int(start)-23)+"\n")

	return browser_lines


def print_pairs(fo, pairs, n, variety, method):

	output = get_output(pairs, int(n), float(variety), method)
	for line in output:
		fo.write(line) 
	return output

def print_browser(fo, pairs, n, variety):

	browser = get_browser(pairs, int(n), float(variety))
	fo.write('track name="pairs" description=" " visibility=3 itemRgb="On"\n')
	for line in browser:
		fo.write(line)

matrix = [
[0,0,0,0],
[0.006447933,-0.005126545,0.01868168,-0.0263186],
[-0.002556713,0.01784427,-0.001228439,-0.02082668],
[0.02067124,0.001810971,0.009265117,-0.03726897],
[-0.01541201,-0.006865736,0.02225981,-0.007179623],
[0.006701975,-0.009327798,0.01266799,-0.01061834],
[0.003990316,0.008476303,0.01163767,-0.03195586],
[-0.005194459,-0.002546783,0.02465747,-0.02884586],
[-0.02852648,0.009015268,0.03048746,-0.02578712],
[-0.009598615,0.002571972,0.01891338,-0.02318647],
[0.01047803,-0.009462136,0.03431422,-0.04550686],
[0.02912254,-0.01410245,0.002359817,-0.01265939],
[0.004315982,-0.006850088,0.02789401,-0.03952176],
[0.01410663,-0.003234736,-0.01502249,0.01281278],
[0.0376519,-0.003153252,-0.02123348,-0.00178908],
[0.01593133,0.025584,-0.03391245,0.000224124],
[0.02561585,-0.02293057,0.01418497,-0.01378385],
[-0.01040574,0.005392964,0.009837631,-0.01119136],
[0.04550361,-0.01997099,0.003354427,-0.0156738],
[0.0248649,-0.04276305,0.03505556,-0.01229142]]

def SeqToIndex(seq):
	index = {'A':0, 'C':1, 'G':2, 'T':3}
	seq_index = []
	for nt in seq:
		seq_index.append(index[nt])
	return seq_index

def ComputeSeqScore(matrix, index):

	score = 0
	for i,j in enumerate(index):
		score += matrix[i][j]
	return score

def calc_score_han_xu(sgRNA, matrix, si, sgrna_info):

	seq = SeqToIndex(sgRNA)
	score = ComputeSeqScore(matrix, seq)
	if score > float(si):
			sgrna_info.append(str(score))
	return score


def print_log(ids, pairs, start_time, options):

	end_time = datetime.datetime.now()
	diff = end_time - start_time 
	minutes, seconds = divmod(diff.seconds, 60)
	complete, incomplete = 0, 0
	dic_pairs = {}

	for pair in pairs:
		id_ = pair.split("(")[0]
		if id_ in dic_pairs.keys():
			dic_pairs[id_] += 1
		else:
			dic_pairs[id_] = 1

	for key in dic_pairs.keys():
		if dic_pairs[key] >= 10:
			complete += 1
		else:
			incomplete += 1

	print '-'*40
	print "Analysis took %d minutes and %d seconds.\n" % (minutes, seconds)
	print "%d Analyzed targets" % (len(ids))
	print "\t%d targets with %d pairs" % (complete, int(options.n))
	print "\t%d targets with less than %d pairs" % (incomplete, int(options.n))
	print "\t%d targets with 0 results" % (len(ids)-complete-incomplete)
	print '-'*40
	print "Designed pairs in: %s.txt" % (options.outfile)
	print "BED file for pairs in: %s.bed" % (options.outfile)

