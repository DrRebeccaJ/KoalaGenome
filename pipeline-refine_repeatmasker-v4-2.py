## Pipeline to take RepeatMasker calls and refine to split interspersed repeats and reconstruct the disrupted insertions

################################################################################
#
#	programme_mit_license
#
#	Copyright (c) 09/12/2016 Amanda Chong
#
#	Permission is hereby granted, free of charge, to any person obtaining a copy
#	of this software and associated documentation files (the "Software"), to deal
#	in the Software without restriction, including without limitation the rights
#	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#	copies of the Software, and to permit persons to whom the Software is
#	furnished to do so, subject to the following conditions:
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#	THE SOFTWARE.
#
#################################################################################


from __future__ import division

import argparse
import cPickle
import csv
import datetime
import fnmatch
import operator
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

from datetime import timedelta
from decimal import *
from operator import itemgetter
from random import randint

## Imports that require aditional modules

from ruffus import * ## ruffus
from bx.intervals.intersection import IntervalTree, Interval ## bx-python

###############################################################################

start = time.time()
log_file = str(start) + ".log"

###############################################################################
# OPTIONS
###############################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", dest="infile", required=True,  
					help="Full path to the RepeatMasker .out file")
parser.add_argument("-d", "--work_directory", dest="work_dir", 
					help ="Name of project directory")
parser.add_argument("-m", "--repeat_model_file", dest="repeat_family_models", required=True, 
					help ="Name of library of repeats used to create the RepeatMasker output. NOTE: assumes single-line fasta file")
parser.add_argument("-j", "--RetroTector_directory", dest="ReTe_dir", required=True,
					help ="RetroTector")
parser.add_argument("-g", "--reference_genome", dest="genome_file", required=True,
					help ="Reference genome")
parser.add_argument("-a", "--annotation_file", dest="gff_file", required=True,
					help ="Gene annotation file")


## Pipeline run mode
					
parser.add_argument("-p", "--print", dest="no_run", action="store_true", 
					default=False, 
					help="Prints tasks to be run to stdout")
parser.add_argument("--graph", dest="graph", action="store_true", 
					default=False, 
					help="Prints graphical output of program to output folder")
parser.add_argument("-r", "--run", dest="run", action="store_true", 
					default=False, 
					help="Runs program")
parser.add_argument("-n", "--nproc", dest="num_proc", default=1, 
					 help="Sets number of cores the program will use. Default is 1")

args = parser.parse_args()

if args.infile:
	infile = args.infile
else:
	print "ERROR: No input file provided! Need RepeatMasker results file"

if args.work_dir:
	work_dir = args.work_dir
	if not os.path.exists(work_dir):
		print "Creating project directory..."
		os.makedirs(work_dir)
else:
	print "Creating project directory..."
	start_time = datetime.datetime.now()
	start_date = start_time.strftime("%y-%m-%d")
	file_dir, file_name = os.path.split(infile)
	proj_name = file_name.split(".")[0]
	work_dir = os.path.join(file_dir, "{0}-{1}".format(start_date, proj_name))
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	print "Project directory is {0}".format(work_dir)

if args.repeat_family_models:
	repeat_family_models = args.repeat_family_models
	if not os.path.isfile(repeat_family_models):
		print "ERROR: Must supply file containing repeat families used for RepeatMasker or have pre-calculated family lengths"
		sys.exit()
else:
	if not os.path.isfile(os.path.join(work_dir, "mean_lengths.table")):
		print "ERROR: Must supply file containing repeat families used for RepeatMasker or have pre-calculated family lengths"
		sys.exit()
		
if args.ReTe_dir:
	ReTe_dir = args.ReTe_dir
	ReTe_workplace = os.path.join(ReTe_dir, "Workplace")
else:
	print "ERROR: Must provide directory containiug the RetroTectos jar files ans associated folders"
	sys.exit()
	
if args.genome_file:
	genome_file = args.genome_file
else:
	print "ERROR: Must provide reference genome"
	sys.exit()

if args.gff_file:
	gff_file = args.gff_file
else:
	print "ERROR: Must provide gene annotation file in gff format"
	sys.exit()
	
if args.num_proc:
	num_proc = int(args.num_proc)

tmp_dir = os.path.join(work_dir, "tmp")
if not os.path.exists(tmp_dir):
	os.mkdir(os.path.join(work_dir, "tmp"))

job_dir = os.path.join(work_dir, "job_stats")
if not os.path.exists(job_dir):
	os.mkdir(job_dir)
	
script_dir = os.path.dirname(os.path.abspath(__file__))
accessory_scripts = os.path.join(script_dir, "refine_repeatmasker_functions.py")

###############################################################################
# EXTRA FUNCTIONS
###############################################################################
def clean_ReTe():
	shutil.rmtree(ReTe_workplace)
	new_dna_dir = os.path.join(ReTe_workplace, "NewDNA")
	print "Creating clean RetroTector workplace: ", new_dna_dir
	os.makedirs(new_dna_dir)

def slice_sequence(sequence_file, start, end):
	if int(start) < int(end):
		coords = int(start) - 1
		coorde = int(end)
	else:
		coords = int(end) - 1
		coorde = int(start)		
	try:
		with open(sequence_file, "r") as infile:
			scaffold_fasta = infile.readlines()
			scaffold = str(scaffold_fasta[1])
			sliced_seq = scaffold[coords:coorde]
			return sliced_seq
	except IOError:
		print "\tError: slice_sequence: Sequence file {0} does not exist".format(sequence_file)
			
def reverse_complement(sequence):
	complement = {
	"A": "T",
	"C": "G",
	"G": "C",
	"T": "A",
	"M": "K",
	"R": "Y",
	"W": "W",
	"S": "S",
	"Y": "R",
	"K": "M",
	"V": "B",
	"H": "D",
	"D": "H",
	"B": "V",
	"X": "X",
	"N": "N",
	}	## LIFTED FROM BIOPYTHON
	bases = list(sequence) 
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	return bases

def check_blast_db(db_file): ## Assumes nucleotide database for RebBase/RepeatModeler libraries
	db_file_types = ["nhr", "nin", "nsq"]
	search_list = []
	db_dir, db_fname = os.path.split(db_file)
	basename = db_fname.rsplit(".", 1)[0]

	for fname in os.listdir(db_dir):
		if basename in fname:
			ftype = fname.rsplit(".", 1)[1]
			if ftype in db_file_types:
				search_list.append(ftype)
	list_comparison = cmp(sorted(db_file_types), sorted(search_list))
	if list_comparison == 0:
		return db_file
	else:
		run_makeblastdb = "makeblastdb -in {0} -dbtype nucl -out {0}".format(db_file)
		print "Creating blast database for {0}".format(db_file)
		subprocess.check_call(run_makeblastdb, shell=True)
		return db_file

def get_mean_class_length():
	print "Creating lookup table for repeat class sizes"
	
	mean_length_dict = {}
	
	## Check if model family lengths have been calculated
	length_file = os.path.join(work_dir, "mean_class_lengths.table")
	if os.path.isfile(length_file):
		with open(length_file, "rb") as dict_file:
			for line in dict_file:
				rep_class, ave_len = line.strip().split()
				if not rep_class in mean_length_dict.iterkeys():
					mean_length_dict[rep_class] = int(ave_len)

	else:	
		length_dict = {}
		with open(repeat_family_models, "r") as rmout:
			for block in re.split("^>", rmout.read(), flags=re.MULTILINE):
				if block:
					seq_block = block.split("\n")
					if seq_block[0] != "":
						header = seq_block[0].strip()
						repeat_info = header.split()[0]
						repeat_group = repeat_info.split("#")[1]
						if "/" in repeat_group:
							rm_class, rm_family = repeat_group.split("/")
						else:
							rm_class = repeat_group
						if not rm_class in length_dict.iterkeys():
							length_dict[rm_class] = []
					sequence = ""
					for line in seq_block[1:]:
						sequence += line.strip()
					seqlen = len(sequence)
					length_dict[rm_class].append(seqlen)
					
		mean_length_dict = {}
		for rm_class, len_list in length_dict.iteritems():
			count = 0
			total_length = 0
			for i in len_list:
				total_length += i
				count += 1
			ave_len = int(total_length/count)
			with open(length_file, "a") as table_file:
				table_file.write("{0}\t{1}\n".format(rm_class, ave_len))
			
			if not rm_class in mean_length_dict.iterkeys():
				mean_length_dict[rm_class] = int(ave_len)
					
	return mean_length_dict
		
def get_mean_family_length():
	print "Creating lookup table for repeat family sizes"
	
	mean_length_dict = {}
	
	## Check if model family lengths have been calculated
	length_file = os.path.join(work_dir, "mean_family_lengths.table")
	if os.path.isfile(length_file):
		with open(length_file, "rb") as dict_file:
			for line in dict_file:
				rep_class, rep_family, ave_len = line.strip().split()
				if not rep_class in mean_length_dict.iterkeys():
					mean_length_dict[rep_class] = {}
				if not rep_family in mean_length_dict[rep_class].iterkeys():
					mean_length_dict[rep_class][rep_family] = int(ave_len)
	else:	
		length_dict = {}
		with open(repeat_family_models, "r") as rmout:
			for block in re.split("^>", rmout.read(), flags=re.MULTILINE):
				if block:
					seq_block = block.split("\n")
					if seq_block[0] != "":
						header = seq_block[0].strip()
						repeat_info = header.split()[0]
						repeat_group = repeat_info.split("#")[1]
						if "/" in repeat_group:
							rm_class, rm_family = repeat_group.split("/")
						else:
							rm_class = repeat_group
							rm_family = "N/A"
						if not rm_class in length_dict.iterkeys():
							length_dict[rm_class] = {}
						if not rm_family in length_dict[rm_class].iterkeys():
							length_dict[rm_class][rm_family] = []
					sequence = ""
					for line in seq_block[1:]:
						sequence += line.strip()
					seqlen = len(sequence)
					length_dict[rm_class][rm_family].append(seqlen)
					
		mean_length_dict = {}
		for rm_class, class_breakdown in length_dict.iteritems():
			for family, len_list in class_breakdown.iteritems():
				count = 0
				total_length = 0
				for i in len_list:
					total_length += i
					count += 1
				ave_len = int(total_length/count)
				with open(length_file, "a") as table_file:
					table_file.write("{0}\t{1}\t{2}\n".format(rm_class, family, ave_len))
				
				if not rm_class in mean_length_dict.iterkeys():
					mean_length_dict[rm_class] = {}
				if not family in mean_length_dict[rm_class].iterkeys():
					mean_length_dict[rm_class][family] = int(ave_len)
					
	return mean_length_dict

def get_repeat_model_lengths(repeat_family_models):
	print "Creating lookup table for repeat element sizes"
	
	model_length_dict = {}
	
	## Check if model family lengths have been calculated
	length_file = os.path.join(work_dir, "model_repeat_lengths.table")
	if os.path.isfile(length_file):
		with open(length_file, "rb") as dict_file:
			for line in dict_file:
				rep_class, rep_family, rep_name, rep_len = line.strip().split()
				if not rep_class in model_length_dict.iterkeys():
					model_length_dict[rep_class] = {}
				if not rep_family in model_length_dict[rep_class].iterkeys():
					model_length_dict[rep_class][rep_family] = {}
				if not rep_name in model_length_dict[rep_class][rep_family].iterkeys():
					model_length_dict[rep_class][rep_family][rep_name] = int(rep_len)
	else:	
		with open(repeat_family_models, "r") as rmout:
			for block in re.split("^>", rmout.read(), flags=re.MULTILINE):
				if block:
					seq_block = block.split("\n")
					if seq_block[0] != "":
						header = seq_block[0].strip()
						repeat_info = header.split()[0]
						rm_name, repeat_group = repeat_info.split("#")
						if "/" in repeat_group:
							rm_class, rm_family = repeat_group.split("/")
						else:
							rm_class = repeat_group
							rm_family = "N/A"
						if not rm_class in model_length_dict.iterkeys():
							model_length_dict[rm_class] = {}
						if not rm_family in model_length_dict[rm_class].iterkeys():
							model_length_dict[rm_class][rm_family] = {}
					sequence = ""
					for line in seq_block[1:]:
						sequence += line.strip()
					seqlen = len(sequence)
					model_length_dict[rm_class][rm_family][rm_name] = int(seqlen)

		for rm_class, class_breakdown in model_length_dict.iteritems():
			for rm_family, repeat_models in class_breakdown.iteritems():
				for rm_name, length in repeat_models.iteritems():
					with open(length_file, "a") as table_file:
						table_file.write("{0}\t{1}\t{2}\t{3}\n".format(rm_class, rm_family, rm_name, length))
	return model_length_dict	

repeat_length_dict = get_repeat_model_lengths(repeat_family_models)
family_length_dict = get_mean_family_length() ## This might break on a clean run
class_length_dict = get_mean_class_length()

def get_adjacent (rep_class, rep_fam, rep_name, start, end, strand, contig, rstart, rend, element_quicksect):
	feature_dict = {}
	element_length = int(end) - int(start) + 1
	
	## search 5kb downstream for same repeats - updates start and end coordinates
	downstream_list = element_quicksect.after(int(end), 100, 5000) 
	if downstream_list:
		intersecter = IntervalTree()
		intersecter.insert_interval(Interval(int(rstart), int(rend)))

		keep_list = []
		for repeat in downstream_list:
			if repeat.value["rep_class"] == "Simple_repeat" or repeat.value["rep_class"] == "Low_complexity" or repeat.value["rep_class"] == "Satellite":
				keep_list.append(repeat)
			else:  ##  filter out fragments
				try:
					model_length = repeat_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]][repeat.value["rep_name"]]
				except KeyError:
					try:
						if repeat.value["rep_class"] == "srpRNA":
							model_length = family_length_dict["SINE"]["Alu"]
						else:
							model_length = family_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]]
					except KeyError:
						model_length = class_length_dict[repeat.value["rep_class"]]
				repeat_length = int(repeat.end) - int(repeat.start) + 1

				if repeat.value["rep_class"] == rep_class and repeat.value["rep_name"] == rep_name:  ## same repeat
					fragment_start = repeat.value["rstart"]
					fragment_end = repeat.value["rend"]
					check_overlap = intersecter.find(int(fragment_start), int(fragment_end))
					if check_overlap: ## repeat is probably a separate instance
						if int(repeat_length) > int(model_length)/10:
							keep_list.append(repeat)
					else: ## repeat is probably part of the same insertion - add to cumulative repeat length to check if still a fragment
						if int(repeat.start) < int(start) and repeat.strand == strand:
							start = repeat.start
							element_length += repeat_length
						elif int(repeat.end) > int(end) and repeat.strand == strand:
							end = repeat.end
							element_length += repeat_length
						elif repeat.strand != strand:
							keep_list.append(repeat) ## Repeat is on a different strand and therefore still a different insertion	
				elif int(repeat_length) > int(model_length)/10: ## Repeat is not likely to be a fragment
					keep_list.append(repeat)
				
		for repeat in downstream_list:
			if repeat in keep_list:
				repeat_status = "keep"
			else:
				repeat_status = "remove"
			
			rm_contig = repeat.chrom
			rm_start = repeat.start
			rm_end = repeat.end
			rm_class = repeat.value["rep_class"]
			rm_fam = repeat.value["rep_fam"]
			rm_name = repeat.value["rep_name"]
			rm_strand = repeat.strand
			rm_id = "{0}-{1}-{2}".format(rm_contig, rm_start, rm_end)
			
			if not rm_class in feature_dict.iterkeys():
				feature_dict[rm_class] = {}
			if not rm_fam in feature_dict[rm_class].iterkeys():
				feature_dict[rm_class][rm_fam] = {}
			if not rm_name in feature_dict[rm_class][rm_fam].iterkeys():
				feature_dict[rm_class][rm_fam][rm_name] = {}
			if not rm_id in feature_dict[rm_class][rm_fam][rm_name].iterkeys():
				feature_dict[rm_class][rm_fam][rm_name][rm_id] = {}
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Contig"] = rm_contig
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Strand"] = rm_strand
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Status"] = repeat_status
		
		## check original element or extended element to make sure it is not a fragment
		if rep_class == "Simple_repeat" or rep_class == "Low_complexity" or rep_class == "Satellite":
			element_status = "keep"
		else:
			try:
				model_length = repeat_length_dict[rep_class][rep_fam][rep_name]	
			except KeyError:
				try:
					if rep_class == "srpRNA":
						model_length = family_length_dict["SINE"]["Alu"]
					else:
						model_length = family_length_dict[rep_class][rep_fam]
				except KeyError:
					model_length = class_length_dict[rep_class]

			if int(element_length) > int(model_length)/10:
				element_status = "keep"
				
				rep_id = "{0}-{1}-{2}".format(contig, start, end)
				
				if not rep_class in feature_dict.iterkeys():
					feature_dict[rep_class] = {}
				if not rep_fam in feature_dict[rep_class].iterkeys():
					feature_dict[rep_class][rep_fam] = {}
				if not rep_name in feature_dict[rep_class][rep_fam].iterkeys():
					feature_dict[rep_class][rep_fam][rep_name] = {}
				if not rep_id in feature_dict[rep_class][rep_fam][rep_name].iterkeys():
					feature_dict[rep_class][rep_fam][rep_name][rep_id] = {}
					feature_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
					feature_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
					feature_dict[rep_class][rep_fam][rep_name][rep_id]["Status"] = element_status
			else:
				element_status = "remove"

		## check overlapping repeats and filter out fragments - only needed if element is not a fragment
		if element_status == "keep":
			overlap_list = element_quicksect.find(int(start), int(end))
			
			for repeat in overlap_list:
				if repeat.value["rep_class"] == "Simple_repeat" or repeat.value["rep_class"] == "Low_complexity" or repeat.value["rep_class"] == "Satellite":
					repeat_status = "nested"
				else:
					repeat_length = int(repeat.end) - int(repeat.start) + 1
					try:
						model_length = repeat_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]][repeat.value["rep_name"]]
					except KeyError:
						try:
							if repeat.value["rep_class"] == "srpRNA":
								model_length = family_length_dict["SINE"]["Alu"]
							else:
								model_length = family_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]]
						except KeyError:
							model_length = class_length_dict[repeat.value["rep_class"]]

					if int(repeat_length) > int(model_length)/10: ## Repeat is not likely to be a fragment
						if int(repeat.start) >= int(start) and int(repeat.end) <= int(end):
							repeat_status = "nested"
						else:
							repeat_status = "keep"
					else:
						repeat_status = "remove"
							
					rm_contig = repeat.chrom
					rm_start = repeat.start
					rm_end = repeat.end
					rm_class = repeat.value["rep_class"]
					rm_fam = repeat.value["rep_fam"]
					rm_name = repeat.value["rep_name"]
					rm_strand = repeat.strand
					rm_id = "{0}-{1}-{2}".format(rm_contig, rm_start, rm_end)
					
					if not rm_class in feature_dict.iterkeys():
						feature_dict[rm_class] = {}
					if not rm_fam in feature_dict[rm_class].iterkeys():
						feature_dict[rm_class][rm_fam] = {}
					if not rm_name in feature_dict[rm_class][rm_fam].iterkeys():
						feature_dict[rm_class][rm_fam][rm_name] = {}
					if not rm_id in feature_dict[rm_class][rm_fam][rm_name].iterkeys():
						feature_dict[rm_class][rm_fam][rm_name][rm_id] = {}
						feature_dict[rm_class][rm_fam][rm_name][rm_id]["Contig"] = rm_contig
						feature_dict[rm_class][rm_fam][rm_name][rm_id]["Strand"] = rm_strand
						feature_dict[rm_class][rm_fam][rm_name][rm_id]["Status"] = repeat_status
		
	return feature_dict

def find_intersection(rep_class, rep_fam, rep_name, rep_id, start, end, contig, strand, element_quicksect):
	feature_dict = {}
	overlap_list = element_quicksect.find(int(start), int(end))
	keep_list = []
	for repeat in overlap_list:
		if int(repeat.start) < int(start) and repeat.strand == strand:
			start = repeat.start
		elif int(repeat.end) > int(end) and repeat.strand == strand:
			end = repeat.end
		elif repeat.strand != strand:
			keep_list.append(repeat) ## Repeat is on a different strand and therefore still a different insertion
	for repeat in overlap_list:
		if repeat in keep_list:
			status = "keep"
		else:
			status = "remove"
		rm_contig = repeat.chrom
		rm_start = repeat.start
		rm_end = repeat.end
		rm_class = repeat.value["rep_class"]
		rm_fam = repeat.value["rep_fam"]
		rm_name = repeat.value["rep_name"]
		rm_strand = repeat.strand
		rm_id = "{0}-{1}-{2}".format(rm_contig, rm_start, rm_end)
		
		if not rm_class in feature_dict.iterkeys():
			feature_dict[rm_class] = {}
		if not rm_fam in feature_dict[rm_class].iterkeys():
			feature_dict[rm_class][rm_fam] = {}
		if not rm_name in feature_dict[rm_class][rm_fam].iterkeys():
			feature_dict[rm_class][rm_fam][rm_name] = {}
		if not rm_id in feature_dict[rm_class][rm_fam][rm_name].iterkeys():
			feature_dict[rm_class][rm_fam][rm_name][rm_id] = {}
			feature_dict[rm_class][rm_fam][rm_name][rm_id]["Contig"] = rm_contig
			feature_dict[rm_class][rm_fam][rm_name][rm_id]["Strand"] = rm_strand
			feature_dict[rm_class][rm_fam][rm_name][rm_id]["Status"] = status
	
	## Update original element
	rep_id = "{0}-{1}-{2}".format(contig, start, end)

	if not rep_class in feature_dict.iterkeys():
		feature_dict[rep_class] = {}
	if not rep_fam in feature_dict[rep_class].iterkeys():
		feature_dict[rep_class][rep_fam] = {}
	if not rep_name in feature_dict[rep_class][rep_fam].iterkeys():
		feature_dict[rep_class][rep_fam][rep_name] = {}
	feature_dict[rep_class][rep_fam][rep_name][rep_id] = {}
	feature_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
	feature_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
	feature_dict[rep_class][rep_fam][rep_name][rep_id]["Status"] = "keep"

	return feature_dict

def merge_overlaps(contig, ERV_id, ERV_start, ERV_end, rep_class, rep_fam, rep_name, break_coords, element_quicksect):
	feature_dict = {}

	fragment_list = []
	if len(break_coords) > 1: ## Fragmented ERV (broken pipe)
		for broken_pipe in break_coords:
			(bstart, bend) = broken_pipe
			if int(bstart) < int(bend):
				start = int(bstart)
				end = int(bend)
				strand = 1
			else:
				start = int(bend)
				end = int(bstart)
				strand = -1
			fragment_list.extend([start, end])
			repeat_list = element_quicksect.find(start, end) ## Evaluate repeatmasker calls
			for repeat in repeat_list:
				status = "keep"
				if repeat.value["rep_class"] == rep_class or repeat.value["rep_class"] == "Unknown": ## feature is an LTR element
					status = "remove" ## Ignore feature - very likely to be an ERV fragment
				else:
					repeat_length = int(repeat.end) - int(repeat.start)
					if repeat.value["rep_class"] == "Simple_repeat" or repeat.value["rep_class"] == "Low_complexity" or repeat.value["rep_class"] == "Satellite":
						status = "nested"
					elif repeat.value["rep_class"] in repeat_length_dict.iterkeys():
						try:
							model_length = repeat_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]][repeat.value["rep_name"]]
						except KeyError:
							try:
								if repeat.value["rep_class"] == "srpRNA":
									model_length = family_length_dict["SINE"]["Alu"]
								else:
									model_length = family_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]]
							except KeyError:
								model_length = class_length_dict[repeat.value["rep_class"]]
						if repeat_length > int(model_length)/4: ## RM repeat is not a fragment
							status = "keep"					
						else: ##RM repeat is a fragment
							status = "remove"
					else:
						print "Merging overlaps in {0}-{1}: Non-model repeat {2} ({3})found.".format(start, end, repeat.value["rep_name"], repeat.value["rep_class"])
						status = "nested"

				rm_contig = repeat.chrom
				rm_start = repeat.start
				rm_end = repeat.end
				rm_class = repeat.value["rep_class"]
				rm_fam = repeat.value["rep_fam"]
				rm_name = repeat.value["rep_name"]
				rm_strand = repeat.strand
				
				rm_id = "{0}-{1}-{2}".format(rm_contig, rm_start, rm_end)
				
				if not rm_class in feature_dict.iterkeys():
					feature_dict[rm_class] = {}
				if not rm_fam in feature_dict[rm_class].iterkeys():
					feature_dict[rm_class][rm_fam] = {}
				if not rm_name in feature_dict[rm_class][rm_fam].iterkeys():
					feature_dict[rm_class][rm_fam][rm_name] = {}
				if not rm_id in feature_dict[rm_class][rm_fam][rm_name].iterkeys():
					feature_dict[rm_class][rm_fam][rm_name][rm_id] = {}
					feature_dict[rm_class][rm_fam][rm_name][rm_id]["Contig"] = rm_contig
					feature_dict[rm_class][rm_fam][rm_name][rm_id]["Strand"] = rm_strand
					feature_dict[rm_class][rm_fam][rm_name][rm_id]["Status"] = status
					feature_dict[rm_class][rm_fam][rm_name][rm_id]["Coords"] = ""	
					feature_dict[rm_class][rm_fam][rm_name][rm_id]["RT_ref"] = ""	
									
	elif len(break_coords) == 1: ## ERV is intact				
		if int(ERV_start) < int(ERV_end):
			start = int(ERV_start)
			end = int(ERV_end)
			strand = 1
		else:
			start = int(ERV_end)
			end = int(ERV_start)
			strand = -1
		fragment_list.extend([start, end])
		repeat_list = element_quicksect.find(start, end)
		for repeat in repeat_list:
			if repeat.value["rep_class"] == rep_class or repeat.value["rep_class"] == "Unknown": ## feature is an LTR element
				status = "remove" ## Ignore feature - very likely to be an ERV fragment
			else:
				repeat_length = int(repeat.end) - int(repeat.start)
				if repeat.value["rep_class"] in repeat_length_dict.iterkeys():
					try:
						model_length = repeat_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]][repeat.value["rep_name"]]
					except KeyError:
						try:
							if repeat.value["rep_class"] == "srpRNA":
								model_length = family_length_dict["SINE"]["Alu"]
							else:
								model_length = family_length_dict[repeat.value["rep_class"]][repeat.value["rep_fam"]]
						except KeyError:
							model_length = class_length_dict[repeat.value["rep_class"]]
				else:
					print "Merging overlaps in {0}-{1}: Non-model repeat {2} ({3})found.".format(start, end, repeat.value["rep_name"], repeat.value["rep_class"])
					model_length = 50 ## Semi arbitrary figure to ensure that short simple repeats and low complexity repeats are included
				
				if repeat_length > int(model_length)/4: ## RM repeat is not a fragment
					status = "keep"
					fragment_end = int(repeat.start) - 1
					next_fragment_start = int(repeat.end) + 1
					fragment_list.extend([fragment_end, next_fragment_start])					
				else: ##RM repeat is a fragment
					status = "remove"

			rm_contig = repeat.chrom
			rm_start = repeat.start
			rm_end = repeat.end
			rm_class = repeat.value["rep_class"]
			rm_fam = repeat.value["rep_fam"]
			rm_name = repeat.value["rep_name"]
			rm_strand = repeat.strand
							
			rm_id = "{0}-{1}-{2}".format(rm_contig, rm_start, rm_end)
			
			if not rm_class in feature_dict.iterkeys():
				feature_dict[rm_class] = {}
			if not rm_fam in feature_dict[rm_class].iterkeys():
				feature_dict[rm_class][rm_fam] = {}
			if not rm_name in feature_dict[rm_class][rm_fam].iterkeys():
				feature_dict[rm_class][rm_fam][rm_name] = {}
			if not rm_id in feature_dict[rm_class][rm_fam][rm_name].iterkeys():
				feature_dict[rm_class][rm_fam][rm_name][rm_id] = {}
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Contig"] = rm_contig
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Strand"] = rm_strand
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Status"] = status
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["Coords"] = ""
				feature_dict[rm_class][rm_fam][rm_name][rm_id]["RT_ref"] = ""
					
	else: ## There is an error somewhere....
		print "Error: filter_RMout: ERV coordinates lost?", contig, ERV_id, ERV_start, ERV_end, rep_class, rep_fam, rep_name, break_coords

	if len(fragment_list) % 2 == 0: ## Check that all coordinates are there
		status = "keep"
		fragment_list = [int(x) for x in fragment_list]
		sorted_list = sorted(fragment_list)
		ERV_start = sorted_list[0]
		ERV_end = sorted_list[-1]
#		print sorted_list
		coord_list = zip(sorted_list[::2], sorted_list[1::2])
#		print coord_list
		rep_id = "{0}-{1}-{2}".format(contig, ERV_start, ERV_end)
		if not rep_class in feature_dict.iterkeys():
			feature_dict[rep_class] = {}
		if not rep_fam in feature_dict[rep_class].iterkeys():
			feature_dict[rep_class][rep_fam] = {}
		if not rep_name in feature_dict[rep_class][rep_fam].iterkeys():
			feature_dict[rep_class][rep_fam][rep_name] = {}
		if not rep_id in feature_dict[rep_class][rep_fam][rep_name].iterkeys():
			feature_dict[rep_class][rep_fam][rep_name][rep_id] = {}
			feature_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
			feature_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
			feature_dict[rep_class][rep_fam][rep_name][rep_id]["Status"] = status
			feature_dict[rep_class][rep_fam][rep_name][rep_id]["Coords"] = coord_list	
			feature_dict[rep_class][rep_fam][rep_name][rep_id]["RT_ref"] = ERV_id	
		
	return feature_dict

def merge_genes(contig, start, end, interval_tree):
	overlap_list = interval_tree.find(int(start), int(end))
	overlap_feature_types = []
	for gene_feature in overlap_list:
		overlap_feature_types.append(gene_feature.value["feature_type"])
	if not overlap_feature_types:
		keep_status = True
		region_type = "between_genes"
	elif "tRNA" in overlap_feature_types or "CDS" in overlap_feature_types:
		keep_status = False
		region_type = "coding"
	elif "gene" in overlap_feature_types or "mRNA" in overlap_feature_types:
		if "five_prime_UTR" in overlap_feature_types or "three_prime_UTR" in overlap_feature_types:
			keep_status = True
			region_type = "UTR"
		else:
			keep_status = True
			region_type = "intron"
	else:
		print "Warning: merge_with_genes: Unknown feature type in {0}:{1}-{2}. Treating as non-coding".format(contig, start, end)
		keep_status = True
		region_type = "unknown"
	###check overlapping interval types: if CDS, tRNA, protein_match, mRNA, then discard. If overlapping exon/UTR/gene, then retain and flag as intergenic
	return keep_status, region_type ## Keep=False, coding; keep=True, intron, UTR, between

def touch_file(infile):
	with open(infile, "a"):
		os.utime(infile, None)

def read_gz(gz_file):
	fq_file = gzip.open(gz_file, "rb")
	fq_data = fq_file.read()	
	return fq_data
	
###############################################################################
# PIPELINE TASKS
###############################################################################

@split(infile, "{0}/By_contig/*.csv".format(work_dir))
def read_RMout(infile, outfiles):
	out_dir = os.path.join(work_dir, "By_contig")
	if os.path.exists(out_dir):
		for contig_file in os.listdir(out_dir):
			print "WARNING: read_RMout: Removing contig file: ", contig_file
			os.unlink(os.path.join(out_dir, contig_file))
	else:
		os.makedirs(out_dir)
	
	with open(infile, "rb") as rm_outfile:
		for a in range(3):
			rm_outfile.readline()
		for line in rm_outfile:
			if len(line.strip().split()) > 14:
				score, divergence, dels, ins, qseq, qstart, qend, qrem, strand, rep_name, rep_group, rcoord1, rcoord2, rcoord3, rubbish = line.strip().split(None, 14)
			elif len(line.strip().split()) == 14:
				score, divergence, dels, ins, qseq, qstart, qend, qrem, strand, rep_name, rep_group, rcoord1, rcoord2, rcoord3 = line.strip().split(None, 13)
			else:
				print "Error: read_RMout: invalid input: ", line
				sys.exit()
			
			outname = qseq + ".csv"
			
			with open(os.path.join(out_dir, outname), "ab") as outfile:
				writer = csv.writer(outfile, delimiter="\t")
				
				if strand == "+":
					strand = "1"
					rstart = rcoord1
					rend = rcoord2
					rrem = rcoord3.strip("(").strip(")")
				elif strand == "C":
					strand = "-1"
					rstart = rcoord3
					rend = rcoord2
					rrem = rcoord1.strip("(").strip(")")
				else:
					print "Error: read_RMout: invalid strand designation in line: ", line
					sys.exit()
				
				if "/" in rep_group:
					rep_class, rep_fam = rep_group.split("/")
				else:
					rep_class = rep_group
					rep_fam = "N/A"
				writer.writerow([qseq, qstart, qend, qrem, strand, rep_name, rep_class, rep_fam, rstart, rend, rrem])


@subdivide(read_RMout, regex(r"(.+)\/(.+).csv"), [r"\1/\2.dict", r"\1/\2.nested.dict"])
def filter_RMout(infile, outfiles):
	element_dict = {}
	element_list = []
	contig_quicksect = IntervalTree()
	
	with open(infile, "rb") as RM_out:
		for line in RM_out:
			qseq, qstart, qend, qrem, qstrand, rep_name, rep_class, rep_fam, rstart, rend, rrem = line.strip().split()
			rep_id = "{0}-{1}-{2}".format(qseq, qstart, qend)
			repeat_info = (rep_class, rep_fam, rep_name, qseq, int(qstart), int(qend), qstrand, int(rstart), int(rend)) 
			element_list.append(repeat_info)
			contig_quicksect.insert_interval(Interval(int(qstart), int(qend),
								chrom = qseq,
								strand = qstrand,
								value = {"rep_class" : rep_class,
										"rep_fam" : rep_fam,
										"rep_name": rep_name,
										"rstart" : int(rstart),
										"rend" : int(rend)}
								))
			
			if not rep_class in element_dict.iterkeys():
				element_dict[rep_class] = {}
			if not rep_fam in element_dict[rep_class].iterkeys():
				element_dict[rep_class][rep_fam] = {}
			if not rep_name in element_dict[rep_class][rep_fam].iterkeys():
				element_dict[rep_class][rep_fam][rep_name] = {}
			if not rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
				element_dict[rep_class][rep_fam][rep_name][rep_id] = {}
				element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = qseq
				element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = qstrand

	nested_dict = {}	
	for repeat in sorted(element_list, key=itemgetter(5)):
		(rep_class, rep_fam, rep_name, qseq, qstart, qend, qstrand, rstart, rend) = repeat
		if int(rstart) > int(rend):
			rep_start = int(rend)
			rep_end = int(rstart)
		else:
			rep_start = int(rstart)
			rep_end = int(rend)
		filtered_block = get_adjacent(rep_class, rep_fam, rep_name, qstart, qend, qstrand, qseq, rep_start, rep_end, contig_quicksect)
		if filtered_block:
			for rep_class, val1 in filtered_block.iteritems():
				for rep_fam, val2 in val1.iteritems():
					for rep_name, val3 in val2.iteritems():
						for rep_id, info in val3.iteritems():
							contig = info["Contig"]
							strand = info["Strand"]
							status = info["Status"]
							if status == "nested":
								if rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
									del element_dict[rep_class][rep_fam][rep_name][rep_id]
								if not rep_class in nested_dict.iterkeys():
									nested_dict[rep_class] = {}
								if not rep_fam in nested_dict[rep_class].iterkeys():
									nested_dict[rep_class][rep_fam] = {}
								if not rep_name in nested_dict[rep_class][rep_fam].iterkeys():
									nested_dict[rep_class][rep_fam][rep_name] = {}
								if not rep_id in nested_dict[rep_class][rep_fam][rep_name].iterkeys():
									nested_dict[rep_class][rep_fam][rep_name][rep_id] = {}
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
							elif status == "remove":
								if rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
									del element_dict[rep_class][rep_fam][rep_name][rep_id]
							elif status == "keep":
								if rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
									element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
									element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
								else:
									element_dict[rep_class][rep_fam][rep_name][rep_id] = {}
									element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
									element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
		
	
	## clean up remaining overlaps and empty categories
	contig_tree = {}
	list_dict = {}
	for rep_class, values in element_dict.iteritems():
		for rep_fam, values in values.iteritems():
			for rep_name, values in values.iteritems():
				if not rep_name in contig_tree.iterkeys():
					rep_name_quicksect = IntervalTree()
					contig_tree[rep_name] = rep_name_quicksect
				rep_name_quicksect = contig_tree[rep_name]
				if not rep_name in list_dict.iterkeys():
					list_dict[rep_name] = []
				for rep_id, info in values.iteritems():
					contig, start, end = rep_id.split("-")

					rep_name_quicksect.insert_interval(Interval(int(start), int(end),
								chrom = contig,
								strand = info["Strand"],
								value = {"rep_class" : rep_class,
										"rep_fam" : rep_fam,
										"rep_name": rep_name}))

					element_info = (rep_class, rep_fam, rep_name, rep_id, int(start), int(end), info["Strand"])
					list_dict[rep_name].append(element_info)
				contig_tree[rep_name] = rep_name_quicksect	
	
	for rep_name, element_list in list_dict.iteritems():
		for repeat in sorted(element_list, key=itemgetter(4)):
			(rep_class, rep_fam, rep_name, rep_id, start, end, strand) = repeat
			contig, start, end = rep_id.split("-")
			if rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
				overlaps = find_intersection(rep_class, rep_fam, rep_name, rep_id, start, end, contig, strand, contig_tree[rep_name])
				if overlaps:
					for rep_class, values in overlaps.iteritems():
						for rep_fam, values in values.iteritems():
							for rep_name, values in values.iteritems():
								for rep_id, info in values.iteritems():
									contig = info["Contig"]
									strand = info["Strand"]
									status = info["Status"]
									if status == "remove":
										if rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
											del element_dict[rep_class][rep_fam][rep_name][rep_id]
									elif status == "keep":
										if rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
											element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
											element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
										else:
											element_dict[rep_class][rep_fam][rep_name][rep_id] = {}
											element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
											element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
									else:
										print "WARNING: filter_RMout: ", rep_id, info

	## Test print
#	with open(os.path.join(work_dir, "test-{0}.txt".format(contig)), "w") as test_file:
#		for rep_class, values in element_dict.iteritems():
#			test_file.write("{0}\n".format(rep_class)) 
#			for rep_fam, values in values.iteritems():
#				test_file.write("\t{0}\n".format(rep_fam))  
#				for rep_name, values in values.iteritems():
#					test_file.write("\t\t{0}\n".format(rep_name))  
#					for rep_ID, info in values.iteritems():
#						test_file.write("\t\t\t{0}\t{1}\n".format(rep_ID, info)) 

	for rep_class, values in nested_dict.iteritems():
		for rep_fam, values in values.iteritems():
			for rep_name, values in values.iteritems():
				for rep_id, info in values.iteritems():
					if not rep_class in element_dict.iterkeys():
						element_dict[rep_class] = {}
					if not rep_fam in element_dict[rep_class].iterkeys():
						element_dict[rep_class][rep_fam] = {}
					if not rep_name in element_dict[rep_class][rep_fam].iterkeys():
						element_dict[rep_class][rep_fam][rep_name] = {}
					if not rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
						element_dict[rep_class][rep_fam][rep_name][rep_id] = {}
						element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = info["Contig"]
						element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = info["Strand"]
	filtered_dict = {}
	for rep_class, values in element_dict.iteritems():
		for rep_fam, values in values.iteritems():
			for rep_name, values in values.iteritems():
				for rep_id, info in values.iteritems():
					if not rep_class in filtered_dict.iterkeys():
						filtered_dict[rep_class] = {}
					if not rep_fam in filtered_dict[rep_class].iterkeys():
						filtered_dict[rep_class][rep_fam] = {}
					if not rep_name in filtered_dict[rep_class][rep_fam].iterkeys():
						filtered_dict[rep_class][rep_fam][rep_name] = {}
					if not rep_id in filtered_dict[rep_class][rep_fam][rep_name].iterkeys():
						filtered_dict[rep_class][rep_fam][rep_name][rep_id] = {}
						filtered_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = info["Contig"]
						filtered_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = info["Strand"]
	
	## Test print
#	with open(os.path.join(work_dir, "test-{0}.txt".format(contig)), "w") as test_file:
#		test_file.write("Nested dict:\n")
#		for rep_class, values in nested_dict.iteritems():
#			test_file.write("{0}\n".format(rep_class)) 
#			for rep_fam, values in values.iteritems():
#				test_file.write("\t{0}\n".format(rep_fam))  
#				for rep_name, values in values.iteritems():
#					test_file.write("\t\t{0}\n".format(rep_name))  
#					for rep_ID, info in values.iteritems():
#						test_file.write("\t\t\t{0}\t{1}\n".format(rep_ID, info)) 
#	## Test print
#	with open(os.path.join(work_dir, "test-{0}.txt".format(contig)), "a") as test_file:
#		test_file.write("F dict:\n")
#		for rep_class, values in filtered_dict.iteritems():
#			test_file.write("{0}\n".format(rep_class)) 
#			for rep_fam, values in values.iteritems():
#				test_file.write("\t{0}\n".format(rep_fam))  
#				for rep_name, values in values.iteritems():
#					test_file.write("\t\t{0}\n".format(rep_name))  
#					for rep_ID, info in values.iteritems():
#						test_file.write("\t\t\t{0}\t{1}\n".format(rep_ID, info))  									

	with open(outfiles[0], "wb") as dict_file:
		cPickle.dump(filtered_dict, dict_file, 1)

	with open(outfiles[1], "wb") as nested_dict_file:
		cPickle.dump(nested_dict, nested_dict_file, 1)
		

## Directories for refinement of ERV calls
ERV_dir = os.path.join(work_dir, "ERV_refinement")
ReTe_ERVs = os.path.join(ERV_dir, "RetroTector_output")
contig_dir = os.path.join(ERV_dir, "contigs")
chain_dir = os.path.join(ERV_dir, "RetroTector_chainfiles")


@split(genome_file, "{0}/*.fas".format(contig_dir))
def get_ERV_contigs(infile, outfiles):
	if not os.path.exists(contig_dir):
		os.makedirs(contig_dir)	

	## Clean up previous runs
	for contig_file in os.listdir(contig_dir):
		print "WARNING: get_ERV_contigs: Removing contig file: ", contig_file
		os.unlink(os.path.join(contig_dir, contig_file))

	with open(genome_file, "rb") as genome:
		for block in re.split("^>", genome.read(), flags=re.MULTILINE):
			if block:
				seq_block = block.split("\n")
				if seq_block[0] != "":
					seq_ID = seq_block[0]
					sequence = ""
					for line in seq_block[1:]:
						sequence += line.strip()
					if len(sequence) > 1000:
						seq_file_name = "{0}.fas".format(seq_ID)
						with open(os.path.join(contig_dir, seq_file_name), "w") as outfile:
							outfile.write(">{0}\n{1}".format(seq_ID, sequence)) ## Writes to single line file
				else:
					break
					
#@follows(clean_ReTe)
@transform(get_ERV_contigs, regex(r"(.+)/(.+).fas"), r"{0}/NewDNA/\2.fas".format(ReTe_workplace))
def move_contigs(infile, outfile):	
	shutil.copy2(infile, '{0}/NewDNA/'.format(ReTe_workplace))

@merge(move_contigs, "{0}/ReTe_input_files.txt".format(ERV_dir))
def ReTe_init(infile, outfile):
	flagfile = open(outfile, "a")
	if not os.path.exists(ERV_dir):
		os.makedirs(ERV_dir)
	
	## Clean up previous runs in project directory
	if os.path.exists(ReTe_ERVs):
		shutil.rmtree(ReTe_ERVs)
	
	for seq_file in infile:
		header = seq_file.rsplit('/', 1)[1]
		flagfile.write(header + '\n')
	home_dir = os.getcwd()
	os.chdir('{0}'.format(ReTe_dir))
			 
	subprocess.call(["java","-Xmx30720m","-classpath", os.path.join(ReTe_dir, "RetroTectorEngine.jar"), "retrotector/RetroTectorEngine", "SweepDNA", "SweepScripts","Quit"])

	paths = "{0}/Workplace/".format(ReTe_dir)
	shutil.copytree(paths, ReTe_ERVs, ignore = shutil.ignore_patterns("NewDNA"))
	os.chdir(home_dir)
	
@follows(ReTe_init)
@split(ReTe_init, "{0}/solo_LTR.csv".format(ERV_dir))
def get_soloLTR (infile, outfile):
	with open(outfile, "wb") as soloLTR_file:
		listfile = csv.writer(soloLTR_file, delimiter="\t")
		for path, dirs, files in os.walk(ERV_dir):
			for fname in fnmatch.filter(files, "*SelectedChains.txt"):
				fdir = path
				for line in open(os.path.join(fdir, fname)):
					if line.startswith("DNAFile: "):
						scaffseq = line.split()[1]
						seqid = scaffseq.replace(".fas", "").rsplit(".", 1)[0]
				with open(os.path.join(fdir, fname), 'r') as sourcefile:
					for block in re.findall("(^[A-Z]SingleLTRs::.*?::)", 
											sourcefile.read(), 
											re.DOTALL|re.MULTILINE):
						ltrset = block.splitlines()[1:-1]
						for solo in ltrset:
							sololoc = solo.split()[1].strip("()")
							solos, soloe = sololoc.split("-")
							coords = int(solos)
							coorde = int(soloe)
							length = abs(coords - coorde) + 1
							listfile.writerow([seqid, solos, soloe, length])

@split(ReTe_init, "{0}/*.chain".format(chain_dir))
def make_chainfiles (infile, outfile):
	if not os.path.exists(chain_dir):
		os.makedirs(chain_dir)
	for path, dirs, files in os.walk(ReTe_ERVs):
		for fname in fnmatch.filter(files, "*SelectedChains.txt"):
			fdir = path
			for line in open(os.path.join(fdir, fname)):
				if line.startswith("DNAFile: "):
					scaffseq = line.split()[1]
					seqid = scaffseq.replace(".fas", "").rsplit(".", 1)[0]
			with open(os.path.join(fdir, fname), 'r') as sourcefile:
				chainno = 1
				for block in re.findall("(^Chain[A-Z}[0-9]+::.*?::)", sourcefile.read(), re.DOTALL|re.MULTILINE):
					header = "{0}_{1}".format(seqid, chainno)
					blockfile = os.path.join(chain_dir, "{0}.chain".format(header))
					if os.path.isfile(blockfile):
						chainno += 1
						header = "{0}_{1}".format(seqid, chainno)
						chainfile = open(os.path.join(chain_dir, "{0}.chain".format(header)), "w")
						chainfile.write("SeqFile: " + header + "\n" + block)
					else:
						chainfile = open(os.path.join(chain_dir, "{0}.chain".format(header)), "w")
						chainfile.write("SeqFile: " + header + "\n" + block)

@merge (make_chainfiles, os.path.join(ERV_dir, "LTR_sequences.fas"))
def get_LTR (infile, outfile):
	for path, dirs, files in os.walk(chain_dir):
		for fname in fnmatch.filter(files, "*.chain"):
			fdir = path
			chainfile = os.path.join(fdir, fname)
			for line in open(chainfile, 'r'):
				if line.startswith("SeqFile: "):
					scaffid = line.strip().split()[1]
				if line.startswith("Chain"):
					chainid = line.strip()
					ctype = chainid[5:6]
					cno = chainid[6:-2]
					chainno = "{0:0>5}".format(cno)
			seqid = "{0}-{1}{2}".format(scaffid, ctype, chainno)
			lines = open(chainfile, 'r').readlines()
			ltr5 = "N/A"
			ltr3 = "N/A"
			ltr5pat = re.compile("SubGene 5LTR")
			ltr3pat = re.compile("SubGene 3LTR")
			for i in range(len(lines)):
				if ltr5pat.search(lines[i]):
					ltr5line = lines[i+3]
					ltr5 = ltr5line.strip().split()[0]
				if ltr3pat.search(lines[i]):
					ltr3line = lines[i+3]
					ltr3 = ltr3line.strip().split()[0]
			ltrfile = open(outfile, 'a')
			if ltr5 != "N/A":
				ltrfile.write(">{0}-5LTR\n{1}\n".format(seqid, ltr5.upper()))
			if ltr3 != "N/A":
				ltrfile.write(">{0}-3LTR\n{1}\n".format(seqid, ltr3.upper()))

@merge (make_chainfiles, os.path.join(ERV_dir, "ERV_coordinates.csv"))
def get_ERV_info(infile, outfile):
	temp_file = tempfile.NamedTemporaryFile(mode="w", dir=tmp_dir, delete=False)
	ERV_coordinates = csv.writer(temp_file, delimiter="\t")
	for path, dirs, files in os.walk(chain_dir):
		for fname in fnmatch.filter(files, "*.chain"):
			fdir = path
			chainfile = os.path.join(fdir, fname)
			site5 = "N/A"
			site3 = "N/A"
			sbreak = "N/A"
			ebreak = "N/A"
			for line in open(chainfile, "r"):
				if line.startswith("SeqFile: "):
					splitid = line.strip().split()[1]
					scaffid = splitid.rsplit("_", 1)[0]
				if line.startswith("Chain"):
					chainid = line.strip()
					ctype = chainid[5:6]
					cno = chainid[6:-2]
					chainno = "{0:0>5}".format(cno)
				if line.strip().startswith("Broken"):
					break_coords = line.strip().split()[-1]
					try:
						sbreak, ebreak = break_coords.split(">")
					except:
						print "WARNING: get_ERV_info: looking at wrong line?", chainfile
				if line.startswith("Integration sites"):
					lset = line.strip()
					site = lset.rsplit(' ', 1)[1]
					site5 = site.split('<>')[0]
					site3 = site.split('<>')[1]
			seqid = "{0}-{1}{2}".format(scaffid, ctype, chainno)
			lines = open(chainfile, 'r').readlines()
			ltr5s = "N/A"
			ltr5e = "N/A"
			ltr3s = "N/A"
			ltr3e = "N/A"
			pbs = "N/A"
			ppt = "N/A"
			ervchain = re.compile('^Chain')
			ltr5pat = re.compile("SubGene 5LTR")
			ltr3pat = re.compile("SubGene 3LTR")
			pbspat = re.compile("SubGene PBS")
			pptpat = re.compile("SubGene PPT")
			for i in range(len(lines)):
				if ervchain.search(lines[i]):
					start_line = lines[i+11]
					start = start_line.strip()
					end_line = lines[i+12]
					end = end_line.strip()
				if ltr5pat.search(lines[i]):
					ltr5coord = lines[i+1].rsplit(None, 1)[1]
					ltr5s = ltr5coord.replace("[", "").split("-")[0]
					ltr5e = ltr5coord.replace("]", "").split("-")[1]
					ltr5line = lines[i+3]
					ltr5 = ltr5line.strip().split()[0]
				if ltr3pat.search(lines[i]):
					ltr3coord = lines[i+1].rsplit(None, 1)[1]
					ltr3s = ltr3coord.replace("[", "").split("-")[0]
					ltr3e = ltr3coord.replace("]", "").split("-")[1]
					ltr3line = lines[i+3]
					ltr3 = ltr3line.strip().split()[0]
				if pbspat.search(lines[i]):
					pbsline = lines[i+3]
					pbs = pbsline.strip().split()[0]
				if pptpat.search(lines[i]):
					pptline = lines[i+3]
					ppt = pptline.strip().split()[0]
			ERV_coordinates.writerow([seqid, start, end, ltr5s, ltr5e, ltr3s, ltr3e, pbs, ppt, site5, site3, sbreak, ebreak])
	temp_file.close()
	reader = csv.reader(open(os.path.join(work_dir, temp_file.name)), delimiter="\t")
	sortedfile = sorted(reader, key=operator.itemgetter(0))
	writer = csv.writer(open(outfile, "wb"), delimiter="\t")
	coords = []
	for row in sortedfile:
		if row not in coords:
			writer.writerow(row)
			coords.append(row)
	del coords[:]

@transform(get_ERV_info, regex(r"(.+)/ERV_coordinates.csv"), r"\1/ERV_sequences.fas")
def get_ERVs (infile, outfile): 
	reader = csv.reader(open(infile, "rb"), delimiter="\t")
	for row in reader:
		header = row[0]
		splitid = header.rsplit("-", 1)[0]
		scaffold = splitid.split("_")[0]
		for path, dirs, files in os.walk(contig_dir):
			for fname in fnmatch.filter(files, scaffold + ".fas"):
				fdir = path
				sequence_file = os.path.join(fdir, fname)
		coords = int(row[1])
		coorde = int(row[2])
		if coords < coorde:
			ERV_seq = slice_sequence(sequence_file, coords, coorde)
		else:
			raw_seq = slice_sequence(sequence_file, coords, coorde)
			ERV_seq = reverse_complement(raw_seq)
		seqfile = open(outfile, "a")
		if ERV_seq != "":
			seqfile.write(">{0}\n{1}\n".format(header, ERV_seq.upper()))

@transform(get_ERV_info, regex(r"(.+)/ERV_coordinates.csv"), r"\1/internal_ERV_sequences.fas")
def get_ERV_int(infile, outfile):
	
	## Clean up previous runs
	if os.path.isfile(outfile):
		os.unlink(outfile)
	
	reader = csv.reader(open(infile, "rb"), delimiter="\t")
	for row in reader:
		header = row[0]
		scaffold = header.rsplit("-", 1)[0]
		if "_" in scaffold:
			scaffold = scaffold.split("_")[0]
		fname = scaffold + ".fas"
		scaffold_file = os.path.join(contig_dir, fname)
		ints = row[4]
		if ints == "N/A":
			coords = int(row[1])
		else:
			coords = int(ints)
		inte = row[5]
		if inte == "N/A":
			coorde = int(row[2])
		else:
			coorde = int(inte)
		if row[11] != "N/A": ## Broken pipe
			sbreak = int(row[11])
			ebreak = int(row[12])
			start1 = coords
			end1 = sbreak
			ERV_seq1 = slice_sequence(scaffold_file, start1, end1)
			start2 = ebreak
			end2 = coorde
			ERV_seq2 = slice_sequence(scaffold_file, start2, end2)
			ERV_seq = ERV_seq1 + ERV_seq2	
		else:
			start = coords
			end = coorde
			ERV_seq = slice_sequence(scaffold_file, start, end)
		if coords < coorde: ## Strand: 1
			sequence = ERV_seq
		else: ## Strand: -1
			sequence = reverse_complement(ERV_seq)
		if sequence != "":
			with open(outfile, "a") as seqfile:
				try:
					seqfile.write(">{0}\n{1}\n".format(header, sequence.upper()))
				except:
					print "ERROR: get_ERV_int: invalid input? ", row
					continue

@transform(get_ERV_int, regex(r"(.+)/internal_ERV_sequences.fas"), r"\1/ERV_blastn_output.table")
def blast_ERVs(infile, outfile):
	blast_db = check_blast_db(repeat_family_models)
	print "Blast database: ", blast_db
	blast_input = "blastn -query {0} -db {1} -out {2} -outfmt \"6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore\" -task dc-megablast".format(infile, blast_db, outfile)
	subprocess.check_call(blast_input, shell=True)

@transform(blast_ERVs, regex(r"(.+)/ERV_blastn_output.table"), r"\1/ERV_blastn_output.filtered.table", r"\1/Manual_checklist.txt")
def get_top_hits(infile, outfile, check_file):
	filtered_dict = {}
	checklist = set()
	ReTe_ERV = ""
	top_bitscore = 0
	top_hit = ""
	aln_length = 0

	with open(infile, "r") as blast_output:
		for line in blast_output:
			qseqid, sseqid, pident, qlen, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split()
			min_length = int(qlen)/10
			if qseqid == ReTe_ERV:
				if Decimal(bitscore) > Decimal(top_bitscore):
					top_bitscore = bitscore
					top_hit = sseqid
					top_match = (top_hit, top_bitscore)
					filtered_dict[ReTe_ERV] = top_match
					aln_length = length
					if qseqid in checklist and int(length) >= int(min_length):
						checklist.discard(qseqid)
				elif Decimal(bitscore) == Decimal(top_bitscore):
					if int(length) > int(aln_length):
						aln_length = length
						top_bitscore = bitscore
						top_hit = sseqid
						top_match = (top_hit, top_bitscore)
						filtered_dict[ReTe_ERV] = top_match
						if qseqid in checklist and int(length) >= int(min_length):
							checklist.discard(qseqid)				
			else:
				if int(length) >= int(min_length):
					if not qseqid in filtered_dict.iterkeys():
						ReTe_ERV = qseqid
						top_bitscore = bitscore
						top_hit = sseqid
						top_match = (top_hit, top_bitscore)
						filtered_dict[ReTe_ERV] = top_match
						aln_length = length
					else:
						print "ERROR: get_top_hits: Something is wrong....", ReTe_ERV, top_hit, top_bitscore
				else:
					checklist.add(qseqid)

	with open(outfile, "w") as filtered_out:
		for ERV_ID, blast_hit in sorted(filtered_dict.iteritems()):
			top_match, score = (blast_hit)
			filtered_out.write("{0}\t{1}\n".format(ERV_ID, top_match))
	
	with open(check_file, "w") as check_list:
		for element in checklist:
			check_list.write("{0}\n".format(element))

@merge([get_ERV_info, get_top_hits], ["{0}/ERV_coordinates.classified.dict".format(ERV_dir), "{0}/ERV_info.classified.csv".format(ERV_dir)])
def reclassify_ERVs(infiles, outfiles):

	##Clean up
	for outfile in outfiles:
		if os.path.isfile(outfile):
			os.unlink(outfile)
			
	new_coord_file = outfiles[0]
	classified = outfiles[1]
	
	old_coord_file = infiles[0]
	blast_file = infiles[1]
	
	coord_dict = {}
	
	with open(blast_file, "rb") as infile1:
		for line in infile1:
			qseqid, RM_seqid = line.strip().split()
			rm_name, rm_group = RM_seqid.split("#")
			if "/" in rm_group:
				rm_class, ERV_type = rm_group.split("/")
			else:
				rm_class = rm_group
				ERV_type = "N/A"
			with open(old_coord_file, "rb") as infile2:
				for line2 in infile2:
					if line2 != "":
						ERV_id, start, end, ltr5s, ltr5e, ltr3s, ltr3e, pbs, ppt, site5, site3, sbreak, ebreak = line2.split()
						if ERV_id == qseqid:
							fragment_coords = []
							tmp_file = tempfile.NamedTemporaryFile(mode="w", dir=tmp_dir, delete=False)
							tmp_outfile = os.path.join(tmp_dir, tmp_file.name)
							with open(tmp_outfile, "wb") as tmp_out:
								writer = csv.writer(tmp_out, delimiter="\t")
								if sbreak != "N/A" and ebreak != "N/A":
									fragment1 = ERV_id + "-1"
									fragment1_start = start
									fragment1_end = sbreak
									fragment1_tuple = (int(fragment1_start), int(fragment1_end))
									fragment_coords.append(fragment1_tuple)
									writer.writerow([fragment1, rm_class, ERV_type, rm_name, fragment1_start, fragment1_end])
									fragment2 = ERV_id + "-2"
									fragment2_start = ebreak
									fragment2_end = end
									fragment2_tuple = (int(fragment2_start), int(fragment2_end))
									fragment_coords.append(fragment2_tuple)
									writer.writerow([fragment2, rm_class, ERV_type, rm_name, fragment2_start, fragment2_end])
								else:
									fragment_coords = (start, end)
									writer.writerow([ERV_id, rm_class, ERV_type, rm_name, start, end])
								rm_class_set = set()
								ERV_type_set = set()
								rm_name_set = set()
							with open(tmp_outfile, "rb") as tmp_in:
								for line in tmp_in:
									fragment, masker_class, masker_ERV_type, masker_name, fragment_start, fragment_end = line.strip().split()
									rm_class_set.add(str(masker_class))
									ERV_type_set.add(str(masker_ERV_type))
									rm_name_set.add(str(masker_name))

							if len(rm_class_set) > 1: 
								rm_class = ""
								for x in rm_class_set:
									rm_class += "{0}_".format(x)
								rm_class = rm_class.strip("_")
							elif len(rm_class_set) == 1: 
								for x in rm_class_set:
									rm_class = x
							else:
								print "Error: reclassify_ERVs: empty set??"							
							
							if len(ERV_type_set) > 1: 
								ERV_type = ""
								for x in ERV_type_set:
									ERV_type += "{0}_".format(x)
								ERV_type = ERV_type.strip("_")
							elif len(ERV_type_set) == 1: 
								for x in ERV_type_set:
									ERV_type = x
							else:
								print "Error: reclassify_ERVs: empty set??"							
							
							if len(rm_name_set) > 1: 
								rm_name = ""
								for x in rm_name_set:
									rm_name += "{0}_".format(x)
								rm_name = rm_name.strip("_")
							elif len(rm_name_set) == 1: 
								for x in rm_name_set:
									rm_name = x
							else:
								print "Error: reclassify_ERVs: empty set??"
								
							coord_dict[ERV_id] = {}
							coord_dict[ERV_id]["Class"] = rm_class
							coord_dict[ERV_id]["Type"] = ERV_type
							coord_dict[ERV_id]["RM_name"] = rm_name
							coord_dict[ERV_id]["Start"] = start
							coord_dict[ERV_id]["End"] = end
							coord_dict[ERV_id]["Fragments"] = []
							for coords in fragment_coords:
								coord_dict[ERV_id]["Fragments"].append(coords)
							
							with open(classified, "a") as outfile:
								writer = csv.writer(outfile, delimiter="\t")
								writer.writerow([ERV_id, rm_class, ERV_type, rm_name, start, end, ltr5s, ltr5e, ltr3s, ltr3e, pbs, ppt, site5, site3])

	with open(new_coord_file, "wb") as dict_file:
		cPickle.dump(coord_dict, dict_file, 1)

@subdivide(reclassify_ERVs, regex(r"(.+)/ERV_coordinates.classified.dict"), r"\1/By_contig/*.dict")
def split_ERV_dict(infile, outfiles):
	contig_dir = os.path.join(ERV_dir, "By_contig")
	if not os.path.exists(contig_dir):
		os.makedirs(contig_dir)	

	## Clean up previous runs
	for contig_file in os.listdir(contig_dir):
		print "WARNING: split_ERV_dict: Removing contig file: ", contig_file
		os.unlink(os.path.join(contig_dir, contig_file))
	
	ERV_dict = cPickle.load(open(infile[0], "rb"))
	
	genome_dict = {}
	break_coords = None
	for ERV_id, ERV_info in ERV_dict.iteritems():
		contig = ERV_id.rsplit("-", 1)[0]
		if "_" in contig:
			contig = contig.split("_")[0]
		rm_class = ERV_info["Class"]
		rm_fam = ERV_info["Type"]
		rm_name = ERV_info["RM_name"]
		ERV_start = ERV_info["Start"]
		ERV_end = ERV_info["End"]
		if ERV_info["Fragments"]:
			break_coords = ERV_info["Fragments"]
		else:
			break_coords = (int(ERV_start), int(ERV_end))
			
		if not contig in genome_dict.iterkeys():
			genome_dict[contig] = {}
		
		genome_dict[contig][ERV_id] = {}
		genome_dict[contig][ERV_id]["Class"] = rm_class
		genome_dict[contig][ERV_id]["Type"] = rm_fam
		genome_dict[contig][ERV_id]["RM_name"] = rm_name
		genome_dict[contig][ERV_id]["Start"] = ERV_start
		genome_dict[contig][ERV_id]["End"] = ERV_end
		genome_dict[contig][ERV_id]["Fragments"] = []
		for coords in break_coords:
			genome_dict[contig][ERV_id]["Fragments"].append(coords)
	
	for contig, ERV_info in genome_dict.iteritems():
		outname = contig + ".dict"
		outfile = os.path.join(contig_dir, outname)
		with open(outfile, "wb") as dict_file:
			cPickle.dump(ERV_info, dict_file, 1)
			
@follows(filter_RMout)
@transform(split_ERV_dict, regex(r"(.+)/ERV_refinement/By_contig/(.+).dict"), r"\1/By_contig/\2.dict")
def check_rm_files(infile, outfile):
	if not os.path.exists(outfile):
		touch_file(outfile)
			
@follows(filter_RMout, check_rm_files)
@transform(check_rm_files, regex(r"(.+)/By_contig/([a-zA-Z0-9]*)\.dict"), r"\1/By_contig/\2.nested.dict")
def check_nested_files(infile, outfile):
	if not os.path.exists(outfile):
		touch_file(outfile)

@follows(split_ERV_dict, check_rm_files)
@collate([filter_RMout, check_rm_files], regex(r"(.+)/By_contig/([a-zA-Z0-9]*)\.dict"), r"\1/ERV_refinement/By_contig/\2.dict")
def check_ERV_files(infile, outfile):
	if not os.path.exists(outfile):
		touch_file(outfile)
		
@follows(check_nested_files, split_ERV_dict, check_ERV_files)
@mkdir(filter_RMout, regex(r"(.+)/By_contig/([a-zA-Z0-9]*).dict"), r"\1/Merged")
@collate([filter_RMout, check_rm_files], regex(r"(.+)/By_contig/([a-zA-Z0-9]*).dict"), add_inputs([r"\1/By_contig/\2.nested.dict", r"\1/ERV_refinement/By_contig/\2.dict"]), r"\1/Merged/\2.csv")
def map_ERVs_back(infiles, outfile):
	infiles = infiles[0]
	infile1 = infiles[0] ## RMout .dict file
	additional_files = infiles[1] 
	infile2 = additional_files[1] ## ERV dict
	infile3 = additional_files[0] ## nested repeats
	
	if os.stat(infile3).st_size > 0:
		nested_dict = cPickle.load(open(infile3, "rb"))
		for rep_class, values in nested_dict.iteritems():
			for rep_fam, values in values.iteritems():
				for rep_name, values in values.iteritems():
					for rep_id, info in values.iteritems():
						nested_dict[rep_class][rep_fam][rep_name][rep_id]["Coords"] = ""
						nested_dict[rep_class][rep_fam][rep_name][rep_id]["RT_ref"] = "nested"
	else:
		nested_dict = {}
	
	if os.stat(infile1).st_size > 0:
		element_dict = cPickle.load(open(infile1, "rb"))

		contig_quicksect = IntervalTree()

		for rep_class, values in element_dict.iteritems():
			for rep_fam, values in values.iteritems():
				for rep_name, values in values.iteritems():
					for rep_id, info in values.iteritems():
						contig, start, end = rep_id.split("-")
						contig_quicksect.insert_interval(Interval(int(start), int(end),
															chrom = contig,
															strand = info["Strand"],
															value = {"rep_class" : rep_class,
																	"rep_fam" : rep_fam,
																	"rep_name": rep_name}))
						info["Coords"] = ""
						info["RT_ref"] = ""
		
		if os.stat(infile2).st_size > 0:
			ERV_dict = cPickle.load(open(infile2, "rb"))
			filtered_repeat_dict = {}
			
			for ERV_id, ERV_info in ERV_dict.iteritems():
				rm_class = ERV_info["Class"]
				rm_fam = ERV_info["Type"]
				rm_name = ERV_info["RM_name"]
				ERV_start = ERV_info["Start"]
				ERV_end = ERV_info["End"]
				break_coords = ERV_info["Fragments"]
				sorted_list = sorted(break_coords)
				coord_list = zip(sorted_list[::2], sorted_list[1::2])
				TE_dict = merge_overlaps(contig, ERV_id, ERV_start, ERV_end, rm_class, rm_fam, rm_name, coord_list, contig_quicksect)
				filtered_repeat_dict.update(TE_dict)
			
			for rep_class, values in filtered_repeat_dict.iteritems(): 
				for rep_fam, values in values.iteritems():
					for rep_name, values in values.iteritems():
						for rep_id, info in values.iteritems():
							contig = info["Contig"]
							strand = info["Strand"]
							status = info["Status"]
							if status == "remove":
								del element_dict[rep_class][rep_fam][rep_name][rep_id]
							elif status == "keep":
								if not rep_name in element_dict[rep_class][rep_fam].iterkeys():
									element_dict[rep_class][rep_fam][rep_name] = {}
								if not rep_id in element_dict[rep_class][rep_fam][rep_name].iterkeys():
									element_dict[rep_class][rep_fam][rep_name][rep_id] = {}
								element_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
								element_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
								if info["Coords"]:
									element_dict[rep_class][rep_fam][rep_name][rep_id]["Coords"] = info["Coords"]
								else:
									element_dict[rep_class][rep_fam][rep_name][rep_id]["Coords"] = ""
								if info["RT_ref"]:
									element_dict[rep_class][rep_fam][rep_name][rep_id]["RT_ref"] = info["RT_ref"]
								else:
									element_dict[rep_class][rep_fam][rep_name][rep_id]["RT_ref"] = ""
							elif status == "nested":
								if not rep_class in nested_dict.iterkeys():
									nested_dict[rep_class] = {}
								if not rep_fam in nested_dict[rep_class].iterkeys():
									nested_dict[rep_class][rep_fam] = {}
								if not rep_name in nested_dict[rep_class][rep_fam].iterkeys():
									nested_dict[rep_class][rep_fam][rep_name] = {}
								if not rep_id in nested_dict[rep_class][rep_fam][rep_name].iterkeys():
									nested_dict[rep_class][rep_fam][rep_name][rep_id] = {}
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["Contig"] = contig
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["Strand"] = strand
								if info["Coords"]:
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["Coords"] = info["Coords"]
								else:
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["Coords"] = ""
								if info["RT_ref"]:
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["RT_ref"] = info["RT_ref"]
								else:
									nested_dict[rep_class][rep_fam][rep_name][rep_id]["RT_ref"] = "nested"
	else:
		element_dict = cPickle.load(open(infile2, "rb"))
			
	contig_dict = {}						
	
	for rep_class, values in element_dict.iteritems():
		for rep_fam, values in values.iteritems():
			for rep_name, values in values.iteritems():
				for rep_id, info in values.iteritems():
					contig = info["Contig"]
					strand = info["Strand"]
					contig, start, end = rep_id.split("-")
					if info["Coords"] != "":
						for coord_set in info["Coords"]:
							(coord_s, coord_e) = coord_set
							repeat_info = {"Start":coord_s, "End":coord_e, "ID":rep_id, "Class":rep_class, "Family":rep_fam, "Element_name":rep_name, "Strand":strand, "RT_ref": info["RT_ref"]}
							if not contig in contig_dict.iterkeys():
								contig_dict[contig] = {}
							contig_dict[contig][int(coord_s)] = repeat_info
					elif info["RT_ref"] != "":
						repeat_info = {"Start":start, "End":end, "ID":rep_id, "Class":rep_class, "Family":rep_fam, "Element_name":rep_name, "Strand":strand, "RT_ref": info["RT_ref"]}
						if not contig in contig_dict.iterkeys():
							contig_dict[contig] = {}
						contig_dict[contig][int(start)] = repeat_info
					else:
						repeat_info = {"Start":start, "End":end, "ID":rep_id, "Class":rep_class, "Family":rep_fam, "Element_name":rep_name, "Strand":strand, "RT_ref": "N/A"}
						if not contig in contig_dict.iterkeys():
							contig_dict[contig] = {}
						contig_dict[contig][int(start)] = repeat_info
	
	for rep_class, values in nested_dict.iteritems():
		for rep_fam, values in values.iteritems():
			for rep_name, values in values.iteritems():
				for rep_id, info in values.iteritems():
					contig = info["Contig"]
					strand = info["Strand"]
					contig, start, end = rep_id.split("-")
					if info["Coords"] != "":
						for coord_set in info["Coords"]:
							(coord_s, coord_e) = coord_set
							repeat_info = {"Start":coord_s, "End":coord_e, "ID":rep_id, "Class":rep_class, "Family":rep_fam, "Element_name":rep_name, "Strand":strand, "RT_ref": info["RT_ref"]}
							if not contig in contig_dict.iterkeys():
								contig_dict[contig] = {}
							contig_dict[contig][int(coord_s)] = repeat_info
					if info["RT_ref"] != "":
						repeat_info = {"Start":start, "End":end, "ID":rep_id, "Class":rep_class, "Family":rep_fam, "Element_name":rep_name, "Strand":strand, "RT_ref": info["RT_ref"]}
						if not contig in contig_dict.iterkeys():
							contig_dict[contig] = {}
						contig_dict[contig][int(start)] = repeat_info
					else:
						repeat_info = {"Start":start, "End":end, "ID":rep_id, "Class":rep_class, "Family":rep_fam, "Element_name":rep_name, "Strand":strand, "RT_ref": "N/A"}
						if not contig in contig_dict.iterkeys():
							contig_dict[contig] = {}
						contig_dict[contig][int(start)] = repeat_info
	
	for contig, elements in sorted(contig_dict.iteritems()):
		with open(outfile, "wb") as contig_csv:
			writer = csv.writer(contig_csv, delimiter="\t")
			for start_coord, element_info in sorted(elements.iteritems()):
				coord_s = element_info["Start"]
				coord_e = element_info["End"]
				rep_id = element_info["ID"]
				rep_class = element_info["Class"]
				rep_fam = element_info["Family"]
				rep_name = element_info["Element_name"]
				strand = element_info["Strand"]
				rete_ref = element_info["RT_ref"]
				writer.writerow([contig, coord_s, coord_e, rep_id, rep_class, rep_fam, rep_name, strand, rete_ref])

@split(gff_file, "{0}/Gene_annotations/*.gff".format(work_dir))
def split_genes_by_contig(infile, outfiles):
	out_dir = os.path.join(work_dir, "Gene_annotations")
	if os.path.exists(out_dir):
		for contig_file in os.listdir(out_dir):
			print "WARNING: split_genes_by_contig: Removing contig file: ", contig_file
			os.unlink(os.path.join(out_dir, contig_file))
	else:
		os.makedirs(os.path.join(work_dir, "Gene_annotations"))
	
	with open(infile, "r") as gene_file:
		for line in gene_file:
			if not line.startswith("#"):
				gene_contig = line.strip().split()[0]
				
				outname = gene_contig + ".gff"
				outfile = os.path.join(out_dir, outname)
				with open(outfile, "ab") as contig_file:
					contig_file.write(line)
					
@follows(split_genes_by_contig)
@collate([filter_RMout, check_rm_files], regex(r"(.+)/By_contig/([a-zA-Z0-9]*)\.dict"), r"\1/Gene_annotations/\2.gff")
def check_gene_files(infile, outfile):
	if not os.path.exists(outfile):
		touch_file(outfile)
		
@follows(split_genes_by_contig, check_gene_files)
@transform(map_ERVs_back, regex(r"(.+)\/Merged\/(.+)\.csv"), add_inputs(r"\1/Gene_annotations/\2.gff"), r"\1/Merged/\2.gff")
def merge_with_genes(infiles, outfile): 
	infile1 = infiles[0]
	infile2 = infiles[1]
	if os.path.exists(outfile):
		print "WARNING: file exists. Overwriting..."
		os.unlink(outfile)
	
	keep_feature_list = ["tRNA", "mRNA", "gene", "five_prime_UTR", "CDS", "three_prime_UTR"] 
	
	contig_quicksect = IntervalTree()
	
	if os.stat(infile2).st_size > 0:
		with open(infile2, "r") as gene_gff:
			reader = csv.reader(gene_gff, delimiter="\t")
			for row in reader:
				if not row[0].startswith("#"):
					gene_contig = row[0]
					gene_feature_type = row[2]
					gene_start = row[3]
					gene_end = row[4]
					gene_strand = row[6]
					if gene_strand == "+":
						gene_strand = "1"
					elif gene_strand == "-":
						gene_strand = "-1"
					if gene_feature_type in keep_feature_list:
						contig_quicksect.insert_interval(Interval(int(gene_start), int(gene_end), 
															chrom = gene_contig, 
															strand = gene_strand,
															value = {"feature_type":gene_feature_type}))

		id_list = []
		with open(infile1, "r") as repeat_file, open(outfile, "a") as out_gff:
#			out_gff.write("##gff-version 3\n")
			writer = csv.writer(out_gff, delimiter="\t")
			for line in repeat_file:
				contig, coord_s, coord_e, rep_id, rep_class, rep_fam, rep_name, strand, info = line.strip().split()
				keep_status, region_type = merge_genes(contig, coord_s, coord_e, contig_quicksect) 
				if strand == "1":
					strand = "+"
				elif strand == "-1":
					strand = "-"
				id_list.append(rep_id)
				id_count = id_list.count(rep_id)
				gff_id = rep_id + "." + str(id_count)
				if info == "N/A":
					attributes = "ID={0};Name={1};Parent={2};Note={3}/{4}, repeat_type={5}".format(gff_id, rep_name, rep_id, rep_class, rep_fam, region_type)
				elif info == "nested":
					attributes = "ID={0};Name={1};Parent={2};Note={3}/{4}, repeat_type={5}, nested=True".format(gff_id, rep_name, rep_id, rep_class, rep_fam, region_type)
				else:
					attributes = "ID={0};Name={1};Parent={2};Note={3}/{4}, repeat_type={5}, RetroTector_reference={6}".format(gff_id, rep_name, rep_id, rep_class, rep_fam, region_type, info)
				if keep_status == True:
					writer.writerow([contig, "custom", "repeat_region", coord_s, coord_e, ".", strand, ".", attributes])
	else:
		id_list = []
		with open(infile1, "r") as repeat_file, open(outfile, "a") as out_gff:
			writer = csv.writer(out_gff, delimiter="\t")
			for line in repeat_file:
				contig, coord_s, coord_e, rep_id, rep_class, rep_fam, rep_name, strand, info = line.strip().split()
				region_type = "between_genes"
				if strand == "1":
					strand = "+"
				elif strand == "-1":
					strand = "-"
				id_list.append(rep_id)
				id_count = id_list.count(rep_id)
				gff_id = rep_id + "." + str(id_count)
				if info == "N/A":
					attributes = "ID={0};Name={1};Parent={2};Note={3}/{4}, repeat_type={5}".format(gff_id, rep_name, rep_id, rep_class, rep_fam, region_type)
				elif info == "nested":
					attributes = "ID={0};Name={1};Parent={2};Note={3}/{4}, repeat_type={5}, nested=True".format(gff_id, rep_name, rep_id, rep_class, rep_fam, region_type)
				else:
					attributes = "ID={0};Name={1};Parent={2};Note={3}/{4}, repeat_type={5}, RetroTector_reference={6}".format(gff_id, rep_name, rep_id, rep_class, rep_fam, region_type, info)
				writer.writerow([contig, "custom", "repeat_region", coord_s, coord_e, ".", strand, ".", attributes])
			
## List for gff file: "tRNA", "gene", "mRNA", "exon", "five_prime_UTR", "CDS", "expressed_sequence_match", "protein_match", "three_prime_UTR"
## Use "repeat_region" in final gft
## gff format: "scaffold", "source", "feature_type", "start", "end", "score", "strand", "phase", "attributes"
## attributes = repeat class, family, name, id. See http://gmod.org/wiki/GFF3 for terms


@split("{0}/Merged".format(work_dir), "{0}/Merged/*.gff".format(work_dir))
def dummy_job(infile, outfiles):
	pass

@mkdir(dummy_job, regex(r"(.+)\/Merged\/(.+)\.gff"), r"\1/Final")
@transform(dummy_job, regex(r"(.+)\/Merged\/(.+)\.gff"), r"\1/Final/\2.gff")	
def check_final_overlaps(infile, outfile):
	contig_quicksect = IntervalTree()
	contig_dict = {}
	with open(infile, "r") as raw_gff:
		for line in raw_gff:
			scaffold, input, input_type, start, end, dot1, strand, dot2, info = line.strip().split("\t")
			id_field = info.split(";")[0]
			repeat_id = id_field.split("=", 1)[1]	
			if strand == "+":
				num_strand = 1
			elif strand == "-":
				num_strand = -1
			else:
				print "[WARNING]: check_final_overlaps: unknown strand designation. Assuming + strand"
				num_strand = 1
			contig_quicksect.insert_interval(Interval(int(start), int(end), chrom = scaffold, strand = num_strand, value = {"repeat_id" : repeat_id}))
			contig_dict[int(start)] = {"ID":repeat_id, "Start":start, "End":end, "info":info}
	
	with open(infile, "r") as repeat_calls, open(outfile, "w") as filtered:
		writer = csv.writer(filtered, delimiter="\t")
		for line in repeat_calls:
			scaffold, input, input_type, start, end, dot1, strand, dot2, info = line.strip().split("\t")
			if "nested" in info:
				##check for overlaps, update info
				overlaps = contig_quicksect.find(int(start), int(end))
				if len(overlaps) == 1: ## not nested after filtering
					info = info.rsplit(",", 1)[0]
				elif len(overlaps) < 1:
					print "[ERROR]: check_final_overlaps: repeat interval does not exist?!", line
			writer.writerow([scaffold, input, input_type, start, end, dot1, strand, dot2, info])

@collate(check_final_overlaps, regex(r"(.+)\/Final\/(.+)\.gff"), r"\1/Repeat_calls.gff")
def summary_gff(infiles, outfile):
	with open(outfile, "w") as gff_file:
		gff_file.write("##gff-version 3\n")
		for contig_file in infiles:
			with open(contig_file, "r") as contig_gff:
				for line in contig_gff:
					gff_file.write(line)
			
@subdivide(check_final_overlaps, regex(r"(.+)\/Final\/(.+)\.gff"), [r"\1/Final/\2.rtype.dict", r"\1/Final/\2.location.dict"])
def summary_repeats_bycontig(infile, outfiles):
	type_summary = {}
	location_summary = {}
	
	non_nucleotides = ["N", "n", "X", "x"]
	
	with open(infile, "rb") as repeats:
		reader = csv.reader(repeats, delimiter="\t")
		line_count = 0
		for row in reader:
			line_count += 1
			contig = row[0]
			coord_s = row[3]
			coord_e = row[4]
			strand = row[6]
			attributes = row[8]
			gff_notes = attributes.split(";")[-1]
			nested = False
			RT_ref = None
			if len(gff_notes.split(",")) == 2:
				higher_classifications, repeat_type = gff_notes.split(",")
			elif len(gff_notes.split(",")) == 3:
				higher_classifications, repeat_type, other = gff_notes.split(",")
				if "nested" in other:
					nested = True
				elif "RetroTector" in other:
					RT_ref = other.split("=")[1]
			if nested == False:
				repeat_classifiction = higher_classifications.split("=")[1]
				rep_class, rep_fam = repeat_classifiction.split("/", 1)
				repeat_location = repeat_type.split("=")[1]
				
				sequence_file_name = contig + ".fas"
				sequence_file = os.path.join(contig_dir, sequence_file_name)
				repeat_sequence = slice_sequence(sequence_file, coord_s, coord_e)
				for non_base in non_nucleotides:
					if non_base in repeat_sequence:
						clean_sequence = repeat_sequence.replace(non_base, "")
						repeat_sequence = clean_sequence
				repeat_length = len(repeat_sequence)
				if rep_class not in type_summary.iterkeys():
					type_summary[rep_class] = {}
				if rep_fam not in type_summary[rep_class].iterkeys():
					type_summary[rep_class][rep_fam] = {"count":0, "length":0}
				type_summary[rep_class][rep_fam]["count"] += 1
				type_summary[rep_class][rep_fam]["length"] += repeat_length
				
				if rep_class not in location_summary.iterkeys():
					location_summary[rep_class] = {}
				if rep_fam not in location_summary[rep_class].iterkeys():
					location_summary[rep_class][rep_fam] = {"between_genes":0, "UTR":0, "intron":0, "unknown":0}
				location_summary[rep_class][rep_fam][repeat_location] += 1
		print line_count
	
	with open(outfiles[0], "wb") as type_file:
		cPickle.dump(type_summary , type_file, 1)
				
	with open(outfiles[1], "wb") as location_file:
		cPickle.dump(location_summary , location_file, 1)

@collate(summary_repeats_bycontig, regex(r"(.+)\/Final\/(.+)\.rtype\.dict"), r"\1/Repeat_fractions.out")
def summary_repeats(infiles, outfile):
	genome_length = 0
	
	with open(genome_file, "r") as genome:
		non_nucleotides = ["N", "n", "X", "x"]
		for line in genome:
			if not line.startswith(">"):
				for non_base in non_nucleotides:
					if non_base in line:
						nucleotides = line.strip().replace(non_base, "")
						line = nucleotides
				line_length = len(line)
				genome_length += int(line_length)
	
	repeat_dict = {}
	
	for dict_file in infiles:
		contig_summary = cPickle.load(open(dict_file, "rb"))
		for rep_class, values in sorted(contig_summary.iteritems()):
			if rep_class not in repeat_dict.iterkeys():
				repeat_dict[rep_class] = {}
			for rep_fam, attributes in sorted(values.iteritems()):
				if rep_fam not in repeat_dict[rep_class].iterkeys():
					repeat_dict[rep_class][rep_fam] = {"count":0, "length":0}
				repeat_dict[rep_class][rep_fam]["count"] += int(attributes["count"])
				repeat_dict[rep_class][rep_fam]["length"] += int(attributes["length"])
					
	with open(outfile, "w") as summary_file:
		summary_file.write("Repeat_class\tRepeat_family\tCount\t%_genome\n")
		for rep_class, values in sorted(repeat_dict.iteritems()):
			for rep_fam, attributes in sorted(values.iteritems()):
				genome_percentage = int(attributes["length"])/genome_length
				summary_file.write("{0}\t{1}\t{2}\t{3}\n".format(rep_class, rep_fam, attributes["count"], genome_percentage))

@collate(summary_repeats_bycontig, regex(r"(.+)\/Final\/(.+)\.location\.dict"), r"\1/Repeat_locations.out")
def summary_locations(infiles, outfile):

	repeat_dict = {}
	
	for dict_file in infiles:
		contig_summary = cPickle.load(open(dict_file, "rb"))
		for rep_class, values in sorted(contig_summary.iteritems()):
			if rep_class not in repeat_dict.iterkeys():
				repeat_dict[rep_class] = {}
			for rep_fam, locations in sorted(values.iteritems()):
				if rep_fam not in repeat_dict[rep_class].iterkeys():
					repeat_dict[rep_class][rep_fam] = {}
				for genomic_location, count in locations.iteritems():
					if genomic_location not in repeat_dict[rep_class][rep_fam].iterkeys():
						repeat_dict[rep_class][rep_fam][genomic_location] = 0
					repeat_dict[rep_class][rep_fam][genomic_location] += int(count)
					
	with open(outfile, "w") as summary_file:
		summary_file.write("Repeat_class\tRepeat_family\tIntergenic\tUTR\tIntronic\tUnknown\n")
		for rep_class, values in sorted(repeat_dict.iteritems()):
			for rep_fam, locations in sorted(values.iteritems()):
				summary_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(rep_class, rep_fam, locations["between_genes"], locations["UTR"], locations["intron"], locations["unknown"]))
		
@follows(summary_gff, summary_repeats, summary_locations)
def clean_up():
	for fname in os.listdir(tmp_dir):
		tmp_file = os.path.join(tmp_dir, fname)
		os.unlink(tmp_file)
	os.rmdir(tmp_dir)
	end = time.time()
	runtime = end - start
	print timedelta(seconds=runtime)
	print "DONE"

###############################################################################
# PIPELINE RUN OPTIONS
###############################################################################

if __name__	 == "__main__":
	if args.no_run:
		pipeline_printout(sys.stdout, [clean_up], forcedtorun_tasks=[], verbose = 5, verbose_abbreviated_path = 0)
	if args.graph:  ## DOES NOT WORK ON HPC RUN IN INTERACTIVE MODE AND COPY OUTPUT INTO NEW DOC ON LOCAL
		pipeline_printout_graph(os.path.join(work_dir, os.path.basename(__file__).rsplit('.', 1)[0] + ".svg"), "svg", [clean_up])
	if args.run:
#		pipeline_run([clean_up], forcedtorun_tasks=[map_ERVs_back], verbose=5, multiprocess=num_proc) ## FOR FUNCTION TESTING
#		pipeline_run([check_rm_files, check_nested_files, check_ERV_files], forcedtorun_tasks=[check_rm_files, check_nested_files, check_ERV_files], verbose=2, multiprocess=num_proc) ## PICK-UP BROKEN RUNS
#		pipeline_run([merge_with_genes], touch_files_only=True, verbose=2, multiprocess=num_proc) ## MARK TASKS UP TO THIS POINT AS UP-TO-DATE
		pipeline_run([clean_up], verbose=2, multiprocess=num_proc)
