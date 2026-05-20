#!/usr/bin/env python

"""
run_scDeepLUCIA.py
------------------------------
Call loops by scDeepLUCIA.

Outputs (per chromosome, per resolution window):
  - scDeepLUCIA_result/sample/chrom/score.xls.gz   : the scores with columns "chrom", "index_one", "index_two", "prob"
  - scDeepLUCIA_result/sample/chrom/loop.bedpe   : called loops in bedpe format

Usage:
    python run_scDeepLUCIA.py \
        sample_name chromosome_name

Example:
  cd ~/work/ ## all required files saved in work, including the souce of scDeepLUCIA, the scDeepLUCIA_inputs, and scDeepLUCIA.sif container
  mkdir tmp
  export APPTAINER_TMPDIR=${PWD}/tmp
  export TMPDIR=${PWD}/tmp
  apptainer shell -B ${PWD} --nv scDeepLUCIA.sif
  cd ~/work/
  run_scDeepLUCIA.py sample1 chr1
"""

import sys
import gzip
import itertools
import os

from pathlib import Path

import numpy
import pandas 
import tensorflow 

import matplotlib.pyplot as plt
from sklearn import metrics

#pwd=os.getcwd()
pwd=Path(__name__).parent.resolve()
sys.path.append(os.path.join(pwd, 'scDeepLUCIA/deeplucia_toolkit/'))
import make_dataset
import misc

resolution = 5000

keras_model_filename = os.path.join(pwd, "scDeepLUCIA/model/trained.h5")
sample = sys.argv[1] #"X0"
chrom = sys.argv[2]#"chr10"

scan_start = 0     # scan_start = 1000 # scan_start = 1000
scan_end   = 26137 # scan_end   = 2000 # scan_end   = 1020

marker_type = "r2_030M"
genome_version = "mm10" 

feature_dirname =  Path(os.path.join(pwd, "scDeepLUCIA_inputs"))
result_dirname =  Path(os.path.join(pwd, "scDeepLUCIA_result", sample, chrom))
result_dirname.mkdir(parents=True,exist_ok=True)

prediction_filename = result_dirname / "score.xls.gz"
loop_bedpe_filename = result_dirname / "loop.bedpe"

model = misc.load_model(keras_model_filename)

seq_array_dirname,epi_array_dirname,con_array_dirname = misc.get_directory(feature_dirname,chrom,sample,genome_version)
chrom_sample = (chrom,sample)
chrom_sample_list = [chrom_sample]

chrom_to_seq_array = make_dataset.load_seq_array_dir(chrom_sample_list, seq_array_dirname)
chrom_sample_to_epi_array = make_dataset.load_epi_array_dir(chrom_sample_list , marker_type , epi_array_dirname)
chrom_sample_to_con_array = make_dataset.load_con_array_dir(chrom_sample_list , con_array_dirname)

con_array = chrom_sample_to_con_array[chrom_sample]
scan_end = con_array.shape[0] - 2
scanning_loop_candidate_gen = make_dataset.gen_scanning_loop_candidate(chrom,sample,scan_start,scan_end)

pair_list = []
prob_list = []
# 32768 = 2^15
for _,chunk in itertools.groupby(enumerate(scanning_loop_candidate_gen) , lambda x : x[0]//32768):
    loop_candidate_list = []
    for _,loop_candidate in chunk:
        pair = loop_candidate[2]
        loop_candidate_list.append(loop_candidate)
        pair_list.append(pair)
    batched_feature,_ = make_dataset.extract_seq_epi_dataset_nonshuffle(loop_candidate_list, chrom_to_seq_array, chrom_sample_to_epi_array)
    output = model.predict(batched_feature)
    batched_prob_pred = numpy.squeeze(output,axis=1)
    for prob in batched_prob_pred:
        prob_list.append(prob)

if len(pair_list) == len(prob_list):
    with gzip.open(prediction_filename,"wt") as prediction_file:
        prediction_file.write("chrom\tindex_one\tindex_two\tprob\n")
        for pair,prob in zip(pair_list,prob_list):
            prediction_file.write("\t".join(map(str,[chrom, pair[0],pair[1], prob])) + "\n")

loop_df = pandas.read_table(prediction_filename)
dist_filtered_loop_df = misc.filter_by_distance(loop_df)
cutoff_filtered_loop_df = misc.filter_by_quantile(dist_filtered_loop_df,con_array)
clustered_loop_df = misc.form_loop_cluster(cutoff_filtered_loop_df)
misc.save_as_bedpe(clustered_loop_df,chrom,resolution,loop_bedpe_filename)

