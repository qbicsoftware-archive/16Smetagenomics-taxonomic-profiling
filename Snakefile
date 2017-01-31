import os
import sys
import subprocess
import json
import re
from os.path import join as pjoin
from os.path import exists as pexists
import glob
import numpy as np
import pandas as pd
from itertools import product

configfile: "config.json"
workdir: config["var"]
SNAKEDIR = config['src']

try:
    VERSION = subprocess.check_output(
        ['git', 'describe', '--tags', '--always', '--dirty'],
        cwd=SNAKEDIR
    ).decode().strip()
except subprocess.CalledProcessError:
    VERSION = 'unknown'

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
ETC = config['etc']

def data(path):
    return os.path.join(DATA, path)

def ref(path):
    return os.path.join(REF, path)

def log(path):
    return os.path.join(LOGS, path)

def result(path):
    return os.path.join(RESULT, path)

def etc(path):
    return os.path.join(ETC, path)

try:
    with open(etc("params.json")) as f:
        parameters = json.load(f)
except OSError as e:
    print("Could not read parameter file: " + str(e), file=sys.stderr)
    sys.exit(1)
except ValueError as e:
    print("Invalid parameter file: " + str(e), file=sys.stderr)
    sys.exit(1)

filename = os.listdir(data(""))[0]


for ending in ['fastq.gz', 'fastq', 'fq', 'fq.gz']:
    if filename.endswith(ending):
        file_ending = ending
        for p in ['_1', '_R1', '.1']:
            if p in filename:
                match = re.search('\%s(.*)%s' % (p, "\.{}".format(ending)), filename)
                seperator = ''
                if match:
                   match2 = re.search(r'[\_\.]', match.group(1))
                   seperator = match2.group(0) if match2 else seperator
                pair_id = p
                SAMPLES, LANES= glob_wildcards(data("{smp}" + p + "{lane}" + ending))


def filter_combinator(combinator, pair_id,  greenlist):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            wc_comb_dict = {}
            wc_comb_dict = dict(wc_comb)
            filename = wc_comb_dict['smp'] + pair_id + wc_comb_dict['lane']
            if any(filename in x for x in greenlist):
                samp = (('smp',wc_comb_dict['smp']),('lane', wc_comb_dict['lane']))
                yield samp
    return filtered_combinator

filter_for = os.listdir(data(""))
filtered_product = filter_combinator(product, pair_id, filter_for)

pair_id_2 = pair_id.replace("1","2")
rule all:
    input: [x[:-1] for x in expand(result("{smp}{lane}"), filtered_product, smp=set(SAMPLES), lane=set(LANES))]


rule clipAndMerge:
    input:
        L = data("{smp}" + pair_id + seperator + "{lane}." + file_ending) if seperator != '' else data("{smp}{lane}" + pair_id + "." + file_ending),
        R = data("{smp}" + pair_id_2 + seperator + "{lane}." + file_ending) if seperator != '' else data("{smp}{lane}" + pair_id_2 + "." + file_ending)
    log:
        log("{smp}{lane}_log.txt")
    output:
        "ClipAndMerge/{smp}" + seperator + "{lane}." + file_ending
    run:
        shell("java -jar " + os.path.join(os.environ['CLIPANDMERGE_BIN_DIR'], 'ClipAndMerge-1.7.5.jar') + " -in1 " + input['L'] + " -in2 " + input['R'] + " -o {output} -log {log}")


rule runMalt:
    input:
        "ClipAndMerge/{smp}" + seperator + "{lane}." + file_ending
    log:
        log("{smp}{lane}_log.txt")
    output:
        result("{smp}" + seperator + "{lane}")
    run:
        shell("mkdir {output}")
        shell("malt-run -m BlastN -at SemiGlobal -t 64 -wlca -mq 25 -d /lustre_cfc/qbic/reference_genomes/16SMicrobial -o {output} -i {input} >> {log}")
