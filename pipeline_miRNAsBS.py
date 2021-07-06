from ruffus import *

import sys
import os
import sqlite3
import shutil
import CGATCore.Experiment as E
from CGATCore import Pipeline as P
import re
import glob
import collections
import CGAT.GTF as GTF
import CGATCore.IOTools as IOTools
import CGAT.BamTools.bamtools as BamTools
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC



# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS["projectsrc"] = os.path.dirname(__file__)
#for key, value in PARAMS.iteritems():
#    print "%s:\t%s" % (key,value)


PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh
######################################################
# Extract 3' UTR from gtf file and from QUASAR output
######################################################

@follows(mkdir("Threeprime_gtf.dir"))
@transform("Quasar_stats_results/*.results",
           regex("Quasar_stats_results/(.+).results"),
           add_inputs(PARAMS["gtf"]),
           r"Threeprime_gtf.dir/\1_three_prime_transcriptome.gtf.gz")

def threeprimegtf(infiles, outfile):

    Quasar, gtf = infiles
    statement = ''' zcat %(gtf)s | awk '$3=="three_prime_utr"' | 
                    bedtools intersect -u -a - -b %(Quasar)s | 
                    bedtools intersect -u -a %(gtf)s -b - | 
                    awk '$3=="transcript"'   | 
                    bedtools    intersect -u -a %(gtf)s -b - |  gzip > %(outfile)s '''
    job_memory = "8G"
    P.run(statement)

#######################################################

@follows(threeprimegtf,mkdir("Transcriptome_sequence.dir"))

@transform("Quasar_stats_results/*.results",
           regex("Quasar_stats_results/(.+).results"),
           add_inputs(r"Threeprime_gtf.dir/\1_three_prime_transcriptome.gtf.gz"),
           [r"Transcriptome_sequence.dir/\1_threeprime.WT.fasta",r"Transcriptome_sequence.dir/\1_threeprime.Mutant.fasta"])

def Sequence(infiles, outfiles):

    Quasar,gtf = infiles
    genome = PARAMS["mapfasta"]
    WT,mutant = outfiles
    out = P.snip(WT,".WT.fasta")
 
    statement = ''' python /data/md1srd/Software/utr_sequences/Transcript_UTR_WT_MUT.py -g %(genome)s -a %(gtf)s --stdin %(Quasar)s --stdout %(out)s'''
    job_memory = "8G"
    P.run(statement)

##############################################################

@follows(Sequence,mkdir("miRNA_BS.dir"))

@transform("Transcriptome_sequence.dir/*.fasta",
           regex("Transcriptome_sequence.dir/(.+).fasta"),
           add_inputs(PARAMS["miRNA_db"]), r"miRNA_BS.dir/\1.hits")

def miRNAScan(infiles, outfile):
    Transcriptome,miRNAs = infiles
    statement = '''miranda %(miRNAs)s %(Transcriptome)s -en -20 -strict -out %(outfile)s &&
                  grep -A 1 "Scores for this hit:" %(outfile)s | sort | grep '>' >  %(outfile)s.hit '''

    job_memory = "8G"
    P.run(statement)




@follows(threeprimegtf,Sequence,miRNAScan)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))




