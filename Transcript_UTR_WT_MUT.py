from ruffus import *
import glob

from CGAT import GTF
from CGAT import Genomics
import pandas
from TranscriptCoordInterconverter import TranscriptCoordInterconverter
import os
import sys
import re
from CGAT import IndexedFasta
import CGATCore.IOTools as IOTools
import CGATCore.Experiment as E


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

	# setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])
    
    parser.add_option("-f", "--input-format", dest="format", type="choice",
                      choices = ["bed", "vcf"],
                      default = "bed",
                      help="Format of varients input file, default is output from")

    parser.add_option("-a", "--gtf_file", dest="gtf", type="string",
      			help="Ensembl gzip gtf file")

    parser.add_option("-g", "--genome-fasta", dest="fasta", type="string",
                        help="Location of indexed fasta sequence for genome")

#    parser.add_option("-f", "--input-format", dest="format", type="choice",
#                      choices = ["bed", "vcf"],
#                      default = "bed",
#                      help="Format of varients input file, default is output from"
#						"intersection quasar input and bed")

#    parser.add_option("-w", "--Wild_type", dest="wt_fasta", action="store_true",
#                      default="fasta",

#                      help="Output Wild type sequence per transcript")
    
#   parser.add_option("-m", "--mutant_sequence", dest="mut_fasta", action="store_true",
#                       default="fasta",
#                       help="Output mutant sequence per transcript")

	# add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)
	
    SNPs_file = pandas.read_csv(options.stdin,sep="\t",
				header=None)

    SNPs_file.columns=["chr", "start", "end", "rsID", "ref", "alt","beta.se","combined_pvalue","Pval","avg_beta","n","BH"]

    #print (SNPs_file)								
    SNPs_file = SNPs_file.set_index(["chr","start"]).sort_index()
	
    gtf_file = IOTools.open_file(options.gtf)
    gtf_iterator = GTF.iterator(gtf_file)

    genome = IndexedFasta.IndexedFasta(options.fasta)
   
   
    stdout1 = options.stdout.name + ".WT.fasta"
    stdout2 = options.stdout.name + ".Mutant.fasta"


    stdout1 = IOTools.open_file(stdout1, "w")  
    stdout2 = IOTools.open_file(stdout2, "w")
   
	
    for transcript in GTF.transcript_iterator(gtf_iterator):

        utr = [e for e in transcript if e.feature == "three_prime_utr"]

        if len(utr) == 0:
            continue
        overlapping_snps = pandas.DataFrame()

        for e in utr:
            e.feature = "exon"
            this_chrom = SNPs_file.loc[e.contig]
            snps = this_chrom.loc[e.start:e.end]
            overlapping_snps = pandas.concat([overlapping_snps, snps])

        if overlapping_snps.shape[0] == 0:
            continue

# get the wt sequence
        wt_sequence = []
        utr = sorted(utr, key=lambda x: x.start)

        for exon in utr: 
             wt_sequence.append(genome.getSequence(contig=exon.contig, 
                                strand="+",
                                start=exon.start,
                                end=exon.end))
        
        if utr[0].strand == "-":
            wt_sequence = [Genomics.reverse_complement(e) for e in wt_sequence[::-1]]
        wt_sequence = "".join(wt_sequence)
        wt_sequence = wt_sequence.upper()
        converter = TranscriptCoordInterconverter(utr)
     #   print (">" + str((utr[0].transcript_id))+ "\n"+ wt_sequence + "\n")
        #print (wt_sequence +"\n")   
# output wildtype
            
        stdout1.write(">" +  str((utr[0].transcript_id)) + "\n" + str(wt_sequence) + "\n")
        #options.stdout1.write(str(wt_sequence) +"\n")

        for _, snp in overlapping_snps.iterrows():
            snp_position = int(converter.genome2transcript([snp["end"]-1])[0])
            if utr[0].strand == "-":
                ref = Genomics.reverse_complement(snp["ref"])
                alt = Genomics.reverse_complement(snp["alt"])
            else:
                ref = snp["ref"]
                alt = snp["alt"]
            assert wt_sequence[snp_position] == ref , "position %s in %s should be %s, but is %s" % (snp_position, utr[0].transcript_id, ref, wt_sequence[snp_position])

            id = [utr[0].transcript_id,
                 snp["rsID"],
                 snp["ref"],
                 snp["alt"],
                 str(snp_position),
                 str(snp["end"]-1)]
			
            mut_seq = list(wt_sequence)
            mut_seq[snp_position] = alt
            mut_seq = "".join(mut_seq)
        
            stdout2.write(">" + "_".join(id) + "\n" + str(mut_seq + "\n"))
            #options.stdout2.write(mut_seq +"\n")        
            E.stop()



if __name__ == "__main__":
    sys.exit(main(sys.argv))

