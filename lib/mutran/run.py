#! /usr/bin/env python

import sys, os, subprocess, shutil
import ConfigParser

import utils

def main(args):

    mutation_file = args.mutation_file
    junction_file = args.junction_file
    output_prefix = args.output_prefix
    # annotation_dir = args.annotation_dir
    is_anno = True if args.f == "anno" else False

    parser = ConfigParser.SafeConfigParser()
    parser.read(args.param)

    reference_genome = parser.get("alignment", "reference_genome")

    output_prefix_dir = os.path.dirname(output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    ##########
    # processing mutation file

    # convert mutation data to vcf (if --anno is on)
    if is_anno == True:
        utils.convert_anno2vcf(mutation_file, output_prefix + ".mutran_tmp.unsorted.vcf", reference_genome)
    else:
        utils.remove_vcf_header(mutation_file, output_prefix + ".mutran_tmp.unsorted.vcf")

    hout = open(output_prefix + ".mutran_tmp.vcf", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", output_prefix + ".mutran_tmp.unsorted.vcf"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in sorting vcf file"
        sys.exit(1)


    s_ret = subprocess.call([parser.get("software", "bgzip"), "-f", output_prefix + ".mutran_tmp.vcf"])
    if s_ret != 0:
        print >> sys.stderr, "Error in bgzip compression"
        sys.exit(1)


    s_ret = subprocess.call([parser.get("software", "tabix"), "-p", "vcf", output_prefix + ".mutran_tmp.vcf.gz"])
    if s_ret != 0:
        print >> sys.stderr, "Error in tabix indexing"
        sys.exit(1)
    ##########


    ##########
    # processing splicing junction file
    utils.proc_star_junction(junction_file, output_prefix + ".mutran_tmp.junction.txt", 
                             parser.getint("star_junction_filt", "read_num_thres"), 
                             parser.getint("star_junction_filt", "overhang_thres"),
                             parser.getboolean("star_junction_filt", "remove_annotated"))

    ##########
     
    """
    s_ret = subprocess.call(["bgzip", "-f", output_prefix + ".vcf"])

    # compress and add index

    utils.filterImproper(input_bam, output_prefix + ".filt.bam", mapq_thres)


    hOUT = open(output_prefix + ".filt.bed12", 'w')
    s_ret = subprocess.call([bedtools_path + "/bedtools", "bamtobed", "-bed12", "-i", output_prefix + ".filt.bam"], stdout = hOUT)
    hOUT.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating filt.bed12"
        sys.exit(1)


    hOUT = open(output_prefix + ".edge.bed", 'w')
    s_ret = subprocess.call([bedtools_path + "/bedtools", "intersect", "-a", output_prefix + ".filt.bed12", "-b", annotation_dir + "/edge.bed", 
                     "-split", "-wao"], stdout = hOUT)
    hOUT.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating edge.bed"
        sys.exit(1)

    hOUT = open(output_prefix + ".edge_broaden.bed", 'w')
    s_ret = subprocess.call([bedtools_path + "/bedtools", "intersect", "-a", output_prefix + ".filt.bed12", "-b", annotation_dir + "/edge_broaden.bed", 
                     "-split", "-wao"], stdout = hOUT)
    hOUT.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating edge_broaden.bed"
        sys.exit(1)

    utils.summarize_edge(output_prefix + ".edge.bed", output_prefix + ".edge_broaden.bed", output_prefix + ".intron.bed", 5)

    subprocess.call(["rm", "-rf", output_prefix + ".filt.bam"])
    subprocess.call(["rm", "-rf", output_prefix + ".filt.bed12"])
    subprocess.call(["rm", "-rf", output_prefix + ".exon.bed"])
    subprocess.call(["rm", "-rf", output_prefix + ".exon2base.txt"])
    """

