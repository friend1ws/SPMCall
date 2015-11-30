#! /usr/bin/env python

import sys, os, subprocess, shutil
import ConfigParser

import utils
import annotate
import associate 

def main(args):

    mutation_file = args.mutation_file
    junction_file = args.junction_file
    output_prefix = args.output_prefix
    # annotation_dir = args.annotation_dir
    is_anno = True if args.f == "anno" else False
    control_file = args.ctrl

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
                             control_file,
                             parser.getint("star_junction_filt", "read_num_thres"), 
                             parser.getint("star_junction_filt", "overhang_thres"),
                             parser.getboolean("star_junction_filt", "remove_annotated"))


    annotate.annot_junction(output_prefix + ".mutran_tmp.junction.txt",
                              output_prefix + ".mutran_tmp.junction.annot.txt",
                              parser.get("annotation", "annotation_dir"))

    ##########

    ##########
    # associate mutation and junction
    associate.get_snv_junction(output_prefix + ".mutran_tmp.junction.annot.txt",
                               output_prefix + ".splicing_mutation.txt",
                               output_prefix + ".mutran_tmp.vcf.gz",
                               parser.get("annotation", "annotation_dir"))
    ##########

    if parser.getboolean("debug", "debug_mode") != True:
        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.unsorted.vcf"])
        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.vcf.gz"])
        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.vcf.gz.tbi"])
        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.junction.txt"])
        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.junction.annot.txt"])

