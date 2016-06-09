#! /usr/bin/env python

import sys, os, subprocess, shutil
# import ConfigParser

import utils
import junc_utils.utils
import junc_utils.annotate
import associate 

def main(args):

    mutation_file = args.mutation_file
    junction_file = args.junction_file
    output_prefix = args.output_prefix
    # annotation_dir = args.annotation_dir
    is_anno = True if args.f == "anno" else False
    control_file = args.ctrl
    is_sv = True if args.sv else False
    is_debug = True if args.debug else False 

    # parser = ConfigParser.SafeConfigParser()
    # parser.read(args.param)

    reference_genome = args.reference_genome

    output_prefix_dir = os.path.dirname(output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    ##########
    # processing mutation file

    if not is_sv: 

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


        s_ret = subprocess.call(["bgzip", "-f", output_prefix + ".mutran_tmp.vcf"])
        if s_ret != 0:
            print >> sys.stderr, "Error in bgzip compression"
            sys.exit(1)


        s_ret = subprocess.call(["tabix", "-p", "vcf", output_prefix + ".mutran_tmp.vcf.gz"])
        if s_ret != 0:
            print >> sys.stderr, "Error in tabix indexing"
            sys.exit(1)

    else:
        
        utils.convert_genosv2bed(mutation_file, output_prefix + ".mutran_tmp.unsorted.bedpe")

        hout = open(output_prefix + ".mutran_tmp.bedpe", 'w') 
        s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", output_prefix + ".mutran_tmp.unsorted.bedpe"], stdout = hout)
        hout.close()
        
        if s_ret != 0:
            print >> sys.stderr, "Error in sorting bedpe file"
            sys.exit(1)
    
        
        s_ret = subprocess.call(["bgzip", "-f", output_prefix + ".mutran_tmp.bedpe"])
        if s_ret != 0:
            print >> sys.stderr, "Error in bgzip compression"
            sys.exit(1)
    
        
        s_ret = subprocess.call(["tabix", "-p", "bed", output_prefix + ".mutran_tmp.bedpe.gz"])
        if s_ret != 0:
            print >> sys.stderr, "Error in tabix indexing"
            sys.exit(1)

    ##########

    ##########
    # processing splicing junction file
    junc_utils.utils.proc_star_junction(junction_file, output_prefix + ".mutran_tmp.junction.txt", 
                             control_file, args.read_num_thres, args.overhang_thres, not args.keep_annotated, False)


    junc_utils.annotate.annot_junction(output_prefix + ".mutran_tmp.junction.txt",
                              output_prefix + ".mutran_tmp.junction.annot.txt",
                              args.annotation_dir, 3, 10)

    ##########

    ##########
    # associate mutation and junction
    if not is_sv:
        associate.get_snv_junction(output_prefix + ".mutran_tmp.junction.annot.txt",
                                   output_prefix + ".splicing_mutation.txt",
                                   output_prefix + ".mutran_tmp.vcf.gz",
                                   args.annotation_dir)
    else:
        associate.get_sv_junction(output_prefix + ".mutran_tmp.junction.annot.txt",
                                  output_prefix + ".splicing_sv.txt",
                                  output_prefix + ".mutran_tmp.bedpe.gz",               
                                  args.annotation_dir)
    ##########

    if is_debug != True:
        if not is_sv:
            subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.unsorted.vcf"])
            subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.vcf.gz"])
            subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.vcf.gz.tbi"])
        else:
            subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.unsorted.bedpe"])
            subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.bedpe.gz"])
            subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.bedpe.gz.tbi"])

        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.junction.txt"])
        subprocess.call(["rm", "-rf", output_prefix + ".mutran_tmp.junction.annot.txt"])

