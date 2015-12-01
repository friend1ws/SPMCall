#! /usr/bin/env python

import sys
import random
import subprocess

input_file = sys.argv[1]
output_dir = sys.argv[2]
repeat_num = sys.argv[3]
additional_params = sys.argv[4]

ID2mutation = {}
ID2junction = {}

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        ID2mutation[F[0]] = F[1]
        ID2junction[F[0]] = F[2]

        # alert for duplication?


##########
# right combination

hout = open(output_dir + ("/splicing_sv_list.txt" if "--sv" in additional_params else "/splicing_mutation_list.txt"), 'w')
for ID in sorted(ID2junction): 
    output_prefix = output_dir + '/' + ID
    params = additional_params.split(' ')
    print ' '.join(["mutrans", ID2mutation[ID], ID2junction[ID], output_prefix] + params)
    subprocess.call(["mutrans", ID2mutation[ID], ID2junction[ID], output_prefix] + params)
    
    output_file = output_prefix + (".splicing_sv.txt" if "--sv" in additional_params else ".splicing_mutation.txt")
    with open(output_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            print >> hout, ID + '\t' + line

hout.close()


##########
# permutation
for i in range(int(repeat_num)):

    print "permutation:" + str(int(i) + 1)

    perm_success_flag = 0
    while perm_success_flag == 0:
        ID_perm = {}
        ID_remain = ID2junction.keys()
        for ID in ID2junction:
            ID_set = list(set(ID_remain) - set([ID]))
            if len(ID_set) == 0: break 
            ID_selected = random.choice(ID_set)
            ID_perm[ID] = ID_selected
            ID_remain.remove(ID_selected)

        if len(ID_perm) == len(ID2junction): perm_success_flag = 1
 

    hout = open(output_dir + "/permtutation_" + str(int(i) + 1) + ".txt", 'w')    
    for ID in sorted(ID2junction):
        output_prefix = output_dir + '/' + ID + '.' + ID_perm[ID]
        params = additional_params.split(' ')
        print ' '.join(["mutrans", ID2mutation[ID], ID2junction[ID_perm[ID]], output_prefix] + params)
        subprocess.call(["mutrans", ID2mutation[ID], ID2junction[ID_perm[ID]], output_prefix] + params)

        output_file = output_prefix + (".splicing_sv.txt" if "--sv" in additional_params else ".splicing_mutation.txt")
        with open(output_file, 'r') as hin:
            for line in hin:
                line = line.rstrip('\n')
                print >> hout, ID + '\t' + ID_perm[ID] + '\t' + line

    hout.close()
##########

