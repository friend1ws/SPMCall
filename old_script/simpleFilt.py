#! /usr/local/bin/python

import sys

inputFile = sys.argv[1]
thres = sys.argv[2]

hIN = open(inputFile, 'r')

for line in hIN:
    F = line.rstrip('\n').split('\t')

    if int(F[22]) >= int(thres):
        print '\t'.join(F[0:3] + F[22:23])

