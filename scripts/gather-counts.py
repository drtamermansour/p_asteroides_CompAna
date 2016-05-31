#! /usr/bin/env python
"""
This script gathers & converts Salmon output counts into something that
edgeR can read ("counts files").

Run it in a directory above all of your Salmon output directories, and
it will create a bunch of '.counts' files that you can load into R.

See https://github.com/ngs-docs/2015-nov-adv-rna/ for background info.

C. Titus Brown, 11/2015
"""
import os, os.path
import sys, getopt
import csv

def process_quant_file(root, filename, outname):
    """
    Convert individual quant.sf files into .counts files (transcripts\tcount).
    """
    print >>sys.stderr, 'Loading counts from:', root, filename
    outfp = open(outname, 'w')
    print >>outfp, "transcript\tcount"

    d = {}
    full_file = os.path.join(root, filename)
    for line in open(full_file):
        if line.startswith('#'):
            continue
        name, length, tpm, count = line.strip().split('\t')

        print >>outfp, "%s\t%s" % (name, float(count))


def main(argv):
    """
    read command line arguments
    """
    start_dir = '.'
    try:
       opts, args = getopt.getopt(argv,"hi:",["start_dir="])
    except getopt.GetoptError:
       print 'gather-counts.py -i <start_dir>'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'gather-counts.py -i <start_dir>'
          sys.exit()
       elif opt in ("-i", "--start_dir"):
          start_dir = arg
    """
    Find all the quant.sf files, convert them into properly named .counts
    files.

    Here, "proper name" means "directory.counts".
    """
    quantlist = []

    print >>sys.stderr, 'Starting in:', os.path.abspath(start_dir)
    for root, dirs, files in os.walk('.'):
        for filename in files:
            if filename.endswith('quant.sf'):
                dirname = os.path.basename(root)
                outname = dirname + '.counts'
                process_quant_file(root, filename, dirname + '.counts')
                quantlist.append(outname)
                
                break

    print ",\n".join([ "\"%s\"" % i for i in sorted(quantlist)])

if __name__ == '__main__':
    main(sys.argv[1:])
