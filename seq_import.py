#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Tool for importing sequence files in Zim.'''

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

from biotools import make_header

#TODO Use SeqRecord as standard format for easy export.
#TODO Print useful information such as length.

def main(seq_dir, root_page):
    # Remove extension.
    root = root_page[:-4]
    # Make sure there is a directory for sequences.
    try:
        os.mkdir(root)
    except:
        pass
    # Get sequences from dir.
    sequences = os.listdir(seq_dir)
    # Only filenames with .seq extension.
    sequences = [seq for seq in sequences if seq.endswith('.seq')]

    # Iterate over each sequence file.
    for seq_filename in sequences:

        sequence = open(os.path.join(seq_dir, seq_filename))

        # Parse filename as ID and ID as filename.
        seq_id = os.path.basename(sequence.name.split('-')[0])
        #seq_name

        # Define locus file.
        seq_file = os.path.join(root, '%s.txt' % seq_id)

        # Verify if page already exists.
        try:
            seq_page = open(seq_file, 'r')
            seq_page.close()
            # If it does, skip page.
            continue
        except:
            # If not, create page.
            pass

        # Create locus_page and write header.
        seq_page = open(seq_file, 'w')
        seq_page.write(make_header(os.path.basename(seq_id)))

        # Write sequence in FASTA format.
        seq_page.write('@sequence \n')
        seq_page.write("'''\n")
        seq_page.write(sequence.read())
        seq_page.write("\n'''\n")

        # Close locus file.
        seq_page.close()

    # Exit program.
    sys.exit(0)

if __name__ == '__main__':
    # arg1: directory with sequence files
    # arg2: root page
    main(sys.argv[1], sys.argv[2])
