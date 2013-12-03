#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Import sequence files to Zim.
    
    This script takes all .seq files of a folder and creates a new Zim page for each individual sequence.

    Why: used to import sequence files returned by our genome sequencing facilities.
'''

import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biolib import *


def main(seq_dir, root_page):

    # Define subfolder by removing file extension from root page..
    root = remove_extension(root_page)

    # Make sure there is a directory for imported sequences.
    make_dir(root)

    # Get sequences from directory.
    sequences = os.listdir(seq_dir)

    # Get only filenames with .seq extension.
    sequences = [seq for seq in sequences if seq.endswith('.seq')]

    # Iterate over each sequence file.
    for seq_filename in sequences:

        # Open sequence file for parsing.
        seq_filepath = os.path.join(seq_dir, seq_filename)
        sequence = open(seq_filepath)

        # New filepath for sequence file.
        seq_basename = os.path.splitext(seq_filename)[0]
        seq_basename = sanitize_zim(seq_basename)
        seq_zimpath = os.path.join(root, '%s.txt' % seq_basename)

        # If file already exists, skip sequence; never overwrite.
        if file_exists(seq_zimpath):
            continue

        # Write new file.

        # Create locus_page and write header.
        seq_page = open(seq_zimpath, 'w')
        seq_page.write(make_header(seq_basename))

        # Create SeqRecord for sequence.
        seq_itself = sequence.read().replace('\r', '').replace('\n', '').replace('\t', '')
        seq_record = SeqRecord(Seq(seq_itself), id=seq_basename, name=seq_basename, description='')

        # Write sequence in fasta format.
        seq_page.write(write_fasta(seq_record, tag='sequence'))

        # Close locus file.
        seq_page.close()

    # Exit program.
    sys.exit(0)

if __name__ == '__main__':
    # arg1: directory with sequence files
    # arg2: root page
    main(sys.argv[1], sys.argv[2])
