#!/usr/bin/python
# -*- coding: utf-8 -*-

'''FASTA importer tool for Zim.

Dependencies:
    - Python 2.7
    - Biopython 1.56

Install:
    Tools > Custom Tools > +

    Name: Import FASTA
    Description: Import FASTA file creating individual subpages for each entry
    command: ~/.local/share/zim/plugins/zim-biotools/custom_tools/import_fasta.py %t %s

    Icon: no default

    UNCHECKED   Command does not modify data
    UNCHECKED   Output should replace current selection
    CHECKED     Show in the toolbar

Usage:
    Paste the path to a fasta file. Select it and run the custom tool with
    Tools > "Import FASTA". One locus per page will be created under a new
    subpage named "Loci".
'''

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from biolib import make_header, make_dir

# Arguments.
fasta_file = sys.argv[1]
root_page = sys.argv[2]

# Get absolute paths.
absroot = os.path.abspath(root_page)
absfasta = os.path.abspath(os.path.expanduser(fasta_file))

# Remove extension.
root = os.path.splitext(absroot)[0]

# Make sure there is a folder to the root.
make_dir(root)

# Define organism tag by reading the root name.
organism = os.path.basename(root)

# Parse FASTA file.
loci = SeqIO.parse(absfasta, 'fasta')

# Create loci page and directory.
loci_file = os.path.join(root, 'Loci.txt')
loci_page = open(loci_file, 'w')
loci_page.write(make_header('Loci'))
loci_page.close()

# Define Loci directory.
loci_dir = os.path.join(root, 'Loci')
make_dir(loci_dir)

# Iterate over each locus.
for locus in loci:

    # Sanitize locus.id.
    locus_id = locus.id.replace('/', '-')

    # Define locus file.
    locus_file = os.path.join(loci_dir, '{}.txt'.format(locus_id))

    # Verify if page already exists.
    try:
        locus_page = open(locus_file, 'r')
        locus_page.close()
        # If it does, skip page.
        continue
    except:
        # If not, create page.
        pass

    # Identify the frame and the strand by parsing the description.
    # Example: >Lrub_5432 | frame: +1 | candidates: Six3-6, Optix
    frame = int(locus.description.split('|')[-2][-2])
    frame_step = frame - 1
    strand = locus.description.split('|')[-2][-3]

    # Always garantee that the sequence is a plus strand.
    if strand == '-':
        locus.seq = locus.seq.reverse_complement()
        locus.description = locus.description.replace('frame: {}{}'.format(strand, frame), 'frame: +{}'.format(frame))

    # Translate using the correct frame.
    translated_seq = locus.seq[frame_step:].translate()

    # Create SeqRecord for protein.
    protein = SeqRecord(translated_seq, id=locus.id, name=locus.name, description=locus.description)

    # Create locus_page and write header.
    locus_page = open(locus_file, 'w')
    locus_page.write("{header}".format(header=make_header(locus_id)))

    # Write organism name.
    locus_page.write("@{}\n\n".format(organism))

    # Write sequence in FASTA format.
    locus_page.write("@locus {} bp \n'''\n".format(len(locus.seq)))
    locus_page.write("'''\n")
    locus_page.write(locus.format('fasta'))
    locus_page.write("\n'''\n")

    # Write protein sequence.
    locus_page.write('@protein \n')
    locus_page.write("'''\n")
    locus_page.write(protein.format('fasta'))
    locus_page.write("\n'''\n")

    # Close locus file.
    locus_page.close()

# Exit program.
sys.exit(0)
