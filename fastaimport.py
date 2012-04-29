#!/usr/bin/python
# -*- coding: utf-8 -*-

'''FASTA importer tool for Zim.


Dependencies:
    Python 2.7
    Biopython 1.56

'''

import os
import sys
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import translate

#FIXME If there are no subpages, create them.

def make_header(title):
    '''Return string with header of a Zim page.'''

    timestamp = datetime.now().replace(microsecond=0)
    header = '''Content-Type: text/x-zim-wiki\nWiki-Format: zim 0.4\nCreation-Date: %s\n\n====== %s ======\nCreated %s\n\n''' % (timestamp.isoformat(), title, timestamp.strftime('%A %d %B %Y'))
    return header

def print_seq(sequence, page):
    '''Print well-formatted sequence.'''
    i = 0
    for char in sequence:
        if i == 60:
            page.write('\n')
            i = 0
        page.write(char)
        i += 1

def main(fasta_file, root_page):
    # Remove extension.
    root = root_page[:-4]

    # Define organism tag by reading the root name.
    organism = os.path.basename(root)

    # Parse FASTA file.
    loci = SeqIO.parse(fasta_file, 'fasta')

    # Create loci page and directory.
    #XXX Write some information in the loci page? Maybe log?
    loci_file = os.path.join(root, 'Loci.txt')
    loci_page = open(loci_file, 'w')
    loci_page.write(make_header('Loci'))
    loci_page.close()

    # Define Loci directory.
    loci_dir = os.path.join(root, 'Loci')
    try:
        os.mkdir(loci_dir)
    except:
        pass

    # Iterate over each locus.
    for locus in loci:

        # Sanitize locus.id.
        locus_id = locus.id.replace('/', '-')

        # Define locus file.
        locus_file = os.path.join(loci_dir, '%s.txt' % locus_id)

        # Verify if page already exists.
        try:
            locus_page = open(locus_file, 'r')
            locus_page.close()
            # If it does, skip page.
            continue
        except:
            # If not, create page.
            pass

        # Make plus strand and protein sequence.
        frame = int(locus.description.split('|')[-2][-2])
        frame_step = frame - 1
        strand = locus.description.split('|')[-2][-3]

        # Always garantee that the sequence is a plus strand.
        if strand == '-':
            locus.seq = locus.seq.reverse_complement()
            locus.description = locus.description.replace('frame: %s%d' % (strand, frame), 'frame: +%d' % frame)

        # Convert to plain string to translate into protein.
        sequence = str(locus.seq)
        sequence = sequence[frame_step:]
        protein = translate(sequence)

        # Create locus_page and write header.
        locus_page = open(locus_file, 'w')
        locus_page.write(make_header(locus_id))

        # Write organism name.
        locus_page.write('@%s ' % organism)
        locus_page.write('\n\n')

        # Write sequence in FASTA format.
        #TODO Format into verbatim.
        locus_page.write('@locus \n')
        locus_page.write("'''\n")
        locus_page.write(locus.format('fasta'))
        locus_page.write("\n'''\n")

        # Write protein sequence.
        locus_page.write('@protein \n')
        locus_page.write("'''\n")
        print_seq(protein, locus_page)
        locus_page.write("\n'''\n")

        # Close locus file.
        locus_page.close()

    # Exit program.
    sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
