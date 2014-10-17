#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Protein fetcher tool for Zim.

Fetches aminoacid sequences using NCBI identifier. Sequence is written to the
notebook page with the corresponding URL for quick access.

Dependencies:
    - Python 2.7
    - Biopython 1.56

Install:
    Tools > Custom Tools > +

    Name: NCBI Protein
    Description: Fetches aminoacid sequence of a gene from NCBI
    command: ~/.local/share/zim/plugins/zim-biotools/custom_tools/ncbi_fetch_protein.py your@email.com %t %s

    Icon: no default

    UNCHECKED   Command does not modify data
    UNCHECKED   Output should replace current selection
    CHECKED     Show in the toolbar

Usage:
    Paste the gene protein ID in the page where you want to have the sequence.
    Select the ID and run the custom tool with Tools > "NCBI Protein".
'''

import sys
from Bio import Entrez, SeqIO

# User email
user_email = sys.argv[1]

# Sequence ID
seq_id = sys.argv[2]

# Output file
page = sys.argv[3]

#Set email address required for Entrez.
Entrez.email = user_email

# Fetch data via Entrez using the sequence id as FASTA.
handle = Entrez.efetch(db='protein', id=seq_id, rettype='fasta')

# Create record from FASTA.
record = SeqIO.read(handle, 'fasta')

# Parse gene id to generate url.
gene_id = record.id.split('|')[1]
url = 'http://www.ncbi.nlm.nih.gov/protein/%s?report=fasta' % gene_id

# Parse organism (dirtily...).
species = record.description.split('[')[1].split(']')[0].replace(' ', '_')

# Open notebook page as file.
f = open(page, 'a')

# Print protein sequence.
f.write("""@{species}\n\n{url}\n\n@gene\n'''\n{sequence}\n'''\n""".format(
    url=url, species=species, sequence=record.format('fasta')))

# Close file.
f.close()

# Exit program.
sys.exit(0)
