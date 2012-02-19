#!/usr/bin/python 
# -*- coding: utf-8 -*-

'''Gene Fetcher tool for Zim.

Fetch aminoacid sequences using NCBI identifier. Sequence is written to the
notebook page with the corresponding URL for quick access.

Although any identifier can be used, there may be unexpected problems
(eg, Zim might hang if sequence is too long; url will not be correct).
Thus, only protein sequences are currently supported.

Dependencies:
    Python 2.7
    Biopython 1.56

'''

import sys
from Bio import Entrez, SeqIO
Entrez.email = 'email@email.com' # Fill with your email address.

def main(seq_id, page):
    # Open notebook page as file.
    f = open(page, 'a')
    # Fetch data via Entrez using the sequence id as FASTA.
    handle = Entrez.efetch(db='protein', id=seq_id, rettype='fasta')
    # Create record from FASTA.
    record = SeqIO.read(handle, 'fasta')
    # Parse gene id to generate url.
    gene_id = record.id.split('|')[1]
    url = 'http://www.ncbi.nlm.nih.gov/protein/%s?report=fasta' % gene_id

    # Write data to file.
    f.write('\n')
    f.write(url)
    f.write('\n\n')
    f.write(record.format('fasta'))
    # Close file.
    f.close()
    # Exit program.
    sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
