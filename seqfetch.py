#!/usr/bin/python 
# -*- coding: utf-8 -*-

'''Gene fetcher tool for Zim.

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

def main(user_email, seq_id, page):
    #Set email address required for Entrez.
    Entrez.email = user_email

    # Open notebook page as file.
    f = open(page, 'a')

    # Fetch data via Entrez using the sequence id as FASTA.
    handle = Entrez.efetch(db='protein', id=seq_id, rettype='fasta')

    # Create record from FASTA.
    record = SeqIO.read(handle, 'fasta')

    # Parse gene id to generate url.
    gene_id = record.id.split('|')[1]
    url = 'http://www.ncbi.nlm.nih.gov/protein/%s?report=fasta' % gene_id

    # Parse organism (dirtily...).
    species = record.description.split('[')[1].split(']')[0].replace(' ', '_')

    # Write data to file.
    f.write('\n')

    # Default tags.
    f.write('@%s ' % species)
    f.write('\n\n')

    # NCBI sequence url.
    f.write(url)
    f.write('\n\n')

    # Sequence itself.
    f.write('@gene \n')
    f.write("'''\n")
    f.write(record.format('fasta'))
    f.write("\n'''\n")

    # Close file.
    f.close()

    # Exit program.
    sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
