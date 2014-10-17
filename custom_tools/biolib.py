#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Common library for zim-biotools.'''

import os
from datetime import datetime


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


def write_fasta(seq_record, tag='locus'):
    '''Use SeqRecord format function to write fasta.'''
    string = "@%s %d bp \n" \
                "'''\n" \
                "%s"\
                "\n'''\n" % (tag, len(seq_record.seq), seq_record.format('fasta'))
    return string


def make_dir(path):
    '''Make dir, if it does not exist.'''
    try:
        os.mkdir(path)
    except:
        pass


def file_exists(path):
    '''Return True if file exists.'''
    try:
        a_file = open(path, 'r')
        a_file.close()
        return True
    except:
        return False


def remove_extension(filename):
    '''Gracefully remove file extension.'''
    root = os.path.splitext(filename)[0]
    return root


def sanitize_zim(string):
    '''Replace illegal characters for Zim pages.

    At the moment it only replaces "/" for "-"
    '''
    return string.replace('/', '-')