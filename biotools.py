#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Common library for zim-biotools.'''

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
