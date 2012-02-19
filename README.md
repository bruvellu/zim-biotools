Simple gene fetcher tool for Zim
================================

This is a simple custom tool to fetch aminoacid sequences from NCBI and append 
it to a notebook page in [Zim Desktop Wiki](http://zim-wiki.org/).

Example of output:

    @Drosophila_melanogaster @gene 
    
    http://www.ncbi.nlm.nih.gov/protein/17647207?report=fasta
    
    >gi|17647207|ref|NP_523989.1| boule, isoform B [Drosophila melanogaster]
    MHKIAAAPPPSATPGGGLETPLAAPKYGTLIPNRIFVGGISGDTTEADLTRVFSAYGTVK
    STKIIVDRAGVSKGYGFVTFETEQEAQRLQADGECVVLRDRKLNIAPAIKKQPNPLQSIV
    ATNGAVYYTTTPPAPISNIPMDQFAAAVYPPAAGVPAIYPPSAMQYQPFYQYYSVPMNVP
    TIWPQNYQENHSPLLHSPTSNPHSPHSQSHPQSPCWSIEDLRDTLPRV

Dependencies
------------

For the script to work you need to have Python and Biopython installed; 
standard repository packages should suffice. Do not forget Zim :P

Installation
------------

1. Copy `seqfetch.py` to any folder.
2. Open **Zim** and go to Tools > Custom Tools.
3. Click the `+` button to add a new tool.
4. Fill in the fields with the following:

Name: Gene Fetcher
Description: Fetches aminoacid sequences from NCBI
Command: `/home/you/seqfetch.py your@email.com %t %s`

5. Activating "Show in the toolbar" checkbox is useful, but not required (as 
   well as an icon).
6. Click OK and you are ready.

Email is a requirement by Entrez, otherwise the script will not work.

Usage
-----

All you need to fetch a sequence is the identifier, it looks like this: **NP_523989.1**.

1. Paste the identifier to a notebook page and select it with the cursor.
2. Run Gene Fetcher tool by clicking in the toolbar icon or via Tools > Gene 
   Fetcher.
3. Data will be appended to the end of the page.
