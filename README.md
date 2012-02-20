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

You need **Python** (which is most likely installed) and **Biopython** packages, so the command below should be enough:

    sudo apt-get install python-biopython

Do not forget Zim :P

Installation
------------

1. Copy `seqfetch.py` to any folder.
2. Open **Zim** and go to **Tools** > **Custom Tools**.
3. Click the `+` button to add a new tool.
4. Fill in the fields with the following:

        Name: Gene Fetcher
        Description: Fetches aminoacid sequences from NCBI
        Command: /home/you/seqfetch.py your@email.com %t %s

    Arguments are **email address**, **sequence identifier**, and **page filepath**.

5. Activating "Show in the toolbar" checkbox is useful, but not required (as 
   well as an icon).
6. Click OK and you are ready.

Remember to change to the correct file location in your system and to include 
a valid email address to the **Command** field. A valid email is required to 
use the Entrez Programming Utilities of NCBI, please refer to the [usage 
guidelines](http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen) for further information.

Usage
-----

All you need is the identifier of the sequence you want to fetch, it looks like this: **NP_523989.1**.

1. Paste the identifier to a notebook page and select it with the cursor.
2. Run Gene Fetcher tool by clicking in the toolbar icon or via Tools > Gene 
   Fetcher.
3. Data will be appended to the end of the page.
