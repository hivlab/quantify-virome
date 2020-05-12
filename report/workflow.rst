Taxonomy was assigned to virus reads from metagenomic sequence library using `quantify-virome workflow`_:
Reads were mapped onto host and bacterial genomes with `BBMap`_, and both optical and PCR duplicates were removed with clumpify_.
Similar nucleotide sequences were clustered using CD-HIT_.
Low complexity and short-period tandem repeats were masked in nucleotide sequences using tantan_.
Interspersed repeats and low complexity DNA sequences were masked in nucleotide sequences using RepeatMasker_.
Taxonomy was assigned to reads using simultaneous BLAST_ queries against virus and nonredundant database.  
Quality control was performed with FastQC_, `FastQ Screen`_ and aggregated into an interactive report via MultiQC_.

.. _BBMap: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. _clumpify: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/
.. _quantify-virome workflow: https://github.com/avilab/sarscov2
.. _MultiQC: http://multiqc.info/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _FastQ Screen: https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _CD-HIT: http://weizhongli-lab.org/cd-hit/
.. _tantan: http://cbrc3.cbrc.jp/~martin/tantan/
.. _RepeatMasker: http://www.repeatmasker.org
.. _BLAST: https://www.ncbi.nlm.nih.gov/books/NBK279690/
