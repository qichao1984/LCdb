# LCdb
Alternative links containing unzipped database files: https://www.alipan.com/s/Gs62Gxax7ii 

# LCdb: A Curated Functional Gene Database for Metagenomic Profiling of Lignin Catabolism Pathways

Chen, J., Lin, L., Tu, Q., Peng, Q., Wang, X., Liang C., Zhou, J.,and Yu, X., (2023), LCdb: A Curated Functional Gene Database for Metagenomic Profiling of Lignin Catabolism Pathways.

Lignin, as an abundant organic carbon, plays a vital role in the global carbon cycle, and the lignin catabolism driven by microorganisms is an important biogeochemical cycling process of the Earth's biosphere. Shotgun metagenome sequencing has opened a new avenue to advance our understanding of lignin catabolism microbial communities. However, accurate metagenomic profiling of lignin catabolism microbial communities remains technically challenging, mainly due to low coverage of lignin catabolism genes/pathways, difficulties in distinguishing homologous genes and a long research time on publicly available orthology databases. It is essential to develop a comprehensive and accurate database for characterizing lignin catabolism microbial communities in metagenomic studies. To solve those problems, we constructed a manually curated lignin catabolism database (LCdb) for metagenome sequencing data analysis of lignin catabolism microbial communities in the environment. 

The developed LCdb contains 474 gene families and 471,705 representative sequences affiliated with 62 phyla of bacteria/archaea, and 379,649 homologous orthology groups were also included to reduce false positive sequence assignments. 

Four files are included in LCdb:

<b>1. LCdb.zip</b>: fasta format representative sequences obtained by clustering curated sequences at 100% sequence identity. This file can be used for "BLAST" searching LCdb genes in shotgun metagenomes.

<b>2. id2genemap</b>: a mapping file that maps sequence IDs to gene names, only sequences belonging to LCdb gene families are included. Sequences for LCdb homologs are not included. This file is used to generate LCdb profiles from BLAST-like results against the LCdb database. 

<b>3. LCdb_FunctionProfiler.PL</b>: a perl script for functional profiling of lignin catabolism genes.

<b>4. LCdb_TaxonomyProfiler.PL</b>: a perl script for taxonomical profiling of lignin catabolism microbial communities.

<b>DOWNLOAD/INSTALLATION</b>

git clone https://github.com/qichao1984/LCdb.git

<b>Dependencies and Tools</b>

<i>Perl modules that can be easily installed via cpan:</i>
<p>List::Util</p>
<p>Getopt::Long</p>
<i>Dependencies for LCdb_FunctionProfiler.PL, currently supported database searching tools are: </i>
<p>usearch: https://www.drive5.com/usearch/download.html
<p>diamond: https://github.com/bbuchfink/diamond/releases
<p>blast: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz</p>
<i>Dependencies for LCdb_TaxonomyProfiler.PL:</i>
<p>seqtk: https://github.com/lh3/seqtk.git
<p>kraken2: https://github.com/DerrickWood/kraken2.git

<b>USAGE</b>

Before getting started, please modify both scripts (LCdb_FunctionProfiler.PL, LCdb_TaxonomyProfiler.PL) at lines 6-18 to specify the locations of third party tools and their parameters. If the tools are already in the system path, no revision is needed. By default, basic parameters are used for these tools. Users are encouraged to make revisions in cases of short reads and/or expecting more strict/relaxed results. We also encourage users to develop useful implementations based on LCdb.

Note: Kraken2 database could be downloaded from https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads, or built locally.

<b>Example for using LCdb_FunctionProfiler.PL:</b>

perl LCdb_FunctionProfiler.PL -d \<workdir\> -m \<diamond|usearch|blast\> -f \<filetype\> -s \<seqtype\> -si \<sample size info file\> -rs \<random sampling size\> -o \<outfile\>
  
Detailed explanations: 

-d : specify the directory where your fasta/fastq (or gzipped) files are located. 

-m : specify the database searching program you plan to use, currently diamond, usearch and blast are supported. 

-f : specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s : sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have, note that no file extensions should be included here.

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

-o : the output file for N cycle gene profiles.   


<b>Example for using LCdb_TaxonomyProfiler.PL:</b>

perl LCdb_TaxonomyProfiler.PL -d \<workdir\> -m \<diamond|usearch|blast\> -f \<filetype\> -s \<seqtype\> -si \<sample size info file\> -rs \<random sampling size\> 
  
Detailed explanations: 

-d : specify the directory where your fasta/fastq (or gzipped) files are located. 

-m : specify the database searching program you plan to use, currently diamond, usearch and blast are supported. 

-f : specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s : sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have, note that no file extensions should be included here.

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

  
