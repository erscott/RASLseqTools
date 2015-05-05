RASLseqTools
==============
RASLseq FASTQ reads to RASLprobe counts.


RASL-seq is a powerful and inexpensive method to assess gene expression without the need for RNA isolation1. We recently published a modified protocol using the RNA ligase Rnl2 which demonstrated dramatically increased ligation efficiency2. This python package offers an alignment method leveraging BLASTn, pandas, py-editdist, and NumPy. Optimizations will follow in the near future.

<BR>
<BR>

<h4>19 April 2015</h4>
<BR>Initial public release RASLseqTools, Default Aligner is now STAR (https://github.com/alexdobin/STAR)<BR>
<BR>
Example IPython Notebook: http://nbviewer.ipython.org/github/erscott/RASLseqTools/blob/master/ipynb/RASLseqTools_STAR_example.ipynb



<BR>
<h4>STAR Usage</h4>
python /path/to/RASLseqAnalysis_STAR.py [required args -f -a -p -w -d -o] [optional args -P -A -n -o5 -o3 -ws -we]

<BLOCKQUOTE>

  -f : str, absolute path to fastq file(s) (accepts gzip files, comma-separated list if multiple fastqs) <BR>
  
  -a : str, absolute path to STAR bin directory <BR>
  
  -p : str, absolute path to probes file <BR>

  -w : str, absolute path to annotations file <BR>

  -d : str, absolute path to output directory <BR>

  -o : str, absolute path of output file <BR>

  <BR>
  
  -P : bool, verbose printing <BR>
  
  -A : bool, Write STAR alignments to disk, will be written in output directory <BR>
  
  -n : int, number of jobs, currently requires 2 processors <BR>
  
  -o5: int, number of bases to clip from 5-prime end of read to isolate probe sequence, default=24 <BR>
  
  -o3: int, number of bases to clip from 3-prime end of read to isolate probe sequence, default=22 <BR>
  
  -ws: int, index position of the wellbarcode start base in read, default=0 <BR>
  
  -we: int, index position of the wellbarcode end base in read, default=8 <BR>
  
  
  <BR>

  example command: <BR>
  python /path/to/RASLseqAnalysis_STAR.py -f /path/to/your.fastq.gz -a /path/to/STAR_binary/ -p /paht/to/RASL.probes -w /path/to/annotations.bc -d /path/to/write_directory/ -o /path/to/blastdb/write_file.txt  -P -A -n 1 -o5 25 -03 20 -ws 0 -we 8 <BR>
  <BR>
  example data: can be found in the data directory <BR>
  



<BR>
<h4>BLASTn Usage</h4>
python /path/to/RASLseqAnalysis_BLAST.py [required args -f -s -p -w -d -b -o] [optional args -P]

<BLOCKQUOTE>

  -f : absolute path to fastq file (accepts gzip files) <BR>
  
  -p : absolute path to probes file <BR>

  -w : absolute path to annotations file <BR>

  -d : absolute path to write directory for blast database <BR>

  -b : absolute path to blast bin directory <BR>

  -o : absolute path of output file <BR>
  
  -s : specifies sequencer id in fastq index line, e.g. @HISEQ <BR>

  <BR>
  
  -P : verbose printing <BR>


  <BR>

  example command: <BR>
  python /path/to/RASLseqAnalysis_BLAST.py -f /path/to/your.fastq.gz -s @HISEQ -p /paht/to/RASL.probes -w /path/to/annotations.bc -d  /path/to/blastdb/write_dir/ -b /path/to/blast/ncbi-blast-2.2.26+/bin/ -P -o /path/to/output.txt <BR>
  
  example data: can be found in the data directory <BR>
  
  

</BLOCKQUOTE>

<h6>NOTE:RASLseqAnalysis_NAR.py is provided for transparency and requires manual parameter settings to run</h6>



<BR>
<BR>
<h4>Input File Formats</h4>

FASTQ: standard FASTQ format (optionally gzipped)<BR>
<BR>
X.probes: tab-separated file describing the RASLseq Probes with the following columns and column headers <BR>
<BLOCKQUOTE>
  AcceptorProbeSequence <BR>
  DonorProbeSequence <BR>
  AcceptorAdaptorSequence <BR>
  DonorAdaptorSequence <BR>
  ProbeName <BR>
  <BR>
  Please see example file in data/ directory <BR>
</BLOCKQUOTE>
<BR>
X.bc: tab-separated file describing each well in the experiment with the following columns and column headers<BR>
<BLOCKQUOTE>
  REQUIRED: <BR>
  PlateBarcode <BR>
  WellBarcode <BR>
  OPTIONAL: additional columns with well metadata, column headers are user defined, e.g. drug_concentration
  <BR>
  Please see example file in data/ directory <BR>
</BLOCKQUOTE>

<BR>
<BR>

<h4>Dependencies</h4>
<ul>
STAR aligner
</ul>

<ul>
BLASTn
</ul>
<ul>
pandas
</ul>
<ul>
Levenshtein editdist
</ul>
<ul>
NumPy
</ul>




<BR>
<BR>
<h4>References</h4>

1. H. Li, J. Qiu, X.-D. Fu, RASL-seq for massive parallel and quantitative analysis of gene expression, 
Curr. Protocol. Mol. Biol., 98 (2012), pp. 4.13.1–4.13.9

2. Larman HB, Scott ER, Wogan M, Oliveira G, Torkamani A, Schultz PG,  Sensitive, multiplex and direct quantification of RNA sequences using a modified RASL assay, Nucleic Acids Res. 2014;42(14):9146-57



