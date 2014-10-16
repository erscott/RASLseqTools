RASLseqAligner
==============

RASLseq FASTQ reads to RASLprobe counts.


RASL-seq is a powerful and inexpensive method to assess gene expression without the need for RNA isolation1. We recently published a modified protocol using the RNA ligase Rnl2 which demonstrated dramatically increased ligation efficiency2. This python package offers an alignment method leveraging BLASTn, pandas, py-editdist, and NumPy. Optimizations will follow in the near future.

<BR>
<BR>
<h4>Usage</h4>
python /path/to/RASLseqAnalysis.py [required args -f -s -p -w -d -b -o] [optional args -P]

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
  python /path/to/RASLseqAnalysis.py -f /path/to/your.fastq.gz -s @HISEQ -p /paht/to/RASL.probes -w /path/to/annotations.bc -d  /path/to/blastdb/write_dir/ -b /path/to/blast/ncbi-blast-2.2.26+/bin/ -P -o /path/to/output.txt <BR>
  
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
  plate_bc <BR>
  well_bc <BR>
  OPTIONAL: additional columns with well metadata, column headers are user defined, e.g. drug_concentration
  <BR>
  Please see example file in data/ directory <BR>
</BLOCKQUOTE>

<BR>
<BR>

<h4>Dependencies</h4>
<ul>
BLASTn
</ul>
<ul>
pandas
</ul>
<ul>
py-editdist
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



