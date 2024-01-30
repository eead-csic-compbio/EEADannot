## EEADannot

This repository holds the scripts and data files used for the manual curation of DNA motifs and cis regulatory sites, 
mostly from plants, that are eventually added as a separate library/collection  
in the database 
[footprintDB](https://floresta.eead.csic.es/footprintdb), also available at [RSAT::Plants](https://rsat.eead.csic.es/plants),
which are part of the [INB/ELIXIR-ES](https://inb-elixir.es) resources portfolio.

### Citation

Contreras-Moreira B, Sebastian A. FootprintDB: Analysis of Plant Cis-Regulatory Elements, 
Transcription Factors, and Binding Interfaces. Methods Mol Biol. 2016; 1482:259-77. 
doi: [10.1007/978-1-4939-6396-6_17](https://doi.org/10.1007/978-1-4939-6396-6_17). PMID: [27557773](https://pubmed.ncbi.nlm.nih.gov/27557773).

### Example motif in TRANSFAC format

    MO  EREBP1
    NA  Os02g54160.1
    P0      A      C      G      T
    01   1.00   0.00   0.00   0.00
    02   0.00   0.00   1.00   0.00
    03   0.00   1.00   0.00   0.00
    04   0.50   0.50   0.00   0.00
    05   0.00   0.00   1.00   0.00
    06   0.00   0.50   0.50   0.00
    07   0.00   1.00   0.00   0.00
    RX  PUBMED:23703395
    RL  Serra TS et al (2013) OsRMC, a negative regulator of salt stress response in rice...Plant Mol Biol 82(4-5): 439-455
    RX  PUBMED:12913152
    RL  Cheong YH et al (2003) BWMK1, a rice mitogen-activated protein kinase... Plant Physiol 132(4):1961-72
    XX
    FA  EREBP1
    NA  Ethylene-responsive transcription factor 1; EREBP1; OsEREBP1; ERF1_ORYSJ; Q6K7E6; Q9SE28
    SQ  MCGGAIIHHLKGHPEGSRRATEGLLWPEKKKPRWGGGGRRHFGGFVEEDDEDFEADFEEFEVDSGDSDLELGEEDDDDVVEI...
    OS  Oryza sativa
    CC  family:AP2/ERF
    XX
    SI  EREBP1_1
    SQ  gAGCCGCCa
    XX
    SI  EREBP1_2
    SQ  gAGCAGGCa
    XX
    //

### Production steps

* Add new TFs to [TFsequences.faa](./TFsequences.faa)
  - Make sure 'FullName' has no blanks.

* Add new motifs to [PWM.tab](./PWM.tab).
  - In name line, 1st word is motif name [no spaces], 2nd is 1+ TF names separated by commas [,].
  - Make sure separators among weights/columns are TABs.
  - To convert MEME/HOMER motifs you can use a one liner such as:
  
      perl -ane 'next if(/^>/ || /^#/); $f++; for $c (1 .. @F){ $data[$f][$c]=$F[$c-1] }; $maxc=@F if(@F>$maxc); END{ for $c (1 .. $maxc){ for $ff (1 .. $f){ printf("%1.3f\t",$data[$ff][$c]) } print "\n"} }' motif.meme

* Add individual sites, if any, to [sites.tab](./sites.tab).

* Add new papers to [references.tab](./references.tab)
  - Make sure 1st field matches a motif name in PWM.tab.
  - Second field is one or more TF name/primary key separated by commas [,].
  - Third field is PubMed id.

* Actually format the library in footprintDB format:

    $ perl create_library4footprintdb.pl
