## EEADannot

This repository holds the scripts and data files used for the manual curation of DNA motifs and cis regulatory sites, 
mostly from plants, that are eventually added as a separate library/collection  
in the database 
[footprintDB](https://floresta.eead.csic.es/footprintdb).
This library of motifs was reported for the first time in protocol 
[https://doi.org/10.1007/978-1-4939-6396-6_17](https://pubmed.ncbi.nlm.nih.gov/27557773).

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
    RL  Serra TS, Figueiredo DD, Cordeiro AM, Almeida DM, Lourenco T, Abreu IA, Sebastian A, Fernandes L, Contreras-Moreira B,Oliveira MM, Saibo NJM (2013) OsRMC, a negative regulator of salt stress response in rice, is regulated by two AP2/ERF transcription factors. Plant Mol Biol 82(4-5): 439-455
    RX  PUBMED:12913152
    RL  Cheong YH, Moon BC, Kim JK, Kim CY, Kim MC, Kim IH, Park CY, Kim JC, Park BO, Koo SC, Yoon HW, Chung WS, Lim CO, Lee SY, Cho MJ (2003) BWMK1, a rice mitogen-activated protein kinase, locates in the nucleus and mediates pathogenesis-related gene expression by activation of a transcription factor. Plant Physiol 132(4):1961-72
    XX
    FA  EREBP1
    NA  Ethylene-responsive transcription factor 1; EREBP1; OsEREBP1; ERF1_ORYSJ; Q6K7E6; Q9SE28
    SQ  MCGGAIIHHLKGHPEGSRRATEGLLWPEKKKPRWGGGGRRHFGGFVEEDDEDFEADFEEFEVDSGDSDLELGEEDDDDVVEIKPAAFKRALSRDNLSTITTAGFDGPAAKSAKRKRKNQFRGIRQRPWGKWAAEIRDPRKGVRVWLGTFNSAEEAARAYDAEARRIRGKKAKVNFPEAPTTAQKRRAGSTTAKAPKSSVEQKPTVKPAFNNLANANAFVYPSANFTSNKPFVQPDNMPFVPAMNSAAPIEDPIINSDQGSNSFGCSDFGWENDTKTPDITSIAPISTIAEVDESAFIKSSTNPMVPPVMENSAVDLPDLEPYMRFLLDDGAGDSIDSLLNLDGSQDVVSNMDLWSFDDMPVSDFY*
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
  - Make sure separators among weights/columns are TABs.
  - To convert MEME/HOMER motifs you can a one liner such as:
  
      $ perl -ane 'next if(/^>/ || /^#/); $f++; for $c (1 .. @F){ $data[$f][$c]=$F[$c-1] }; $maxc=@F if(@F>$maxc); END{ for $c (1 .. $maxc){ for $ff (1 .. $f){ printf("%1.3f\t",$data[$ff][$c]) } print "\n"} }' motif.meme

* Add individual sites, if any, to [sites.tab](./sites.tab).

* Add new papers to [references.tab](./references.tab)
  - Make sure first field matched a PWM name.
  - Second field is a TF name/primary key.
  - Third field is PubMed id.

* Actually format the library in footprintDB format:

    $ perl create_library4footprintdb.pl

