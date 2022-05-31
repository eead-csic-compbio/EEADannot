## EEADannot

This repository holds the scripts and data files used for the manual curation of DNA motifs and cis regulatory sites, 
mostly from plants, that are eventually added as a separate library 
([EEADannot](https://floresta.eead.csic.es/footprintdb/index.php?database=28&type=motif&page=1)) 
in the database 
[footprintDB](https://floresta.eead.csic.es/footprintdb/index.php).
This library of motifs was reported for the first time in protocol 
[https://doi.org/10.1007/978-1-4939-6396-6_17](https://pubmed.ncbi.nlm.nih.gov/27557773).

### Steps

* Add new TFs to [./TFsequences.faa](TFsequences.faa)
  - Make sure 'FullName' has no blanks.
* Add new motifs to [./PWM.tab](PWM.tab). 
  - To convert MEME/HOMER motifs you can a one liner such as:
  
      perl -ane 'next if(/^>/ || /^#/); $f++; for $c (1 .. @F){ $data[$f][$c]=$F[$c-1] }; $maxc=@F if(@F>$maxc); END{ for $c (1 .. $maxc){ for $ff (1 .. $f){ printf("%1.3f\t",$data[$ff][$c]) } print "\n"} }' MdDAM4-MdSVPa.meme

* Add new papers to [./references.tab](references.tab)
  - Make sure first field matched a PWM name.
  - Second field is a TF name/primary key.
  - Third field is PubMed id.

* Actually format the library in footprintDB format:

    perl create_library4footprintdb.pl

