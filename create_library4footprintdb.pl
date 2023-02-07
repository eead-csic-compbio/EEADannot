#!/usr/bin/perl -w
use strict;

# Script to make a footprintDB library of motifs & transcription factors (TF) from local files
#
# $TFSEQFILE : FASTA file with peptide TF sequences with metadata in header; 1st field is primary key
# $SITEFILE  : TSV file with 1/cis regulatory site per line; linked to TF through primary key, there might be several lines for the same TF
# $REFSFILE  : CSV file with one PubMed entry per line; it might be repeated for several motifs
# $PWMFILE   : consensus TSV file where each matrix has two metadata: motif name\tprimary key; there might be several for same TF

my $PWMFILE     = 'PWM.tab';
my $REFSFILE    = 'references.tab';
my $TFSEQFILE   = 'TFsequences.faa';
my $SITEFILE    = 'sites.tab';

my $WLIBRARY    = 'EEADannot.fdb'; 

my ($sec,$min,$hour,$day,$mon,$year) = localtime(time);
my $timestamp = sprintf("%s%02d%02d",$year+1900,$mon+1,$day);
my $date = sprintf("%s-%02d-%02d",$year+1900,$mon+1,$day);

my $VERSION = $date; #'1.03';

my $LIBRARYHEADER = <<EOH;
VV  Name: EEADannot; Version: $VERSION; Date: $date; 
VV  Authors: Contreras-Moreira B, Sebastian A.
VV  Url: http://floresta.eead.csic.es/footprintdb; Email: compbio\@eead.csic.es
XX
EOH

###################################################

my (%TFsequence,@acc,%synonyms,%species,%family,%fullname,%sites,%pubmeds,%references,%seen);
my ($tfaccession,$syns,$site,$pubmed,$reference,$nt,$bs,$p,$weights,$motifname,$pwm_width);

my ($n_of_references,$n_of_sites,$n_of_pwms) = (0,0,0);

my (%bindingdata,%motif);

## 1) read TF sequences
open(TFSEQ,$TFSEQFILE) || die "# $0 : cannot read $TFSEQFILE\n";
while(my $line = <TFSEQ>)
{
  #>Os01g64790.1 | Species: Oryza sativa | Symbols: EREBP2; OsEREBP2; ERF99 | Family: AP2/ERF | FullName: EREBP2

  next if($line =~ /^$/);
  chomp($line); 
  if($line =~ /^>(\S+)/) {
    $tfaccession = $1;
    push(@acc,$tfaccession);
    if($line =~ /Species: (.*?) \|/) {
      $species{$tfaccession} = $1; 
    }
    if($line =~ /Symbols: (.*?) \|/) {
      $syns = $1; 
      $syns =~ s/,/;/g;
      $synonyms{$tfaccession} = $syns;
    }      
    if($line =~ /Family: (.*?) \|/) {
      $family{$tfaccession} = $1;
    }
    if($line =~ /FullName: (.*)/) {
      $fullname{$tfaccession} = $1;
    }
  } else { 
    $line =~ s/ //g;
    $TFsequence{$tfaccession} .= $line; 
  }
}	
close(TFSEQ); 
 
print "# read ".scalar(keys(%TFsequence))." sequences\n"; 

## 2) read sites
open(SITES,$SITEFILE) || die "# $0 : cannot read $SITEFILE\n";
while(my $line = <SITES>) {
  #motif;object;site
  #EREBP1;Os02g54160.1;gAGCCGCCa
  
  next if($line =~ /^$/ || $line =~ /^#/);
  chomp($line); 
  
  ($motifname,$tfaccession,$site) = split(/;/,$line,3);
  
  push(@{$sites{$tfaccession}{$motifname}},$site);
  $n_of_sites++;
}	
close(TFSEQ); 

print "# read $n_of_sites sites\n"; 

## 3) read references
open(REFS,$REFSFILE) || die "# $0 : cannot read $REFSFILE\n"; 
while(my $line = <REFS>) {

  #Motif;Object;PubMed;FullReference
  #EREBP1;Os02g54160.1  23703395  Serra TS, .. Saibo NJM (2013) OsRMC, a negative regulator .. Plant Mol Biol 82(4-5): 439-455

  next if($line =~ /^$/ || $line =~ /^#/);
  chomp($line); 
  
  ($motifname,$tfaccession,$pubmed,$reference) = split(/;/,$line,4);
  
  push(@{$pubmeds{$tfaccession}{$motifname}},$pubmed);
  push(@{$references{$tfaccession}{$motifname}},$reference);
  $n_of_references++;
}
close(REFS);  

print "# read $n_of_references references\n"; 

## 4) parse PWMS and create library
print "# creating $WLIBRARY...\n";
open(LIB,">$WLIBRARY") || die "# $0 : cannot create $WLIBRARY\n";
print LIB $LIBRARYHEADER;

open(PWM,$PWMFILE) || die "# $0 : cannot read $PWMFILE\n"; 
while(<PWM>) {
  next if(/^$/);
	
  #EREBP1  Os02g54160.1
  #A:	1.00	0.00	0.00	0.50	0.00	0.00	0.00
  #C:	0.00	0.00	1.00	0.50	0.00	0.50	1.00
  #G:	0.00	1.00	0.00	0.00	1.00	0.50	0.00
  #T:	0.00	0.00	0.00	0.00	0.00	0.00	0.00
  
  if(/^([ACGT]):[\t|\s](.*)/) {
    ($bs,$weights) = ($1,$2); 
	
    push(@{$motif{$bs}},split(/\s+/,$weights)); 
    $pwm_width = scalar(@{$motif{$bs}});
		
    if($bs eq 'T') {
      if(scalar(keys(%motif) != 4)){ die "# wrong format ($motifname)\n"; }
		
      $n_of_pwms++;
			
      print LIB "MO  $motifname\n"; 
      print LIB "NA  $tfaccession\n";
      print LIB "P0      A      C      G      T\n";
			
      foreach $p (0 .. $pwm_width-1){
        printf(LIB "%02d",$p+1);
        foreach $nt ('A','C','G','T'){ printf(LIB "%7s",$motif{$nt}[$p]) }
        print LIB "\n"; 
      }	
      
      foreach $p (0 .. scalar(@{$pubmeds{$tfaccession}{$motifname}})-1) {
        printf(LIB "RX  PUBMED:%s\n",$pubmeds{$tfaccession}{$motifname}[$p]);
        printf(LIB "RL  %s\n",$references{$tfaccession}{$motifname}[$p]);
      }	
      print LIB "XX\n";
	
      print LIB "FA  $fullname{$tfaccession}\n" if($fullname{$tfaccession}); 
      print LIB "NA  $synonyms{$tfaccession}\n" if($synonyms{$tfaccession});
      print LIB "SQ  $TFsequence{$tfaccession}\n";
      print LIB "OS  $species{$tfaccession}\n";
      print LIB "CC  family:$family{$tfaccession}\n" if($family{$tfaccession});
      print LIB "XX\n";
      
      if($sites{$tfaccession}{$motifname}) {
        my $site_count = 1;
        foreach $p (0 .. scalar(@{$sites{$tfaccession}{$motifname}})-1)	{
  				printf(LIB "SI  %s_%d\n",$motifname,$site_count);
          printf(LIB "SQ  %s\n",$sites{$tfaccession}{$motifname}[$p]);
          print LIB "XX\n";
          $site_count++;
  	}	
      }
      
      # close this repository entry  
      print LIB "//\n";
    }	
  }
  else	{
    %motif = ();
    ($motifname,$tfaccession) = (split)[0,1];
    if($seen{$motifname}) {
      die "# $0 : motif $motifname found twice in $PWMFILE, exit\n";
    }
    $seen{$motifname}=1;
  }
}	
close(PWM); 
close(LIB);
	
print "# created $WLIBRARY with $n_of_pwms motifs\n";

