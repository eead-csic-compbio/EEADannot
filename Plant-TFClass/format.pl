use strict;

my (%fam,$id,$family);

open(FAM,$ARGV[0]) || die "# cannot read $ARGV[0]\n";
while(<FAM>) {
  chomp;	
  ($id, $family) = split(/\t/,$_); #print "$id -> $family\n";
  $family =~ s/\s+$//;  
  $fam{$id} = $family;
  $fam{$id} =~ s/\s+/___/g; # keep full name in BLAST tabular output 
}
close(FAM);

open(RAWFASTA,$ARGV[1]) || die "# cannot read $ARGV[1]\n";
while(<RAWFASTA>) {
  if(/^>[^|]+\|([^|]+)/) {
    $id = $1;
    print ">$fam{$id}\n";
  } else {
    print
  }
}
close(RAWFASTA);
