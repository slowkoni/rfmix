#!/usr/bin/perl -w
use strict;

while(<STDIN>) {
  chomp;
  next if m/^#/;
  next if m/^\s*$/;
  
  my ($chm, $spos, $epos, $sgpos, $egpos, $n_snps, @d) = split/\t/;

  # Take this out if we start using 1-indexed populations in rfmix v2
  for(my $j=0; $j < @d; $j++) {
    $d[$j]++;
  }

  # currently mimicking the output style of the old rfmix, no headers
  # or leading columns, space delimited
  for(my $i=0; $i < $n_snps; $i++) {
    print join(" ",@d),"\n";
  }
}
