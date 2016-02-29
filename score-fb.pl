#!/usr/bin/perl -w
use strict;
use cmdline_args;

if (!$ARGV[0] || $ARGV[0] eq "help") {
  print STDERR <<EOF;

  Script to score the forward-backward results of rfmix (v2) for simulated data.

  Usage: -f <rfmix forward-backward output file> -t <true result file> -k <number of subpops>

  The "true result file" is expected to simply list the index number of correct
  subpopulation at each SNP for each sample. Note that the forward-backward output
  file gives results per window rather than per SNP, so this script will
  automatically expand the output internally to compare SNP by SNP with the
  expected correct result.

  
EOF
  exit -1;
}

my ($fb_fname, $truth_fname, $n_subpops, $window_size) = ("", "", 0, 0);
my %args;
$args{'-f'} = [ \$fb_fname,     1 ];
$args{'-t'} = [ \$truth_fname,  1 ];
$args{'-k'} = [ \$n_subpops,    1 ];
$args{'-w'} = [ \$window_size,  1 ];
cmdline_args::get_options(\%args, \@ARGV);

if ($fb_fname eq "") {
  print STDERR "\nSpecify forward-backward results from rfmix (v2) with -f option.\n";
  exit -1;
}
if ($truth_fname eq "") {
  print STDERR "\nSpecify expected correct results file with -t option.\n";
  exit -1;
}
if ($n_subpops < 2) {
  print STDERR "\nNumber of subpops (-k) must be 2 or more.\n";
  exit -1;
}
if ($window_size < 1) {
  print STDERR "\nSpecify number of SNPs per window with -w option.\n";
  exit -1;
}

my @m;
my @d;
for(my $j=0; $j < $n_subpops; $j++) {
  $m[$j] = [];
  for(my $k=0; $k < $n_subpops; $k++) {
    $m[$j]->[$k] = 0.;
  }
  $d[$j] = 1e-7;
}

open F, "<$fb_fname"
  or die "Can't open forward-backward file $fb_fname ($!)";

open T, "<$truth_fname"
  or die "Can't open expected results file $truth_fname ($!)";

# pull off first line header, we don't need it
$_ = <T>;

while(<F>) {
  chomp;
  next if m/^#/;
  next if m/^\s*$/;

  my ($chm, $pos, $gpos, $snp_idx, @rest) = split/\t/;
  my $i = 0;
  my $n_snps = 0;
  while(defined($_ = <T>) && $i < $window_size) {
    chomp;
    next if m/^#/;
    next if m/^\s*$/;

    my ($tchm, $tpos, @truth) = split/\t/;

    for(my $j=0; $j < @rest && $j < @truth; $j++) {
      my @p = split/\ /,$rest[$j];
      my $t = $truth[$j] - 1;
      for(my $k=0; $k < $n_subpops; $k++) {
	$m[$t]->[$k] += $p[$k];
      }
      $d[$t]++;
    }
    
    $i++;
    $n_snps++;
  }
  last if eof(T);
}

close F;
close T;

for(my $j=0; $j < $n_subpops; $j++) {
  for(my $k=0; $k < $n_subpops; $k++) {
    $m[$j]->[$k] /= $d[$j];
  }
}

for(my $j=0; $j < $n_subpops; $j++) {
  printf "%3.1f", $m[0]->[$j] * 100.;
  for(my $k=1; $k < $n_subpops; $k++) {
    printf "\t%3.1f", $m[$k]->[$j] * 100.;
  }
  print "\n";
}
