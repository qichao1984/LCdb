#!/usr/bin/perl
use strict;
use List::Util qw(shuffle);

my %id2gene;
open( FILE, "id2genemap.txt" ) || die "#1\n";
while (<FILE>) {
  chomp;
  my @items = split( "\t", $_ );
  $id2gene{ $items[0] } = $items[1];
}
close FILE;

my %abundance;
my %samples;
foreach my $file ( glob("*.diamond") ) {
  $file =~ /(.*?)\.assembled.diamond/;
  my $sample = $1;
  $samples{$sample} = 1;
  my %hit;
  open( FILE, "$file" ) || die "#2\n";
  while (<FILE>) {
    chomp;
    my @items = split( "\t", $_ );
    my $gene = $id2gene{ $items[1] };
    if ( !$hit{ $items[0] } ) {
      $abundance{$sample}{$gene}++ if $gene;
      $hit{ $items[0] } = 1;
    }
  }
  close FILE;
}

my %size;
my @sizes;
open( FILE, "sampleinfo.txt" ) || die "#3\n";
while (<FILE>) {
  chomp;
  my @items = split( "\t", $_ );
  $size{ $items[0] } = $items[1];
  push( @sizes, $items[1] );
}
close FILE;
foreach my $sample(keys %samples){
  die "$sample was not found in sampleinfo, please check!\n" if !$size{$sample};
}

my $randomsampling;
@sizes = sort { $a <=> $b } @sizes;
my $rs = $sizes[0] if !defined $randomsampling;
$rs = $randomsampling if defined $randomsampling;
my %abundance_rs = &RandomSampling( \%abundance, \%size, $rs );

my @samples = keys %samples;
open( OUT, ">Lignin degradation _function.txt" ) || die "#2\n";
print OUT "#random sampling: $rs\n";
print OUT "Gene\t", join( "\t", @samples ), "\n";
foreach my $gene ( sort keys %abundance_rs ) {
  print OUT "$gene";
  foreach my $sample (@samples) {
    print OUT "\t$abundance_rs{$gene}{$sample}" if $abundance_rs{$gene}{$sample};
    print OUT "\t0"                             if !$abundance_rs{$gene}{$sample};
  }
  print OUT "\n";
}
close OUT;

sub RandomSampling() {
  my ( $abundance, $size, $rs ) = @_;
  my %abundance = %$abundance;
  my %size      = %$size;
  my %sum;
  foreach my $sample ( keys %abundance ) {
    foreach my $gene ( keys %{ $abundance{$sample} } ) {
      $sum{$sample} += $abundance{$sample}{$gene};
    }
  }
  my %abundance_rs;
  foreach my $sample ( keys %size ) {
    my @array = shuffle( 1 .. $size{$sample} );
    @array = @array[ 0 .. $rs - 1 ];
    @array = grep { $_ <= $sum{$sample} } @array;
    my $i = 1;
    foreach my $gene ( keys %{ $abundance{$sample} } ) {
      my @tmp = grep { $_ >= $i && $_ <= ( $abundance{$sample}{$gene} + $i - 1 ) } @array;
      $abundance_rs{$gene}{$sample} = @tmp;
      $i += $abundance{$sample}{$gene};
    }
  }
  return %abundance_rs;
}
