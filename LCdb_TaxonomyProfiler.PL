#!/usr/bin/perl
use strict;
use List::Util qw(shuffle);

my %size;
my @sizes;
open(FILE,"/clusterfs/node2/sampleinfo.txt")||die"#3\n";
while(<FILE>){
  chomp;
  my @items = split( "\t", $_ );
  $size{$items[0]} += $items[1];
}
close FILE;
foreach my $sample(keys %size){
  push(@sizes,$size{$sample});
}

my %abundance;
foreach my $sample(keys %size){
  open(KRAKEN,"$sample\.lignin.kraken2.report")||die"#4\n";
  while(<KRAKEN>){
	chomp;
	my @items=split("\t",$_);
	if($items[0]=~/d__(Bacteria|Archaea)(.*?)\|p__([^\|]+)$/){
	  my $phylum=$3;
	  $abundance{"phylum"}{$sample}{$phylum}=$items[1];
	}
	if($items[0]=~/d__(Bacteria|Archaea)(.*?)\|c__([^\|]+)$/){
	  my $class=$3;
	  $abundance{"class"}{$sample}{$class}=$items[1];
	}
	if($items[0]=~/d__(Bacteria|Archaea)(.*?)\|o__([^\|]+)$/){
	  my $order=$3;
	  $abundance{"order"}{$sample}{$order}=$items[1];
	}
	if($items[0]=~/d__(Bacteria|Archaea)(.*?)\|f__([^\|]+)$/){
	  my $family=$3;
	  $abundance{"family"}{$sample}{$family}=$items[1];
	}
	if($items[0]=~/d__(Bacteria|Archaea)(.*?)\|g__([^\|]+)$/){
	  my $genus=$3;
	  $abundance{"genus"}{$sample}{$genus}=$items[1];
	}
	if($items[0]=~/d__(Bacteria|Archaea)(.*?)\|s__([^\|]+)$/){
	  my $species=$3;
	  $abundance{"species"}{$sample}{$species}=$items[1];
	}
  }
  close KRAKEN;
}

@sizes=sort{$a<=>$b}@sizes;
my $rs=$sizes[0];
my %abundance_rs;
foreach my $taxlevel(keys %abundance){
  %{$abundance_rs{$taxlevel}}=&RandomSampling(\%{$abundance{$taxlevel}},\%size,$rs);
}

my @samples = keys %size;
foreach my $taxlevel(keys %abundance_rs){
  open(OUT,">Sulfur_tax_$taxlevel.txt")||die"#5\n";
  print OUT "#random sampling: $rs\n";
  print OUT "Tax\t", join( "\t", @samples ), "\n";
  foreach my $tax ( sort keys %{$abundance_rs{$taxlevel}} ) {
	print OUT "$tax";
	foreach my $sample (@samples) {
	  print OUT "\t$abundance_rs{$taxlevel}{$tax}{$sample}" if $abundance_rs{$taxlevel}{$tax}{$sample};
	  print OUT "\t0"                          if !$abundance_rs{$taxlevel}{$tax}{$sample};
	}
	print OUT "\n";
  }
  close OUT;
}

sub RandomSampling() {
  my ( $abundance, $size, $rs ) = @_;
  my %abundance = %$abundance;
  my %size      = %$size;
  my %sum;
  foreach my $sample(keys %abundance){
	foreach my $gene(keys %{$abundance{$sample}}){
	  $sum{$sample}+=$abundance{$sample}{$gene};
	}
  }
  my %abundance_rs;
  foreach my $sample(keys %size){
	my @array=shuffle (1..$size{$sample});
	@array=@array[0..$rs-1];
	@array=grep{$_<=$sum{$sample}}@array;
	my $i=1;
	foreach my $gene(keys %{$abundance{$sample}}){
	  my @tmp=grep{$_>=$i && $_<=($abundance{$sample}{$gene}+$i-1)} @array;
	  $abundance_rs{$gene}{$sample}=@tmp;
	  $i=$abundance{$sample}{$gene}+1;
	}
  }   
  return %abundance_rs;
}

