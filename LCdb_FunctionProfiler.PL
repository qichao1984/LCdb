#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(shuffle);
use File::Spec;
use File::Basename;

# Define variables for command-line arguments
my ($workdir, $method, $filetype, $seqtype, $sample_info_file, $randomsampling, $output_file);

# Parse command-line arguments
GetOptions(
    'd=s'   => \$workdir,          
    'm=s'   => \$method,           
    'f=s'   => \$filetype,         
    's=s'   => \$seqtype,          
    'si=s'  => \$sample_info_file, 
    'rs=i'  => \$randomsampling, 
    'o=s'   => \$output_file       
) or die usage();

# Check for mandatory arguments
foreach my $arg ($workdir, $method, $filetype, $seqtype, $sample_info_file, $output_file) {
    unless (defined $arg) {
        die "Error: Missing required arguments.\n", usage();
    }
}

# Validate method
my %valid_methods = map { $_ => 1 } qw(diamond usearch blast);
unless (exists $valid_methods{lc $method}) {
    die "Error: Invalid method '$method'. Valid options are: diamond, usearch, blast.\n", usage();
}

# Change to working directory
chdir $workdir or die "Error: Cannot change to directory '$workdir': $!\n";

if (lc $method eq 'diamond') {
    my $lcdb_faa = File::Spec->catfile($workdir, 'LCdb.faa');
    my $diamond_db = 'LCdb'; 
    my $diamond_output_pattern = '%s.assembled.diamond'; 

    # Read sample names 
    my @samples;
    open(my $si_fh, '<', $sample_info_file) || die "Error: Cannot open sample info file '$sample_info_file': $!\n";
    while (my $line = <$si_fh>) {
        chomp $line;
        my @items = split(/\t/, $line);
        next unless @items >= 2;
        push @samples, $items[0];
    }
    close $si_fh;

    unless (@samples) {
        die "Error: No samples found in sample info file '$sample_info_file'.\n";
    }

    print "Found " . scalar(@samples) . " samples in '$sample_info_file'.\n";

    # Check if LCdb.faa exists
    unless (-e $lcdb_faa) {
        die "Error: LCdb.faa file not found in '$workdir'.\n";
    }

    if (-e "$diamond_db.dmnd") {
        print "DIAMOND database already exists.\n";
    } else {
        # Run diamond makedb (only if database doesn't exist)
        print "Running diamond makedb...\n";
        my $makedb_cmd = "diamond makedb --in $lcdb_faa --db $diamond_db";
        system($makedb_cmd) == 0
            or die "Error: Failed to execute '$makedb_cmd'.\n";

        print "diamond makedb completed successfully.\n";
    }

    # Iterate over each sample to run diamond blastp
    foreach my $sample (@samples) {
        my $query_fa = File::Spec->catfile($workdir, "$sample.$filetype");
        my $diamond_output = sprintf($diamond_output_pattern, $sample);

        if (-e $diamond_output) {
            print "Assembled diamond file '$diamond_output' exists for sample '$sample'. Skipping diamond blastp for this sample.\n";
            next;
        }

        unless (-e $query_fa) {
            warn "Warning: Query file '$query_fa' for sample '$sample' does not exist. Skipping this sample.\n";
            next;
        }

        print "Running diamond blastp for sample '$sample'...\n";
        my $blastp_cmd = "diamond blastp --db $diamond_db -q $query_fa -p 48 -o $diamond_output -f 6 -e 1e-5";
        system($blastp_cmd) == 0
            or die "Error: Failed to execute '$blastp_cmd'.\n";

        print "diamond blastp for sample '$sample' completed successfully. Output saved to '$diamond_output'.\n";
    }
}

# Read id2genemap.txt
my %id2gene;
open( FILE, "id2genemap.txt" ) || die "Error: Cannot open id2genemap.txt: $!\n";
while (<FILE>) {
  chomp;
  my @items = split( "\t", $_ );
  $id2gene{ $items[0] } = $items[1];
}
close FILE;
print "Loaded gene mapping from id2genemap.txt.\n";

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
open( FILE, $sample_info_file ) || die "Error: Cannot open sample info file '$sample_info_file': $!\n";
while (<FILE>) {
  chomp;
  my @items = split( "\t", $_ );
  $size{ $items[0] } = $items[1];
  push( @sizes, $items[1] );
}
close FILE;
print "Loaded sample size information from '$sample_info_file'.\n";
foreach my $sample(keys %samples){
  die "$sample was not found in sampleinfo, please check!\n" if !$size{$sample};
}

@sizes = sort { $a <=> $b } @sizes;
my $rs = $sizes[0] if !defined $randomsampling;
$rs = $randomsampling if defined $randomsampling;
my %abundance_rs = &RandomSampling( \%abundance, \%size, $rs );

my @samples = keys %samples;
open( OUT,'>', $output_file ) || die "Error: Cannot write to output file '$output_file': $!\n";
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

sub usage {
    return <<"END_USAGE";
Usage: perl $0 -d <workdir> -m <diamond|usearch|blast> -f <filetype> -s <seqtype> -si <sample size info file> -rs <random sampling size> -o <outfile>

Options:
    -d      <workdir>           Working directory containing input files
    -m      <diamond|usearch|blast>  Method to process (diamond, usearch, or blast)
    -f      <filetype>          File type extension (e.g., fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz)
    -s      <seqtype>           Sequence type (nucl or prot)
    -si     <sample info file>  Path to sample size info file (e.g., sampleinfo.txt)
    -rs     <random sampling size>  Integer specifying random sampling size (optional)
    -o      <outfile>           Path to output file

Example:
    perl $0 -d ./data -m diamond -f fa -s prot -si sampleinfo.txt -rs 1000 -o output.txt

END_USAGE
}
