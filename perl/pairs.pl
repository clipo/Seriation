#!/usr/bin/perl
use strict;
use Data::Dumper;
use Graph::Directed;
use Graph::Writer::VCG;
use Graph::Writer::Dot;
use Math::Combinatorics;
use Excel::Writer::XLSX;
use Getopt::Long qw(HelpMessage);
use Pod::Usage;
use Time::HiRes;
use Array::Utils qw(:all);
use Statistics::Descriptive;
use Statistics::PointEstimation;
require Term::Screen;
use List::MoreUtils qw/ uniq /;
use GD::Graph::histogram;

my $xyfile;

GetOptions(
    'xyfile=s'                  => \$xyfile,
) or pod2usage(2);

my @xyAssemblages;
my %xAssemblage;
my %yAssemblage;

   open( INFILE, $xyfile ) or die "Cannot open XY File: $xyfile.\n";

   my $typecount;
   while (<INFILE>) {
      chomp;
      my @line = split /\s+/, $_;
      my $label = shift @line;
    if ($label ne "ID") {
        push @xyAssemblages, $label;
        ### note: this is made for UTMs -- which are listed as Northing/Easting. So Y is first -- X is second...
        $yAssemblage{ $label } = $line[0];
        $xAssemblage{ $label } = $line[1];
      }
   }
   
   my @distances;
   
    ## We use the Math::Combinatorics to get all the combinations of 2
   my $assemblagePairs = Math::Combinatorics->new( count => 2, data => [@xyAssemblages] );
   ## Go through all of the combinations
   while ( my @combo = $assemblagePairs->next_combination ) {
      my $distance = sqrt( ($xAssemblage{ $combo[0] } -$xAssemblage{$combo[1]})**2 + ($yAssemblage{$combo[0]} -$yAssemblage{$combo[1]})**2) ;
         push @distances, $distance;
      #print "pairname: $pairname: ", $distance, "\n\r";
   }

my $graph = new GD::Graph::histogram(400,600);
    $graph->set( 
                x_label         => 'Distance',
                y_label         => 'Count',
                title           => "Counts of assemblage pairs",
                x_labels_vertical => 1,
                bar_spacing     => 0,
                shadow_depth    => 1,
                shadowclr       => 'dred',
                histogram_bins => 30,
                transparent     => 0,
            ) 
            or warn $graph->error;
        
    my $gd = $graph->plot(\@distances) or die $graph->error;
    open(IMG, ">$xyfile-dist-histogram.png") or die $!;
    binmode IMG;
    print IMG $gd->png;