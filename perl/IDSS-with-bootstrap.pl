#!/usr/bin/perl
use strict;
use Data::Dumper;
use Graph::Directed;
use Graph::Undirected;
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
use Devel::Size qw(size total_size);


my $debug                   = 0;
my $filterflag              = 0; ## do you want to try to output all of the solutions (filtered for non trivial)
my $largestonly             = 0; #  only output the largest set of solutions
my $individualfileoutput    = 0; ## create files for all the indivdual networks
my $bootstrapCI             = 0; ## flag for the CI bootstrap
my $bootstrapSignificance   = 99.5;
my $man                     = 0;
my $help                    = 0;
my $inputfile;
my $threshold               = 0.5;
my $noscreen                = 0;     ## flag for screen output
my $excel                   = 0;       ## flag for excel file output (not implemented yet)
my $xyfile                  = "";
my $mst                     = 0; ## minimum spanning tree
my $pairwiseFile            = "";
my $stats                   = 0; ## output stats, histograms of counts, etc
my $nosum                   = 0;
my $allSolutions            = 0;
my $memusage                = 0; ## memory usage flag
## find the largest valuein a hash
sub largest_value_mem (\%) {
    my $hash   = shift;
    my ($key, @keys) = keys   %$hash;
    my ($big, @vals) = values %$hash;

    for (0 .. $#keys) {
        if ($vals[$_] > $big) {
            $big = $vals[$_];
            $key = $keys[$_];
        }
    }
    $key
}
# process command line options; if none, print a usage/help message.
# note - manual page for explaining the options, what they do, and how to use them
# is at the bottom of this file in simple POD format.

GetOptions(
    'debug'                     => \$debug,
    'bootstrapCI'               => \$bootstrapCI,
    'bootstrapSignificance=f'   => \$bootstrapSignificance, 
    'filtered'                  => \$filterflag,
    'largestonly'               => \$largestonly,
    'indivfiles'                => \$individualfileoutput,
    'help'                      => sub { HelpMessage() },
    'input=s'                   => \$inputfile,
    'excel'                     => \$excel,
    'threshold=f'               => \$threshold,
    'noscreen'                  => \$noscreen,
    'xyfile=s'                  => \$xyfile,
    'pairwise=s'                => \$pairwiseFile,
    'mst'                       => \$mst,
    'stats'                     => \$stats,
    'nosum'                     => \$nosum,
    'allsolutions'              => \$allSolutions,
    'memusage'                  => \$memusage,
    man                         => \$man
) or pod2usage(2);

my $DEBUG = $debug;    # our "$debug level"

if ($DEBUG) {
    print "Verbose debugging output is on!!!\n";
    print "Dont keep track of ALL solutions: $nosum\n";
    print "Processing input file: $inputfile\n";
    print "filterflag: ", $filterflag, "\n";
    print "bootstrapCI: ", $bootstrapCI, "\n";
    print "bootstrapSignificance: ", $bootstrapSignificance, "\n";
    print "largestonly: ", $largestonly, "\n";
    print "individualfileouput: ", $individualfileoutput, "\n";
    print "threshold is currently set to: $threshold\n";
    print "noscreen: ", $noscreen, "\n";
    print "excel:  ", $excel, "\n";
    print "xyfile: ", $xyfile,"\n";
    print "memory usage: ", $memusage, "\n";
    print "pairwise:", $pairwiseFile,"\n";
    print "mst:  ", $mst, "\n";
    print "stats: ", $stats, "\n";
}

# start the clock to track how long this run takes
my $start = Time::HiRes::gettimeofday();

########################################### READ IN THE DATA  ####################################
## The expected data format is <assemblage name><tab><type 1 number><tab><type 2 number>...<return>
## Tab delimited, no headers, raw counts for types (not percentages)
## Example data can be found in the testdata directory

# define some key lists and vars
my @assemblages = ();
my @rowtotals   = ();
my @labels      = ();
my $numrows     = 0;
my $cols        = 0;
my %collections;
my @assemblageNumber      = ();
my @arrayOfSeriations     = ();
my %assemblageFrequencies = ();
my %assemblageValues =      ();
my %assemblageSize = ();
my @allNetworks           = ();
my $maxnumber;
my $count = 0;
my $screen = 1;    ## if the noscreen option is 1, then set the screen flag to 0. Otherwise use the screen
if ($noscreen==1) { $screen=0;}
if ($DEBUG) { $screen = 0; }
my $scr;

## Set up the screen display (default).
## the debug option does not use this since it gets messay
$screen and $scr = new Term::Screen;
$screen and $scr->clrscr();  # clear the screen

## open the data file
my $useOutputFile = substr($inputfile,0,-4);
open( INFILE, $inputfile ) or die "Cannot open $inputfile.\n";
open( OUTFILE, ">$useOutputFile.vna" ) or die "Can't open file $useOutputFile.vna to write.\n";
open( OUTDOTFILE, ">$useOutputFile.dot") or die "Can't open file $useOutputFile.dot to write.\n";
open( OUTPAIRSFILE, ">$useOutputFile-pairs.vna") or die "Can't open file $useOutputFile-pairs.vna to write.\n";
if ($mst) {
    open(OUTBOOTSTRAPFILE, ">$useOutputFile-mst.vna") or die "Can't open file $useOutputFile-mst.vna to write\n";
    open(OUTDISTANCEFILE, ">$useOutputFile-mst-distance.vna") or die "Can't open file $useOutputFile-mst-distance.vna to write\n";
}
my %pairwise = {};
my %pairwiseError = {};

if ($pairwiseFile) {
    open (PAIRWISE,$pairwiseFile ) or die "Cannot open $pairwiseFile.\n";
    while (<PAIRWISE>) {
        my @line = split /\s+/, $_;
        my $pair = $line[0]."#".$line[1];
        $pairwise{ $pair } = $line[2];
        $pairwiseError{ $pair } = $line[3];
    }
}
## some output so we know what we are doing 
$screen and $scr->at(1,1)->puts("Filename:  $inputfile");
$screen and $scr->at(2,1)->puts("Threshold: $threshold");

## Read in the data
# the input of the classes -- note these are absolute counts, not frequencies
# might want to change that...\
$screen and $scr->at(1,40)->puts("STEP: Read in data...");
my $typecount;

while (<INFILE>) {
    #print;
    chomp;
    my @line = split /\s+/, $_;
    my $label = shift @line;
    if ($label) {
        push @labels,           $label;
        push @assemblageNumber, $count;
        $cols = scalar(@line);
        my $rowtotal = 0;
        for (@line) {
            $rowtotal += $_;
        }
        push @rowtotals, $rowtotal;
        

        #print "rowtotal: $rowtotal\n";
        my $freq = [];
        my $values = [];
        $typecount=0;
        for (@line) {
            # push @freq, $_;
            my $in = $_;
            my $f = $in / $rowtotal;
            my $ff = sprintf( "%.4f", $f );
            push @$freq, $ff;
            push @$values, $in;
            $typecount++;
        }
        push @assemblages, [@$freq];
        $assemblageFrequencies{$label} = [@$freq];
        $assemblageValues{ $label } = [@$values];
        $assemblageSize{ $label }= $rowtotal;
        $count++;
    }
    #print "---- row end ----\n";
}
$numrows = scalar(@assemblages);
$DEBUG and print Dumper(\@assemblages),"\n";
$DEBUG and print Dumper( \%assemblageFrequencies ),"\n";

my $maxSeriations = $count;
$screen and $scr->at(3,1);
$screen and $scr->puts("Maximum possible seriation solution length: $count");

#######################################XY File####################################################################
## if you want to characterize or sort by cluster distance, provide a file with X Y coordiantes for each assemblage
##
##
my %xAssemblage;
my %yAssemblage;
my @xyAssemblages;
my $largestX=0;
my $largestY=0;
my %distanceBetweenAssemblages;

if ($xyfile) {
   ## open the xy file
   open( INFILE, $xyfile ) or die "Cannot open XY File: $xyfile.\n";
   $screen and $scr->at(1,40)->puts("STEP: Read in XY data...");
   my $typecount;
   while (<INFILE>) {
      chomp;
    my @line = split /\s+/, $_;
    my $label = shift @line;
    if ($label) {
        push @xyAssemblages, $label;
        ### note: this is made for UTMs -- which are listed as Northing/Easting. So Y is first -- X is second...
        $yAssemblage{ $label } = $line[0];
        $xAssemblage{ $label } = $line[1];
      }
   }
   ## We use the Math::Combinatorics to get all the combinations of 2
   my $assemblagePairs = Math::Combinatorics->new( count => 2, data => [@xyAssemblages] );
   ## Go through all of the combinations
   while ( my @combo = $assemblagePairs->next_combination ) {
      my $pairname      = $combo[0]. "*" .$combo[1];
      my $pairname2      = $combo[1]. "*" .$combo[0];
      my $distance = sqrt( ($xAssemblage{ $combo[0] } -$xAssemblage{$combo[1]})**2 + ($yAssemblage{$combo[0]} -$yAssemblage{$combo[1]})**2) ;
      $distanceBetweenAssemblages{ $pairname }= $distance;
      $distanceBetweenAssemblages{ $pairname2 }= $distance;
      #print "pairname: $pairname: ", $distance, "\n\r";
   }
   $largestX = $xAssemblage{ largest_value_mem(%xAssemblage) };
   $largestY = $yAssemblage{ largest_value_mem(%yAssemblage)};
}



#############################################  THRESHOLD DETERMINATION ####################################
## first compare each assemblage and see if the threshold is exceeded.
## By this I mean the maximum difference between frequencies of any type is greater than what is specified.
## When threshold = 0, all combinations are used. Later the "ends" of solutions are not evaluated if the
## difference between the last assemblage and any other free assemblage is > the threshold.
## This arbitrary setting is to keep from arbitrary solutions being stuck on that come from the "ends"
## of solutions.
##
## Precalculate all of the max differences between types in assembalge pairs. 
$screen and $scr->at(1,40)->puts("STEP: Calculate thresholds....");
my %assemblageComparison = ( "pairname" => 0 );

## We use the Math::Combinatorics to get all the combinations of 2
my $pairs = Math::Combinatorics->new( count => 2, data => [@assemblageNumber] );
## Go through all of the combinations
while ( my @permu = $pairs->next_combination ) {
    #$DEBUG and print $labels[ $permu[0] ] . " * " . $labels[ $permu[1] ] . "\n";
    my $pairname      = $labels[ $permu[0] ] . " * " . $labels[ $permu[1] ];
    my $pairname2     = $labels[ $permu[1] ] . " * " . $labels[ $permu[0] ];
    my $maxDifference = 0;
    for ( my $i = 0 ; $i < $cols ; $i++ ) {
        my $ass1 = $assemblages[ $permu[0] ][$i];
        my $ass2 = $assemblages[ $permu[1] ][$i];
        my $diff = abs( $ass1 - $ass2 ) ;
        if ( $diff > $maxDifference ) {
            $maxDifference = $diff;
        }
    }
    $assemblageComparison{ $pairname } = $maxDifference;
    $assemblageComparison{ $pairname2 } = $maxDifference;
}

$DEBUG and print Dumper( \%assemblageComparison ), "\n";

########################################### BOOTSTRAP SECTION ####################################
#bootrap stuff

srand($start);
my %perror  = ();
my %pvalue  = ();
my $results = 0;
my $ptr1;
my $classes = 0;
my %typeFrequencyLowerCI = ();
my %typeFrequencyUpperCI = ();

# now do ALL the pairwise assemblage comparisons
# go to sleep and come back later.

if ($bootstrapCI) {
    $screen and $scr->at(1,40)->puts("STEP: Bootstrap CIs...        ");
    my $countup=0;
    ## for each assemblage
    $DEBUG and print "Calculating bootstrap confidence intervals\n\r";
    # for each assemblage
    
    foreach my $currentLabel ( sort keys %assemblageFrequencies) {
        my $label = $labels[$countup];
        my @a =  $assemblageFrequencies{ $currentLabel };    
        my $currentAssemblageSize = $assemblageSize{ $currentLabel };
            
        ## create an array of stats objects - one for each type
        my @arrayOfStats = ();
        my $stat;   ## this will be the stat objects
  
        for ( $count = 0 ; $count < $typecount ; $count++ ) {
            push @arrayOfStats, Statistics::PointEstimation->new();
            $arrayOfStats[$count]->set_significance($bootstrapSignificance);
        }

        ## size of bootstrapping (how many assemblages to create)
        my $loop = 1000;
        ## bootstrap the assemblages -- $loop times for the statistics
       
        while ($loop) {
            my $assemsize = $currentAssemblageSize;
            my @new_assemblage   = ();
            my $class;
            
            my ( @cumulate, $count );
            $classes = $typecount;
            my $index = 0;
            my $total = 0.0;
            $count   = 0;
            
            ## now count through the classes and set up the frequencies  
            for ( $count = 0 ; $count < $classes ; $count++ ) {
                $cumulate[$index] = $a[0][$count];  ## this is the index of the frequency for this class
                $total += $a[0][$count];            ## should ultimate equal 100
                $index++;                        ## index should be total # of types at end
            }
            ## now build the assemblages of the same size.
            my $rand;
            while ($assemsize) {
                $rand  = rand;              ## random number from 0-1
                $class = 0;
                #print "Got a $rand\n\r";
                ### continue while we 
                while (( $class < $index ) && ( $rand > $cumulate[$class] ) ) {
                    $rand -= $cumulate[$class];
                    $class++;
                }
                #print "This goes in class $class\n\r";
                push @new_assemblage, $class;
                $assemsize--;
            }
            ## this should result in a new assemblage of the same size
            my ( @ahat, %aholder, %bholder );
            %aholder = ();
            
            ## initialize arrauy
            my $indexN = 0;
            for ( $indexN = 0 ; $indexN < $classes ; $indexN++ ) {
                $ahat[$indexN] = 0;
            }
            ## count up the classes
            for (@new_assemblage) {
                $aholder{$_}++;
            }
            
            my $classCount=0;
            foreach my $stat (sort @arrayOfStats) {
                my $results = $aholder{ $classCount };
                $stat->add_data($results / $currentAssemblageSize);
                $classCount++;
             }
            $loop--;
        }
        my @lowerCI;
        my @upperCI;
        foreach my $stat (sort @arrayOfStats) {
            push @lowerCI, $stat->lower_clm();
            push @upperCI, $stat->upper_clm();
        }
        $typeFrequencyLowerCI{ $label } = \@lowerCI;
        $typeFrequencyUpperCI{ $label } = \@upperCI;
        $results = 0;
        $countup++;
    }    
}

############## pre calcuate the valid pairs of assemblages and stuff into hash ############################
my %validComparisonsArray;
my @validComparisonAssemblages;

foreach my $label ( @labels ) {
    my @cAssemblages;
    foreach my $comparativeLabel (@labels ) {
        if ($assemblageComparison{ $label . " * " . $comparativeLabel }  <= $threshold && $comparativeLabel ne $label) {
            push @cAssemblages, $comparativeLabel;
        }
    }
    $validComparisonsArray{ $label } = [@cAssemblages];
}    

#print Dumper(\%validComparisonAssemblages);
#print Dumper([keys %validComparisonAssemblages]);
    

########################################### FIND ALL THE VALID TRIPLES  ####################################
my $directionstate;
$numrows = scalar(@assemblages);
my @triples;
my @triplettype;
my @tripletNames = ();
my @tripletArray = ();
my $net          = Graph::Undirected->new;
my $comparison12;
my $comparison23;
my $error;
my $numberOfTriplets = 0;

## This uses the Math::Combinatorics to create all the permutations of 3. This is the simplest solution programmatically
## Why recreate the wheel? Use Perl!
$screen and $scr->at(1,40)->puts("STEP: Find valid triples....      ");
my $permutations =Math::Combinatorics->new( count => 3, data => [@assemblageNumber] );

while ( my @permu = $permutations->next_combination ) {
    my $tripletname = $labels[ $permu[0] ] . " * ". $labels[ $permu[1] ] . " * ". $labels[ $permu[2] ];
    $comparison12 = "";
    $comparison23 = "";
    $error        = 0;
    for ( my $i = 0 ; $i < $cols ; $i++ ) {
        my $difscore;
        my $difscore2;

        my $ass1 = $assemblages[ $permu[0] ][$i];
        my $ass2 = $assemblages[ $permu[1] ][$i];
        my $ass3 = $assemblages[ $permu[2] ][$i];
         
        ## first compare assemblages 1 and 2
        if ($bootstrapCI ) {
            my $upperCI_1 = $typeFrequencyUpperCI{ $labels[ $permu[0] ] }->[$i];
            my $lowerCI_1 = $typeFrequencyUpperCI{ $labels[ $permu[0] ] }->[$i];
            my $upperCI_2 = $typeFrequencyUpperCI{ $labels[ $permu[1] ] }->[$i];
            my $lowerCI_2 = $typeFrequencyUpperCI{ $labels[ $permu[1] ] }->[$i];
            my $dif1 = $ass1 - $ass2;
            if ( $upperCI_1 < $lowerCI_2 )  {
                $difscore = -1;
            } elsif ($lowerCI_1 > $upperCI_2 ) {
                $difscore = 1;
            } else {
                $difscore = 0;
            }
        } else {   ### if the bootstrapCI is not being used
            my $dif1 = $ass1 - $ass2;
            if ( $dif1 < 0 )  { $difscore = -1; }
            if ( $dif1 > 0 )  { $difscore = 1;  }
            if ( $dif1 == 0 ) { $difscore = 0;  }
        }

        ## now compare assemblages 2 and 3
        if ($bootstrapCI ) {  ## boostrap confidence intervals
            my $upperCI_2 = $typeFrequencyUpperCI{ $labels[ $permu[1] ] }->[$i];
            my $lowerCI_2 = $typeFrequencyUpperCI{ $labels[ $permu[1] ] }->[$i];
            my $upperCI_3 = $typeFrequencyUpperCI{ $labels[ $permu[2] ] }->[$i];
            my $lowerCI_3 = $typeFrequencyUpperCI{ $labels[ $permu[2] ] }->[$i];
            if ( $upperCI_2 < $lowerCI_3 )  {
                $difscore = -1;
            } elsif ($lowerCI_2 > $upperCI_3 ) {
                $difscore = 1;
            } else {
                $difscore = 0;
            }
        } else {            ### if the bootstrapCI is not being used
            my $dif2 = $ass3 - $ass2;
            if ( $dif2 < 0 )  { $difscore2 = -1; }
            if ( $dif2 > 0 )  { $difscore2 = 1;  }
            if ( $dif2 == 0 ) { $difscore2 = 0;  }
        }
        
        
        if ( ( $difscore == 1 ) && ( $difscore2 == 1 ) ) {      ## F1 > F2 < F3 ## criteria not met
            $error++;
        }
        elsif ( ( $difscore == 1 ) && ( $difscore2 == -1 ) ) {   ## F1 > F2 > F3 OK
            $comparison12 .= "U";
            $comparison23 .= "D";
        }
        elsif ( ( $difscore == -1 ) && ( $difscore2 == -1 ) ) {   #  F1 < F2 >F3 OK
            $comparison12 .= "X";
            $comparison23 .= "X";
        }
        elsif ( ( $difscore == -1 ) && ( $difscore2 == 1 ) ) {    # F1 < F2 < F3
            $comparison12 .= "D";
            $comparison23 .= "U";
        }
        elsif ( ( $difscore == 0 ) && ( $difscore2 == 1 ) ) {
            $comparison12 .= "M";
            $comparison23 .= "U";
        }
        elsif ( ( $difscore2 == 0 ) && ( $difscore == 1 ) ) {
            $comparison12 .= "U";
            $comparison23 .= "M";
        }
        elsif ( ( $difscore == 0 ) && ( $difscore2 == -1 ) ) {
            $comparison12 .= "M";
            $comparison23 .= "D";
        }
        elsif ( ( $difscore == -1 ) && ( $difscore2 == 0 ) ) {
            $comparison12 .= "D";
            $comparison23 .= "M";
        }
        elsif ( ( $difscore == 0 ) && ( $difscore2 == 0 ) ) {
            $comparison12 .= "M";
            $comparison23 .= "M";
        }
        else {
            $error++;
            print "\n\rNo match to our possibility of combinations. Difscore 1: $difscore Difscore 2: $difscore2 \n\r";
            print "I must quit. Debugging required.\n\r";
            exit();
        }
    }

    if ( $error == 0 ) {
        undef $net;
        $net = Graph::Undirected->new();
        $net->set_graph_attribute("GraphID", $numberOfTriplets);
        $net->add_vertex( $labels[ $permu[0] ] );
        $net->set_vertex_attribute($labels[ $permu[0] ] ,"End",1);
        $net->add_vertex( $labels[ $permu[1] ] );
        $net->set_vertex_attribute($labels[ $permu[1] ] ,"End",0);
        $net->add_vertex( $labels[ $permu[2] ] );
        $net->set_vertex_attribute($labels[ $permu[0] ] ,"End",1);
        
        $net->add_weighted_edge(
            $labels[ $permu[1] ],
            $labels[ $permu[0] ],
            $comparison12
        );
        $net->set_edge_attribute($labels[ $permu[1]], $labels[ $permu[0] ] , "GraphID", $numberOfTriplets);
        $net->set_edge_attribute($labels[ $permu[1]], $labels[ $permu[0] ] , "End", 1);
        $net->add_weighted_edge(
            $labels[ $permu[1] ],
            $labels[ $permu[2] ],
            $comparison23
        );
        $net->set_edge_attribute($labels[ $permu[1]], $labels[ $permu[2] ] , "GraphID", $numberOfTriplets);
        $net->set_edge_attribute($labels[ $permu[1]], $labels[ $permu[2] ]  , "End", 1);
        ## set the names to the end...
        $net->set_graph_attribute("End_1",$labels[ $permu[0]]);
        $net->set_graph_attribute("End_2",$labels[ $permu[2]]);
        $DEBUG and print "VALID SOLUTION: " . $labels[ $permu[0] ] . " * " . $labels[ $permu[1] ] . " * " . $labels[ $permu[2] ] . "\n";
        $DEBUG and print "VALID SOLUTION: \t  $comparison12\t  ---   $comparison23\n";
        push @triples, $net;
        $numberOfTriplets++;
    }
    $error = 0;

}
$DEBUG and print "------- ALL VALID TRIPLES ---------\n";
$DEBUG and print Dumper(\@triples),"\n";
$DEBUG and print "------- ALL VALID TRIPLES ---------\n";

########################################### THIS IS THE MAIN SERIATION SORTING SECTION ####################################
## now go through the combination of triplets that are valid (i.e., allthose that have no Z in them)
my @array = ();
## now go through all of the assemblages and each one to the top and bottom of the triplets to see if they fit.
## if they do, add them to the network

$screen and $scr->at(1,40)->puts("STEP: Main seriation sorting... ");

#my $currentMaxSeriations = 3;
my $currentMaxSeriationSize = 4;

my $maxEdges  = 2;      ## keeps track of the current largest size of solution
my $stepcount = 0;     ## keeps track of the current step (from 0, adds one each step)
my $match     = 0;      ## keeps track of whether a solution is found for each step. Once no new solution is found, the search is over!

my @networks;  ## array of solutions from previous step (starts out with the 3s) (directed and deep graphs)
my @newnets;   ## array of new solutions (directed and deep graphs)
my @solutions;  ## Array of all solutions (undirected and shallow copies of graphs)
my $solutionCount=scalar(@triples);   ## an index of the current number of solutions

### This is just preparing the intial set of threes into the array that will be used as the seed for all future
### seriation sets.

while ( $currentMaxSeriationSize <= $maxSeriations ) {
    #my $index    = 0;
    ### first time through copy the triples...
    if ($currentMaxSeriationSize==4) {
        @solutions = @{ \@triples };
        @networks = @{ \@triples };
    }
    $stepcount++;
    $DEBUG and print "__________________________________________________________________________________________\n";
    $DEBUG and print "Step number:  $currentMaxSeriationSize\n";
    $DEBUG and print "__________________________________________________________________________________________\n";
    $screen and $scr->at(4,1)->puts("Step number:                     ");
    $screen and $scr->at(4,1)->puts("Step number:  $currentMaxSeriationSize ");
    my $netnum=scalar(@networks);
    $DEBUG and print "Number of solutions from previous step: $netnum\n";
    $screen and $scr->at(5,1)->puts("Number of solutions from previous step:         ");
    $screen and $scr->at(5,1)->puts("Number of solutions from previous step: $netnum");
    $match = 0;      ## set the current match to zero for this step (sees if there are any new solutions for step)
    ## look through the set of existing valid networks.
    foreach my $nnetwork (@networks) {
        $DEBUG and print "-----------------------------------------------------------------------------------\n";
        $DEBUG and print "Network: ", $nnetwork, "\n";
        $DEBUG and print "-----------------------------------------------------------------------------------\n";
        ## find the ends
        ## given the ends, find the valid set of assemblages that can be potentially added
        ## this list is all assemblages meet the threshold requirements
        #print Dumper($nnetwork);
        my $whichEnd=0;

        foreach my $endAssemblage ($nnetwork->get_graph_attribute("End_1"), $nnetwork->get_graph_attribute("End_2")) {
            $whichEnd++; ## either a 1 or a 2
            ##print "The end of assemblages of $nnetwork are: ". $nnetwork->get_graph_attribute("End_1") . "  and ". $nnetwork->get_graph_attribute("End_2"), "\n";
            foreach my $testAssemblage ( @ { $validComparisonsArray{ $endAssemblage } }) {
                
                $DEBUG  and print "\t\tChecking assemblage: ", $testAssemblage, " to see if it fits on the end of the current solution.\n";
                $DEBUG  and print "\t\tFirst check to see if it is included already. If it has, move on.\n";
                my @vertices = $nnetwork->vertices;
                
                if ( (! grep { $_ eq $testAssemblage} @vertices) ) {   ## if the assemblage is NOT in the list of existing vertices.
                
                    # get the exterior vertices (should be 2)
                    $DEBUG  and print "\t\tFind the ends of the network. Do this by getting all the vertices \n";
                    $DEBUG  and print " \t\tand looking for the ones with only 1 connection. There should be just 2 here.\n";
                    ## loop through all of the edges to see if they can be stuck on the ends of the networks.
                    
                    my @newassemblage = ();
                    my @oldassemblage = ();
                    my $comparisonMap;

                    $DEBUG and print "\t\t", $endAssemblage, " is on the edge since it only has one vertice.\n";
                    
                    #################### THRESHOLDING STUFF #################################### 
                    #first determine if the pairs are within the threshold value (0 = all assemblages)
                    #my $pairname = $endAssemblage. " * " . $testAssemblage;
                    #my $diff     = $assemblageComparison{$pairname};
                    #$DEBUG and print "\t\t\tFor $pairname the max frequency difference is $diff.\n";
                    #if ($diff eq undef) {
                    #       print "\n\rError: pairname: $pairname not found in hash lookup!\n\r";
                    #       print "Probably a problem with the names of the assemblages. Check for weird characters. Exiting.\n\r";
                    #       exit();
                    #}
                    ## go through process only if threshold is 0 or difference value is below threshold
                    ## this should mean that network will not grow unless the above conditions are met.
                    my $error = 0;
                    ##if (  ($threshold>0 ) and ($diff > $threshold ) )  {
                    ##   $DEBUG and print "\t\t\tThreshold = $threshold and Diff = $diff. Since $diff < $threshold, continue.\n";
                    ##   $error++; # this should ensure future failure....
                    ##    next;
                    ##}
                    ############################################################################
                    
                    @newassemblage = @{ $assemblageFrequencies{ $testAssemblage } };
                    @oldassemblage = @{ $assemblageFrequencies{ $endAssemblage } };
                    
                    #### FIND INNER EDGE RELATIVE TO THE EXISTING END ASSEMBLAGE ##############
                    my @neighbors = $nnetwork->neighbors($endAssemblage);
                    if (scalar(@neighbors) > 1 || scalar(@neighbors)==0){
                        print "\r\n\r\n\r\nThere are too many or two few neighbors (should only be 1!). Error!\n\r";
                        print "\r\nWe are testing $endAssemblage and got ", Dumper(\@neighbors);
                        print Dumper($nnetwork);
                        print $nnetwork;
                        exit();
                    }
                    $DEBUG and print "\t\t\t The number of neighbors at $endAssemblage is ", scalar(@neighbors), " (should be just one).\n";
                    my $g = $nnetwork->get_edge_weight( $neighbors[0], $endAssemblage );
                    $DEBUG and print "\t\t\tThere should be just 1 neighbor to $endAssemblage and that is: $neighbors[0]\n";
                    $DEBUG and print "\t\t\t\t it has a relation of $g\n";
                    my $outerEdge= $endAssemblage;
                    my $innerEdge= $neighbors[0];
                    ##########################################################################
                    
                    
                    my @comparison = split //, $g;
                    my $cols = scalar(@newassemblage);
                    $comparisonMap ="";

                    # go through the columns
                    if (!$error) {
                           for ( my $i = 0 ; $i < $cols ; $i++ ) {
                                my ( $difscore, $difscore2 );
                                my $val1 = $newassemblage[$i];
                                my $val2 = $oldassemblage[$i];
                                $DEBUG and print "\t\tComparing Assemblage: $testAssemblage    and    Assemblage: $endAssemblage            ########\n ";
                                $DEBUG and print "\t\t\t\tType $i - Type $i - Type $i - Type $i - Type $i - Type $i - Type $i  ########  \n";
                                $DEBUG and print "\t\t\t\tType $i:  $testAssemblage 1: ", $newassemblage[$i], " $endAssemblage 2: ", $oldassemblage[$i], "\n";
                                 ## ALL COMBINAATIONS
                                   #   dif	comparison	result	comparisonMap
                                   #   1	U	      okay	U
                                   #   1	M	      okay	U
                                   #   1	X	      bad	--
                                   #   1	D	      bad	--
                                   #   0	U	      okay	U
                                   #   0	M	      okay	M
                                   #   0	D	      okay	D
                                   #   0	X	      okay	X
                                   #   -1	U	      bad	--
                                   #   -1	M	      okay	M
                                   #   -1	D	      okay	D
                                   #   -1	X	      okay	D
                                if ($bootstrapCI ) {
                                    my $upperCI_test = $typeFrequencyUpperCI{ $testAssemblage  }->[$i];
                                    my $lowerCI_test = $typeFrequencyUpperCI{ $testAssemblage  }->[$i];
                                    my $upperCI_end = $typeFrequencyUpperCI{ $endAssemblage }->[$i];
                                    my $lowerCI_end = $typeFrequencyUpperCI{ $endAssemblage }->[$i];
                                    if ( $upperCI_test < $lowerCI_end )  {
                                        $difscore = -1;
                                    } elsif ($lowerCI_test > $upperCI_end ) {
                                        $difscore = 1;
                                    } else {
                                        $difscore = 0;
                                    }
                                } else {
                                    my $dif1 = $newassemblage[$i] - $oldassemblage[$i];
                                    if ( $dif1 < 0 )  { $difscore = -1; }
                                    if ( $dif1 > 0 )  { $difscore = 1; }
                                    if ( $dif1 == 0 ) { $difscore = 0; }
                                }
                                $DEBUG and print "\t\t\t\tType $i: - comparison is:  ", $comparison[$i], " a score of: ", $difscore, "\n";
                                 if (   ( $difscore == 1 ) && ( $comparison[$i] =~ "U" ) ) {                                                 #### 1 U
                                    $comparisonMap .= "U";
                                    $DEBUG and print  "\t\t\t\tType $i: Got a difscore of 1 and a comparison of a U. This works. \n";
                                    $DEBUG and print " \t\t\t\tAdding $testAssemblage to vertices $endAssemblage\n";
                                    
                                } elsif (( $difscore == 1 ) && ( $comparison[$i] =~ "M" ) ) {                                                ### 1 M
                                    # this is okay - its a match and the new value is greater. New value shoudl be U
                                    # need to find what was happening on the previous comparison to know whether this needs
                                    # to be up or down.
                                    $DEBUG and print "\t\t\t\tType $i: Got a difscore of 1 and a comparison of a M. This could be okay.\n";
                                    my $xerror        =     0;
                                    $DEBUG and print "\t\t\t\tType $i:   Matching case A (1, M)  \n";
                                    $DEBUG and print "\t\t\t\t\ This will only work if there no Xs anywhere previously OR if the opposite end doesnt ALSO go up!\n";
                                    my $stopFlag      =     0;   ## use this flag to determine if one needs to keep checking through the pairs.
                                    my ($outwardEdge, $inwardEdge);
                                    my $old_inner= $outerEdge;
                                    $DEBUG and print "\t\t\t\t\t ", $nnetwork;
                                    my @currentEdges  = $nnetwork->edges05($innerEdge);
                                    my $ccount;
                                    foreach my $checkEdge (@currentEdges) {     ### no need to go in order -- jsut look at all the other edges to see if there is an X
                                       $DEBUG and print "\t\t\t\t\ Now on @$checkEdge[0] @$checkEdge[1] \n";
                                       $inwardEdge= @$checkEdge[0];
                                       $outwardEdge = @$checkEdge[1];
                                       my $comparison = $nnetwork->get_edge_weight( $outwardEdge, $inwardEdge) || $nnetwork->get_edge_weight( $inwardEdge, $outwardEdge );
                                       my @compArray = split //, $comparison;
                                       if (!$compArray[$i]) {
                                          print "Comparison is empty. Error! Stopping.\n\r\n\r";
                                          exit();
                                       }
                                       $DEBUG and print "\t\t\t\tType $i: Here is what we get for comparison # $ccount \n";  ## note that i is the current type
                                       $DEBUG and print " \t\t\t\t\t $inwardEdge - $outwardEdge: ", $comparison, "->", $compArray[$i], "\n";
                                       if ($compArray[$i] =~ "X" || $compArray[$i] =~ "U" )  {
                                             $xerror = 1;  ### BLARGH a previous X or an UP ! This will not be tolerated!
                                             $stopFlag = 1; ## We can stop
                                       }
                                       $DEBUG and print "\t\t\t\t\t Since I got $compArray[$i] my potential new value is still X.\n";
                                       $DEBUG and print "\t\t\t\t\t Now going to get the next pair of assembalges to examine in the chain\n";
                                       $ccount++;
                                    }
                                    if ( $xerror ) {
                                        $error++;
                                        last;
                                    } else {
                                        $comparisonMap .= "U";
                                        $DEBUG and print "\t\t\t\t Type $i: For this type, OK to add $testAssemblage to vertices $endAssemblage \n";
                                        $DEBUG and print "\t\t\t\t\t No previous X values anywhere. \n";
                                        $DEBUG and print "\t\t\t\t Type $i: Adding an U to the comparisons for type $i. \n";
                                        $DEBUG and print "\t\t\t\t\t Comparison map is now $comparisonMap\n";
                                    }
                                }
                                elsif
                                  ## error the new value is greater but shoudl be less. Error!
                                  ( ( $difscore == 1 ) &  ( $comparison[$i] =~ "D" ) ) {                                                  ## 1 D
                                    #print "mismatch!\n";
                                    $error++;
                                    next;
                                    $DEBUG and print "\t\t\t\tType $i: Value 1:  ", $newassemblage[$i], " value 2: ", $oldassemblage[$i], "\n";
                                    $DEBUG and print "\t\t\t\tType $i: Comparison is:  ", $comparison[$i], " a score of: ", $difscore, "\n";
                                    $DEBUG and print "\t\t\t\tType $i: Rejecting $testAssemblage from $endAssemblage \n";
                                    $DEBUG and print "\t\t\t\t\t because value is 1 and comparison is D.\n";
                                }
                                elsif
                                  ## new score is less and the Comparison is up .. Error!
                                  ## note -- need to check to see if there is a previous change in direction because
                                  ## its possible that the this is a mode change.
                                  ## do this by logging all modes in the original triplet -- keep track of modes
                                  ## per type
                                  (    ( $difscore == -1 ) && ( $comparison[$i] =~ "U" ) )                                               # -1 U
                                {
                                    ## first check to see if there is already and X in this column somewhere else.
                                    my $xerror   = 0;
                                    $DEBUG and print "\t\t\t\tType $i:  Case B (-1, U). Potentially can add $testAssemblage and vert $endAssemblage\n ";
                                    $DEBUG and print "\t\t\t\tType $i:  But need to check the rest of the chain for X's (can't already be an X). \n";
                                    my $stopFlag      =   0;   ## use this flag to determine if one needs to keep checking through the pairs.
                                    my ($outwardEdge, $inwardEdge);
                                    my @currentEdges  = $nnetwork->edges05($innerEdge);
                                    my $ccount;
                                    foreach my $checkEdge (@currentEdges) {     ### no need to go in order -- just look at all the other edges to see if there is an X
                                       $DEBUG and print "\t\t\t\t\ Now on @$checkEdge[0] @$checkEdge[1] \n";
                                       $inwardEdge= @$checkEdge[0];
                                       $outwardEdge = @$checkEdge[1];
                                       my $comparison = $nnetwork->get_edge_weight( $outwardEdge, $inwardEdge) || $nnetwork->get_edge_weight( $inwardEdge, $outwardEdge );
                                       my @compArray = split //, $comparison;
                                       if (!$compArray[$i]) {
                                          print "Comparison is empty. Error! Stopping.\n\r\n\r";
                                          exit();
                                       }
                                       $DEBUG and print "\t\t\t\tType $i: Here is what we get for comparison # $ccount \n";  ## note that i is the current type
                                       $DEBUG and print " \t\t\t\t\t $inwardEdge - $outwardEdge: ", $comparison, "->", $compArray[$i], "\n";
                                       if ($compArray[$i] =~ "X")  {
                                             $xerror += 1;  ### BLARGH a previous X! This will not be tolerated!
                                             $stopFlag = 1; ## We can stop
                                       }
                                    }
                                    if ($xerror > 0) {
                                        $error += 1;
                                        last;
                                        $DEBUG and print "\t\t\t\tType $i: Rejecting $testAssemblage from $endAssemblage) because there was an X \n";
                                        $DEBUG and print "\t\t\t\t\t  This would make it multimodal - so error.\n";
                                    } else {
                                        $comparisonMap .= "X";   ## this is an X unless there is an error....
                                        $DEBUG and print "\t\t\t\tType $i:Definitely OK to add $testAssemblage to vertices $endAssemblage because score \n";
                                        $DEBUG and print "\t\t\t\t\tis -1 and the comparison is U but no other Xs in the previous linkages.\n";
                                        $DEBUG and print "\t\t\t\tType $i: Adding an X to the comparisons for type $i. \n";
                                        $DEBUG and print "\t\t\t\t\tComparison map is now $comparisonMap\n";
                                    }    #end if if check error (xerror)
                                }
                                elsif ## new score is less and the comparison is down .. Yes!
                                  (    ( $difscore == -1 ) && ( $comparison[$i] =~ "D" ) )                                          ## -1   D
                                {
                                    $comparisonMap .= "D";
                                    $DEBUG and print "\t\t\t\tType $i: Adding a D to the comparisons for type $i. Comparison map is now $comparisonMap\n";
                                
                                } elsif  (    ( $difscore == -1 )  || ( $comparison[$i] =~ "M" ) ) {                               ## -1 M
                                
                                    # new score is less but comparison is Match. Okay
                                    #my $vHere    = $endAssemblage;
                                    my $xerror   = 0;  ## count for errors
                                    $DEBUG and print "\t\t\t\t #### For type $i we have a matching Case C (-1, M) \n";
                                    $DEBUG and print "\t\t\t\t\t We can potentially add $testAssemblage and vert $endAssemblage but need to check further\n";
                                    $DEBUG and print "\t\t\t\t\t because score is -1 and the comparison is M.\n";
                                    $DEBUG and print" \t\t\t\t\t Could be X or U or M or D\n";
                                    my $change        =   "M";  ## the new comparison variable to use
                                    my $stopFlag      =   0;   ## use this flag to determine if one needs to keep checking through the pairs.
                                    ## now get the next set of comparisons
                                    $DEBUG and print "\t\t\t\t ", $nnetwork;
                                    my @currentEdges  = $nnetwork->edges05($innerEdge);
                                    my $ccount;
                                    my ($outwardEdge, $inwardEdge);
                                    foreach my $checkEdge (@currentEdges) {     ### no need to go in order -- jsut look at all the other edges to see if there is an X
                                       $DEBUG and print "\t\t\t\t\ Now on @$checkEdge[0] @$checkEdge[1] \n";
                                       $inwardEdge= @$checkEdge[0];
                                       $outwardEdge = @$checkEdge[1];
                                       my $comparison = $nnetwork->get_edge_weight( $outwardEdge, $inwardEdge) || $nnetwork->get_edge_weight( $inwardEdge, $outwardEdge );
                                       my @compArray = split //, $comparison;
                                       if (!$compArray[$i]) {
                                          print "Comparison is empty. Error! Stopping.\n\r\n\r";
                                          exit();
                                       }
                                       $DEBUG and print "\t\t\t\tType $i: Here is what we get for comparison # $ccount \n";  ## note that i is the current type
                                       $DEBUG and print " \t\t\t\t\t $inwardEdge - $outwardEdge: ", $comparison, "->", $compArray[$i], "\n";
                                       if ( $compArray[$i] =~ "U" ) {
                                             $change = "X";
                                             $stopFlag=1; ## we can stop
                                       } elsif (($compArray[$i] =~ "X") or ($compArray[$i] =~ "D"))  {
                                             $change = "D";
                                             $stopFlag = 1; ## We can stop
                                       } elsif ($compArray[$i] =~ "M") {
                                             $change = "M";
                                             ## in this case we have to keep going
                                       } else {
                                          print "ERROR: missing value -- comparison is $compArray[$i]. Must have a value.\n\r\n\r";
                                          exit();
                                       }
                                       $DEBUG and print "\t\t\t\t\t Since I got $compArray[$i] my potential new value is $change.\n";
                                       $DEBUG and print "\t\t\t\t\t Now going to get the next pair of assembalges to examine in the chain\n";
                                       
                                       $ccount++;
                                    }
                                    ## in this case I dont think there are any errors possible. types can always go down from any other value
                                    $comparisonMap = $comparisonMap . $change;      ## use the value from above.         
                                    $DEBUG and print "\t\t\t\tType $i: OK to add $testAssemblage to vertices $endAssemblage because \n";
                                    $DEBUG and print "\t\t\t\t score is -1 and the comparison is D. ComparisonMap is now $comparisonMap\n";
                                    if (!$comparisonMap ) {
                                       print "\n\rERROR: comparisonMap can't be empty. Bug here. \n\r\n\r";
                                       exit();
                                    }
                                } elsif  # new score is match but comparison is Match. Okay
                                  ( ( $difscore == 0 ) && ($comparison[$i] =~ "U") ) {                                              ## 0  U
                                    $comparisonMap .= "U";
                                    $DEBUG and print "\t\t\t\tType $i:  Ok to add  $testAssemblage to vertices $endAssemblage because its a match. \n";
                                    $DEBUG and print "\t\t\t\tType $i: ComparisonMap is now $comparisonMap\n";
                                } elsif  # new score is match but comparison is Match. Okay
                                  ( ( $difscore == 0 ) && ($comparison[$i] =~ "D") ) {                                              ## 0 D
                                    $comparisonMap .= "D";
                                    $DEBUG and print "\t\t\t\tType $i: Ok to add $testAssemblage to vertices $endAssemblage because its a match. \n";
                                    $DEBUG and print "\t\t\t\tType $i: ComparisonMap is now $comparisonMap\n";
                                } elsif  # new score is match but comparison is Match. Okay
                                  ( ( $difscore == 0 ) && ($comparison[$i] =~ "M") ) {                                              ## 0 M
                                    $comparisonMap .= "M";
                                    $DEBUG and print "\t\t\t\tType $i:  Ok to add  $testAssemblage to vertices $endAssemblage because its a match. \n";
                                    $DEBUG and print "\t\t\t\tType $i: ComparisonMap is now $comparisonMap\n";
                                } elsif # newscore is down but comparison is X. This means that there was already a peak
                                  (    ( $difscore == -1 ) && ( $comparison[$i] =~ "X" ) )                                          ## -1 X
                                {
                                    ## this is okay since it is down from a mode peak
                                    $DEBUG and print "\t\t\t\tType $i:  Ok to add  $testAssemblage to vertices $endAssemblage because \n";
                                    $DEBUG and print " \t\t\t\tscore is -1 and the comparison is D. ComparisonMap is now $comparisonMap\n";
                                    $comparisonMap .= "D";
                                } elsif (( $difscore == 1 ) && ( $comparison[$i] =~ "X" ) )                                         ## 1  X
                                {
                                    ## new score is up but comparison is X.. no cant work because past peak
                                    $error++;
                                    last;
                                    $DEBUG and print "\t\t\t\tType $i: Rejecting $testAssemblage from $endAssemblage]. We can't go up \n";
                                    $DEBUG and print " \t\t\t\tafter a peak. so error. Error now $error\n";
                                } elsif # newscore is down but comparison is X. This means that there was already a peak
                                  (    ( $difscore == 0 ) && ( $comparison[$i] =~ "X" ) )                                           ## 0  X
                                {
                                    ## this is okay since it is down from a mode peak
                                    $comparisonMap .= "X";
                                    $DEBUG and print "\t\t\t\tType $i:  Ok to add  $testAssemblage to vertices $endAssemblage because score \n";
                                    $DEBUG and print "\t\t\t\t is 0 and the comparison is X. ComparisonMap is now $comparisonMap\n";
                                } else {
                                    print "\n\r\t\t\t\tERROR!!!! Not found match combination! MUST FIX! Some combination\n\r ";
                                    print "\t\t\t\t is not being caught correctly... Exiting.\n\r";
                                    print "\t\t\t\tHere is the score of the differences in  for Type $i: $difscore\n\r";
                                    print "\t\t\t\tHere is the comparison value: ", $comparison[$i], "\n\r";
                                    $error++;
                                    exit();
                                }
                                $DEBUG
                                  and print "\t\t\t\tType $i:  Error so far $error\n\r\n\r";
                            }
                        }
                        if ( $error==0 )  {
                            $DEBUG and print "--------------------------------------------------\n\r\n\r";
                            $DEBUG and print "Original network: ", $nnetwork, "\n\r";
                            $DEBUG and print "New comparison map is: $comparisonMap\n\r";
                            
                            ## no errors so add vertice added to the new network
                            #my $oldedgenum = $nnetwork->edges;
                            my @vertices = $nnetwork->vertices;
                            #if ( ! grep { $_ eq $testAssemblage} @vertices) {
                                #first make a copy
                                #print Dumper($nnetwork);
                                my $new_network = $nnetwork->deep_copy_graph;
                                #print Dumper($new_network);
                                $solutionCount++;    ## increment the # of solutions
                                $new_network->set_graph_attribute("GraphID", $solutionCount);
                                ## now add to this *new one*
                                $new_network->add_vertex($testAssemblage);
                                ## mark this vertice as the new "END"
                                $new_network->set_vertex_attribute($testAssemblage,"End", "1" );
                                ## mark the interior vertice as not and "END"
                                $new_network->set_vertex_attribute($endAssemblage,"End", "0" );
                                #### This adds the comparison to the new edge that has been added.
                                $new_network->add_weighted_edge( $endAssemblage, $testAssemblage,  $comparisonMap );
                                ## Mark this edge and a new end edge
                                $new_network->set_edge_attribute($endAssemblage , $testAssemblage, "End", "1");
                                $new_network->set_edge_attribute($endAssemblage , $testAssemblage, "GraphID", $solutionCount);
                                ## mark the previous edge as no longer the end 
                                my @n = $new_network->neighbours($endAssemblage);  ## assume 0 is the only neighbor (should be only one!)
                                $new_network->set_edge_attribute($n[0], $endAssemblage , "End", "0");
                                
                                #print "Which end: ", $whichEnd, "\n";
                                if ($whichEnd==1) {
                                    $new_network->set_graph_attribute("End_1", $testAssemblage);
                                } else {
                                    $new_network->set_graph_attribute("End_2", $testAssemblage);
                                }
                                 ## increment the ends (we start with 0, then 1)
                                
                                $DEBUG and print "New network (with addition): ", $new_network, "\n\r";
                                ## copy this solution to the new array of networks
                                #print Dumper($new_network);
                                #print $new_network, "\n";
                                push @newnets, $new_network;   ## contains solutions for just this step - add to list
                                #push @networks, $new_network;
                                if ($nosum==0 ) {
                                    push @solutions, $new_network; ## this is a list of all the solutions (shallow a)
                                } 
                                
                                my $currentTotal =  scalar(@newnets);
                                if (($new_network->unique_edges) > $maxEdges) {
                                    $maxEdges = $new_network->unique_edges;
                                    $screen and $scr->at(6,1)->puts("Current Max Edges: $maxEdges   ");
                                }
                                $screen and $scr->at(7,1)->puts("Sum of all solutions up to this step: $solutionCount");
                                $screen and $scr->at(8,43)->puts("                                           ");
                                $screen and $scr->at(8,1)->puts("Current number of seriation linkages at this step: $currentTotal");
                                if ($memusage) {
                                my $text = "Memory used for array @ newnets: ". total_size(\@newnets);
                                    $screen and $scr->at(9,1)->puts($text);
                                }
                                $DEBUG and print "-------------------------------------------------\n\r";
                            ##}
                            
                        }
                }    # end of iterate through the existing network link
            }    # end of if assemblage not already in netowrk check
        }    #end of assemblage loop
    }
    #print Dumper(\@networks);
    ########################################### AFTER ALL THE SORTING AND WINNOWING ######################
    $DEBUG and print "Number of current solutions now: ", scalar(@newnets), "\n\r";

    ## now push the current step solutions to the list that serves as the basis of the next step. 
    #print Dumper(\@solutions);
    #print "num of new nets: ", scalar(@newnets),"\n";
    ## no match at this point so no point in going forward.
    if ( scalar(@newnets) == 0 ) {
        my $text = "Max seriation size reached - Largest solution set: ". scalar(@networks). " out of ". $solutionCount; 
        $screen and $scr->at(9,1)->puts( $text );
        $screen and $scr->at(10,1)->puts("Maximum # edges in largest solution is: $maxEdges (note # of edges = # of assemblages - 1)");
        $DEBUG and print "Maximum # edges in largest solution is: $maxEdges (note # of edges = # assemblages -1) \n\r\n";
        $currentMaxSeriationSize = $maxSeriations+1;
    } else {
        ## copy the array of new solutions back to the working set to continue. over time this should get smaller and smaller...
        ### release the memory
        undef @networks;
        @networks= @{\@newnets};
        undef @newnets;
        @newnets=();  ## clear the array for the next new set of assemblages.
        $currentMaxSeriationSize++;
    }    #end of network loop
}    #end of master loop through iterations

if (scalar(@networks)==0 || $networks[0] eq undef ) {
    print "\n\r\n\r\n\r\n\r\n\rNo solutions Found!!\n\r";
    exit();
}

########################################### DO SOME FILTERING (OR NOT) ####################################
# now do some weeding. Basically start with the first network that is the largest, and work backwards. Ignore any
# network that is already represented in the smaller ones since these are trivial (e.g., A->B->C->D alreacy covers
# A->C->D.) That should then leave all the unique maximum solutions (or so it seems)

## first need to sort the networks by size
my @filteredarray = ();
if ( $filterflag == 1 ) {
    $screen and $scr->at(1,40)->puts("STEP: Filter to get uniques... ");
    $DEBUG and print "---Filtering solutions so we only end up with the unique ones.\n";
    $DEBUG and print "----Start with ", scalar(@solutions), " solutions. \n";
   for (my $i=scalar(@solutions)-1;$i>-1;$i--) {
      my $exists=0;
      foreach my $tnetwork (@filteredarray) {
        #print "tnetwork: ", $stepSeriationList{ $tnetwork} , "\n\r";
        my @fnetworkArray = $solutions[$i]->vertices;
        my @tnetworkArray = $tnetwork->vertices;
            my @minus = array_minus( @fnetworkArray, @tnetworkArray );
            if (scalar(@minus)== 0) {
                  $exists++;
            }
         }
      if (!$exists) {
         ##print "pushing $fnetwork to list\n\r";
         push @filteredarray, $solutions[$i] ;
      } 
   }
    $DEBUG and print "End with ", scalar(@filteredarray), " solutions.\n";
    my $filterCount= scalar(@filteredarray);
    $screen and $scr->at(11,1)->puts("End with $filterCount solutions.\n");
} else {
    #print "\n\rNow going to print just the largest network out of a pool of ", scalar(@networks), "\n\r";
    #@filteredarray = @{\@networks}; ; ## just the largest one
    @filteredarray = @networks ; ## just the largest one
    }

if ($allSolutions>0) {
    @filteredarray = @solutions; ### all of the solutions as a default
}

###########################################  OUTPUT INDIVIDUAL .vcg and .dot FILES ###################
if ($individualfileoutput) {
    $screen and $scr->at(1,40)->puts("STEP: Individual file output  ");
    my $writer = Graph::Writer::VCG->new();
    $count = 0;
    my $name;
    foreach my $network (@filteredarray) {
        my $E = $network->unique_edges;
        if ( $E > $maxnumber - 1 ) {
            $count++;
            $name = $count . '.vcg';
            $writer->write_graph( $network, $name );
        }

        #print Dumper($network);
    }

    my $writer2 = Graph::Writer::Dot->new();
    $count = 0;
    foreach my $network (@filteredarray) {
        my $V = $network->vertices;
        if ( $V == $maxEdges ) {
            $count++;
            $name = $count . '.dot';
            $writer2->write_graph( $net, $name );
        }

        #print Dumper($network);
    }
}


########################################### OUTPUT SECTION ####################################
$screen and $scr->at(13,1)->puts( "Now printing output file... ");
    $screen and $scr->at(1,40)->puts("STEP: Output files...         ");
print OUTFILE "*Node data\n";
print OUTFILE "ID AssemblageSize X Y Easting Northing\n";
print OUTPAIRSFILE "*Node data\n";
print OUTPAIRSFILE "ID AssemblageSize X Y Easting Northing\n";
print OUTDOTFILE "graph seriation \n{\n";
print OUTDOTFILE "\n/* list of nodes */\n";
if ($mst) {
    print OUTBOOTSTRAPFILE "*Node data\nID AssemblageSize X Y Easting Northing\n";
    print OUTDISTANCEFILE "*Node data\nID AssemblageSize X Y Easting Northing\n";
}
$count = 0;
$screen and $scr->at(1,40)->puts("STEP: Printing list of nodes....     ");
## note this assumes the use of UTM coordinates (northing and easting)
foreach my $l (@labels) {
    #print OUTFILE $l, "\n";
    my $x = $xAssemblage{ $l }/1000000 || 0;
    my $y = ($largestY-$yAssemblage{ $l })/100000 || 0;
    my $easting = $xAssemblage{ $l } || 0;
    my $northing = $yAssemblage{ $l } || 0;
    print OUTFILE $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
    print OUTPAIRSFILE $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
    print OUTDOTFILE "\"".$l."\";\n";
    if ($mst){
        print OUTBOOTSTRAPFILE  $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
        print OUTDISTANCEFILE  $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
    }
}
print OUTFILE "*Node properties\nID AssemblageSize X Y Easting Northing\n";
print OUTPAIRSFILE "*Node properties\nID AssemblageSize X Y Easting Northing\n";
if ($mst) {
    print OUTBOOTSTRAPFILE "*Node properties\nID AssemblageSize X Y Easting Northing\n";
    print OUTDISTANCEFILE "*Node properties\nID AssemblageSize X Y Easting Northing\n";
}
$screen and $scr->at(1,40)->puts("STEP: Printing list of nodes attributes... ");
foreach my $l (@labels) {
   my $x = $xAssemblage{ $l }/1000000 || 0;
    my $y = ($largestY-$yAssemblage{ $l })/100000 || 0;
    my $easting = $xAssemblage{ $l } || 0;
    my $northing = $yAssemblage{ $l } || 0;
    print OUTFILE $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
    print OUTPAIRSFILE $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
    print OUTDOTFILE "\"".$l."\";\n";
    if ($mst) {
        print OUTBOOTSTRAPFILE $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
        print OUTDISTANCEFILE $l . " ". $assemblageSize{ $l }." ".$x." ".$y." ".$easting." ".$northing."\n";
    }
}

## This prints out counts of the edges as they appear in ALL of the solutions
$screen and $scr->at(1,40)->puts("STEP: Going through and counting pairs...     ");
print OUTPAIRSFILE "*Tie data\nFrom To Edge Count\n";
if ($mst) {
    print OUTBOOTSTRAPFILE "*Tie data\nFrom To Edge End Weight ID\n";
    print OUTDISTANCEFILE "*Tie data\nFrom To Edge End Weight ID\n";
}
## first count up all of the edges by going through the solutions and each edge
## put the edge count in a hash of edges
my %edgeHash=();

foreach my $network (@filteredarray) {
    #print "\n\r netowrk type ", ref($network),"\n\r";
    my @Edges;
    if (ref($network) eq 'REF') {
        @Edges = $$network->unique_edges;
    } elsif ($network eq undef ) {
        next;
    } else {
        @Edges = $network->unique_edges
    }
    my $eCount=0;
    foreach my $e (@Edges) {   
        my $edge0 = @$e[0];
        my $edge1 = @$e[1];
        my $pairname= $edge0." ".$edge1;
        $edgeHash{ $pairname }++;
    }
}
## now go through the edgeHash and print out the edges
## do this is sorted order of the counts. For fun.
$screen and $scr->at(1,40)->puts("STEP: Doing the pair output...                ");
foreach (sort { ($edgeHash{$b} cmp $edgeHash{$a}) || ($b cmp $a) } keys %edgeHash)  {
    print OUTPAIRSFILE $_, " 1 ", $edgeHash{$_}, "\n";
}

print OUTFILE "*Tie data\nFrom To Edge Weight Network End pValue pError meanSolutionDistance\n";
$screen and $scr->at(1,40)->puts("STEP: Eliminating duplicates...     ");
my @uniqueArray = uniq @filteredarray;

$screen and $scr->at(1,40)->puts("STEP: Printing edges...     ");
print OUTDOTFILE "\n/* list of edges */\n";

my %distanceHash=();
my %seriationHash;
## only print unique ones...

foreach my $network (@uniqueArray) {
    if (ref($network) eq 'REF') {
        $network = $$network;
    } elsif ($network eq undef ) {
        next;
    }
    $screen and $scr->at(14,1)->puts( "Now on solution: ");
    $screen and $scr->at(14,18)->puts($network->get_graph_attribute("GraphID") );
    my $eCount;   
    
    if ($largestonly>0) {
        if ($network->unique_edges == $maxEdges) {
            
            my $groupDistance=0;
            my @Edges = $network->unique_edges;
            my $meanDistance=0.0;
            my $eCount=0;
            if ($xyfile) {
               foreach my $e (@Edges) {
                  my $edge0 = @$e[0];
                  my $edge1 = @$e[1];
                  my $pairname= $edge0."*".$edge1;
                  $groupDistance += $distanceBetweenAssemblages{ $pairname };
                  $eCount++;
               }
               $meanDistance = $groupDistance/$eCount;      ## use the average for the group for now
               #print "\n\rMean distance for this group is: ", $meanDistance, "\n\r";
               $network->set_graph_attribute("meanDistance", $meanDistance);
            } else {
                $meanDistance="0";
                $network->set_graph_attribute("meanDistance", "0");
            }
            my $text;
            foreach my $e (@Edges) {
               my $edge0 = @$e[0];
               my $edge1 = @$e[1];
               my ($pVal, $pErr);
                if ( $pairwiseFile ) {
                    my $pairname= $edge0."#".$edge1;
                    $pVal = $pairwise{ $pairname };
                    $pErr = $pairwiseError{ $pairname };
                } else {
                    $pVal = 0.0;
                    $pErr = 0.0;
                }
                print OUTDOTFILE "\"",@$e[0], "\""," -- ", "\"", @$e[1], "\"", " [weight = \"", $network->get_edge_weight(@$e[0], @$e[1]),"\" ];\n";
                $text = @$e[0]. " ". @$e[1]." 1 ".scalar(@Edges). " ". $network->get_graph_attribute("GraphID"). " ". $network->get_edge_attribute(@$e[0], @$e[1], "End")." ". $pVal." ". $pErr;
                print OUTFILE $text, " ", $meanDistance, "\n";
            }
            $network->set_graph_attribute("meanDistance", $meanDistance);
            $distanceHash{ $network->get_graph_attribute("GraphID") }= $meanDistance;
            $seriationHash{ $network->get_graph_attribute("GraphID") }->{'meanDistance'}= $meanDistance;
            $seriationHash{ $network->get_graph_attribute("GraphID") }->{'ID'}=$network->get_graph_attribute("GraphID");
            $seriationHash{ $network->get_graph_attribute("GraphID") }->{'size'}=scalar(@Edges);
        }
    } else {  ## not just the largest, but ALL seriation solutions
            
        my @Edges = $network->unique_edges;
        my $groupDistance=0;
        my $meanDistance=0.0;
        my $eCount=0;
        my $text;
        if ($xyfile) {
           foreach my $e (@Edges) {
              my $edge0 = @$e[0];
              my $edge1 = @$e[1];
              my $pairname= $edge0."*".$edge1;
              $groupDistance += $distanceBetweenAssemblages{ $pairname };
              $eCount++;
           }
           $meanDistance = $groupDistance/$eCount;         ##use the average distance as the metric
        } else {
           $meanDistance = "0";
        }
        foreach my $e (@Edges) {
           my $edge0 = @$e[0];
           my $edge1 = @$e[1];
           my ($pVal, $pErr);
            if ( $pairwiseFile ) {
                my $pairname= $edge0."#".$edge1;
                $pVal = $pairwise{ $pairname };
                $pErr = $pairwiseError{ $pairname };
            } else {
                $pVal = 0.0;
                $pErr = 0.0;
            }
            print OUTDOTFILE "\"", @$e[0],"\"", " -- ", "\"", @$e[1], "\"", " [weight = \"", $network->get_edge_weight($edge0, $edge1),"\" ];\n";
            $text = @$e[0]. " ". @$e[1]." 1 ".scalar(@Edges). " ". $network->get_graph_attribute("GraphID"). " ". $network->get_edge_attribute(@$e[0], @$e[1], "End")." ". $pVal." ". $pErr;
            print OUTFILE $text, " ", $meanDistance, "\n";        
        }

        $network->set_graph_attribute("meanDistance", $meanDistance);
        $distanceHash{ $text }= $meanDistance;
        $seriationHash{ $network->get_graph_attribute("GraphID") }->{'meanDistance'}= $meanDistance;
        $seriationHash{ $network->get_graph_attribute("GraphID") }->{'ID'}=$network->get_graph_attribute("GraphID");
        $seriationHash{ $network->get_graph_attribute("GraphID") }->{'size'}=scalar(@Edges);
    }
}

## Do a MST analysis

if ($mst) {
    my $megaNetwork = Graph::Undirected->new;
    my $distanceNetwork = Graph::Undirected->new;
    
    foreach my $network (@uniqueArray) {
        if (ref($network) eq 'REF') {
            $network = $$network;
        } elsif ($network eq undef ) {
            next;
        }
        $screen and $scr->at(14,1)->puts( "MST creation - solution: ");
        $screen and $scr->at(14,26)->puts($network->get_graph_attribute("GraphID"));
        my $eCount;   
        my @Edges = $network->unique_edges;
        my $graphID = $network->get_graph_attribute("GraphID");
        my $groupDistance=0;
        my $meanDistance=0.0;    
        # make a gigantic nework
      
        foreach my $e (@Edges) {
            my $edge0 = @$e[0];
            my $edge1 = @$e[1];
            if ( !$pairwiseFile ) {
                $perror{ @$e[0] . "#" . @$e[1] } = 0.0;
                $pvalue { @$e[0] . "#" . @$e[1] } = 0.0;
            }
            if ( $pairwiseFile ) {
                my $pairname= $edge0."#".$edge1;
                $perror{ @$e[0] . "#" . @$e[1] } = $pairwiseError{ $pairname };
                $pvalue{ @$e[0] . "#" . @$e[1] } = $pairwise{ $pairname };
            }
            my $pairname= $edge0."#".$edge1;
            $megaNetwork->add_vertex( $edge0 );
            $megaNetwork->add_vertex( $edge1 );
            $megaNetwork->add_edge($edge0, $edge1);
            $megaNetwork->set_edge_weight($edge0, $edge1, $pairwise{ $pairname } );
            $megaNetwork->set_edge_attribute($edge0, $edge1, "GraphID", $graphID);
            my $distpairname= $edge0."*".$edge1;
            $distanceNetwork->add_vertex( $edge0 );
            $distanceNetwork->add_vertex( $edge1 );
            my $distance = $distanceBetweenAssemblages{ $distpairname };
            $distanceNetwork->add_edge($edge0,$edge1);
            $distanceNetwork->set_edge_weight($edge0, $edge1, $distance);
            $distanceNetwork->set_edge_attribute($edge0, $edge1, "GraphID", $graphID);
        }
    }
        
    my $mstgBootstrap = $megaNetwork->minimum_spanning_tree;
    $count=1;
    my @Edges = $mstgBootstrap->unique_edges;
    foreach my $e (@Edges) {
        my $edge0 = @$e[0];
        my $edge1 = @$e[1];
        my $pairname= $edge0."#".$edge1;
        print OUTBOOTSTRAPFILE $edge0, " ", $edge1, " ", $count, " ", $megaNetwork->get_edge_attribute( $edge0, $edge1, "End"), " ", $pairwise{ $pairname }, " ", $megaNetwork->get_edge_attribute( $edge0, $edge1, "GraphID"), "\n";
        $count++;
    }
    my $mstgDistance = $distanceNetwork->minimum_spanning_tree;
    $count=1;
    @Edges = $mstgDistance->unique_edges;
    foreach my $e (@Edges) {
        my $edge0 = @$e[0];
        my $edge1 = @$e[1];
        my $pairname= $edge0."*".$edge1;
        print OUTDISTANCEFILE $edge0, " ", $edge1, " ", $count, " ", $megaNetwork->get_edge_attribute( $edge0, $edge1, "End"), " ", $distanceBetweenAssemblages{ $pairname }, " ", $distanceNetwork->get_edge_attribute( $edge0, $edge1, "GraphID"),"\n";
        $count++;
    }
}

if ($xyfile ) {
    #print "SeriationHash: ", Dumper(\%seriationHash);
    my @sortedKeys = SortHashByMultipleColumns(\%seriationHash,["meanDistance:asc","size:dsc"]);
    my @data = [];
    $count=0;
    open(OUTSIZEFILE, ">$inputfile-distances.txt") or die $!;
    print OUTSIZEFILE "Seriation_Solution Mean_Distance Solution_Size\n";
    foreach my $sortedKey(@sortedKeys){
        print OUTSIZEFILE $sortedKey. " ". $seriationHash{$sortedKey}->{'meanDistance'} . " " . $seriationHash{$sortedKey}->{'size'} . "\n";
        $data[ $count ] = int($seriationHash{$sortedKey}->{'meanDistance'});
        $count++;
    }
    my $graph = new GD::Graph::histogram(400,600);
    $graph->set( 
                x_label         => 'Mean Distance Between Assemblages',
                y_label         => 'Count',
                title           => "Seriation Solutions for $inputfile",
                x_labels_vertical => 1,
                bar_spacing     => 0,
                shadow_depth    => 1,
                shadowclr       => 'dred',
                transparent     => 0,
                 histogram_bins => 15,
            ) 
            or warn $graph->error;
        
    my $gd = $graph->plot(\@data) or die $graph->error;
    open(IMG, ">$useOutputFile-histogram.png") or die $!;
    binmode IMG;
    print IMG $gd->png;
}

if ($stats>0 ) {
    my %vertHash;
    foreach my $network (@uniqueArray) {
        if (ref($network) eq "REF") {
            $network = $$network;
        }  elsif ($network eq undef ) {
            next;
        }
        my @verts = $network->vertices;      
        foreach my $vert (@verts) {
            $vertHash{ $vert }++;
        }
    }
    my @data;
    foreach my $key (keys %vertHash){
        push @{$data[0]},$key;
        push @{$data[1]}, $vertHash{ $key };
    }
    my $graph = new GD::Graph::bars(400,600);
    $graph->set( 
                x_label         => 'Assemblages',
                y_label         => 'Count',
                title           => "Counts of assemblages in solutions for $inputfile",
                x_labels_vertical => 1,
                bar_spacing     => 0,
                shadow_depth    => 1,
                shadowclr       => 'dred',
                transparent     => 0,
            ) 
            or warn $graph->error;
    
    my $gd = $graph->plot(\@data) or die $graph->error;
    open(IMG, ">$useOutputFile-count-histogram.png") or die $!;
    binmode IMG;
    print IMG $gd->png;
}
print OUTFILE "\n";
print OUTDOTFILE "}\n";



###########################################
if ($excel) {
    print "\n\rNow printing excel output file... \n\r";
    # Create a new Excel workbook
    my $filename = $useOutputFile."-seriationList.xlsx";
    my $workbook = Excel::Writer::XLSX->new( $filename );
    
    # Add a worksheet
    my $worksheet = $workbook->add_worksheet();
    my $percentWorksheet = $workbook->add_worksheet();
    my $row = 1;
    my $solutionNumber=1;
    foreach my $network (@uniqueArray) {    
        if (ref($network) eq 'REF') {
            $network = $$network;
        } elsif ($network eq undef ) {
            next;
        }
        #my $undirectedNetwork = $network->undirected_copy;
        my $col=0;
        if ($largestonly>0) {
            if ($network->unique_edges == $maxEdges) {
                $worksheet->write( $row, $col, "Seriation Solution");
                $percentWorksheet->write($row,$col, "Seriation Solution");
                $col++;
                
                my $endAssemblage;
                my @verts = $network->longest_path;
                
                my $vs= scalar(@verts);
                for (my $j; $j<$vs; $j++) {
                    $worksheet->write( $row, $col, $network->get_graph_attribute("GraphID")  );
                    $percentWorksheet->write($row, $col, $network->get_graph_attribute("GraphID")  );
                    $col++;
                    $worksheet->write( $row, $col, $verts[$j] );
                    $percentWorksheet->write( $row, $col, $verts[$j] );
                    my @assemblage= @{ $assemblageValues{ $verts[$j] } };
                    my @assemblageFreq = @{ $assemblageFrequencies{ $verts[$j] } };
                    my $cols = scalar(@assemblage);
                    for ( my $i = 0 ; $i < $cols ; $i++ ) {
                        $col++;
                        $worksheet->write($row, $col, $assemblage[$i]);
                        $percentWorksheet->write($row,$col, $assemblageFreq[$i]);
                    }
                    $col=1;
                    $row++;
                }
                $row = $row+2;
                $solutionNumber++;
            }
        } else {
            $worksheet->write( $row, $col, "Seriation Solution");
            $percentWorksheet->write( $row, $col, "Seriation Solution");
            $col++;
            my $endAssemblage;
            my @verts = $network->longest_path;
            my $vs= scalar(@verts);
            for (my $j; $j<$vs; $j++) {
                $worksheet->write( $row, $col, $network->get_graph_attribute("GraphID") );
                $percentWorksheet->write( $row, $col, $network->get_graph_attribute("GraphID") );
                $col++;
                $worksheet->write( $row, $col, $verts[$j] );
                $percentWorksheet->write($row, $col, $verts[$j] );
                my @assemblage= @{ $assemblageValues{ $verts[$j] } };
                my @assemblageFreq = @{ $assemblageFrequencies{ $verts[$j] } };
                my $cols = scalar(@assemblage);
                for ( my $i = 0 ; $i < $cols ; $i++ ) {
                    $col++;
                    $worksheet->write($row, $col, $assemblage[$i]);
                    $percentWorksheet->write($row,$col, $assemblageFreq[$i]);
                }
                $col=1;
                $row++;
            }
            $row = $row+2;
            $solutionNumber++;
        }
    }       
}
my $timediff= Time::HiRes::gettimeofday() - $start;
$screen and $scr->at(14,1)->puts( "Time for processing: $timediff seconds");
print "\n\r\n\r\n\r\n\r";


sub SortHashByMultipleColumns{
    my($hashRef,$sortInfoAR) = @_;
    my $sortCriteria;
    
    foreach my $sortInfo(@$sortInfoAR){
        my($sortColumn,$sortDirection) = split(/\:/,$sortInfo);
        
        my $sortType;
        # VERIFY THAT THIS IS A VALID COLUMN
        if(! defined $hashRef->{((keys %$hashRef))[0]}->{$sortColumn}){
            print "SortHashByMultipleColumns Error: $sortColumn is not a column in this hash.\n";
            exit(0);
        # VALID COLUMN, FIGURE OUT IF IT IS A NUMERIC COLUMN OR ALPHA
        } else {
            if($hashRef->{((keys %$hashRef))[0]}->{$sortColumn} =~ /^\d+/){
                $sortType = "numeric";
            } else {
                $sortType = "alpha";
            }
        }
            
        # want to sort a number 
        if($sortType eq "numeric"){
            # sort it asc
            if($sortDirection =~ /asc/i){
                # add an or if we already have something in the sort criteria
                if($sortCriteria){
                    $sortCriteria .= qq| or |;
                }
                $sortCriteria .= '$hashRef->{$a}->{\'' . $sortColumn . '\'} <=> $hashRef->{$b}->{\'' . $sortColumn . '\'}';
            # sort it desc
            } else {
                # add an or if we already have something in the sort criteria
                if($sortCriteria){
                    $sortCriteria .= qq| or |;
                }
                $sortCriteria .= '$hashRef->{$b}->{\'' . $sortColumn . '\'} <=> $hashRef->{$a}->{\'' . $sortColumn . '\'}';
            }
        # want to sort it by alpha 
        } else {
            # sort it asc
            if($sortDirection =~ /asc/i){
                # add an or if we already have something in the sort criteria
                if($sortCriteria){
                    $sortCriteria .= qq| or |;
                }
                $sortCriteria .= '$hashRef->{$a}->{\'' . $sortColumn . '\'} cmp $hashRef->{$b}->{\'' . $sortColumn . '\'}';
            # sort it desc
            } else {
                # add an or if we already have something in the sort criteria
                if($sortCriteria){
                    $sortCriteria .= qq| or |;
                }
                $sortCriteria .= '$hashRef->{$b}->{\'' . $sortColumn . '\'} cmp $hashRef->{$a}->{\'' . $sortColumn . '\'}';
            }
        }
    }
    
    my $sortfunc = eval "sub { $sortCriteria }";
       my @sortedKeys = sort $sortfunc keys %$hashRef;

    return @sortedKeys;
}

__END__

    =head1 IDSS-with-bootstrap.pl

    Iterative Determinisitic Seriation Solutions with bootstrap significance testing

    =head1 SYNOPSIS

    IDSS-with-bootstrap.pl [options] 

     Options:
       -help                brief help message
       -man                 full documentation
       -input=<filename.txt>    filename of data to seriate
       -xyfile=<filename>   filename of XY data for assemblages
       -threshold=<value>   value specified as the maximum % difference between assemblages examined
       -largestonly         only the largest seriation solutions are printed in output
       -filtered            filter the output set to get just the unique solutions (no repeats)
       -indivfiles          create individual .vna files for each solution
       -bootstrapCI         use bootstrap confidence intervals for comparison
       -bootstrapSignficance    specify the significance level of the bootstrap confidence intervals (default=95)
       -debug               print debugging output
       -noscreen            don't use terminal output - just standard out
       -excel               output excel files for creating graphical seriation (not working yet)
       -mst                 minimum spanning tree
       -pairwise=<filename> pairwise comparisons
       -stats               statistics output
       
    Output will be:
        <filename>.vna      Netdraw file
        <filename>.dot      DOT file

    Example:
        perl IDSS-with-bootstrap.pl -input=../testdata/pfg-cpl.txt -xyfile=../testdata/pfgXY.txt -bootstrapCI -bootstrapsignificance=99.5 -pairwise=../testdata/pfg-cpl-bootstrap.txt -mst
        
    =head1 OPTIONS

    =over 8

    =item B<-help>

    Print a brief help message and exits.

    =item B<-man>

    Prints the manual page and exits.

    =back

    =head1 DESCRIPTION

    B<This program> reads given input file(s) and create seriations.

    =cut
