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
require Term::Screen;

my $debug      = 0;
my $filterflag = 0; ## do you want to try to output all of the solutions (filtered for non trivial)
my $largestOnly          = 0; #       # only output the largest set of solutions
my $individualfileoutput = 0; ## create files for all the indivdual networks
my $bootstrap            = 0; ## flag for bootstrap
my $man                  = 0;
my $help                 = 0;
my $inputfile;
my $bootstrapdebug = 0;
my $threshold      = 0;
my $noscreen        = 0;     ## flag for screen output
my $excel          = 0;       ## flag for excel file output (not implemented yet)

# process command line options; if none, print a usage/help message.
# note - manual page for explaining the options, what they do, and how to use them
# is at the bottom of this file in simple POD format.

GetOptions(
    'debug'          => \$debug,
    'bootstrap'      => \$bootstrap,
    'bootstrapdebug' => \$bootstrapdebug,
    'filtered'       => \$filterflag,
    'largestonly'    => \$largestOnly,
    'indivfiles'     => \$individualfileoutput,
    'help'           => sub { HelpMessage() },
    'input=s'        => \$inputfile,
    'excel'          => \$excel,
    'threshold=f'      => \$threshold,
    'noscreen'         => \$noscreen,
    man              => \$man
) or pod2usage(2);

my $DEBUG = $debug;    # our "$debug level"

if ($DEBUG) {
    print "Verbose debugging output is on!!!\n";
    print "Processing input file: $inputfile\n";
    $filterflag           or print "filterflag is off\n";
    $bootstrap            or print "bootstrap is off\n";
    $largestOnly          or print "output largest solutions is off\n";
    $individualfileoutput or print "individual network file output is off\n";
    print "threshold is currently set to: $threshold\n";
    $noscreen or print "noscreen output is currently off (which means that there is output to the screen)\n";
    $excel or print "excel output is off\n";
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
my %assemblageFrequencies = {};
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
open( INFILE, $inputfile ) or die "Cannot open $inputfile.\n";
open( OUTFILE, ">$inputfile.vna" ) or die "Can't open file $inputfile.vna to write.\n";
open( OUTDOTFILE, ">$inputfile.dot") or die "Can't open file $inputfile.dot to write.\n";

## some output so we know what we are doing 
$screen and $scr->at(1,1)->puts("Filename:  $inputfile");
$screen and $scr->at(2,1)->puts("Threshold: $threshold");

## Read in the data
# the input of the classes -- note these are absolute counts, not frequencies
# might want to change that...
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
        for (@line) {
            # push @freq, $_;
            my $f = $_ / $rowtotal;
            my $ff = sprintf( "%.4f", $f );
            push @$freq, $ff;
        }
        push @assemblages, [@$freq];
        $assemblageFrequencies{$label} = $freq;
        $count++;
    }
    #print "---- row end ----\n";
}
$numrows = scalar(@assemblages);

my $maxSeriations = $count;
$screen and $scr->at(3,1);
$screen and $scr->puts("Maximum possible seriation solution length: $count");


#############################################  THRESHOLD DETERMINATION ####################################
## first compare each assemblage and see if the threshold is exceeded.
## By this I mean the maximum difference between frequencies of any type is greater than what is specified.
## When threshold = 0, all combinations are used. Later the "ends" of solutions are not evaluated if the
## difference between the last assemblage and any other free assemblage is > the threshold.
## This arbitrary setting is to keep from arbitrary solutions being stuck on that come from the "ends"
## of solutions.
##
## Precalculate all of the max differences between types in assembalge pairs. 

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
        my $diff = abs( $ass1 - $ass2 );
        if ( $diff > $maxDifference ) {
            $maxDifference = $diff;
        }
    }
    $assemblageComparison{$pairname} = $maxDifference;
    $assemblageComparison{$pairname2} = $maxDifference;
}

$DEBUG and print Dumper( \%assemblageComparison ), "\n";


########################################### FIND ALL THE VALID TRIPLES  ####################################
my $directionstate;
$numrows = scalar(@assemblages);
my @nets;
my @triplettype;
my @tripletNames = ();
my @tripletArray = ();
my $net          = Graph::Directed->new;
my $comparison12;
my $comparison23;
my $error;
my $numberOfTriplets;

## This uses the Math::Combinatorics to create all the permutations of 3. This is the simplest solution programmatically
## Why recreate the wheel? Use Perl!

my $permutations =
  Math::Combinatorics->new( count => 3, data => [@assemblageNumber] );

while ( my @permu = $permutations->next_combination ) {
    #$DEBUG and print $labels[ $permu[0] ] . " * ". $labels[ $permu[1] ] . " * ". $labels[ $permu[2] ] . "\n";
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
        #$DEBUG and print "TESTING:  ". $ass1. "-" . $ass2 . "-" . $ass3 . "\n";

        my $dif1 = $ass1 - $ass2;

        #my $dif1 = $permu[0][$i] - $permu[1][$i];
        if ( $dif1 < 0 )  { $difscore = -1; }
        if ( $dif1 > 0 )  { $difscore = 1; }
        if ( $dif1 == 0 ) { $difscore = 0; }

        my $dif2 = $ass3 - $ass2;

        #my $dif2 = $permu[1][$i] - $permu[2][$i];
        if ( $dif2 < 0 )  { $difscore2 = -1; }
        if ( $dif2 > 0 )  { $difscore2 = 1; }
        if ( $dif2 == 0 ) { $difscore2 = 0; }

        ## F1 > F2 < F3 ## criteria not met
        if ( ( $difscore == 1 ) && ( $difscore2 == 1 ) ) {
            $error++;
        }
        elsif ( ( $difscore == 1 ) && ( $difscore2 == -1 ) ) {   ## F1 > F2 > F3
            $comparison12 .= "U";
            $comparison23 .= "D";
        }
        elsif ( ( $difscore == -1 ) && ( $difscore2 == -1 ) ) {   #  F1 < F2 >F3
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
            print "Not Match! (other)  Difscore 1: $difscore Difscore 2: $difscore2 \n";
            exit();
        }
    }

    if ( $error == 0 ) {
        undef $net;
        $net = Graph::Directed->new;
        $net->add_vertex( $labels[ $permu[0] ] );
        $net->add_vertex( $labels[ $permu[1] ] );
        $net->add_vertex( $labels[ $permu[2] ] );
        $net->add_weighted_edge(
            $labels[ $permu[0] ],
            $labels[ $permu[1] ],
            $comparison12
        );
        $net->add_weighted_edge(
            $labels[ $permu[1] ],
            $labels[ $permu[2] ],
            $comparison23
        );
        $DEBUG and print "VALID SOLUTION: " . $labels[ $permu[0] ] . " * " . $labels[ $permu[1] ] . " * " . $labels[ $permu[2] ] . "\n";
        $DEBUG and print "VALID SOLUTION: \t  $comparison12\t  ---   $comparison23\n";
        push @nets, $net;
        $numberOfTriplets++;
    }
    $error = 0;

}

########################################### THIS IS THE MAIN SERIATION SORTING SECTION ####################################
## now go through the combination of triplets that are valid (i.e., allthose that have no Z in them)
my @array = ();
## now go through all of the assemblages and each one to the top and bottom of the triplets to see if they fit.
## if they do, add them to the network

#my $currentMaxSeriations = 3;
my $currentMaxSeriationSize = 3;

#print "Number of Triplets:  ",$numberOfTriplets, "\n";
my $maxEdges  = 0;      ## keeps track of the current largest size of solution
my $stepcount = 0;     ## keeps track of the current step (from 0, adds one each step)
my $match     = 0;      ## keeps track of whether a solution is found for each step. Once no new solution is found, the search is over!

my @networks;  ## array of solutions from previous step (starts out with the 3s) (directed and deep graphs)
my @newnets;   ## array of new solutions (directed and deep graphs)
my @stepSeriationList = [];  ## array of solutions at this step (undirected and shallow copies of graphs)
my @solutions;  ## Array of all solutions (undirected and shallow copies of graphs)
my $solutionCount=0;   ## an index of the current number of solutions

### This is just preparing the intial set of threes into the array that will be used as the seed for all future
### seriation sets.
foreach my $n (@nets) {
      push @networks, \$n;        ## this is the array of the current successful set
      push @solutions, \$n;      ## This is an array of ALL solutions to date
      $stepSeriationList[ $solutionCount ] = \$n;     ## this is an ordered list of all solutions (order in which found for each step)
      $solutionCount++;          ## counts the current set of solutions
}
my $solutionSum = scalar(@networks);  ## number of solutions to this point (which is equal to all 3s);
my %seriationStep={};        ## hash of the array of solutions for this step

$seriationStep{$currentMaxSeriationSize}=\@stepSeriationList;  ## add to the list of solutions at this step
@stepSeriationList=();  ## clear out the step list. 

while ( $currentMaxSeriationSize < $maxSeriations ) {
    $currentMaxSeriationSize++;
    #my $index    = 0;
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
    foreach my $network (@networks) {
        #my $index = 0;
        my $nnetwork = $$network;
        $DEBUG and print "-----------------------------------------------------------------------------------\n";
        $DEBUG and print "Network: ", $nnetwork, "\n";
        $DEBUG and print "-----------------------------------------------------------------------------------\n";
        ## find the ends
        ## given the ends, find the valid set of assemblages that can be potentially added
        ## this list is all assemblages meet the threshold requirements
        foreach my $testAssemblage (@labels) {
            $DEBUG  and print "\t\tChecking assemblage: ", $testAssemblage, " to see if it fits on the end of the current solution.\n";
            $DEBUG  and print "\t\tFirst check to see if it is included already. If it has, move on.\n";
            if ( ! $nnetwork->has_vertex($testAssemblage) ) {
                # get the exterior vertices (should be 2)
                my @V = $nnetwork->vertices;    ## list of all the vertices
                $DEBUG  and print "\t\tFind the ends of the network. Do this by getting all the vertices \n";
                $DEBUG  and print " \t\tand looking for the ones with only 1 connection. There should be just 2 here.\n";
                ## loop through all of the edges to see if they can be stuck on the ends of the networks.
                foreach my $endAssemblage (@V) {
                    $DEBUG and print "\t\tChecking vertice: ", $endAssemblage, "\n";
                    my @Edges = $nnetwork->edges_at($endAssemblage);
                    my $edges = scalar(@Edges);
                    $DEBUG and print "\t\t\tThis vertice: $endAssemblage has this number of edges:  ", $edges, "\n";
                    my @newassemblage = ();
                    my @oldassemblage = ();
                    my $comparisonMap;
                    ## only if it has one edge. We only want to consider the ends of the netowrk -- so skip the others.
                    ## Just the ends.
                    if ( $edges == 1 ) {
                        $DEBUG and print "\t\t", $endAssemblage, " is on the edge since it only has one vertice.\n";
                        @newassemblage = @{ $assemblageFrequencies{ $testAssemblage } };
                        @oldassemblage = @{ $assemblageFrequencies{ $endAssemblage } };
                        my @edge = $nnetwork->edges_at($endAssemblage);
                        $DEBUG and print "\t\t\t The number of edges at $endAssemblage:  ", scalar(@edge), " (should be just one).\n";
                        my $connectedAssemblage = $edge[0][1];
                        my $g = $$network->get_edge_weight( $edge[0][0], $edge[0][1] );

                        #print Dumper $g;
                        $DEBUG and print "\t\t\tThere should be just 2 vertices here 0: $edge[0][0] and 1: $edge[0][1]\n";
                        $DEBUG and print "\t\t\t\t with a relation of $g\n";
                        my $numedges1 = $nnetwork->edges_at($edge[0][0]);
                        my $numedges2 = $nnetwork->edges_at($edge[0][1]);
                        my ($innerEdge, $outerEdge);
                        $DEBUG and print "\t\t\t\tNUMBER OF EDGES - $edge[0][0] - $numedges1  - $edge[0][1] - $numedges2\n";
                        if ( $numedges1 == 1 ) {
                           $outerEdge= $edge[0][0];
                           $innerEdge= $edge[0][1];
                        } else {
                           $outerEdge= $edge[0][1];
                           $innerEdge= $edge[0][0];                          
                        }
                        $DEBUG and print "\t\t\t Outer edge is $outerEdge\n";
                        $DEBUG and print "\t\t\t Inner edge is $innerEdge\n";
                        $DEBUG and print "\t\t\t Comparison is $g or -- ", $nnetwork->get_edge_weight($outerEdge,$innerEdge), "\n";
               
                        #first determine if the pairs are within the threshold value (0 = all assemblages)
                        my $pairname = $testAssemblage . " * " . $endAssemblage;
                        my $diff     = $assemblageComparison{$pairname};
                        $DEBUG and print "\t\t\tFor $pairname the max frequency difference is $diff.\n";
                        if ($diff=="") {
                           print "pairname not found in hash lookup!\n";
                           exit();
                        }
                        ## go through process only if threshold is 0 or difference value is below threshold
                        ## this should mean that network will not grow unless the above conditions are met.
                        my $error = 0;
                        if (  ($threshold>0 ) and ($diff > $threshold ))  {
                           $DEBUG and print "\t\t\tThreshold = $threshold and Diff = $diff. Since $diff < $threshold, continue.\n";
                           $error++; # this should ensure future failure....
                        }
                        my @comparison = split //, $g;
                        my $cols = scalar(@newassemblage);
                        my $comparisonMap;

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
                                
                                my $dif1 = $newassemblage[$i] - $oldassemblage[$i];
                                if ( $dif1 < 0 )  { $difscore = -1; }
                                if ( $dif1 > 0 )  { $difscore = 1; }
                                if ( $dif1 == 0 ) { $difscore = 0; }

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
                                    $DEBUG and print "\t\t\t\t\ This will only work if there no Xs anywhere previously.\n";
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
                                       if ($compArray[$i] =~ "X")  {
                                             $xerror = 1;  ### BLARGH a previous X! This will not be tolerated!
                                             $stopFlag = 1; ## We can stop
                                       }
                                       $DEBUG and print "\t\t\t\t\t Since I got $compArray[$i] my potential new value is still X.\n";
                                       $DEBUG and print "\t\t\t\t\t Now going to get the next pair of assembalges to examine in the chain\n";
                                       $ccount++;
                                       if ($ccount>$numrows) {
                                          print "PROBLEM. Too many checks. Quitting\n\r\n\r";
                                          exit();
                                       }
                                    }
                                    if ( $xerror ) {
                                        $error++;
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
                        
                        if ( !$error)  {
                            $DEBUG and print "--------------------------------------------------\n\r\n\r";
                            $DEBUG and print "Original network: ", $$network, "\n\r";
                            $DEBUG and print "New comparison map is: $comparisonMap\n\r";

                            ## no errors so add vertice added to the new network
                            #my $oldedgenum = $nnetwork->edges;
                            $nnetwork->add_weighted_edge( $testAssemblage, $endAssemblage, $comparisonMap );
                            $DEBUG and print "New network (with addition): ", $nnetwork, "\n\r";

                            ## copy this solution to the new array of networks
                            push @newnets, $nnetwork;   ## contains solutions for just this step
                            my $currentTotal =  scalar(@newnets);
                            if ($nnetwork->edges > $maxEdges) {
                              $maxEdges = $nnetwork->edges;
                           }
                            #print "MAX EDGES!!! ", $maxEdges ,"\n";
                            $screen and $scr->at(7,1)->puts("Sum of all solutions up to this step: $solutionSum");
                            $screen and $scr->at(8,43)->puts("                   ");
                            $screen and $scr->at(8,1)->puts("Current nunmber of seriation linkages at this step: $currentTotal");
                            $DEBUG and print "-------------------------------------------------\n\r";
                        }
                    }
                    else {
                        $DEBUG and print "\t\t$endAssemblage has too many edges ( $edges is more than 1) so skipping\n\r";
                    }    # end of if check for the end edges
                }    # end of iterate through the existing network link
            }    # end of if assemblage not already in netowrk check
        }    #end of assemblage loop
    }    #end of network loop
   ########################################### AFTER ALL THE SORTING AND WINNOWING ######################
   $DEBUG and print "Number of current solutions now: ", scalar(@newnets), "\n\r";

   ## now push the current step solutions to the list that serves as the basis of the next step. 
   foreach my $n (@newnets) {
      push @networks,$n;  ## use this for the next iteration -- its all the new solutions. 
      push @solutions, $n; ## this is a list of all the solutions (shallow a)
      $stepSeriationList[$solutionCount]= $n;  #This is supposed to be an array that keeps track of the solutions by the order they were created (low = small)
      $solutionCount++;
   }
   $solutionSum  =  scalar(@solutions);
   @stepSeriationList=();  ## clear the step seriation list so that we can get the next set for the next set of solutions. 
   
   @newnets=();  ## clear the array for the next new set of assemblages. 
   ## no match at this point so no point in going forward.
    if ( $match == 0 ) {
        $screen and $scr->at(9,1)->puts( "Maximum seriation size reached - no more assemblages added that iteration. ");
        $screen and $scr->at(10,1)->puts("Maximum # assemblages in largest solution is: $maxEdges");
        $DEBUG and print "Maximum # assemblages in largest solution is: $maxEdges\n\r\n";
        $maxnumber = $currentMaxSeriationSize - $maxSeriations;
        ## to break now...
        $currentMaxSeriationSize = $maxSeriations;
        ## copy the old network set to the new array (just in this case -- otherwise its empty
    }
}    #end of master loop through iterations


########################################### DO SOME FILTERING (OR NOT) ####################################
# now do some weeding. Basically start with the first network that is the largest, and work backwards. Ignore any
# network that is already represented in the smaller ones since these are trivial (e.g., A->B->C->D alreacy covers
# A->C->D.) That should then leave all the unique maximum solutions (or so it seems)

## first need to sort the networks by size
my @filteredarray = ();
if ( $filterflag == 1 ) {
    $DEBUG and print "---Filtering solutions so we only end up with the unique ones.\n";
    $DEBUG and print "----Start with ", scalar(@solutions), " solutions. \n";
    my $count     = scalar(@solutions)-1;
    while ($count>-1) {
         my $exists;
         my $netcheck = $solutions[$count];
         #print ref($netcheck),"\n";
         foreach my $fnetwork (@filteredarray) {
            #print "loop: ", ref($fnetwork),"\n";
            if ($netcheck ne $fnetwork) {
               $exists++; 
            }
         }
         unless ($exists) {
             push @filteredarray, $netcheck;
         }
         $count--;
    }
    $DEBUG and print "End with ", scalar(@filteredarray), " solutions.\n";
    my $filterCount= scalar(@filteredarray);
    $screen and $scr->at(11,1)->puts("End with $filterCount solutions.\n");
} elsif ($largestOnly) {
    @filteredarray = @networks; ## just the largest one
} else {
    @filteredarray = @solutions; ### all of the solutions as a default
}

###########################################  OUTPUT INDIVIDUAL .vcg and .dot FILES ###################
if ($individualfileoutput) {
    my $writer = Graph::Writer::VCG->new();
    $count = 0;
    my $name;
    foreach my $network (@filteredarray) {
        my $E = $network->edges;
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

########################################### BOOTSTRAP SECTION ####################################
#bootrap stuff

$numrows = scalar(@assemblages);
srand($start);
my %perror  = ();
my %pvalue  = ();
my $results = 0;
my $loop    = 100;
my $bigloop = 5;
my $loop2   = $loop;
my $ptr1    = 0;
my $ptr2    = 1;
my $classes = 0;

# now do ALL the pairwise assemblage comparisons
# go to sleep and come back later.

if ($bootstrap) {
    while ( $ptr1 < $numrows ) {
        while ( $ptr2 < $numrows ) {

            my $stat = Statistics::Descriptive::Full->new();

            my @a    = @{ $assemblages[$ptr1] };
            my @b    = @{ $assemblages[$ptr2] };
            my $numa = $rowtotals[$ptr1];
            my $numb = $rowtotals[$ptr2];

            # calculate the directionality for later
            my @direction = ();
            my $num       = scalar(@a);
            my $c         = 0;
            for ( $c = 0 ; $c < $num ; $c++ ) {
                if    ( $a[$c] < $b[$c] ) { push @direction, -1; }
                elsif ( $a[$c] > $b[$c] ) { push @direction, +1; }
                else                      { push @direction, 0; }
            }

            my ( @cum_a, @cum_b, $count );
            $classes = scalar(@a);
            my $index_a = 0;
            my $total_a = 0.0;
            $count   = 0;
            for ( $count = 0 ; $count < $classes ; $count++ ) {
                $cum_a[$index_a] = $a[$count];
                $total_a += $a[$count];
                $index_a++;
            }
            $classes = scalar(@b);
            my $index_b = 0;
            my $total_b = 0.0;
            $count = 0;
            for ( $count = 0 ; $count < $classes ; $count++ ) {
                $cum_b[$index_b] = $b[$count];
                $total_b += $b[$count];
                $index_b++;
            }

            $index_a--;
            $index_b--;

            # now we loop 100 times and keep track
            my $cycle = $bigloop;
            while ($cycle) {

                #print "($debug) cycle value: $cycle\n";

                # now we loop 1000 times and resample
                $loop = $loop2;
                while ($loop) {
                    my $assemsize = $numa;
                    my @assem_a   = ();
                    my $class;
                    my $total = scalar(@a);
                    my $rand;
                    while ($assemsize) {

                        #$rand = ( truly_random_value() % 10000 ) / 10000 ;
                        $rand  = rand;
                        $class = 0;
                        while (( $class < $index_a )
                            && ( $rand > $cum_a[$class] ) )
                        {
                            $rand -= $cum_a[$class];
                            $class++;
                        }
                        push @assem_a, $class;
                        $assemsize--;
                    }

                    $assemsize = $numb;
                    my @assem_b   = ();
                    $total = scalar(@b);
                    while ($assemsize) {

                        #$rand = ( truly_random_value() % 10000 ) / 10000 ;
                        $rand  = rand;
                        $class = 0;
                        while (( $class < $index_b )
                            && ( $rand > $cum_b[$class] ) )
                        {
                            $rand -= $cum_b[$class];
                            $class++;
                        }
                        push @assem_b, $class;
                        $assemsize--;
                    }

                    my ( @ahat, @bhat, %aholder, %bholder );
                    %aholder = ();
                    %bholder = ();
                    my $index = 0;

                    for ( $index = 0 ; $index < $cols ; $index++ ) {
                        $ahat[$index] = 0;
                        $bhat[$index] = 0;
                    }

                    for (@assem_a) {
                        $aholder{$_}++;
                    }

                    for (@assem_b) {
                        $bholder{$_}++;
                    }

                    for ( keys %aholder ) {
                        $ahat[$_] = ( $aholder{$_} / $numa );
                    }

                    for ( keys %bholder ) {
                        $bhat[$_] = ( $bholder{$_} / $numb );
                    }

                    # calculate the directionality for later
                    my @dir = ();
                    my $num = scalar(@ahat);
                    my $c   = 0;
                    for ( $c = 0 ; $c < $num ; $c++ ) {
                        $bootstrapdebug and print "loop $loop ", $ahat[$c] - $bhat[$c], "\t";
                        if    ( $ahat[$c] < $bhat[$c] ) { push @dir, -1; }
                        elsif ( $ahat[$c] > $bhat[$c] ) { push @dir, +1; }
                        else                            { push @dir, 0; }
                    }
                    $bootstrapdebug and print "\n";

                    # compare the two sets of results
                    $num = scalar(@dir);
                    $c   = 0;
                    my $diff = 0;
                    for ( $c = 0 ; $c < $num ; $c++ ) {
                        $bootstrapdebug and print "loop $loop ", $direction[$c], "/",
                          $dir[$c], "\t";
                        if ( $direction[$c] == $dir[$c] ) { next; }
                        $diff++;
                    }
                    $bootstrapdebug and print "\n";
                    if ( $diff == 0 ) { $results++; }

                    $loop--;
                }

                #print "Results:  $results matches of $loop2 trials\n";
                #print "Probability: ", $results / $loop2, "\n";

                $stat->add_data( $results / $loop2 );
                $cycle--;

                $results = 0;
            }
            my $label1 = $labels[$ptr1] . "-" . $labels[$ptr2];
            my $label2 = $labels[$ptr2] . "-" . $labels[$ptr1];
            $pvalue{$label1} = $stat->mean();
            $perror{$label1} = $stat->standard_deviation();
            $pvalue{$label2} = $stat->mean();
            $perror{$label2} = $stat->standard_deviation();
            $DEBUG and print $labels[$ptr1], "\t", $labels[$ptr2], "\t";
            $DEBUG and print $stat->mean(), "\t";
            $DEBUG and print $stat->standard_deviation(), "\n";

            undef $stat;
            $ptr2++;
        }
        $ptr1++;
        $ptr2 = $ptr1 + 1;

    }
}

########################################### OUTPUT SECTION ####################################
$screen and $scr->at(13,1)->puts( "Now printing output file... ");

print OUTFILE "*Node data\n";
print OUTDOTFILE "graph seriation \n{\n";
print OUTDOTFILE "\n/* list of nodes */\n";
$count = 0;
foreach my $l (@labels) {
    print OUTFILE $l, "\n";
    print OUTDOTFILE "\"".$l."\";\n";
}
print OUTFILE "*Tie data\n";
print OUTFILE "From To Edge Weight Network pValue pError\n";


my @uniqueArray;
foreach my $compareNetwork (@filteredarray) {
    my $exists = 0;

    if (ref($compareNetwork) eq "REF") {
      $compareNetwork = $$compareNetwork;
    }
    foreach my $uarray (@uniqueArray) {
      if (ref($uarray) eq "REF"){
         $uarray = $$uarray;
      }
        if ( $compareNetwork eq $uarray ) {
            $exists++;
        }
    }
    if ( !$exists ) {
        push @uniqueArray, $compareNetwork;
    }
}

print OUTDOTFILE "\n/* list of edges */\n";

## only print unique ones...
foreach my $network (@uniqueArray) {
   if (ref($network) eq "REF") {
      $network = $$network;
   }
    $count++;
    my $E = $network->edges;
    if ($largestOnly) {
        if ( $E == $maxEdges ) {
            my @Edges = $network->unique_edges;
            foreach my $e (@Edges) {
               my $edge0 = @$e[0];
               my $edge1 = @$e[1];
                if ( !$bootstrap ) {
                    $perror{ @$e[0] . "-" . @$e[1] } = 0.0;
                    $pvalue{ @$e[0] . "-" . @$e[1] } = 0.0;
                }
                print OUTFILE @$e[0], " ", @$e[1], ", 1, ", scalar(@Edges), ", ", $count, ", ";
                print OUTFILE $pvalue{ @$e[0] . "-" . @$e[1] }, ", ", $perror{ @$e[0] . "-" . @$e[1] }, "\n";
                print OUTDOTFILE "\"",@$e[0], "\""," -- ", "\"", @$e[1], "\"", " [weight = \"", $network->get_edge_weight(@$e[0], @$e[1]),"\" ];\n";
            }
            print OUTFILE "---------------------------\n";

        }
    } else {
        my @Edges = $network->unique_edges;
        foreach my $e (@Edges) {
               my $edge0 = @$e[0];
               my $edge1 = @$e[1];
            if ( !$bootstrap ) {
                $perror{ @$e[0] . "-" . @$e[1] } = 0.0;
                $pvalue { @$e[0] . "-" . @$e[1] } = 0.0;
            }
            print OUTFILE @$e[0], " ", @$e[1], ", 1, ", scalar(@Edges), ", ", $count, ", ";
            print OUTFILE $pvalue{ @$e[0] . "-" . @$e[1] }, ", ", $perror{ @$e[0] . "-" . @$e[1] }, "\n";
            print OUTDOTFILE "\"", @$e[0],"\"", " -- ", "\"", @$e[1], "\"", " [weight = \"", $network->get_edge_weight($edge0, $edge1),"\" ];\n";
        }
          print OUTFILE "---------------------------\n";
    }
   
}

print OUTFILE "\n";
print OUTDOTFILE "}\n";

###########################################
if ($excel) {
    print "\n\rNow printing excel output file... \n\r";
    ## nothing yet
}
my $timediff= Time::HiRes::gettimeofday() - $start;
$screen and $scr->at(14,1)->puts( "Time for processing: $timediff seconds");
print "\n\r\n\r\n\r\n\r";

__END__

    =head1 IDSS-with-bootstrap.pl

    Iterative Determinisitic Seriation Solutions with bootstrap significance testing

    =head1 SYNOPSIS

    IDSS-with-bootstrap.pl [options] [file ...]

     Options:
       -help            brief help message
       -man             full documentation
       --input         input file name

    =head1 OPTIONS

    =over 8

    =item B<-help>

    Print a brief help message and exits.

    =item B<-man>

    Prints the manual page and exits.

    =back

    =head1 DESCRIPTION

    B<This program> will read the given input file(s) and do something
    useful with the contents thereof.

    =cut
