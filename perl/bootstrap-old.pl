#!/usr/bin/perl

#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long qw(HelpMessage);
use Pod::Usage;
use Time::HiRes;
use Array::Utils qw(:all);
use Statistics::Descriptive;
use Statistics::PointEstimation;
use Math::Combinatorics;

# process command line options; if none, print a usage/help message.
# note - manual page for explaining the options, what they do, and how to use them
# is at the bottom of this file in simple POD format.

my ($debug, $inputfile, $man);
my $loop    = 100;

GetOptions(
    'debug'                     => \$debug,
    'help'                      => sub { HelpMessage() },
    'input=s'                   => \$inputfile,
    'samples=s'                 => \$loop,
    man                         => \$man
) or pod2usage(2);

my $DEBUG = $debug;    # our "$debug level"

if ($DEBUG) {
    print "Verbose debugging output is on!!!\n";
    print "Processing input file: $inputfile\n";
}

my @labels;
my @assemblages;
my %assemblageFrequencies;
my %assemblageSize;
my @assemblageNumber;
my @rowtotals;

## open the data file
my $useOutputFile = substr($inputfile,0,-4);
open( INFILE, $inputfile ) or die "Cannot open $inputfile.\n";
open( OUTFILE, ">$useOutputFile-bootstrap.txt" ) or die "Can't open file $useOutputFile.vna to write.\n";


my $numrows = 0;
my $cols = 0;

while( <INFILE> ) {
	print;
	chomp;
	my @line = split /\s+/, $_;
	my $label = shift @line;
	push @labels, $label;
	$cols = scalar( @line );
	my $rowtotal = 0;
	for ( @line ) {
		$rowtotal += $_;				
	}
	push @rowtotals, $rowtotal;
	#print "rowtotal: $rowtotal\n";
	my @freq = ();
	for ( @line ) {
		push @freq, ( $_ / $rowtotal );
	}
	push @assemblages, [ @freq ];
	#print "---- row end ----\n";
}
$numrows = scalar( @assemblages );




my $results = 0;
my $bigloop = 5;
my $loop2 = $loop;
my $ptr1 = 0;
my $ptr2 = 1;


# now do ALL the pairwise assemblage comparisons
# go to sleep and come back later.

while( $ptr1 < $numrows ) {
	while( $ptr2 < $numrows ) {
		my $stat = new Statistics::PointEstimation;
		my @a =  @{ $assemblages[ $ptr1 ] } ;
		my @b =  @{ $assemblages[ $ptr2 ] } ;
		my $numa = $rowtotals[ $ptr1 ];
		my $numb = $rowtotals[ $ptr2 ];
		
		# calculate the directionality for later
		my @direction = ();
		my $num = scalar( @a );
		my $c = 0;
		for ( $c = 0; $c < $num; $c++ ) {
			if 	( $a[ $c ] < $b[ $c ] ) { push @direction, -1; }
			elsif 	( $a[ $c ] > $b[ $c ] ) { push @direction, +1; }
			else 	{ push @direction, 0;  }
		}
		
		my ( @cum_a, @cum_b, $count );
		my $classes = scalar( @a ); 
		my $index_a = 0;
		my $total_a = 0.0;
		$count = 0;
		for( $count = 0; $count < $classes; $count++ ) {
			$cum_a[ $index_a ] = $a[ $count ];
			$total_a += $a[ $count ];
			$index_a++;
		}
		$classes = scalar( @b ); 
		my $index_b = 0;
		my $total_b = 0.0;
		$count = 0;
		for( $count = 0; $count < $classes; $count++ ) {
			$cum_b[ $index_b ] = $b[ $count ];
			$total_b += $b[ $count ];
			$index_b++;
		}
		
		$index_a--;
		$index_b--;
			
		# now we loop 100 times and keep track
		my $cycle = $bigloop;
		while( $cycle ) {
			# now we loop 1000 times and resample
			$loop = $loop2;
			while( $loop ) {
				my $assemsize = $numa;
				my @assem_a = ();
				my $class;
				my $total = scalar( @a );
				my $rand;
			
				# start the clock to track how long this run takes
				my $start = Time::HiRes::gettimeofday();
				srand($start);
				
				while( $assemsize ) {
				
					$rand = rand;
					$class = 0;
					while(( $class < $index_a ) && ( $rand > $cum_a[ $class ] )) {
						$rand -= $cum_a[ $class ];
						$class++;
					}
					push @assem_a, $class;
					$assemsize--; 
				}
				
				$assemsize = $numb;
				my @assem_b = ();
				$total = scalar( @b );
				while( $assemsize ) {
				
					$rand = rand;
					$class = 0;
					while(( $class < $index_b ) && ( $rand > $cum_b[ $class ] )) {
						$rand -= $cum_b[ $class ];
						$class++;
					}
					push @assem_b, $class;
					$assemsize--; 
				}
				
				
				my ( @ahat, @bhat, %aholder, %bholder );
				%aholder = ();
				%bholder = ();
				my $index = 0;
				
				for ( $index = 0; $index < $cols; $index++ ) {
					$ahat[ $index ] = 0;
					$bhat[ $index ] = 0;
				}
				
				
				for ( @assem_a ) {
					$aholder{ $_ }++;
				}  
				
				for ( @assem_b ) {
					$bholder{ $_ }++;
				}  
				
				for ( keys %aholder ) {
					$ahat[ $_ ] = ( $aholder{ $_ } / $numa );
				}
				
				for ( keys %bholder ) {
					$bhat[ $_ ] = ( $bholder{ $_ } / $numb );
				}
				
				# calculate differences in % between assembalges
				my @dir = ();
				my $num = scalar( @ahat );
				my $c = 0;
				my $cumulativeDiff = 0;
				for ( $c = 0; $c < $num; $c++ ) {
					$debug && print "loop $loop ",$ahat[ $c ] - $bhat[ $c ],"\t";
					$cumulativeDiff += abs($ahat[ $c ] - $bhat[ $c ])
				}
				$debug && print "\n";
			
	
				#print "Results:  $results matches of $loop2 trials\n";
				#print "Probability: ", $results / $loop2, "\n"; 
			
				$stat->add_data( $cumulativeDiff );
				$cycle--;
			
				$results = 0;
			}
			print OUTFILE $labels[ $ptr1 ] , "\t", $labels[ $ptr2 ], "\t";
			print OUTFILE $stat->mean(),"\t"; 
			print OUTFILE $stat->standard_deviation(),"\n"; 
			undef $stat;
			$ptr2++;
		}
	}
	$ptr1++;
	$ptr2 = $ptr1 + 1;
}


