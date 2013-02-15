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

my ($debug, $inputfile, $man, $bootloop);
$bootloop    = 1000;

GetOptions(
    'debug'                     => \$debug,
    'help'                      => sub { HelpMessage() },
    'input=s'                   => \$inputfile,
    'samples=s'                 => \$bootloop,
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
print "Number of assemblage: ", $numrows, "\n";
my $results = 0;
my $ptr1 = 0;
my $ptr2 = 1;

# now do ALL the pairwise assemblage comparisons
# go to sleep and come back later.


my $pairSet = Math::Combinatorics->new( count => 2, data => \@labels);
while ( my @pairs = $pairSet->next_combination ) {
	
#while( $ptr1 < $numrows ) {
#	while( $ptr2 < $numrows ) {
	#	my $stat = new Statistics::PointEstimation;
	#	my @a =  @{ $assemblages[ $ptr1 ] } ;
	#	my @b =  @{ $assemblages[ $ptr2 ] } ;
	#	my $numa = $rowtotals[ $ptr1 ];
	#	my $numb = $rowtotals[ $ptr2 ];
		
	my $stat = new Statistics::PointEstimation;
	my @a    = $assemblageFrequencies{ $pairs[0]  };
	my @b    = $assemblageFrequencies{ $pairs[1]  };
	my $numa = $assemblageSize{ $pairs[0] };
	my $numb = $assemblageSize{ $pairs[1] };
    
    
		#print "comparing: ", $labels[ $ptr1 ], " with: ", $labels[$ptr2], "\n";
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
			
		# now we loop a bunch of times times and keep track

		my $cycle = $bootloop;
		while( $cycle ) {
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
				$debug && print "loop $cycle ",$ahat[ $c ] - $bhat[ $c ],"\t";
				$cumulativeDiff += abs($ahat[ $c ] - $bhat[ $c ])
			}
			$debug && print "\n";
		

			#print "Results:  $results matches of $loop2 trials\n";
			#print "Probability: ", $results / $loop2, "\n"; 
		
			$stat->add_data( $cumulativeDiff );
			$cycle--;
		
			$results = 0;
		}
		#print $labels[ $ptr1 ] , "\t", $labels[ $ptr2 ], "\t";
		print OUTFILE $pairs[0], "\t", $pairs[1], "\t";
		print  $stat->mean(),"\t";
		print  $stat->standard_deviation(),"\n"; 

		#print OUTFILE $labels[ $ptr1 ] , "\t", $labels[ $ptr2 ], "\t";		
		print OUTFILE $pairs[0], "\t", $pairs[1], "\t";
	#	print OUTFILE $stat->mean(), "\t";
	#	print OUTFILE $stat->standard_deviation(), "\n";
		print OUTFILE $stat->mean(),"\t"; 
		print OUTFILE $stat->standard_deviation(),"\n";

		## opposite combination too	
	#	print OUTFILE $labels[ $ptr2 ] , "\t", $labels[ $ptr1 ], "\t";
		print OUTFILE $pairs[0], "\t", $pairs[1], "\t";
		print OUTFILE $stat->mean(),"\t"; 
		print OUTFILE $stat->standard_deviation(),"\n"; 
		undef $stat;
	#	$ptr2++;
	#}
	#$ptr1++;
	#$ptr2 = $ptr1 + 1;
}


