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
my $loop    = 10;

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

my ($typecount, $count, $cols);
while (<INFILE>) {
    #print;
    chomp;
    my @line = split /\s+/, $_;
    my $label = shift @line;
    if ($label) {
        push @labels, $label;
        push @assemblageNumber, $count;
        $cols = scalar(@line);
        my $rowtotal = 0;
        for (@line) {
            $rowtotal += $_;
        }
        push @rowtotals, $rowtotal;

        #print "rowtotal: $rowtotal\n";
        my $freq = [];
        $typecount=0;
        for (@line) {
            # push @freq, $_;
            my $f = $_ / $rowtotal;
            my $ff = sprintf( "%.4f", $f );
            push @$freq, $ff;
            $typecount++;
        }
        push @assemblages, [@$freq];
        $assemblageFrequencies{$label} = [@$freq];
        $assemblageSize{ $label }= $rowtotal;

        $count++;
    }
    #print "---- row end ----\n";
}

print Dumper(\%assemblageSize);
my $numrows = scalar(@assemblages);
$DEBUG and print Dumper(\@assemblages),"\n";
$DEBUG and print Dumper( \%assemblageFrequencies ),"\n";
$numrows = scalar(@assemblages);
# start the clock to track how long this run takes
my $start = Time::HiRes::gettimeofday();
srand($start);
my %perror  = ();
my %pvalue  = ();  
my $results = 0;

my $bigloop = 5;
my $loop2   = $loop;

my $classes = 0;

my $pairSet = Math::Combinatorics->new( count => 2, data => \@labels);
while ( my @pairs = $pairSet->next_combination ) {
    #print Dumper(\@pairs);
    my $stat = new Statistics::PointEstimation;
    my @a    = $assemblageFrequencies{ $pairs[0]  };
    my @b    = $assemblageFrequencies{ $pairs[1]  };
    my $numa = $assemblageSize{ $pairs[0] };
    my $numb = $assemblageSize{ $pairs[1] };
    print Dumper(@a);
    print Dumper(@b);
    #exit();
    print "numa: $numa\n";
    print "numb: $numb\n";
    
    # calculate the directionality for later
    my @direction = ();
    my $num       =$typecount;
    print "Num: ", $num, "\n";

    my $c         = 0;
    for ( $c = 0 ; $c < $num ; $c++ ) {
        if    ( $a[$c] < $b[$c] ) { push @direction, -1; }
        elsif ( $a[$c] > $b[$c] ) { push @direction, +1; }
        else                      { push @direction, 0;  }
    }

    my ( @cum_a, @cum_b, $count );
    $classes = $typecount;
    my $index_a = 0;
    my $total_a = 0.0;
    $count   = 0;
    for ( $count = 0 ; $count < $classes ; $count++ ) {
        $cum_a[$index_a] = $a[$count];
        $total_a += $a[$count];
        $index_a++;
    }
    $classes = $typecount;
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
    
    print Dumper(@cum_a);
    print Dumper(@cum_b);
    exit();
    
    # now we loop 100 times and keep track
    my $cycle = $bigloop;
    while ($cycle) {
        
     # now we loop $loop (default 1000) times and resampl
        $loop = $loop2;
        while ($loop) {
            my $assemsize = $numa;
            my @assem_a   = ();
            my $class;
            my $total = scalar(@a);
            my $rand;
            while ($assemsize) {
                $rand  = rand;
                $class = 0;
                while (( $class < $index_a ) && ( $rand > $cum_a[$class] ) )
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
            print "Assemblage size: ", $assemsize, "\n";
            while ($assemsize) {
                $rand  = rand;
                $class = 0;
                while (( $class < $index_b ) && ( $rand > $cum_b[$class] ) )
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

            for (\@assem_a) {
                $aholder{ $_ }++;
            }

            for (\@assem_b) {
                $bholder{ $_ }++;
            }

            for ( keys %aholder ) {
                $ahat[$_] = ( $aholder{$_} / $numa );
            }

            for ( keys %bholder ) {
                $bhat[$_] = ( $bholder{$_} / $numb );
            }
            print Dumper(\@ahat);
            print Dumper(\@bhat);
            # calculate the directionality for later
            my @dir = ();
            my $num = scalar(@ahat);
            my $c   = 0;
            for ( $c = 0 ; $c < $num ; $c++ ) {
                if    ( $ahat[$c] < $bhat[$c] ) { push @dir, -1; }
                elsif ( $ahat[$c] > $bhat[$c] ) { push @dir, +1; }
                else                            { push @dir, 0; }
            }

            # compare the two sets of results
            $num = scalar(@dir);
            $c   = 0;
            my $diff = 0;
            for ( $c = 0 ; $c < $num ; $c++ ) {
                if ( $direction[$c] == $dir[$c] ) { next; }
                $diff++;
            }
            if ( $diff == 0 ) {
                $results++;
            }
            $loop--;
        }
        $stat->add_data( $results / $loop2 );
        $cycle--;
        $results = 0;
    }
    $pvalue{ $pairs[0] } = $stat->mean();
    $perror{ $pairs[0] } = $stat->standard_deviation();
    $pvalue{ $pairs[1] } = $stat->mean();
    $perror{ $pairs[1] } = $stat->standard_deviation();
    print OUTFILE $pairs[0], "\t", $pairs[1], "\t";
    print OUTFILE $stat->mean(), "\t";
    print OUTFILE $stat->standard_deviation(), "\n";
    undef $stat;
    exit();
}



