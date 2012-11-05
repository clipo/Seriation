#!/usr/bin/perl


use Statistics::Descriptive;
use Time::HiRes;

my $debug = 0;
my $start = Time::HiRes::gettimeofday();

my @assemblages = ();
my @rowtotals = ();
my @labels = ();
my $numrows = 0;
my $cols = 0;

while( <> ) {
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


srand( $start );

my $results = 0;
my $loop = 100;
my $bigloop = 5;
my $loop2 = $loop;
my $ptr1 = 0;
my $ptr2 = 1;


# now do ALL the pairwise assemblage comparisons
# go to sleep and come back later.

while( $ptr1 < $numrows ) {

while( $ptr2 < $numrows ) {

my $stat = Statistics::Descriptive::Full->new();

my @a =  @{ $assemblages[ $ptr1 ] } ;
my @b =  @{ $assemblages[ $ptr2 ] } ;
my $numa = $rowtotals[ $ptr1 ];
my $numb = $rowtotals[ $ptr2 ];


# calculate the directionality for later
my @direction = ();
my $num = scalar( @a );
my $c = 0;
for ( $c = 0; $c < $num; $c++ ) {
	if 		( $a[ $c ] < $b[ $c ] ) { push @direction, -1; }
	elsif 	( $a[ $c ] > $b[ $c ] ) { push @direction, +1; }
	else 							{ push @direction, 0;  }
}

 

my ( @cum_a, @cum_b, $count );
$classes = scalar( @a ); 
my $index_a = 0;
my $total_a = 0.0;
my $count = 0;
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

#print "(debug) cycle value: $cycle\n";

# now we loop 1000 times and resample
$loop = $loop2;
while( $loop ) {
my $assemsize = $numa;
my @assem_a = ();
my $class;
my $total = scalar( @a );
my $rand;
while( $assemsize ) {
	#$rand = ( truly_random_value() % 10000 ) / 10000 ; 
	$rand = rand;
	$class = 0;
	while(( $class < $index_a ) && ( $rand > $cum_a[ $class ] )) {
		$rand -= $cum_a[ $class ];
		$class++;
	}
	push @assem_a, $class;
 	$assemsize--; 
}


my $assemsize = $numb;
my @assem_b = ();
$total = scalar( @b );
while( $assemsize ) {
	#$rand = ( truly_random_value() % 10000 ) / 10000 ; 
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




# calculate the directionality for later
my @dir = ();
my $num = scalar( @ahat );
my $c = 0;
for ( $c = 0; $c < $num; $c++ ) {
	$debug && print "loop $loop ",$ahat[ $c ] - $bhat[ $c ],"\t"; 
	if 		( $ahat[ $c ] < $bhat[ $c ] ) { push @dir, -1; }
	elsif 	( $ahat[ $c ] > $bhat[ $c ] ) { push @dir, +1; }
	else 							{ push @dir, 0;  }
}
$debug && print "\n";


# compare the two sets of results
$num = scalar( @dir );
$c = 0;
my $diff = 0;
for ( $c = 0; $c < $num; $c++ ) {
	$debug && print "loop $loop ",$direction[ $c ],"/",$dir[ $c ],"\t"; 
	if ( $direction[ $c ] == $dir[ $c ] ) { next; }
	$diff++;
}
$debug && print "\n";
if ( $diff == 0 ) { $results++; }

$loop--;
}

#print "Results:  $results matches of $loop2 trials\n";
#print "Probability: ", $results / $loop2, "\n"; 

$stat->add_data( $results/ $loop2 );
$cycle--;

$results = 0;

}

print $labels[ $ptr1], "\t", $labels[ $ptr2 ], "\t";
print $stat->mean(),"\t"; 
print $stat->standard_deviation(),"\n"; 

undef $stat;
$ptr2++;
}
$ptr1++;
$ptr2 = $ptr1 + 1;

}


