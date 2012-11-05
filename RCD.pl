#!/usr/bin/perl
use strict;
use Data::Dumper;
use Graph;
use Graph::Writer::VCG;
use Graph::Writer::Dot;
use Math::Combinatorics;
use Time::HiRes;
use Array::Utils qw(:all);

use Statistics::Descriptive;

my $start = Time::HiRes::gettimeofday();

# define some key lists and vars
my $debug = 0;
my @assemblages = ();
my @rowtotals = ();
my @labels = ();
my $numrows = 0;
my $cols = 0;
my %collections;
my @assemblageNumber= ();
my @arrayOfSeriations=();
my %assemblageFrequencies={};
my $maxnumber;

print "Input file name (e.g., input.txt) : ";
my $file = <STDIN>;
chomp($file);
open (INFILE, $file) or die "Cannot open $file.\n";

open (OUTFILE, ">$file.vna") or die "Can't open file $file.vna.\n";
my $count;



#################################################
# the input of the classes -- note these are absolute counts, not frequencies 
# might want to change that...
my $count=0;

while( <INFILE> ) {
	#print;
	chomp;
	my @line = split /\s+/, $_;
	my $label = shift @line;
	push @labels, $label;
	push @assemblageNumber, $count;
	$cols = scalar( @line );
	my $rowtotal = 0;
	for ( @line ) {
		$rowtotal += $_;				
	}
	push @rowtotals, $rowtotal;
	#print "rowtotal: $rowtotal\n";
	my $freq = [];
	for ( @line ) {	
#		push @freq, $_;
		my $f = $_ / $rowtotal;
		my $ff = sprintf("%.2f", $f);
		push @$freq, $ff;
	}
	push @assemblages, [ @$freq ];
	$assemblageFrequencies{ $label }= $freq ; 
	$count++;
	#print "---- row end ----\n";
}


====================================


srand( $time);
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

$stat->AddData( $results/ $loop2 );
$cycle--;

$results = 0;

}

print $labels[ $ptr1], "\t", $labels[ $ptr2 ], "\t";
print $stat->Mean(),"\t"; 
print $stat->StandardDeviation(),"\n"; 

undef $stat;
$ptr2++;
}
$ptr1++;
$ptr2 = $ptr1 + 1;

}




my $maxSeriations=$count;
my $directionstate;
$numrows = scalar( @assemblages );
my @nets;
my @triplettype;
my @tripletNames = ();
my @tripletArray =();
my $net = Graph::Undirected->new;
my $comparison12;
my $comparison23;
my $error;

my $numberOfTriplets;
my $currentMaxSeriationSize = 3;

my $permutations = Math::Combinatorics->new(count => 3,
                                            data => [@assemblageNumber],
                                                             );
	
  	#print "combinations of 3 from: ".join(" ",@assemblages)."\n";
  	#print "------------------------".("--" x scalar(@assemblages))."\n";
  	

	while(my @permu = $permutations->next_combination) {
	
  	my $directionrank;
  	my $newstate;
  	
	#print $labels[$permu[0]]." * ".$labels[$permu[1]]." * ".$labels[$permu[2]]."\n";
	my $tripletname = $labels[$permu[0]]." * ".$labels[$permu[1]]." * ".$labels[$permu[2]];
	$comparison12="";
	$comparison23="";
	$error=0;
	for (my $i =0; $i < $cols; $i++) {
  		my $difscore;
  		my $difscore2;	
		
		my $ass1 = $assemblages[$permu[0]][$i];
		my $ass2 = $assemblages[$permu[1]][$i];
		my $ass3 = $assemblages[$permu[2]][$i];
		#print $ass1."-".$ass2."-".$ass3."\n";

    	my $dif1 = $ass1 - $ass2;
		#my $dif1 = $permu[0][$i] - $permu[1][$i];
    	if ($dif1 < 0 ) { $difscore = -1; }
    	if ($dif1 >0 ) { $difscore = 1; }
    	if ($dif1 == 0) { $difscore = 0; }
    
    	my $dif2 = $ass2 - $ass3;
    	#my $dif2 = $permu[1][$i] - $permu[2][$i];
    	if ($dif2 < 0 ) { $difscore2 = -1; }
    	if ($dif2 >0 ) { $difscore2 = 1; }
    	if ($dif2 == 0) { $difscore2 = 0; }
  	
   		## F1 > F2 < F3 ## criteria not met
   		if (($difscore == 1) && ($difscore2 == -1)) {
  			$comparison12 .= "Z";
  			$comparison23 .= "Z";
			$newstate="Z";
			$error++; 
    	} elsif 
    		## F1 > F2 > F3
    		(($difscore ==1 ) && ($difscore2 == 1 )) {
    		$directionrank++;
    		$newstate ="D";
			$comparison12 .= "U";
  			$comparison23 .= "D";

    	} elsif  
    		# F1 < F2 >F3
      		(($difscore ==-1) && ($difscore2 == 1 )) {
    		$directionrank++;
    		$newstate ="X";
			$comparison12 .= "D";
  			$comparison23 .= "D";
    	} elsif 
    		# F1 < F2 < F3
        	(($difscore == -1 ) && ($difscore2 == -1)) {
    		$directionrank++;
    		$newstate ="U";
			$comparison12 .= "D";
  			$comparison23 .= "U";
    	} elsif 
    		## F1 > F2 = F3 or F1= F2 > F3 or F1 < 
    		(($difscore == 0 ) || ($difscore2 ==0)) { 
    		$newstate="M";
			$comparison12 .= "U";
  			$comparison23 .= "D";
			$directionrank++;
    	} else {
			$error++;
    		#print "Not Match! (other)\n";
    	}
		$directionstate .= $newstate;
  	} 	
  	  
  	if ($error==0) 	{		
		#print "comp 12: ", $comparison12, "\n";
		#print "comp 23: ", $comparison23, "\n";

		undef $net;

		$net = Graph::Undirected->new;
		$net->add_vertex($labels[$permu[0]]);
		$net->add_vertex($labels[$permu[1]]);
		$net->add_vertex($labels[$permu[2]]);
		$net->add_weighted_edge( $labels[$permu[0]], $labels[$permu[1]],$comparison12);
		$net->add_weighted_edge( $labels[$permu[1]], $labels[$permu[2]],$comparison23);
		
  		push @nets, $net;
		$numberOfTriplets++;
  		#print "SCORE:  $directionrank\n";
  		#print $directionstate, "\n";
		#print "new triplet.\n";
		#print $net, "\n\n";
  	}
  	$directionrank==0;
  	$directionstate="";
	$error=0;
  	
  }

## now go through the combination of triplets that are valid (i.e., allthose that have no Z in them)

my @array =();

## now go through all of the assemblages and each one to the top and bottom of the triplets to see if they fit.
## if they do, add them to the network 

#foreach  my $network ( @nets ) {
#	print Dumper($network);
#}


my $currentMaxSeriations=3;
#print "Number of Triplets:  ",$numberOfTriplets, "\n";
my $maxEdges =0;
my $stepcount=0;
my @newnets=();
while ($currentMaxSeriationSize < $maxSeriations) {
	$currentMaxSeriationSize++;
	print "__________________________\n";
	print "Minimum size of seriation solutions:  $currentMaxSeriationSize \n";
	print "Number of current solutions: ", scalar(@nets), "\n";
	my $match=0;
	$stepcount++;
	foreach my $label ( @labels ){
		
		#print "now checking: ", $label,"\n";
		my $index=0;
		my @networks =();
		#if ($stepcount==1) { 
		#		@networks = @nets;
		#} else {
		#		@networks = @newnets;
		#}
		foreach  my $network ( @nets ) {
			#print $network;
			## if its not already there.
			#print "has vertex?", $network->has_vertex( $label), "\n";
			if (!$network->has_vertex( $label )){	
				# get the exterior vertices (should be 2)
				my @V = $network->vertices;
				foreach my $v ( @V) {
					#print "vertice: ", $v, "\n";
					my @Edges = $network->edges_at($v);
					my $edges=scalar(@Edges);
					#print "nubmer of edges:  ", $edges,"\n";
					my @newassemblage=();
					my @oldassemblage=();
					if ($edges==1) {
						#print $label, " is on the edge!\n";
						@newassemblage = @{$assemblageFrequencies{$label}};
						@oldassemblage = @{$assemblageFrequencies{ $v }};
						my @edge = $network->edges_at($v);
						#print Dumper @edge;
						#print "$edge[0][0] -> $edge[0][1],\n";
						#print "_____________________________________\n";
						my $g = $network->get_edge_weight($edge[0][0], $edge[0][1]);
						#print Dumper $g;
						my @comparison=split //, $g;
						my $cols = scalar(@newassemblage), "\n";
						my $comparisonMap;
						my $error=0;
						# go through the columns
						for (my $i =0; $i < $cols; $i++) {
  							my ($difscore, $difscore2);	
							my $val1 = $newassemblage[$i];
							my $val2 = $oldassemblage[$i];
							#print "value 1:  ", $newassemblage[$i], " value 2: ", $oldassemblage[$i], "\n";
    						my $dif1 = $newassemblage[$i] - $oldassemblage[$i];
							#my $dif1 = $permu[0][$i] - $permu[1][$i];
					    	if ($dif1 < 0 ) { $difscore = -1; }
					    	if ($dif1 >0 ) { $difscore = 1; }
					   		if ($dif1 == 0) { $difscore = 0; }
							#print "comparison is:  ", $comparison[$i], " a score of: ", $difscore, "\n";
					  	
					   		## F1 > F2 < F3 ## criteria not met
					   		if (($difscore == 1) && ($comparison[$i] =~ "U")) {
					  			$comparisonMap .= "U";
							} elsif
								(($difscore ==1) && ($comparison[$i] =~ "M")) {
									# this is okay - its a match and the new value is greater. New value shoudl be U
								$comparisonMap .= "U";
					    	} elsif 
					    		## F1 > F2 > F3
								## error the new value is greater but shoudl be less. Error!
					    		(($difscore ==1 ) && ($comparison[$i] =~ "D" )) {
								#print "mismatch!\n";
								$error++;
					    	} elsif  
					    		## new score is less and the Comparison is up .. Error!
					      		(($difscore ==-1) && ($comparison[$i] =~ "U" )) {
								#print "mismatch!\n";
								$error++;
					    	} elsif 
					    		## new score is less and the comparison is down .. Yes!
					      		(($difscore ==-1) && ($comparison[$i] =~ "D" )) {
								$comparisonMap .= "D";
					    	} elsif 
					    		# new score is less but comparison is Match. Okay
					    		(($difscore == -1 ) || ($comparison[$i] =~ "M" )) { 
								$comparisonMap .= "D";
					    	} elsif 
					    		# new score is match but comparison is Match. Okay
					    		(($difscore == 0 ) || ($comparison[$i] =~ "M" )) { 
								$comparisonMap .= "M";
					    	} else {
								#print "lost mismatch!\n";
								$error++;
					    		#print "Not Match! (other)\n";
					    	}
  						} 						
						if ($error==0){
							## no errors so it can be added to the network
							$network->add_weighted_edge( $label, $v, $comparisonMap);
							my $e=$network->edges;
							if ($e>$maxEdges) {
								$maxEdges =$e;
								print "New maximum edges: $e\n";
								print $network, "\n";
							}
							#print "match!\n";
							$match++;
							#print $network, "\n";
							## copy this solution to the new array of networks
							push @nets, $network;
							#splice(@nets, $index, $index, $network);
							#push @nets, $network;
						}
					}
				}
			}

		}
	}
	## no match at this point so no point in going forward.
	if ($match==0) {
		print "Maximum seriation size reached - no more assemblages added that iteration. Max # assemblages is: ", $maxEdges, "\n";
		$maxnumber = $currentMaxSeriationSize - $maxSeriations;
		## to break now...
		$currentMaxSeriationSize = $maxSeriations;
	}
	print "Number of current solutions now: ", scalar(@nets), "\n";
}

# now do some weeding. Basically start with the first network that is the largest, and work backwards. Ignore any
# network that is already represented in the smaller ones since these are trivial (e.g., A->B->C->D alreacy covers
# A->C->D.) That should then leave all the unique maximum solutions (or so it seems)

## first need to sort the networks by size
my @filteredarray=();

print "Filtering solutions so we only end up with the unique ones.\n";
print "Start with ", scalar(@nets), " solutions. \n";
my $count=0;
my $edgecount =$maxEdges;
while ($edgecount>1) {
	print "### not on sets of $edgecount \n";
	foreach  my $network ( @nets ) {
		my $E = $network->edges;
		if ($E==$edgecount) {
			if ($edgecount==$maxEdges) { ## include all the largest ones... 
				push @filteredarray, $network;	
				$count++;
				#print "Up to $count solutions (still in initial set)\n";
			} else {  ## only filter the ones less than the max
					## check that the elements of this network are not all represented in the larger sets
					## first get the larger set

				my @V=$network->vertices; ##array of elements
				my $filteredNumber=scalar(@filteredarray);
				my $extracount=0;
				foreach my $fnetwork (@filteredarray) {
					## now check for subsets
					my @fV=$fnetwork->vertices;
					## is @V included in @fV if no, then push @
					my %diff ={};
					#print "network: ", Dumper(@V), "\n";
					#print "filtered: ", Dumper(@fV), "\n";
					@diff{ @V } = @V;
					delete @diff{ @fV };
					my @k = (keys %diff);
					my $extraAssemblages = scalar(@k);
					#print Dumper(@k), "\n";
					#print "extraAssemblages: $extraAssemblages\n";
					if ($extraAssemblages>2) {
						$extracount++;
						print "Not found in the network: $extracount\n";
					}
				}
				if ($extracount==$filteredNumber) {
					print "solution doesnt exist in any of the existing networks - so add to list. \n";
					push @filteredarray, $network;
				}
			}
						
		}
	}
$edgecount--;
}
print "End with ", scalar(@filteredarray), " solutions.\n";


my $writer=Graph::Writer::VCG->new();
$count=0;
my $name;
foreach  my $network (@newnets ) {
	my $E = $network->edges;

	if ($E>$maxnumber-1) {
		$count++;
		$name= $count.'.vcg';
		#$writer->write_graph($network, $name);
	}
	#print Dumper($network);
}

my $writer2=Graph::Writer::Dot->new();
$count=0;
foreach  my $network ( @newnets ) {
	my $V = $network->vertices;
	if ($V==$maxEdges) {
		$count++;
		$name= $count.'.dot';
		#$writer2->write_graph($net, $name);
	}
	#print Dumper($network);
}


print "Now printing output file... \n";

print OUTFILE "*Node data\n";
$count=0;
foreach my $l (@labels){
	print OUTFILE $l, "\n";
}
print OUTFILE "*Tie data\n";
print OUTFILE "From To Edge Weight\n";
foreach  my $network ( @filteredarray ) {
	$count++;
	print "$count: $network\n";
	#my $E = $network->edges;
	#if ($E==$maxEdges) {
		my @Edges=$network->edges;
		foreach my $e (@Edges) {
			print OUTFILE @$e[0], " ", @$e[1], ", 1, ", scalar(@Edges), "\n";
			#print @$e[0], " ", @$e[1], "\n";
		} 
	#}
	#print Dumper($network);
}
print OUTFILE "\n";
my $end = Time::HiRes::gettimeofday();
print "Time for processing:  ";
printf("%.2f", $end - $start);
print " seconds.\n"; 