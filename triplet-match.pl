#!/usr/bin/perl
use strict;
use Data::Dumper;
use Graph;
use Graph::Writer::VCG;
use Graph::Writer::Dot;
use Math::Combinatorics;
use Time::HiRes;

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


# the input of the classes -- note these are absolute counts, not frequencies 
# might want to change that...
my $count=0;

while( <> ) {
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
while ($currentMaxSeriationSize < $maxSeriations) {
	$currentMaxSeriationSize++;
	print "__________________________\n";
	print "MAX SIZE:  $currentMaxSeriationSize \n";
	my $match=0;
	foreach my $label ( @labels ){
		
		#print "now checking: ", $label,"\n";
		my $index=0;
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
							splice(@nets, $index, $index, $network);
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
}

my $writer=Graph::Writer::VCG->new();
my $count=0;
my $name;
foreach  my $network ( @nets ) {

	my $V = $network->vertices;

	if ($V>$maxnumber-1) {
		$count++;
		$name= $count.'.vcg';
		#$writer->write_graph($network, $name);
	}
	#print Dumper($network);
}

my $writer2=Graph::Writer::Dot->new();
$count=0;
foreach  my $network ( @nets ) {
	my $V = $network->vertices;
	if ($V>$maxnumber-1) {
		$count++;
		$name= $count.'.dot';
		#$writer2->write_graph($net, $name);
	}
	#print Dumper($network);
}

open (OUTFILE, ">>seriations-out.vna") or die "Can't open file seriations-out.txt.\n";
print OUTFILE "*Node data\n";
$count=0;
foreach my $l (@labels){
	print OUTFILE $l, "\n";
}
print OUTFILE "*Tie data\n";
foreach  my $network ( @nets ) {
	my $E = $network->edges;
	if ($E==$maxEdges) {
		my @Edges=$network->edges;
		foreach my $e (@Edges) {
			print OUTFILE @$e[0], " ", @$e[1], ", 1, 1\n";
			#print @$e[0], " ", @$e[1], "\n";
		} 
	}
	#print Dumper($network);
}

my $end = Time::HiRes::gettimeofday();
print "Time for processing:  ";
printf("%.2f", $end - $start);
print " seconds.\n"; 