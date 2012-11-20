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


my $debug = 0;
my $filterflag = 0;           ## do you want to try to output all of the solutions (filtered for non trivial)
my $largestOnly = 0;  #       # only output the largest set of solutions
my $individualfileoutput=0;   ## create files for all the indivdual networks
my $bootstrap = 0;            ## flag for bootstrap 
my $man = 0;
my $help = 0;
my $inputfile;
my $bootstrapdebug = 0;
my $excel = 0;                ## flag for excel file output

# process command line options; if none, print a usage/help message.  
# note - manual page for explaining the options, what they do, and how to use them
# is at the bottom of this file in simple POD format.  

GetOptions(
  'debug' => \$debug, 
  'bootstrap' => \$bootstrap,
  'bootstrapdebug' => \$bootstrapdebug,
  'filtered' => \$filterflag,
  'largestonly' => \$largestOnly,
  'indivfiles' => \$individualfileoutput,
  'help' => sub { HelpMessage() },
  'input=s' => \$inputfile,
  'excel' => \$excel,
  man => \$man) or pod2usage(2);



my $DEBUG = $debug;   # our "$debug level"


if ($DEBUG) {
   print "Verbose debugging output is on!!!\n";
   print "Processing input file: $inputfile\n";
   $filterflag or print "filterflag is off\n";
   $bootstrap or print "bootstrap is off\n";
   $largestOnly or print "output largest solutions is off\n";
   $individualfileoutput or print "individual network file output is off\n";
   $excel or print "excel output is off"
}

# start the clock
my $start = Time::HiRes::gettimeofday();


# define some key lists and vars
my @assemblages = ();
my @rowtotals = ();
my @labels = ();
my $numrows = 0;
my $cols = 0;
my %collections;
my @assemblageNumber= ();
my @arrayOfSeriations=();
my %assemblageFrequencies={};
my @allNetworks=();
my $maxnumber;

#print "Input file name (e.g., input.txt) : ";
#my $file = <STDIN>;
#chomp($file);
open (INFILE, $inputfile) or die "Cannot open $inputfile.\n";

open (OUTFILE, ">$inputfile.vna") or die "Can't open file $inputfile.vna.\n";
my $count;

# the input of the classes -- note these are absolute counts, not frequencies 
# might want to change that...
my $count=0;

while( <INFILE> ) {
	#print;
	chomp;
	my @line = split /\s+/, $_;
	my $label = shift @line;
   if ($label) {
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
#		   push @freq, $_;
		   my $f = $_ / $rowtotal;
		   my $ff = sprintf("%.4f", $f);
		   push @$freq, $ff;
	   }
	   push @assemblages, [ @$freq ];
	   $assemblageFrequencies{ $label }= $freq ; 
	   $count++;
	}
   #print "---- row end ----\n";
}
my $numrows = scalar( @assemblages );

my $maxSeriations=$count;
print "Maximum possible seriation solution length: ", $count, "\n";
my $directionstate;
$numrows = scalar( @assemblages );
my @nets;
my @triplettype;
my @tripletNames = ();
my @tripletArray =();
my $net = Graph::Directed->new;
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
  	

while (my @permu = $permutations->next_combination) {
	
  	
	$DEBUG and print $labels[$permu[0]]." * ".$labels[$permu[1]]." * ".$labels[$permu[2]]."\n";
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
		$DEBUG and print $ass1."-".$ass2."-".$ass3."\n";

    	my $dif1 = $ass1 - $ass2;
		#my $dif1 = $permu[0][$i] - $permu[1][$i];
    	if ($dif1 < 0 ) { $difscore = -1; }
    	if ($dif1 >0 ) { $difscore = 1; }
    	if ($dif1 == 0) { $difscore = 0; }
    
    	my $dif2 = $ass3 - $ass2;
    	#my $dif2 = $permu[1][$i] - $permu[2][$i];
    	if ($dif2 < 0 ) { $difscore2 = -1; }
    	if ($dif2 >0 ) { $difscore2 = 1; }
    	if ($dif2 == 0) { $difscore2 = 0; }
  	
   	
   	## F1 > F2 < F3 ## criteria not met
      if (($difscore == 1) && ($difscore2 == 1)) {
			  $error++;
    	} elsif (($difscore ==1 ) && ($difscore2 == -1 )) { ## F1 > F2 > F3
			   $comparison12 .= "U";
  			   $comparison23 .= "D";
    	} elsif  
    		(($difscore ==-1) && ($difscore2 == -1 )) {   #  F1 < F2 >F3
    			$comparison12 .= "X";
  			   $comparison23 .= "X";
    	} elsif 
    			(($difscore == -1 ) && ($difscore2 == 1)) {  # F1 < F2 < F3
			   $comparison12 .= "D";
  			   $comparison23 .= "U";
    	} elsif (($difscore == 0 ) && ($difscore2==1)) {
			   $comparison12 .= "M";
            $comparison23 .= "U";
      } elsif (($difscore2 ==0) && ($difscore==1)) { 
            $comparison12 .= "U";
  			   $comparison23 .= "M";
    	} elsif (($difscore == 0 ) && ($difscore2==-1)) {
			   $comparison12 .= "M";
            $comparison23 .= "D";
      } elsif (($difscore==-1) && ($difscore2 ==0)) { 
            $comparison12 .= "D";
  			   $comparison23 .= "M";
      } elsif (($difscore==0) && ($difscore2 ==0)) { 
            $comparison12 .= "M";
  			   $comparison23 .= "M";   
    	} else {
			$error++;
    		print "Not Match! (other)  Difscore 1: $difscore Difscore 2: $difscore2 \n";
         exit();
    	}
  	} 	
  	  
  	if ($error==0) 	{	
		undef $net;
		$net = Graph::Directed->new;
		$net->add_vertex($labels[$permu[0]]);
		$net->add_vertex($labels[$permu[1]]);
		$net->add_vertex($labels[$permu[2]]);
		$net->add_weighted_edge( $labels[$permu[0]], $labels[$permu[1]],$comparison12);
		$net->add_weighted_edge( $labels[$permu[1]], $labels[$permu[2]],$comparison23);
		$DEBUG and print $labels[$permu[0]]." * ".$labels[$permu[1]]."-".$comparison12." * ".$labels[$permu[2]]."-".$comparison23."\n";
  		push @nets, $net;
		$numberOfTriplets++;
  	}
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
my $match=0;
while ($currentMaxSeriationSize < $maxSeriations) {

	$currentMaxSeriationSize++;
	#print "now checking: ", $label,"\n";
	my $index=0;
   my @networks =();
	$stepcount++;
   
   if ($stepcount==1) { 
         foreach my $n (@nets) {
            my $ne=$n->deep_copy_graph;
            push @networks, $ne;
         }
	} else {
         foreach my $n (@newnets) {
            my $ne=$n->deep_copy_graph;
            push @networks, $ne;
			}
	}
   @newnets=();
	print "__________________________\n";
	print "Iteriation Number:  $currentMaxSeriationSize \n";
	print "Number of current solutions: ", scalar(@networks), "\n";
	
   $match=0;
   
   ## look through the set of existing valid networks. 
   foreach  my $network ( @networks ) {
      my $index=0;
      $DEBUG and print "----------------------------------------------------------------------\n";
      $DEBUG and print "Network: ", $network, "\n";   
      $DEBUG and print "----------------------------------------------------------------------\n";

      foreach my $label ( @labels ){
            
         $DEBUG and print "Now checking assemblage: ", $label," to see if it fits on the end of the current solution.\n";
         $DEBUG and print "First check to see if it is included already. If it has, move on.\n";
         
         if (!$network->has_vertex( $label )){	
				
            # get the exterior vertices (should be 2)
				my @V = $network->vertices;  ## list of all the vertices
            
            $DEBUG and print "Find the ends of the network. Do this by getting all the vertices ";
            $DEBUG and print " and looking for the ones with only 1 connection. There should be just 2 here.\n";
         
            ## loop through all of the edges to see if they can be stuck on the ends of the networks.
			   foreach my $v ( @V ) {
               
					$DEBUG and print "Checking vertice: ", $v, "\n";
					my @Edges = $network->edges_at($v);
					my $edges=scalar(@Edges);
					$DEBUG and print "This vertice: $v has this number of edges:  ", $edges,"\n";
					my @newassemblage=();
					my @oldassemblage=();
               my $comparisonMap;
					## only if it has one edge. We only want to consider the ends of the netowrk -- so skip the others.
               ## Just the ends.
               if ($edges==1) {
						$DEBUG and print $v, " is on the edge since it only has one vertice.\n";
						@newassemblage = @{$assemblageFrequencies{$label}};
						@oldassemblage = @{$assemblageFrequencies{ $v }};
						my @edge = $network->edges_at($v);
						$DEBUG and print "Number of edges here at $v:  ", scalar(@edge), " (should be just one).\n";
                  my $connectedAssemblage = $edge[0][1];
						my $g = $network->get_edge_weight($edge[0][0], $edge[0][1]);

						#print Dumper $g;
                  $DEBUG and print  "There should be just 2 vertices here 0: $edge[0][0] and 1: $edge[0][1], with a relation of $g\n";

						my @comparison=split //, $g;

						my $cols = scalar(@newassemblage), "\n";
						my $comparisonMap;
						my $error=0;
						# go through the columns

						for (my $i =0; $i < $cols; $i++) {
  							my ($difscore, $difscore2);	
							my $val1 = $newassemblage[$i];
							my $val2 = $oldassemblage[$i];
   						$DEBUG and print "########  Type $i - Type $i - Type $i - Type $i - Type $i - Type $i - Type $i  ########  \n";	
                     $DEBUG and print "Type $i:  $label 1:", $newassemblage[$i], " $v 2: ", $oldassemblage[$i], "\n";
    						
                     my $dif1 = $newassemblage[$i] - $oldassemblage[$i];							
                     if ($dif1 <0 ) { $difscore = -1; }
					    	if ($dif1 >0 ) { $difscore = 1; }
					   	if ($dif1 == 0) { $difscore = 0; }
                     
                     $DEBUG and print "Type $i: - comparison is:  ", $comparison[$i], " a score of: ", $difscore, "\n";
					  
					   	# if up and previous direction was up, all is good.
					   	if (($difscore == 1) && ($comparison[$i] =~ "U")) {
					  			$comparisonMap .= "U";
                        $DEBUG and print " Type $i: Got a difscore of 1 and a ";
                        $DEBUG and print " comparison of a U. This works. Adding $label to vertices $edge[0][1]\n";

							} elsif
								(($difscore == 1) && ($comparison[$i] =~ "M")) {
									# this is okay - its a match and the new value is greater. New value shoudl be U
						         # need to find what was happening on the previous comparison to know whether this needs
                           # to be up or down. 
                              $DEBUG and print "Type $i: Got a difscore of 1 and a comparison of a M. This could be okay.\n";
                              my $xerror=0;
                              my $stop=0;
                              my $ccount=1;
                              my @EE = $network->unique_edges;
                              my $numEdges = scalar(@EE);
                              my $change;
                              $DEBUG and print " Type $i:   Case A (1, M) : Potentially can add $label ";
                              $DEBUG and print "to vert $v because score is -1 and the comparison is M.\n";
                              $DEBUG and print "###Network is currently: $network\n";
                              $DEBUG and print "Type $i: But need to check the other $numEdges ";
                              $DEBUG and print " comparisons because this will only work if there no X somewhere (or more Ms)\n";
                             
                              $DEBUG and print "Type $i: These combos are already evaluated: ";
                              $DEBUG and print  $edge[0][0]."-".$edge[0][1]. " and ". $edge[0][1]."-".$edge[0][0], "\n";
                              my %checkHash = {};
                              $checkHash { $edge[0][0]."-".$edge[0][1]} = 1;
                              $checkHash { $edge[0][1]."-".$edge[0][0]} =1;
                             ## check all combinations but the one with self!
                             foreach my $ee (@EE) {
                                
                                 if (!$checkHash{ @$ee[0]."-".@$ee[1] } && !$checkHash{ @$ee[1]."-".@$ee[0] })  {
                                       
                                       my $ge = $network->get_edge_weight(@$ee[0],@$ee[1]);
                                       ## turn comparison into an array to find the right element
                                       my @compArray = ();
                                       @compArray = split //, $ge;
                                       
                                       $DEBUG and print "Type $i:Here is what we get for comparisons # $ccount of ($numEdges) ";
                                       $DEBUG and print " between $label and  and $v: ",  $ge, "->", $compArray[$i], "\n";
                                       $DEBUG and print "Type $i:$ccount comparison (of $numEdges) is a $compArray[$i]. \n";

						                     if ($compArray[$i]=~ "X" ) { 
                                          $DEBUG and print "Type $i: $ccount (of $numEdges) comparison is $compArray[$i]  ";
                                          $DEBUG and print " X found -- so this would make it multimodal. Error.\n";
								                  $xerror++;
                                          $ccount=$numEdges;
                                          } else {
                                          
                                          if ($compArray[$i]=~ "D") { $change="X";
                                             #This would be a mode shift so change is X
                                            
                                          } else {
                                             $xerror++;
                                            #$change="U";                                          
                                            #$DEBUG and print "Type $i: Okay since  X in previous edges (just M or U). ";
                                            #$DEBUG and print "Type $i: If D then this is an X. Adding assemblage. \n";
                                            #$DEBUG and print "Type $i: New comparisonMap:  $change\n";
                                          } ## end of comparison array
                                       }
                                } #end of if check
                           } ## end of foreach my $ee
                           if ($xerror>0) {
                                 $error++;
                           } else {
                  		      $comparisonMap .= $change;
                              $DEBUG and print " Type $i: Adding $label to vertices $v because score is 1 ";
                              $DEBUG and print " and the comparison is M but no other Ds or Xs in the previous linkages.\n";
                           }

					    	} elsif 
								## error the new value is greater but shoudl be less. Error!
					    		(($difscore == 1 ) & ($comparison[$i] =~ "D" )) {
								#print "mismatch!\n";
								$error++;
                        $DEBUG and print "Type $i: Value 1:  ", $newassemblage[$i], " value 2: ", $oldassemblage[$i], "\n";
                        $DEBUG and print "Type $i: Comparison is:  ", $comparison[$i], " a score of: ", $difscore, "\n";
                        $DEBUG and print "Type $i: Rejecting $label from $edge[0][1] because value is up and comparison is D.\n";
					    	} elsif  
					    		## new score is less and the Comparison is up .. Error!
                        ## note -- need to check to see if there is a previous change in direction because
                        ## its possible that the this is a mode change.
                        ## do this by logging all modes in the original triplet -- keep track of modes
                        ## per type
					      		(($difscore ==-1) && ($comparison[$i] =~ "U" )) {
                              $DEBUG and print "Got a difscore of -1 and a comparison of a U. This could be an error.\n";
                              ## first check to see if there is already and X in this column somewhere else.
								      #print "mismatch!\n";
                              my @edgelist=();
                              my @EE = $network->unique_edges;
                              my $numEdges = scalar(@EE);
                              my $xerror=0;
                              my $ccount=1;
                              $DEBUG and print "Type $i:  Case B (-1, U). Potentially can add $label and vert $v because ";
                              $DEBUG and print "score is -1 and the comparison is M.\n";
                              $DEBUG and print "Type $i: ###Network is currently: $network\n";
                              $DEBUG and print "Type $i: But need to check the other $numEdges comparisons because ";
                              $DEBUG and print " this will only work if there no X somewhere (or more Ms)\n";
                              
                              my %checkHash = {};
                              $checkHash { $edge[0][0]."-".$edge[0][1]} = 1;
                              $checkHash { $edge[0][1]."-".$edge[0][0]} = 1;
                              ## check all combinations but the one with self!
                              foreach my $ee (@EE) {
                                
                                 if (!$checkHash{ @$ee[1]."-".@$ee[0] } && !$checkHash{ @$ee[0]."-".@$ee[1]})  {  
                                       my $ge = $network->get_edge_weight(@$ee[0], @$ee[1]);
                                       my @compArray = split //, $ge;
                                       $DEBUG and print "Type $i: Here is what we get for # $ccount (of $numEdges) ";        
                                       $DEBUG and print " comparisons between @$ee[0] @$ee[1]: ", $ge, "->", $compArray[$i], "\n";
                                       $DEBUG and print "Type $i: $ccount comparison (of $numEdges)  is a $compArray[$i]. \n";
						                  
                                       if ($compArray[$i]=~ "X") { 
                                          $DEBUG and print " Type $i: I found an X -- so this would make it multimodal. Error.\n";
								                  $xerror++;
                                       
                                       } else {
                                          $DEBUG and print "Okay since not an M, U, or D.\n";
                                    }
                                 $ccount++;
                                 }
                            }
                           if ($xerror) {
                                 $error++;
                                 $DEBUG and print "Type $i: Rejecting $label from $v) because value is down after a ";
                                 $DEBUG and print " previous U (but already has a mode). Multimodal - so error.\n";
                              } else {
                                $comparisonMap .= "X";
                                $DEBUG and print "Type $i:Definitely adding $label to vertices $v because score ";
                                $DEBUG and print " is -1 and the comparison is U but no other Xs in the previous linkages.\n";
                                $DEBUG and print "Type $i: Adding an X to the comparisons for type $i. \n";
                                $DEBUG and print "Comparison map is now $comparisonMap\n";
                             } #end if if check error (xerror)
                        
					    	} elsif 
					    		## new score is less and the comparison is down .. Yes!
					      	(($difscore == -1) && ($comparison[$i] =~ "D" )) {
								$comparisonMap .= "D";
                        $DEBUG and print "Type $i: Adding a D to the comparisons for type $i. Comparison map is now $comparisonMap\n";
					    	} elsif 
					    		# new score is less but comparison is Match. Okay
					    		(($difscore == -1 ) || ($comparison[$i] =~ "M" )) { 
                              my $vHere=$v;
                              my $xerror=0;
                              my $stop=0;
                              my $ccount=1;
                              my @EE = $network->unique_edges;
                              my $numEdges = scalar(@EE);
                              $DEBUG and print "Type $i: Case C (-1, M) Potentially can add $label and vert $v ";
                              $DEBUG and print "because score is -1 and the comparison is M.\n";   
                              $DEBUG and print "Type $i: ##Network is currently: $network\n";
                              $DEBUG and print "Type $i: But need to check the other $numEdges ";
                              $DEBUG and print " comparisons because this will only work if there is a D or an X somewhere (or more Ms)\n";

                              my %checkHash = ();
                              $checkHash { $edge[0][0]."-".$edge[0][1]} = 1;
                              $checkHash { $edge[0][1]."-".$edge[0][0]} = 1;
                              $DEBUG and print Dumper(%checkHash);
                              ## check all combinations but the one with self!
                              foreach my $ee (@EE) {
                                 my $comparisonName= @$ee[0]."-". @$ee[1];

                                 if ($checkHash{ $comparisonName }==0  )  {  
                                       my $ge = $network->get_edge_weight( @$ee[0], @$ee[1] );
                                       my @compArray = split //, $ge;
                                       $DEBUG and print "Type $i: Here is what we get for comparison # $ccount ";
                                       $DEBUG and print " (of $numEdges) between @$ee[0], @$ee[1]: ", $ge, "->", $compArray[$i], "\n";
                                       $DEBUG and print "$ccount comparison (of $numEdges) is a $compArray[$i]. \n";
						                  
                                    if ($compArray[$i]=~ "X") { 
                                       #$DEBUG and print "Type $i: Since I found U found -- so this would be an error. Error.\n";
								               #$xerror++;
                                       #$xcount++;   ## check and see if there is an X yet.  If not, this might 
                                    } else {
                                       $DEBUG and print "Type $i: Okay since an M, X, or D.\n";
                                    } ## end of if checking of the valid condition (U)
                                    $ccount++;
                                    ## keep going back until we see if we get a U -- if not its okay
                                 }
                            }
                           if ($xerror>0) {
                                 $error++;
                                    $DEBUG and print "Type $i: Check error detected (xerror=$xerror) so cant add assemblage.\n";
                           } else {
                  		      $comparisonMap .= "D";
                              $DEBUG and print "Type $i:Adding $label to vertices $v because "; 
                              $DEBUG and print " score is -1 and the comparison is D. ComparisonMap is now $comparisonMap\n";
                           }
					    	} elsif 
					    		# new score is match but comparison is Match. Okay
					    		(($difscore == 0 )) { 
								$comparisonMap .= "M";
                        $DEBUG and print "Type $i: Adding $label to vertices $v because its a match. \n";
                        $DEBUG and print "Type $i: ComparisonMap is now $comparisonMap\n";
                     } elsif
                        # newscore is down but comparison is X. This means that there was already a peak
                        (($difscore== -1) && ($comparison[$i]=~ "X")) {
                        ## this is okay since it is down from a mode peak
                           $DEBUG and print "Type $i:Adding $label to vertices $v because ";
                           $DEBUG and print " score is -1 and the comparison is D. ComparisonMap is now $comparisonMap\n";
                           $comparisonMap .= "D";
                     } elsif 
                        (($difscore == 1) && ($comparison[$i] =~ "X")) {
                        ## new score is up but comparison is X.. no cant work because past peak
                        $error++;
                        $DEBUG and print "Type $i: Rejecting $label from $v]. We can't go up ";
                        $DEBUG and print " after a peak. so error. Error now $error\n";
                     } elsif
                        # newscore is down but comparison is X. This means that there was already a peak
                        (($difscore== 0) && ($comparison[$i]=~ "X")) {
                        ## this is okay since it is down from a mode peak
                        $comparisonMap .= "X";
                        $DEBUG and print "Type $i: Adding $label to vertices $edge[0][1] because score ";
                        $DEBUG and print " is 0 and the comparison is X. ComparisonMap is now $comparisonMap\n";
                        
					    	} else {
             
								print "ERROR!!!! Not found match combination! MUST FIX! Some combination ";
                        print " is not being caught correctly... Exiting.\n";
                        print "Here is the score of the differences in %s for Type $i: $difscore\n";
                        print "Here is the comparison value: ", $comparison[$i],"\n";
								$error++;
                        exit();
					    	}
  						$DEBUG and print "Type $i:  Error so far $error\n";
                  }
                  if ($error==0){
                     $DEBUG and print "----------^^^^^^^^^^^----------------------------------\n";
           
                     $DEBUG and print "Original network: ", $network, "\n";

   					   ## no errors so add vertice added to the new network
                     my $newnet = $network->deep_copy_graph;
                     my $oldedgenum= $network->edges;
			   		   $newnet->add_weighted_edge( $label, $v, $comparisonMap);
                     $DEBUG and print "New network: ", $newnet, "\n";
					      my $e=$newnet->edges;
					      if ($e>$maxEdges) {
							   	$maxEdges=$e;
                           $DEBUG and print "Old max edges: $oldedgenum\n";
								   $DEBUG and print "New maximum edges: $e\n";
						   }
						   $match++;
						   ## copy this solution to the new array of networks
						   push @newnets, $newnet;
                     push @allNetworks, $newnet;
                     $DEBUG and print "----------^^^^^^^^^^^----------------------------------\n";
                  }
                 } else {
                     $DEBUG and print "$v has too many edges ( $edges is more than 1) so skipping\n";
                 } # end of if check for the end edges
	 				}	# end of iterate through the existing network link
				} # end of if assemblage not already in netowrk check
			} #end of assemblage loop

		} #end of network loop
	## no match at this point so no point in going forward.
	if ($match==0) {
		print "Maximum seriation size reached - no more assemblages added that iteration. Max # assemblages is: ", $maxEdges, "\n";
		$maxnumber = $currentMaxSeriationSize - $maxSeriations;
		## to break now...
		$currentMaxSeriationSize = $maxSeriations;
         
      ## copy the old network set to the new array (just in this case -- otherwise its empty
      foreach my $n (@networks) {
            my $ne=$n->deep_copy_graph;
            push @newnets, $ne;
	   }

	}
	print "Number of current solutions now: ", scalar(@newnets), "\n";
} #end of master loop through iterations

# now do some weeding. Basically start with the first network that is the largest, and work backwards. Ignore any
# network that is already represented in the smaller ones since these are trivial (e.g., A->B->C->D alreacy covers
# A->C->D.) That should then leave all the unique maximum solutions (or so it seems)

## first need to sort the networks by size
my @filteredarray=();

if ($filterflag ==1) {

print "Filtering solutions so we only end up with the unique ones.\n";
print "Start with ", scalar(@newnets), " solutions. \n";
my $count=0;
my $edgecount =$maxEdges;
while ($edgecount>1) {
	$DEBUG and print "### not on sets of $edgecount \n";
	foreach  my $network ( @newnets ) {
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
						$DEBUG and print "Not found in the network: $extracount\n";
					}
				}
				if ($extracount==$filteredNumber) {
					$DEBUG and print "solution doesnt exist in any of the existing networks - so add to list. \n";
					push @filteredarray, $network;
				}
			}
						
		}
	}
$edgecount--;
}
print "End with ", scalar(@filteredarray), " solutions.\n";
}  else {
   @filteredarray = @newnets;
}

if ($individualfileoutput) {
   my $writer=Graph::Writer::VCG->new();
   $count=0;
   my $name;
   foreach  my $network (@filteredarray ) {
   	my $E = $network->edges;
	   if ($E>$maxnumber-1) {
		   $count++;
		   $name= $count.'.vcg';
		   $writer->write_graph($network, $name);
	   }
	   #print Dumper($network);
   }

   my $writer2=Graph::Writer::Dot->new();
   $count=0;
   foreach  my $network ( @filteredarray ) {
	   my $V = $network->vertices;
	   if ($V==$maxEdges) {
		   $count++;
		   $name= $count.'.dot';
		   $writer2->write_graph($net, $name);
	   }
	   #print Dumper($network);
   }
}

############################################
#bootrap stuff


$numrows = scalar( @assemblages );
srand( $start );
my %perror = ();
my %pvalue = ();
my $results = 0;
my $loop = 100;
my $bigloop = 5;
my $loop2 = $loop;
my $ptr1 = 0;
my $ptr2 = 1;
my $classes=0;

# now do ALL the pairwise assemblage comparisons
# go to sleep and come back later.

if ($bootstrap) {
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

         #print "($debug) cycle value: $cycle\n";

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
	               $bootstrapdebug and print "loop $loop ",$ahat[ $c ] - $bhat[ $c ],"\t"; 
	               if 		( $ahat[ $c ] < $bhat[ $c ] ) { push @dir, -1; }
      	         elsif 	( $ahat[ $c ] > $bhat[ $c ] ) { push @dir, +1; }
	               else 							{ push @dir, 0;  }
            }
            $bootstrapdebug and  print "\n";


            # compare the two sets of results
            $num = scalar( @dir );
            $c = 0;
            my $diff = 0;
            for ( $c = 0; $c < $num; $c++ ) {
	            $bootstrapdebug and print "loop $loop ",$direction[ $c ],"/",$dir[ $c ],"\t"; 
	            if ( $direction[ $c ] == $dir[ $c ] ) { next; }
	            $diff++;
            }
            $bootstrapdebug and  print "\n";
            if ( $diff == 0 ) { $results++; }

            $loop--;
         }

         #print "Results:  $results matches of $loop2 trials\n";
         #print "Probability: ", $results / $loop2, "\n"; 

         $stat->add_data( $results/ $loop2 );   
         $cycle--;
   
         $results = 0;
      }
      my $label1= $labels[ $ptr1]."-".$labels[ $ptr2 ];
      my $label2= $labels[ $ptr2]."-".$labels[ $ptr1 ];
      $pvalue{ $label1 } = $stat->mean();
      $perror{ $label1 } =$stat->standard_deviation();
      $pvalue{ $label2 } = $stat->mean();
      $perror{ $label2 } =$stat->standard_deviation();
      $DEBUG and print $labels[ $ptr1], "\t", $labels[ $ptr2 ], "\t";
      $DEBUG and print $stat->mean(),"\t"; 
      $DEBUG and print $stat->standard_deviation(),"\n"; 

      undef $stat;
      $ptr2++;
   }
   $ptr1++;
   $ptr2 = $ptr1 + 1;

}
}



###########################################

print "Now printing output file... \n";

print OUTFILE "*Node data\n";
$count=0;
foreach my $l (@labels){
	print OUTFILE $l, "\n";
}
print OUTFILE "*Tie data\n";
print OUTFILE "From To Edge Weight Network pValue pError\n";
my $compareNetwork;
my @uniqueArray;
foreach my $compareNetwork (@allNetworks) {
   my $exists=0;
   foreach my $uarray (@uniqueArray) {
      if ($compareNetwork eq $uarray) { 
         $exists++;
       }
   }
   if (!$exists) {
      push @uniqueArray, $compareNetwork;
   }
}

## only print unique ones...
foreach  my $network ( @uniqueArray ) {
   #print Dumper($network);
   $count++;
   #print "$count: $network\n";
	my $E = $network->edges;
	if ($largestOnly) {
     if ($E==$maxEdges) {
	     my @Edges=$network->unique_edges;
		  foreach my $e (@Edges) {
                     if (!$bootstrap) { $perror{ @$e[0]."-".@$e[1] }=0.0; 
                           $pvalue{ @$e[0]."-".@$e[1] }=0.0; 
                        } 
                     
		     print OUTFILE @$e[0], " ", @$e[1], ", 1, ", scalar(@Edges), ", ", $count, ", ";
                     print OUTFILE $pvalue{ @$e[0]."-".@$e[1]}, ", ", $perror{ @$e[0]."-".@$e[1]},"\n";
			  #print @$e[0], " ", @$e[1], "\n";
		  } 
	  }
    } else {
		   my @Edges=$network->unique_edges;
		   foreach my $e (@Edges) {
		     if (!$bootstrap) { $perror{ @$e[0]."-".@$e[1] }=0.0; 
                        pvalue{ @$e[0]."-".@$e[1] }=0.0; 
		     }
		     print OUTFILE @$e[0], " ", @$e[1], ", 1, ", scalar(@Edges), ", ", $count, ", ";
                     print OUTFILE $pvalue{ @$e[0]."-".@$e[1]}, ", ", $perror{ @$e[0]."-".@$e[1]},"\n";
			   #print @$e[0], " ", @$e[1], "\n";
         }
    }
}

print OUTFILE "\n";



###########################################

if ($excel) {
   
   print "Now printing excel output file... \n";

}





printf ("Time for processing: %.2f seconds\n", Time::HiRes::gettimeofday() - $start);

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
