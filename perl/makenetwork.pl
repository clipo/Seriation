use strict;
use Data::Dumper;
use Graph;
use Graph::Writer::VCG;
use Graph::Writer::Dot;
use Getopt::Long qw(HelpMessage);
use Pod::Usage;

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


my ($debug, $inputfile, $man,$xyfile);
GetOptions(
    'debug'                     => \$debug,
    'help'                      => sub { HelpMessage() },
    'input=s'                   => \$inputfile,
    'xyfile=s'			=> \$xyfile,
    man                         => \$man
) or pod2usage(2);

## open the data file
my $useOutputFile = substr($inputfile,0,-4);
open( FILE, $inputfile ) or die "Cannot open $inputfile.\n";
open OUTFILE, ">$useOutputFile.vna";

my $count;
my $length;
my $net = Graph::Undirected->new;
my @arrayOfNodes;
my $numberOfNodes;
my %hashOfNodes;
my @xyAssemblages;
my %xAssemblage;
my %yAssemblage;

if ($xyfile) {
   ## open the xy file
   open( XYFILE, $xyfile ) or die "Cannot open XY File: $xyfile.\n";
   my $typecount;
   while (<XYFILE>) {
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
}
my $largestX = $xAssemblage{ largest_value_mem(%xAssemblage) };
my $largestY = $yAssemblage{ largest_value_mem(%yAssemblage)};

## the first part of the input file must contain the list of node names.
## the second part contains the ties
## the first line must contain the # of nodes 
my $lineskip;
my $line;

while (<FILE>)
{
	if ($lineskip<1) {
		$line =$_;
		#first line is number of nodes
		$numberOfNodes = $line;
		$lineskip=1;
	}
}
close FILE;

open( FILE, $inputfile ) or die "Cannot open $inputfile.\n";
my $lcount=0;
$count=0;
print OUTFILE "*Node Data\n";
print OUTFILE "ID X Y Easting Northing\n";
while (<FILE>) {
	## skip first line
	if ($lcount == 0){
		# do nothing
		#print "first line - nothing to do\n";
	} elsif ($count < $numberOfNodes) {
		$line= $_;
		chomp($line);
		$arrayOfNodes[ $count ] = $line;
		$hashOfNodes{ $line } = $count;
		my $x = $xAssemblage{ $line }/1000000 || 0;
		my $y = ($largestY-$yAssemblage{ $line })/100000 || 0;
		print OUTFILE $line," ", $x, " ", $y, " ", $xAssemblage{ $line }, " ", $yAssemblage{ $line }, "\n";		
		#print "Node: $count  Name: $line\n";
		#$net->add_vertex($count);
		$count++;
	}
	$lcount++;
}
$lineskip = $lineskip+$count;

my $linecount=0;
open( FILE, $inputfile ) or die "Cannot open $inputfile.\n";
my %hashOfEdges;
my %newHash;
while( <FILE> ) {
	if ($linecount < $lineskip) {
		#do nothing
	} else {
		# now evaluate the lines.
		my $line=$_;
		chomp($line);
		my ($node1, $node2, $edge) = split("\t", $line);
		#print "Node1: $node1 --> Node2: $node2 --> Edge: $edge \n";
		my $edgename=$node1."%".$node2;
		my $test = $net->is_reachable($hashOfNodes{ $node1}, $hashOfNodes{ $node2 });
		my $newvert;
		if (!$test) {  #first test if NOT connected -- always connect
			$newvert = $net->has_vertex( $hashOfNodes{ $node1 });
			if (!$newvert) { 
				$newHash  { $hashOfNodes{ $node1 } } = $edge;
			}
			$newvert = $net->has_vertex( $hashOfNodes{ $node2 });
			if (!$newvert) { 
				$newHash  { $hashOfNodes{ $node2 } } = $edge;
			}
			#print "Now adding edge:", $hashOfNodes{ $node1},"-", $hashOfNodes{ $node2},"-", $edge,"\n";
			$net->add_weighted_edge($hashOfNodes{ $node1}, $hashOfNodes{ $node2}, $edge);
			$hashOfEdges{ $edgename } = $edge;
		} elsif ( $newHash{ $hashOfNodes{ $node1 } } == $edge ) {
			$newvert = $net->add_weighted_edge($hashOfNodes{ $node1}, $hashOfNodes{ $node2}, $edge);
			$hashOfEdges{ $edgename } = $edge;
		} elsif ( $newHash{ $hashOfNodes{ $node2 } } == $edge ) {
			$newvert = $net->add_weighted_edge($hashOfNodes{ $node1}, $hashOfNodes{ $node2}, $edge);	
			$hashOfEdges{ $edgename } = $edge;
		} elsif ($newHash{ $node1 } > $edge) { 
			delete( $newHash{ $node1});
		} elsif ($newHash{ $node2 } > $edge) {
			delete( $newHash{ $node2} );
		}
	}
	$linecount++;

}

print OUTFILE "*tie data\n";
print OUTFILE "from to weight\n";
my $edgeValue;
my $edgeName;
while ( ( $edgeName, $edgeValue) = each %hashOfEdges ) {
	my ($node1, $node2) = split("%", $edgeName);
	print OUTFILE $node1, "\t", $node2, "\t", $edgeValue, "\n";
}

close OUTFILE;
my $writer=Graph::Writer::VCG->new();
$writer->write_graph($net, 'output.vcg');
my $writer2=Graph::Writer::Dot->new();
$writer2->write_graph($net, 'output.dot');



__END__

    =head1 makenetwork.pl

    Create file for pruning network to the minimum set

    =head1 SYNOPSIS

    matchandsort.pl [options] 

     Options:
       -help                brief help message
       -man                 full documentation
       -input=<filename>    filename of data to seriate
	-xyfile=<filename>  filename of the xy data
	
    =head1 OPTIONS

    =over 8

    =item B<-help>

    Print a brief help message and exits.

    =item B<-man>

    Prints the manual page and exits.

    =back

    =head1 DESCRIPTION

    B<This program> reads given input file(s) and create networks.

    =cut
