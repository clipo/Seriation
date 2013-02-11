use strict;
use Spreadsheet::WriteExcel;
use Data::Dumper;
use Getopt::Long qw(HelpMessage);
use Pod::Usage;

my ($debug,$inputfile,$man,$separator,$halfmatchflag,$dissimilarity, $xyfile);

# Are the values separated by (0) -- tabs and frequency data (1) nothing, (2) tabs, or (3) spaces (1,2,3) : ";
$separator = 1;	   #Default is 1 ( tabs and frequency data )
# "Do you want to use (1) 0.5 matching for ? or (0) strict matching where ? are ignored (1, 0) : ";
$halfmatchflag=0;   # Default is 0
# "Do you want to calculate similarity (0) or dissimilarity (1) :";
$dissimilarity=0;   # Default is 0

GetOptions(
    'debug'                     => \$debug,
    'help'                      => sub { HelpMessage() },
    'input=s'                   => \$inputfile,
    'separator=s'		=> \$separator,
    'probabilitymatching'	=> \$halfmatchflag,
    'dissimilarity'		=> \$dissimilarity,
    'xyfile=s'			=> \$xyfile,
    man                         => \$man
) or pod2usage(2);

my $DEBUG = $debug;    # our "$debug level"

if ($DEBUG) {
    print "Verbose debugging output is on!!!\n";
    print "Processing input file: $inputfile\n";
}

## open the data file
my $useOutputFile = substr($inputfile,0,-4);
open( FILE, $inputfile ) or die "Cannot open $inputfile.\n";
open OUTFILE, ">$useOutputFile-occur.out";
my $workbook = Spreadsheet::WriteExcel->new("$useOutputFile-occur.xls");
my $worksheet=$workbook->addworksheet();

my @array;
my $count;
my $length;
my @namearray;
my %hashOfNodes;

while (<FILE>)
{
	my $line =$_;
	chomp($line);
	my ($name, $val);
	$count++;
	if ($separator==0) {
		my @vals= split("\t", $line);
		my $numClasses= scalar(@vals)-1;
		push @namearray,$vals[0];
		#print $vals[0]," - ",$vals[1], " - ", $vals[2], "\n";
		$worksheet->write($count,0, $vals[0]);
		$worksheet->write(0,$count,$vals[0]);
		my $descriptor;
		for (my $i=1; $i<$numClasses; $i++) {
			#print vals[$i],"\n";
			if ($vals[$i]>0) {
				$descriptor .= "1";
			} else {
				$descriptor .= "0";
			}
		}	
		push(@array, $descriptor);
	}
	if ($separator==1) {
		($name, $val) = split("\t", $line);
		push(@namearray, $name);	
		push(@array, $val);
		$length = scalar($val);
		$worksheet->write($count,0, $name);
		$worksheet->write(0,$count,$name);
	}
	if ($separator == 2 || $separator ==3) {
		my @valarray;
		if ($separator==2) {	
			@valarray = split("\t", $line);
		} else {
			@valarray = split("\s", $line);
		}
		push(@namearray, $valarray[0]);
		#print Dumper(@valarray);
		my @newarray;
		for (my $i=1; $i<scalar(@valarray)-1; $i++) {
			$newarray[$i]=$valarray[$i];
		}	

		push(@array, [ @newarray ]);
		$length=scalar(@newarray);
		$worksheet->write($count,0, $valarray[0]);
		$worksheet->write(0,$count,$valarray[0]);
	}
	chomp($name);
}
close FILE;

#print Dumper(\@array);
#outer loop
for (my $m=0; $m<$count; $m++)
{
	my $testname = $namearray[$m];
	for (my $n=$m; $n<$count; $n++)
	{
		my $comparename = $namearray[$n];
#		print "Comparing :-->$testinstance-->with->$compareinstance\n";
		my $score=0;
		if ($separator == 1 || $separator==0) { ## case where nosparator
			my $testinstance = $array[$m];
			my $compareinstance = $array[$n];
			for (my $o=0; $o<length($testinstance); $o++)
			{
				my $test = substr($testinstance, $o, 1);
				my $compare =substr($compareinstance, $o, 1);
				$score=$score + comparevalues($test, $compare); 
			}
		}
		if ($separator == 2    ) { ## case with tab separator
			my @testinstance = @{$array[$m]} ;	
			my @compareinstance = @{$array[$n]};		
		#	print "Comparing :-->Dumper(@testinstance)-->with->Dumper(@compareinstance)\n";
			for (my $o=0; $o< $length; $o++)
			{
				my $test = $testinstance[$o];
				my $compare = $compareinstance[$o];
				#print "compare $test with $compare\n";
				$score= $score + comparevalues($test,$compare);
			}
		
		}		
		if ($dissimilarity) { $score = ($length)-$score-1; }
		$hashOfNodes {	"$testname#$comparename" } = $score;
		$worksheet->write($n+1,$m+1,$score);
		$score = 0;
	}

}

print OUTFILE $count,"\n";
foreach my $value ( @namearray ) {
   print OUTFILE $value,"\n";
}

foreach my $value (sort {$hashOfNodes{$b} <=> $hashOfNodes{$a} } 
           keys %hashOfNodes)
{
	my ($node1, $node2) = split("\#", $value);
	if (!($node1 eq $node2)) {
		print OUTFILE "$node1\t$node2\t$hashOfNodes{ $value }\n"; 
	}
}

close OUTFILE;
if ($xyfile){
	system("perl makenetwork.pl -input=$useOutputFile-occur.out -xyfile=$xyfile");
} else {
	system("perl makenetwork.pl -input=$useOutputFile-occur.out ");
}
print "Done!\n";
exit();

sub comparevalues {
	my ($ctest, $ccompare) =  @_;
	my $tempscore = 0;	
	#print "comparing $ctest with $ccompare\n";
	if ($halfmatchflag ==1 && ($ctest eq "?" || $ccompare eq "?")) {
		$tempscore = .5;
	} elsif ($ctest == $ccompare)
	{
		$tempscore = 1;
	}
return $tempscore;
}


__END__

    =head1 matchandsort.pl

    Create file for generating occurrence network

    =head1 SYNOPSIS

    matchandsort.pl [options] 

     Options:
       -help                brief help message
       -man                 full documentation
       -input=<filename>    filename of data to seriate
       -dissimilarity	    use dissimilarity for comparison
       -probabilitymatching	use 0.5 as chance for matching missing data (?)
       -separator=<0, 1,2,3>	0=frequency data, 1= no separator, 2=tabs, 3=space
       -xyfile=<filename>	filename with XY data

    Example: perl matchandsort.pl -input=../testdata/pfg-cpl.txt -xyfile=../testdata/pfgXY.txt -separator=0
    
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
