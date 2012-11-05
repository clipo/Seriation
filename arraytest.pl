use Data::Dumper;
 use Array::Utils qw(:all);

my @array1 = ('A', 'B', 'C', 'D');
my @array2 = ('A', 'B', 'C', 'D'); #, 'E', 'F');
my @array3 = ('A', 'B', 'C',  'D', 'E', 'F');

print Dumper(@array1), "\n";
print Dumper(@array2), "\n";
print Dumper(@array3), "\n";
my @arraydiff = array_diff( @array1, @array2);
print "Array1-Array 2", scalar(@arraydiff), "\n";
my $ad1=0;
$ad1=scalar(@arraydiff);
print "Array1-Array 2", $ad1, "\n";
print Dumper(@arraydiff), "\n";
print "--------------\n";

my @arraydiff2 = array_diff( @array3, @array1);
my $ad2=0;
$ad2=scalar(@arraydiff2);
print "Array3-Array 1:  ", $ad2, "\n";
print Dumper(@arraydiff2), "\n";

print "Difference 1- 2\n";
my %diff1;
@diff1{ @array1 } = @array1;
delete @diff1{ @array2 };

my @k = (keys %diff1);
print Dumper @k,"\n";
print "difference: ", scalar(@k), "\n";

print "Difference 1 -3 \n";
my %diff2;
@diff2{ @array1 } = @array1;
delete @diff2{ @array3 };

my @m = (keys %diff2);
print Dumper @m,"\n";
print "difference: ", scalar(@m), "\n";
# %diff1 contains elements from '@a1' that are not in '@a2'
