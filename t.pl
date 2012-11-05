use Graph;
use Data::Dumper;

my $net = Graph::Undirected->new;
		$net->add_vertex("A");
		$net->add_vertex("B");
		$net->add_vertex("C");
      $net->add_vertex("D");
      $net->add_vertex("E");
		$net->add_weighted_edge( "A", "B", 1);
$net->add_weighted_edge( "B", "C", 2);
$net->add_weighted_edge( "C", "D", 3);
$net->add_weighted_edge( "D", "E", 4);

print "network :", $net, "\n";

my @v = $net->vertices;

foreach my $v (@v) {
      print "vert: ", $v, "\n";
      @n = $net->neighbours($v);
      #my @edges = $net->edges_at($v);
      my $step ==0;
      my %visited={};
      if (scalar(@n)==1) {
         while($end==0) {
         $step++;
         foreach my $vv (@n) {
            if (!%visited{ $vv } ) {
               $visited{ $vv } ==1;
               print "$step: ", $vv, "\n"; #Dumper($vv), "\n";
            } 
            $next = $net->neighbours($n[0]);
         }
      }
}
