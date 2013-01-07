$latex = "xelatex  --synctex=1 -file-line-error -output-driver=\"xdvipdfmx -q -E\" %O %S";
#$pdflatex = 'pdflatex %O %S';
$pdflatex = "xelatex --synctex=1 -file-line-error -output-driver=\"xdvipdfmx -q -E\" %O %S";

$pdf_previewer = "open -a /Applications/Skim.app"; 
$clean_ext = "paux lox pdfsync out glg gls glo brf ist synctex.gz bbl thm spl tex";

add_cus_dep('glo', 'gls', 0, 'makeglossaries');

sub makeglossaries {
	system "makeglossaries $_[0]";
	if ( -z "$_[0].glo" ) {
		print "Latexmk: Empty glo file, I am making dummy gls file\n";
		open GLS, ">$_[0].gls";
		close GLS;
	}
	return 0;
}

