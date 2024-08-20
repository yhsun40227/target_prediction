#!/usr/bin/perl
use strict;
use Getopt::Std;

# 3rd party program directory
my $RNAHYBRID_IR = "/Users/yhsun/RNAHybrid"; # directory for RNAHyrbid
my $SCAN_FOR_MATCHES_DIR = "/Applications/ScanForMatches"; # directory for scan_for_matches


# command options
our ($opt_i,$opt_d,$opt_s,$opt_p,$opt_f,$opt_o,$opt_r);
getopt ('idspfor');

my $miRNAfile = ($opt_i||"miRNA.fa");                 # -i miRNA sequence fasta file; default value is "miRNA.fa"
my $seqsource = ($opt_d||"db.fa");                    # -d target sequences fasta file; default value is "db.fa"
my $score_test = ($opt_s||3.5);                       # -s penalty score cutoff; default value is 3.5
my $pd_score_test = ($opt_p||4.0);                    # -p position weighted penalty score cut off; default value is 4.0
my $mfe_test = ($opt_f||0.73);                        # -f MFE ratio cut off; default value is 0.73
my $predition_result = ($opt_o||"prediction_output"); # -o prediction output; default value is "prediction_output"
my $predition_dir = ($opt_r||"html");                 # -r directory for predcition results in html format; default value is "html"

print "miRNA sequence file is $miRNAfile\n";
print "Target sequence file is $seqsource\n";
print "Penalty score cutoff is $score_test\n";
print "Position weighted penalty score cut off is $pd_score_test\n";
print "MFE ratio cut off is $mfe_test\n";
print "Output file is $predition_result\n";
print "HTML file directory is $predition_dir\n";

mkdir "$predition_dir";
mkdir "TMP";

open OUT,">$predition_result"||die;


open IN,"$miRNAfile"||die;
$/=">";
while (<IN>){
	chomp;
	my ($miRNAid, $miRNAseq) = split /\n/,$_;
	my @MIRNAID = split / /,$miRNAid;
	$miRNAid=shift @MIRNAID;
	unless ($miRNAid eq ""){
		print "miRNAid = $miRNAid\nmiRNA sequence = $miRNAseq\n";
		$miRNAid=~s/ +//g;
		$miRNAseq=~s/ +//g;
		&run_scan_for_matches ($miRNAid, $miRNAseq, $seqsource);
	}
}
$/="\n";
close IN;
close OUT;
system "rm -rf ./TMP";



sub run_scan_for_matches{

	my ($miRNAid, $miRNAseq, $seqsource)=@_;
	my %P_result;

	$miRNAseq =~ s/U/T/g;
	my $rcmiRNAseq = reverse $miRNAseq;
	$rcmiRNAseq =~ tr/[A,T,C,G]/[T,A,G,C]/;
	my %redundantcheck;

	my @PAT = ("700","310","301");

	foreach my $pat (@PAT){

		my @splitpat = split //, $pat;
		my $mirnaseq = $rcmiRNAseq;
		open (PAT,">./TMP/temp.pat");
		print PAT "$rcmiRNAseq\[$splitpat[0]\,$splitpat[1]\,$splitpat[2]\]\n";
		close (PAT);
		print "scanning and evaluating pattern $miRNAid\ pattern \[$splitpat[0]\,$splitpat[1]\,$splitpat[2]\]\n";
		system "$SCAN_FOR_MATCHES_DIR/scan_for_matches ./TMP/temp.pat < $seqsource > ./TMP/$miRNAid\_$pat\.match";
			
		open INmatch,"./TMP/$miRNAid\_$pat\.match"||die;
		while (<INmatch>){
				
			chomp;
			
			my ($header, $targetseq) = split /\n/,$_;
			
			unless ($header eq ""){
				$targetseq =~ s/ +//g;
				my ($targetid,$hitinfo) = split /\:/, $header;
				$hitinfo =~ s/\[|\]//g;
				my ($startloc, $endloc) = split /\,/, $hitinfo;
				
				unless ($redundantcheck{$targetid}{$startloc}){
					
					$redundantcheck{$targetid}{$startloc}=1;
					$redundantcheck{$targetid}{$startloc-1}=1;
					$redundantcheck{$targetid}{$startloc+1}=1;
					
					my ($alntargetseq, $alnpat, $alnmirseq, $score, $pw_score) = &penalty_score($rcmiRNAseq,$targetseq);
						
					$alntargetseq =~ s/T/U/g;
					$alnmirseq =~ s/T/U/g;

					if ((($score <= $score_test)&&($pw_score <= $pd_score_test))){
					#if ((($score <=3.5)||($pw_score<=4))){
						my $MFEratio = &mferatio($miRNAseq, $targetseq);
						if ($MFEratio >= $mfe_test){
							
							my $P_record = join "\t",($miRNAid,$targetid,$startloc,$endloc,$alntargetseq,$alnpat,$alnmirseq,$score,$pw_score,$MFEratio);
						
							$P_result{$miRNAid}{$pw_score}{$score}{targetid}{$startloc}=$P_record;
								
							print OUT "$miRNAid\t$targetid\t$alntargetseq\t$alnpat\t$alnmirseq\t$score\t$pw_score\t$MFEratio\n";
						}
					}
				}
			}
		}
		close INmatch;
		unlink ("./TMP/$miRNAid\_$pat\.match");
	}
	&print_result(%P_result);
}


sub penalty_score{

		my ($rcmiRNAseq,$targetseq)=@_;
		my ($score, $pw_score, $alignpat, $i, $finalmirnaseq, $finaltargetseq);
		
		if (length $rcmiRNAseq == length $targetseq){
			($score, $pw_score, $alignpat) = &compareseq($rcmiRNAseq, $targetseq);
			$finalmirnaseq=$rcmiRNAseq;
			$finaltargetseq=$targetseq;
		} else {
			my $testtmpscore=40;
			my $longseq;
			my $shortseq;
			if (length $targetseq > length $rcmiRNAseq){
				for ( $i=0; $i <= length($rcmiRNAseq); $i++) {
					my @seq = split //,$rcmiRNAseq;
					splice (@seq, $i, 0, ("-"));
					my $newrcmiRNAseq = join "",@seq;
					my ($tmpscore, $tmppwscore, $tmpalignpat) = &compareseq($newrcmiRNAseq, $targetseq);
					
					if ($tmpscore < $testtmpscore) {
						$testtmpscore = $tmpscore;
						($score, $pw_score, $alignpat, $finalmirnaseq, $finaltargetseq) = ($tmpscore, $tmppwscore, $tmpalignpat, $newrcmiRNAseq, $targetseq);	
					}
				}
			} else {
				for ( $i=0; $i <= length($targetseq); $i++) {
					my @seq = split //,$targetseq;
					splice (@seq, $i, 0, ("-"));
					my $newtargetseq = join "",@seq;
					my ($tmpscore, $tmppwscore, $tmpalignpat) = &compareseq($rcmiRNAseq, $newtargetseq);

					if ($tmpscore < $testtmpscore) {
						$testtmpscore = $tmpscore;
						($score, $pw_score, $alignpat, $finalmirnaseq, $finaltargetseq) = ($tmpscore, $tmppwscore, $tmpalignpat, $rcmiRNAseq, $newtargetseq);	
					}
				}
			}
		}
		
		$finalmirnaseq =~ tr/[A,T,C,G]/[T,A,G,C]/;
		return ($finaltargetseq, $alignpat, $finalmirnaseq, $score, $pw_score);
}


sub compareseq {

	my ($mseq, $tseq) = @_;
	
	my $ltest;
	if (length $mseq >= 20){$ltest = 20} else {$ltest = length $mseq}

	my @MSEQ = split //,$mseq;
	my @TSEQ = split //,$tseq;
	my $startpoint;
	
	if ((scalar @MSEQ) <= $ltest) {$startpoint = 1} else {$startpoint = (scalar @MSEQ) - $ltest + 1}
	
	my $initscore = 100;
	my $comseqscore=0; 
	my $pw_comseqscore=0; 
	my $comseqalignpat="";
	
	for (my $n = 0; $n < $startpoint; $n++) {
		my $initpat = ""; 
		my $pw_tmpcomseqscore=0; 
		my $tmpcomseqscore=0; 
		my $tmpcomseqalignpat="";
		my @pattern;
		my $match=0;
		my $wobble=0;
		my $mismatch=0;
		my $buldge=0;
		my $pw_wobble=0;
		my $pw_mismatch=0;
		my $pw_buldge=0;
		for (my $k = $n; $k < $ltest+$n ; $k++) {
			my $mirloc = length ($mseq)-$k;
		
			if ($MSEQ[$k] eq $TSEQ[$k]){
				$pattern[$k] = "|";
				$match++;
			} elsif (($MSEQ[$k] eq "A")&&($TSEQ[$k] eq "G")||($MSEQ[$k] eq "C")&&($TSEQ[$k] eq "T")) {
				$pattern[$k] = "o";
				if (($mirloc <= 13)&&($mirloc >= 2)){$pw_wobble++} else {$wobble++};
			} elsif (($MSEQ[$k] eq "-")||($TSEQ[$k] eq "-")) {
				$pattern[$k] = "-";
				if (($mirloc <= 13)&&($mirloc >= 2)){$pw_buldge++} else {$buldge++};
			} else {
				$pattern[$k] = "x";
				if (($mirloc <= 13)&&($mirloc >= 2)){$pw_mismatch++} else {$mismatch++};
			}
		}
		$initpat = join "", @pattern;
		my $tailspace = $startpoint-$n-1;
		 $tmpcomseqalignpat = "~"x$n.$initpat."~"x$tailspace;
		 $tmpcomseqscore = (0.5*($wobble+$pw_wobble))+(1*($mismatch+$pw_mismatch))+(2*($buldge+$pw_buldge));
		 $pw_tmpcomseqscore = ((0.5*$wobble)+(1*$mismatch)+(2*$buldge))+(2*((0.5*$pw_wobble)+(1*$pw_mismatch)+(2*$pw_buldge)));
		
		if ($tmpcomseqscore < $initscore) {
			$initscore = $tmpcomseqscore;
			
			($comseqscore, $pw_comseqscore, $comseqalignpat) = ($tmpcomseqscore, $pw_tmpcomseqscore, $tmpcomseqalignpat);
		}
	}
	return ($comseqscore, $pw_comseqscore, $comseqalignpat);
}

sub mferatio {

	my ($miRNAseq, $targetseq) = @_;
	my $dgmfe;
	my $dgtarget;
	
	$targetseq =~ s/T/U/g;
	$miRNAseq =~ s/T/U/g;
	my $rctseq = reverse $targetseq;
	$rctseq =~ tr /[A,U,C,G]/[U,A,G,C]/;

	open TSEQ,">./TMP/T.seq"||die;
	print TSEQ ">T\n",$targetseq,"\n";
	close TSEQ;

	open TrcSEQ,">./TMP/Trc.seq"||die;
	print TrcSEQ ">Trc\n",$rctseq,"\n";
	close TrcSEQ;

	open MSEQ,">./TMP/M.seq"||die;
	print MSEQ ">M\n",$miRNAseq,"\n";
	close MSEQ;
	

	system "$RNAHYBRID_IR/RNAHybrid -t ./TMP/T.seq -q ./TMP/M.seq -d theta > ./TMP/T-M-theta.dG";
	system "$RNAHYBRID_IR/RNAHybrid -t ./TMP/T.seq -q ./TMP/Trc.seq -d theta > ./TMP/T-Trc-theta.dG";
	
	open TM,"./TMP/T-M-theta.dG"||die;
	while (<TM>){
		chomp;
		my @LINE = split /\n/,$_;
		foreach my $line (@LINE){
			if ($line =~ /^mfe/){
				my @bb = split / /,$line;
				$dgtarget = $bb[1];
			}
		}
	}
	close TM;
	unlink "./TMP/T-M-theta.dG";


	open TT,"./TMP/T-Trc-theta.dG"||die;
	while (<TT>){
		chomp;
		my @LINE = split /\n/,$_;
		foreach my $line (@LINE){
			if ($line =~ /^mfe/){
				my @bb = split / /,$line;
				$dgmfe = $bb[1];
			}
		}
	}
	close TT;
	unlink "./TMP/T-Trc-theta.dG";

	unlink "./TMP/T.seq";
	unlink "./TMP/Trc.seq";
	unlink "./TMP/M.seq";
	
	
	my $MFEratio = $dgtarget/$dgmfe;
	
	return ($MFEratio);
	
}

sub print_result {
	my (%RESULT) = (@_);

	foreach my $miRNAID (keys %RESULT){
	
	open OUThtml,">./$predition_dir/$miRNAID\.html"||die;
	print OUThtml "<html>\n<head>\n<title>$miRNAID target prediction result</title>\n</head>\n";
	print OUThtml  "<table border=\"5\" align=\"Center\">\n";
	print OUThtml  "<tr>\n";
	print OUThtml  "<td align\=\"center\"><font face=\"Courier\">Target ID<br><br>miRNA ID</td>\n";
	print OUThtml  "<td align\=\"center\"><font face=\"Courier\">Target region:<br>start - end<br></td>\n";
	print OUThtml  "<td align\=\"center\"><font face=\"Courier\">5\'\-target sequence\-3\'<br>  complementary pattern  <br>3\'\-miRNA sequence\-5\'</font></td>\n";
	print OUThtml  "<td align\=\"center\"><font face=\"Courier\">score</td>\n";
	print OUThtml  "<td align\=\"center\"><font face=\"Courier\">pd score</td>\n";
	print OUThtml  "<td align\=\"center\"><font face=\"Courier\">MFE ratio</td>\n";
	print OUThtml  "</tr>\n";

	foreach my $PWSCORE (sort {$a <=> $b} keys %{$RESULT{$miRNAID}}){
		foreach my $SCORE (sort {$a <=> $b} keys %{$RESULT{$miRNAID}{$PWSCORE}}){
			foreach my $ID (sort {$a cmp $b} keys %{$RESULT{$miRNAID}{$PWSCORE}{$SCORE}}){
				foreach my $SLOC (sort {$a <=> $b} keys %{$RESULT{$miRNAID}{$PWSCORE}{$SCORE}{$ID}}){
							
					my ($miRNAid,$targetid,$startloc,$endloc,$alntargetseq,$alnpat,$alnmirseq,$score,$pw_score,$MFEratio) = split /\t/, $RESULT{$miRNAID}{$PWSCORE}{$SCORE}{$ID}{$SLOC};
	
					print OUThtml  "<tr>\n";
					print OUThtml  "<td align\=\"center\"><font face=\"Courier\">$targetid<br><br>$miRNAid</td>\n";
					print OUThtml  "<td align\=\"center\"><font face=\"Courier\">$startloc\-$endloc<br><br></td>\n";
					print OUThtml  "<td align\=\"left\"><font face=\"Courier\">5\'\-$alntargetseq\-3\'<br>&nbsp;&nbsp;&nbsp;$alnpat<br>3\'\-$alnmirseq\-5\'</font></td>\n";
					print OUThtml  "<td align\=\"center\"><font face=\"Courier\">$score</td>\n";
					print OUThtml  "<td align\=\"center\"><font face=\"Courier\">$pw_score</td>\n";
					print OUThtml  "<td align\=\"center\"><font face=\"Courier\">";
					printf OUThtml "%.2f", $MFEratio;
					print OUThtml "</td>\n";
					print OUThtml  "</tr>\n";
			
				}
			}
		}
	}
	print OUThtml "</table>\n";
	print OUThtml "</html>\n";
	close OUThtml;
	}
}


