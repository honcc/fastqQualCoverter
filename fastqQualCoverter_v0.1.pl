#!/usr/bin/perl -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle);
use threads;
use Statistics::Descriptive;
use Cwd 'abs_path';

######################################################################################################################################################
#
#	Description
#		This is a perl script to convert the quality score of a fastq (in GZ format) from 33 to 64 scale (and vice versa).
#
#	Input
#		--fastqPath=			file path; [compulsory]; path of the fastq file in .gz of .fastq format;
#		--outScale=				score scale; [64]; score scale in 64 or 33;
#		--outDir=				dir path; [./fastqQualCoverter/]; output directory;
#
#	Output
#
#	Usage
#		perl fastqQualCoverter_v0.1.pl --fastqPath=/Volumes/A_MPro2TB/NGS/fastq/tmpTestPolyALen/1301_herbSeq_DMSO_T0H_2.clean.fastq.gz --outScale=100
#
#	Assumption
#
#	Version history
#
#		v0.1
#			-debut;
#
#####################################################################################################################################################

#==========================================================Main body starts==========================================================================#
printCMDLogOrFinishMessage('CMDLog');

my ($fastqPath, $outScale, $outDir) = readParameters();
my ($qualScale) = checkQualityScale($fastqPath);

if ($qualScale == $outScale) {
	my ($fastqName,$fastqDir,$fastqSuffix) = fileparse($fastqPath, qr/\.[^.]*/);
	my $newPath = "$outDir/$fastqName.qual$outScale.fastq.gz";
	print "QualScale of the fastq is already $outScale.\nCopying the original to outDir.\n";
	system "cp $fastqPath $newPath;";
	die "quitting.\n";
}

convertFastqOnTheFly($fastqPath, $outScale, $outDir);
printCMDLogOrFinishMessage('finishMessage');

########################################################################## readParameters
sub readParameters {
	
	my ($fastqPath, $outScale, $outDir);

	$fastqPath = undef;
	$outScale = 64;
	$outDir= "./fastqQualCoverter";

	&GetOptions(	
		'fastqPath=s' => \$fastqPath, 
		'outScale=i' => \$outScale,
		'outDir:s' => \$outDir,
	);

	system ("mkdir -p -m 777 $outDir");
	
	return ($fastqPath, $outScale, $outDir);
}
########################################################################## printEditedFileOnTheFly
sub convertFastqOnTheFly {
	
	my ($fastqPath, $outScale, $outDir) = @_;
	
	my $convertVal = 31;
	$convertVal = -31 if $outScale == 33;
	
	#----choose the zipper
	my $zipper = 'pigz';
	$zipper = 'gzip' if (`pigz --version 2>&1` =~ m/command not found/);

	my ($fastqName,$fastqDir,$fastqSuffix) = fileparse($fastqPath, qr/\.[^.]*/);

	if ($fastqSuffix eq '.gz') {
		open (FASTQ, "$zipper -d -c $fastqPath |");
	} elsif ($fastqSuffix eq '.fastq') {
		open (FASTQ, "<", "$fastqPath");
	} else {
		die "quitting: $fastqPath is in $fastqSuffix. It has to be in .fastq or .gz format\n";
	}
	
	open (OUTFASTQ, "| $zipper -c >$outDir/$fastqName.qual$outScale.fastq.gz");
	my $readProc = 0;
	
	while (chomp(my $theLine = <FASTQ>)) {
		
		$readProc++;
		if (not $readProc % 10000) {
			print "$readProc\treads processed                         \r";
		}
		
		chomp (my $seqHeader = $theLine);
		chomp (my $seq = <FASTQ>);
		chomp (my $qualHeader = <FASTQ>);
		chomp (my $qual = <FASTQ>);

		my @ASCII = unpack("C*", $qual);
		my @newASCII = ();
		foreach my $val (@ASCII) {
			push(@newASCII, $val+$convertVal);
		}
		my $convertedQual = pack("C*", @newASCII);
		
		print OUTFASTQ $seqHeader."\n";
		print OUTFASTQ $seq."\n";
		print OUTFASTQ "+\n";
		print OUTFASTQ $convertedQual."\n";
		
		last if (eof FASTQ);
	}
	
	print "\n";

	close OUTFASTQ;
	close FASTQ;

}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLog {

	#---open a log file if it doesnt exists
	my $scriptNameXext = $0;
	$scriptNameXext =~ s/\.\w+$//;
	open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
	close CMDLOG;
	
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {

		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## checkQualityScale
sub checkQualityScale {

	my ($fastqPath) = @_;

	print "Checking quality scale of $fastqPath\n";

	my ($fastqName,$fastqDir,$fastqSuffix) = fileparse($fastqPath, qr/\.[^.]*/);
	print "Checking format of $fastqName$fastqSuffix\n";

	if ($fastqSuffix eq '.gz') {
		open (FASTQ, "gzip -d -c $fastqPath |");
	} elsif ($fastqSuffix eq '.fastq') {
		open (FASTQ, "<", "$fastqPath");
	} else {
		die "quitting: $fastqPath is in $fastqSuffix. It has to be in .fastq or .gz format\n";
	}

	my $lineCount = 3;
	my @qualCharAry = ();

	my $qualScale = undef;
	my $linesToSample = 1000;
	while (not defined $qualScale) {

		while (my $curntLine = <FASTQ>) {
			chomp $curntLine;
			$lineCount++;
			if ($lineCount % 4 == 0) {#---multiple of 4
				die "fastq have to be in 4 lines per read foramt\n" if $curntLine !~ m/^@/;
				for my $i (1..3) {
					chomp (my $nextLine = <FASTQ>);
					$lineCount++;
					
					push @qualCharAry, split //, $nextLine if $i == 3; 
				}
			}
			last if $lineCount > $linesToSample;
		}
		
		my %qScoreRngHsh = ();
		foreach my $scale64Or33 ((64, 33)) {#---just record the max and min
			${$qScoreRngHsh{$scale64Or33}}{"max"} = -99999;
			${$qScoreRngHsh{$scale64Or33}}{"min"} = 99999;
		}
		
		foreach my $qualChar (@qualCharAry) {
			my $ASCII = ord($qualChar);
			foreach my $scale64Or33 ((64, 33)) {#---just record the max and min
				my $score = $ASCII - $scale64Or33;
				${$qScoreRngHsh{$scale64Or33}}{"max"} = $score if $score > ${$qScoreRngHsh{$scale64Or33}}{"max"};
				${$qScoreRngHsh{$scale64Or33}}{"min"} = $score if $score < ${$qScoreRngHsh{$scale64Or33}}{"min"};
			}
		}
		
		my %validScaleHsh = ();
		foreach my $scale64Or33 ((64, 33)) {#---just record the max and min
			my $max = ${$qScoreRngHsh{$scale64Or33}}{"max"};
			my $min = ${$qScoreRngHsh{$scale64Or33}}{"min"};
			print "$linesToSample reads sampled: Scale $scale64Or33: max=$max min=$min        \r";
			if ((($min>=0) and ($min<=10)) and (($max>=30) and ($max<=41))) {
				print "$linesToSample reads sample in scale $scale64Or33: max=$max min=$min        \n";
				$validScaleHsh{$scale64Or33}++;
				$qualScale = $scale64Or33;
			}
		}
		
		print "WARNING: more than 1 quality scale fits\n" if keys %validScaleHsh > 1;
		$linesToSample *= 2;
		
		die "Cannot determine quality scale after sampling 10 million reads\n" if (($linesToSample > 10000000) and (not defined $qualScale));
	}
	close FASTQ;
	print "Scale $qualScale is detected                                 \n\n";
	return $qualScale;
}
