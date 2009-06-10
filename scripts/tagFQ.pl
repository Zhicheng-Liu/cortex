#!/usr/local/bin/perl5.6.0 -w
use strict;

my ($i, $rp, $count, $name, $oname
	);
 
$name = "";
$oname = "";
my (%seq) = ();
my ($mindiff) = 0;
foreach $name (@ARGV) {
	if($name =~ /(\S+)\.gz$/ or $name =~ /(\S+)\.Z$/) {
		open F,"gunzip -c $name |";
		$name = $1;
	} else {
		open F,$name;
	}
	my ($tag) = "";
	my ($tail) = "";
	my ($length) = 0;
	my ($pos) = -1;
	my ($tellF) = 0.;
	while (<F>) {
		$tellF += length($_);
		chomp;
		my ($line) = $_;
		if ($line =~ /^\@(\S+)\s*(.*)/) {
			if($pos >= 0) {
				printf "$tag %.0f $length $name $tail\n",$pos;
				$length = 0;
			}
			$tag = $1;
			$tail = "";
			$tail = $2 if defined $2;
			$pos = $tellF - length($line) - 1;
			next;
		}
		if ($line =~ /^[acgtACGTnN]/) {
			$length += length($line);
			next;
		}
		if ($line =~ /^\+$tag/) {
			$line = <F>;
			$tellF += length($line);
			next;
		}
	}
	close F;
	printf "$tag %.0f $length $name $tail\n",$pos;
}
