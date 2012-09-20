#!/usr/bin/perl

use strict;
use warnings;

my $result;
my @results;
my @times;
my $original_file_buffer;

open(FILE, "original_results.txt") or die $!;
$original_file_buffer = <FILE>;
close(FILE);

system("make clean");
system("make");

for (1 .. 5) {
	my $new_file_buffer;
	push(@results, `./project1`);
	open(FILE, "results.txt") or die $!;
	$new_file_buffer = <FILE>;
	close(FILE);
	if ($original_file_buffer ne $new_file_buffer) {
		die "Results do not match!";
	}
}

foreach my $time (@results) {
	#Execution Time: 7 sec, 804779251 nsec
	$time =~ /Execution Time: (\d+) sec, (\d+) nsec/;
	my $sec = $1;
	my $nsec = $2;
	print "$sec\.$nsec\n";
	
}

exit 0;