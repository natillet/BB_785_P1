#!/usr/bin/perl

use strict;
use warnings;

my $result;
my @results;
my $times = 0;
my $original_file_buffer;

open(FILE, "original_results.txt") or die $!;
$original_file_buffer = <FILE>;
close(FILE);

for (1 .. 5) {
	print "Iteration $_\n";
	my $new_file_buffer;
	my $cmd_out = `./project1`;
#	print "Result: $cmd_out";
	push(@results, $cmd_out);
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
	my $full_secs = "$sec\.$nsec";
	$times += $full_secs;
}

print "Average: " . ($times/5) . "\n";

exit 0;
