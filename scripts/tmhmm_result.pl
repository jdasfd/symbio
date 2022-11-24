#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

GetOptions(
    "i|in=s"        => \(my $input),
    "h|help"        => \(my $help),
);

die usage() if defined $help;
die ("Input a file please.") if not defined $input;

sub usage{
    my $help_str = <<"EOF";
    perl tmhmm_result.pl -i <input_file>
    A script for tmhmm results extraction and output to STDOUT.

Options:
    -i,--in     input file name
    -h,--help   help information
EOF
    return $help_str;
}

open my $tsv_in, '<', $input;
while (<$tsv_in>){
    chomp;
    my @array = split/\t/,$_;
    $array[4] =~ /^PredHel=(\d+)$/;
    my $hel_num = $1;
    $array[5] =~ /^Topology=(.+)$/;
    my $range = $1;
    while ($hel_num != 0){
        my $range_print = &RANGE_OUT($range);
        print "$array[0]\t$range_print\tTMD\n";
        $range =~ s/^[i|o]\d+?-\d+?([i|o].+$)/$1/;
        $hel_num = $hel_num - 1;
        if ($hel_num != 0){
            redo;
        }
        else{
            next;
        }
    }
}

sub RANGE_OUT{
    my $range_list = shift;
    $range_list =~ /^[i|o](\d+?)-(\d+?)[i|o]/;
    my $start = $1;
    my $end = $2;
    my $out = "$start\t$end";
    return $out;
}
