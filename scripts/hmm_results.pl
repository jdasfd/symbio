#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;

=head1 NAME
hmm_result.pl - Extract info from the HMM results.
=head 1 SYNOPSIS
    perl hmm_result.pl -i <input_file>
Options:
    -i,--in     input file name
    -h,--help   help information
=cut

GetOptions(
    "i|in=s"        => \(my $input),
    "h|help"        => \(my $help),
);

die usage() if defined $help;
die ("Input a file please.") if not defined $input;

sub usage{
    my $help_str = <<"EOF";
    perl hmm_result.pl -i <input_file>
    A script for hmm results extraction and output to STDOUT.

Options:
    -i,--in     input file name
    -h,--help   help information
EOF
    return $help_str;
}

print "QUERY\tDomain\tE_value\tStart\tEnd\tQ_seq\tHit_seq\n";

my $searchio = Bio::SearchIO -> new ( -format => 'hmmer',
                                      -file => $input,
                                      -best => "true", );

while ( my $result = $searchio -> next_result() ) {
    while ( my $hit = $result -> next_hit ) {
        while ( my $hsp = $hit -> next_hsp ){
            my $query = $result -> query_name;
            my $hit_name = $hit -> name();
            my $evalue = $hsp -> evalue();
            my $start = $hsp -> start('query');
            my $end = $hsp -> end('query');
            my $qseq = $hsp -> query_string;
            $qseq = uc $qseq;
            my $hseq = $hsp -> hit_string;
            $hseq = uc $hseq;
            print "$query\t$hit_name\t$evalue\t$start\t$end\t$qseq\t$hseq\n";
        }
    }
}
