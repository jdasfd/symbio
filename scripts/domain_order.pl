#!/usr/bin/perl -w

use strict;
use warnings;

our @domains;
our %RLK;
my $gene;
my $domain_list;

while(<>){
    chomp;
    my @array = split/\t/, $_;
    if ( ! defined $gene ){
        $gene = $array[0];
        push @domains, $array[1];
    }
    else{
        if ( $array[0] eq $gene ){
            push @domains, $array[1];
        }
        else{
            $domain_list = join ("#", @domains);
            if ( $domain_list =~ /.+?#TMD.+?pk.*/i ){
                $RLK{$gene}++;
            }
            $gene = $array[0];
            @domains = ();
            push @domains, $array[1];
        }
    }
}

END{
    $domain_list = join ("#", @domains);
    if ( $domain_list =~ /.+?#TMD.+?pk.*/i ){
        $RLK{$gene}++;
    }
    for my $keys (keys %RLK){
        if ( $RLK{$keys} == 1 ){
            print "$keys\n";
        }
        else{
            next;
        }
    }
}
