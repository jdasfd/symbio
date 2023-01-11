#!/usr/bin/perl -w

use strict;
use warnings;

our @domains;
our %RLK;
our %RLCK;
my $gene;
my $domain_list;

while(<>){
    chomp;
    my @array = split/\t/, $_;
    if ( ! defined $gene ){
        $gene = $array[0];
        push @domains, $array[3];
    }
    else{
        if ( $array[0] eq $gene ){
            push @domains, $array[3];
        }
        else{
            $domain_list = join ("#", @domains);
            if ( $domain_list =~ /.+?#TMD.+?pk.*/i ) {
                $RLK{$gene}++;
            }
            elsif ( $domain_list =~ /^TMD.*pk.*/i ) {
                $RLCK{$gene}++;
            }
            $gene = $array[0];
            @domains = ();
            push @domains, $array[3];
        }
    }
}

END{
    $domain_list = join ("#", @domains);
    if ( $domain_list =~ /.+?#TMD.+?pk.*/i ) {
        $RLK{$gene}++;
    }
    elsif ( $domain_list =~ /^TMD.*pk.*/i ) {
        $RLCK{$gene}++;
    }
    for my $keys (keys %RLK){
        if ( $RLK{$keys} == 1 ){
            print "$keys\tRLK\n";
        }
        else{
            next;
        }
    }
    for my $keys (keys %RLCK){
        if ( $RLCK{$keys} == 1 ){
            print "$keys\tRLCK_with_TMD\n";
        }
    }
}
