#!/usr/bin/perl -w

use strict;
use warnings;
use AlignDB::IntSpan;

my $gene;
our @print;
my $domain_set = AlignDB::IntSpan -> new;

while(<>){
    chomp;
    my $print_list = $_;
    my @array = split/\t/, $_;
    if ( ! defined $gene ){
        $gene = $array[0];
        push @print, $print_list;
        $domain_set -> add_range($array[3], $array[4]);
    }
    else{
        if ( $array[0] eq $gene ){
            my $test_set = AlignDB::IntSpan -> new;
            $test_set -> add_range($array[3], $array[4]);
            my $result = $domain_set -> intersect($test_set);
            if ( $result -> is_empty  ){
                $domain_set -> add_range($array[3], $array[4]);
                push @print, $print_list;
            }
        }
        else{
            print join "\n", @print;
            print "\n";
            $domain_set -> clear;
            $gene = $array[0];
            @print = ();
            push @print, $print_list;
            $domain_set -> add_range($array[3], $array[4])
        }
    }
}

END{
    print join "\n", @print;
    print "\n";
}
