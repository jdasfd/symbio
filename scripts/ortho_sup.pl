#!/usr/bin/perl -w

use strict;
use warnings;

our %OGNUM;
my @group;

while(<>){
    chomp;
    next if /^Ortho/;
    @group = split/\t/, $_;
    $OGNUM{$group[0]} = $group[1];
}

for (my $i = 0; $i < 150270; $i++) {
    my $i_num = sprintf("%07d", $i);
    $i_num = "OG".$i_num;
    if (defined $OGNUM{$i_num}){
        print "$i_num\t$OGNUM{$i_num}\n";
    }
    else{
        print "$i_num\t0\n";
    }
}
