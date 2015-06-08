#!/usr/bin/env perl
use strict;
use warnings;

#Index your reference with the blast formatdb command
#Run megablast with: megablast -i scaffolds.fa -d GRCm38_68.fa -e 0.01 -m 9 -p 90 -W 250
#requires the megablast output to be put through: sort -k1,1d -k7,7n (sort by scaff position)
#e.g. cat *.megablast.out | sort -k1,1d -k7,7n | perl megablast_error_rate.pl

my %hits;
my $lastendchr=-1;my $lastendscaf=-1;my $lastchr='';my $lastscaf='';
while(my $l=<>)
{
    next if($l=~/^#/);
    chomp($l);
    my @s = split(/\t/, $l);
    
    $hits{$s[0]}{$s[1]}+=($s[7]-$s[6]);
    
    $lastchr=$s[1];$lastendchr=$s[9];$lastscaf=$s[0];$lastendscaf=$s[7];
}

my $totalerror=0;my $total=0;
foreach my $scaf(keys(%hits))
{
    my $biggesthit=0;my $biggestchr='';
    foreach my $hitchr(keys(%{$hits{$scaf}}))
    {
        if($hits{$scaf}{$hitchr}>$biggesthit){$biggesthit=$hits{$scaf}{$hitchr};$biggestchr=$hitchr;}
    }
    
    my $errorbp=0;
    foreach my $hitchr(keys(%{$hits{$scaf}}))
    {
        $errorbp+=$hits{$scaf}{$hitchr} if($hitchr ne $biggestchr);
    }
    
    $totalerror+=$errorbp;$total+=$errorbp+$biggesthit;
    print qq[SCAF: $scaf $biggestchr $biggesthit $errorbp\n];
}
my $e = $totalerror/$total;
print qq[TOTAL: $total $totalerror $e\n];

