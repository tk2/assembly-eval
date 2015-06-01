#!/usr/bin/env perl

use strict;
use warnings;

use GFF;

my $BAMFLAGS = 
{
    'paired_tech'    => 0x0001,
    'read_paired'    => 0x0002,
    'unmapped'       => 0x0004,
    'mate_unmapped'  => 0x0008,
    'reverse_strand' => 0x0010,
    'mate_reverse'   => 0x0020,
    '1st_in_pair'    => 0x0040,
    '2nd_in_pair'    => 0x0080,
    'secondary'    => 0x0100,
    'failed_qc'      => 0x0200,
    'duplicate'      => 0x0400,
};

die qq[usage: bam] unless @ARGV==1;
my $input=shift;

=pod
1. Are all the exons for a gene present?

2. Are all the exons for a gene on same chr?

3. Are all the exons for a gene in the same orientation?

4. Are all the exons for a gene present?

5. Are they in the correct order?
=cut

#open(my $ifh, $cigar=~/\.gz$/?qq[zcat $cigar|] : qq[cat $cigar|] ) or die $!;
open(my $ifh, qq[samtools view $input|] ) or die $!;
my %trans_chr_exon; #transcript -> chr -> exon = 1
while(my $l=<$ifh>)
{
    chomp($l);
    my @s=split(/\s+/,$l);
    next if($s[2] eq '*');
    my @trans = @{GFF::getTranscriptsFromExon($s[0])};
    die qq[no transcript $s[1]\n] if(!@trans);
    foreach my $t(@trans)
    {
        die qq[undefined transcript name for $s[1]\n] if !$t;
        $trans_chr_exon{$t}{$s[2]}{$s[0]}=$s[1]&$$BAMFLAGS{'reverse_strand'} ? '-':'+';
    }
}
close($ifh);

print qq[#Gene num_chrs num_exons num_exons_found num_fwd num_rev missing_ex1,missing_ex2,... chr:exon1,exon2,...\n];

foreach my $trans(keys(%trans_chr_exon))
{
    my @exons=@{GFF::getExonsFromTranscript($trans)};
    
    my $chrshit=0;
    my %exonsF;my $chrs=0;my $chrExStr;my $numfwd=0;my $numrev=0;
    foreach my $chr(keys(%{$trans_chr_exon{$trans}}))
    {
        $chrshit++;
        $chrExStr.=qq[ $chr:];
        my @e=keys(%{$trans_chr_exon{$trans}{$chr}});
        foreach my $t(@e)
        {
            $exonsF{$t}=1;$chrExStr.=qq[,$t];
            if( $trans_chr_exon{$trans}{$chr}{$t} eq '+' ){$numfwd++;}else{$numrev++;}
        }
        $chrs++;
    }
    
    #compare how many exons found vs expected
    print qq[$trans $chrs ].scalar(@exons).qq[ ].scalar(keys(%exonsF)).qq[ $numfwd $numrev ];
    
    if( scalar(keys(%exonsF)) < scalar(@exons) )
    {
        #print missing
        foreach my $e(@exons){my $f=0;foreach my $e1(keys(%exonsF)){if($e eq $e1){$f=1;last;}}print qq[$e,]if($f==0);}
    }
    print qq[$chrExStr\n];
}
