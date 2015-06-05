package GFF;
use strict;
use warnings;

my $annotation=qq[gencode.vM4.annotation.gff3.gz];
my $MIN_EXON_SIZE = 50;

my %gene2exons;
my %exon2gene;
my %transcript2exons;
my %exon2transcripts;
my %transcriptorientations;
my %transcripts2gene;
my %exons;
my %gene2name;

sub getExonsFromGene
{
    my $gene = shift;
    
    if(! %gene2exons)
    {
        print STDERR qq[loading gene2exon\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        my %exons;my %genes;my %temp;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            if($$entry{type} eq 'exon'&&$$entry{gene_type} eq 'protein_coding'&&$$entry{exon_id}&&$$entry{exon_id}=~/^ENS/)
            {
                $temp{$$entry{gene_id}}{$$entry{orientation} eq '+' ? qq[$$entry{start}-$$entry{end}] : qq[$$entry{end}-$$entry{start}]}=$$entry{exon_id};
                $exons{$$entry{exon_id}}=1;$genes{$$entry{gene_id}}=1;
            }
        }
        
        #now convert each entry into top strand ordered list of exons
        foreach my $gene(keys(%temp))
        {
            my @ex;
            foreach my $pos(sort {($a=~/^(\d+)-/)[0]<=>($b=~/^(\d+)-/)[0]}keys(%{$temp{$gene}}))
            {
                push(@ex, $temp{$gene}{$pos});
            }
            $gene2exons{$gene} = \@ex;
        }
        
        print STDERR qq[Found ].scalar(keys(%exons)).qq[ unique exons\n];
        print STDERR qq[Found ].scalar(keys(%genes)).qq[ unique genes\n];
        print STDERR qq[loaded gene2exon\n];
    }
    
    #want to return an ordered list as appear on chromosome of exons
    if( ! $gene2exons{$gene}){return undef;}else{return $gene2exons{$gene};}
}

sub getExonsFromTranscript
{
    my $transcript = shift;
    
    if(! %transcript2exons )
    {
        print STDERR qq[loading transcript2exons\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            if($$entry{type} eq 'exon'&&$$entry{gene_type} eq 'protein_coding'&&$$entry{exon_id}&&$$entry{exon_id}=~/^ENS/&&$$entry{transcript_id}&&$$entry{transcript_id}=~/^ENS/)
            {
                $transcript2exons{$$entry{transcript_id}}{$$entry{orientation} eq '+' ? $$entry{start} : $$entry{end}}=$$entry{exon_id};
                $transcriptorientations{$$entry{transcript_id}}=$$entry{orientation};
            }
        }
        print STDERR qq[loaded transcript2exons\n];
        print STDERR qq[Found ].scalar(keys(%transcript2exons)).qq[ transcripts\n];
    }
    
    if( !$transcript2exons{$transcript}){print STDERR qq[Failed to find exons for: $transcript\n]; return undef;}
    else
    {
        my @ex;
        foreach my $pos(sort {$a<=>$b}keys(%{$transcript2exons{$transcript}})){push(@ex,$transcript2exons{$transcript}{$pos});}
        return \@ex;
    }
}

sub getGeneFromExon
{
    my $exon = shift;
    
    if(! %exon2gene )
    {
        print STDERR qq[loading exon2gene\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        my %genes;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            
            if($$entry{type} eq 'exon'&&$$entry{gene_type} eq 'protein_coding'&&$$entry{exon_id}&&$$entry{exon_id}=~/^ENS/)
            {
                $exon2gene{$$entry{exon_id}}=$$entry{gene_id};
            }
        }
        print STDERR qq[loaded exon2gene\n];
    }
    
    if(!$exon2gene{$exon}){return undef}else{return $exon2gene{$exon}}
}

sub getGeneFromTranscript
{
    my $transcript = shift;
    
    if( ! %transcripts2gene)
    {
        print STDERR qq[loading transcripts2gene\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        my %genes;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            
            if($$entry{type} eq 'exon'&&$$entry{gene_type} eq 'protein_coding'&&$$entry{exon_id}&&$$entry{exon_id}=~/^ENS/)
            {
                $transcripts2gene{$$entry{transcript_id}}=$$entry{gene_id};
            }
        }
        print STDERR qq[loaded transcripts2gene\n];
    }
    
    if(!$transcripts2gene{$transcript}){print STDERR qq[cant find gene for $transcript\n];return undef;}else{return $transcripts2gene{$transcript};}
}

sub getTranscriptsFromExon
{
    my $exon = shift;
    
    if( !%exon2transcripts)
    {
        print STDERR qq[loading exon2transcripts\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            
            if($$entry{type} eq 'exon'&&$$entry{gene_type} eq 'protein_coding'&&$$entry{exon_id}&&$$entry{exon_id}=~/^ENS/&&$$entry{transcript_id}&&$$entry{transcript_id}=~/^ENS/)
            {
                if($exon2transcripts{$$entry{exon_id}}){push($exon2transcripts{$$entry{exon_id}},$$entry{transcript_id});}
                else{$exon2transcripts{$$entry{exon_id}}=[$$entry{transcript_id}];}
            }
        }
        print STDERR qq[loaded exon2transcripts\n];
    }
    if($exon2transcripts{$exon}){return $exon2transcripts{$exon};}{print STDERR qq[Cant find transcripts for: $exon\n];return undef;}
}

sub getExons
{
    my $minsize = shift;
    
    if( !%exons)
    {
        print STDERR qq[loading exons\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            if($$entry{type} eq 'exon'&&$$entry{gene_type} eq 'protein_coding'&&$$entry{exon_id}&&$$entry{exon_id}=~/^ENS/&&$$entry{transcript_id}&&$$entry{transcript_id}=~/^ENS/)
            {
                $exons{$$entry{exon_id}} = [$$entry{chr}, $$entry{start}, $$entry{end}, $$entry{orientation}] if( $$entry{end}-$$entry{start}>=$minsize);
            }
        }
        print STDERR qq[loaded exons\n];
    }
    
    return \%exons;
}

sub getNameFromExon
{
    my $exon = shift;
    
    if(!%gene2name)
    {
        print STDERR qq[loading gene2name\n];
        open(my $tfh, qq[zcat $annotation | ] ) or die $!;
        while(my $l=<$tfh>)
        {
            next if($l=~/^#/);
            my $entry=_parseline($l);
            next if(!$entry);
            if($$entry{type} eq 'gene'&&$$entry{gene_type} eq 'protein_coding')
            {
                $gene2name{$$entry{gene_id}}=$$entry{gene_name};
            }
        }
        print STDERR qq[loaded gene2name ].(scalar(keys(%gene2name))).qq[\n];
    }
    
    my $gene = getGeneFromExon($exon);
    if($gene&&$gene2name{$gene}){return $gene2name{$gene};}else{return undef;}
}

sub getFeaturesInfo
{
    my $featuresfofn = shift;
    my $featuretype = shift;
    my $fieldstr = shift; #string of fields to return in tab format with the id
    my $outfile = shift;
    
    my @fields = split(/\s/, $fieldstr);
    my %ids;open(my $ifh,$featuresfofn)or die $!;while(my $l=<$ifh>){chomp($l);$ids{$l}=1;}close($ifh);
    
    open( my $ofh, qq[>$outfile]) or die $!;
    print STDERR qq[Reading annotation file\n];
    open(my $tfh, qq[zcat $annotation | ] ) or die $!;
    while(my $l=<$tfh>)
    {
        next if($l=~/^#/);
        my $entry=_parseline($l);
        next if(!$entry);
        if($$entry{$featuretype}&&$ids{$$entry{$featuretype}})
        {
            foreach my $field(@fields)
            {
                if($$entry{$field}){print $ofh qq[$$entry{$field}\t];}else{print qq[unknown\t];}
            }
            print $ofh qq[\n];
            delete $ids{$$entry{$featuretype}};
        }
    }
    close($ofh);
}

sub _parseline
{
    my $line = shift;
    my %entry;
    my @s = split(/\t/,$line);
    return undef if(($s[4]-$s[3])<$MIN_EXON_SIZE);
    
    $entry{chr}=$s[0]=~/^chr/?substr($s[0],3):$s[0];$entry{chr}='MT' if($entry{chr}=~/^M$/);
    $entry{source}=$s[1];$entry{type}=$s[2];$entry{start}=$s[3];$entry{end}=$s[4];$entry{orientation}=$s[6];
    my @s1 = split(/;/,$s[8]);
    foreach my $e(@s1)
    {
        my @s2=split(/=/,$e);
        if($s2[0] eq 'gene_type'){$entry{gene_type}=$s2[1];}
        if($s2[0] eq 'exon_id'){$entry{exon_id}=$s2[1];}
        if($s2[0] eq 'gene_id'){$entry{gene_id}=$s2[1];}
        if($s2[0] eq 'transcript_id'){$entry{transcript_id}=$s2[1];}
        if($s2[0] eq 'gene_name' ){$entry{gene_name}=$s2[1];}
    }
    return \%entry;
}

sub revCompDNA
{
	die "Usage: revCompDNA string\n" unless @_ == 1;
	my $seq = shift;
	$seq= uc( $seq );
	$seq=reverse( $seq );
	
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

1;
