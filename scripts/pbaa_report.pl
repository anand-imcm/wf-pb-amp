use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case);
use FindBin qw($Bin);
use Sys::Hostname;
use File::Basename;

my $fastq_seqkit_stats = '';
my $passed_cluster_seqkit_stats = '';
my $failed_cluster_seqkit_stats = '';
my $cluster_qc = '';
my $cluster_variants_summ = '';
my $cluster_ontarget_variants_summ = '';
my $prefix = '';

my $help = '';

# GetOptions ('verbose' => \$verbose, 'all' => \$all);
GetOptions ('fastq=s' => \$fastq_seqkit_stats,
            'passSeq=s' => \$passed_cluster_seqkit_stats,
            'failSeq=s' => \$failed_cluster_seqkit_stats,
            'clusterQC=s' => \$cluster_qc,
            'clusterVariants=s' => \$cluster_variants_summ,
            'clusterOnTargetVariants=s' => \$cluster_ontarget_variants_summ,
            'prefix=s' => \$prefix,
			'help'  => \$help
			);

if (!$fastq_seqkit_stats || !$passed_cluster_seqkit_stats || !$failed_cluster_seqkit_stats || !$cluster_qc || !$cluster_variants_summ || !$cluster_ontarget_variants_summ || !$prefix){
    $help = 1;
}

if ($help) {
    print "Usage:\nperl scripts/pbaa_report.pl --fastq hifi_reads_fastq_seq_stats.tab --passSeq hifi_reads_pbaa_passed_cluster_sequences_stats.tab --failSeq hifi_reads_pbaa_failed_cluster_sequences_stats.tab --clusterQC hifi_reads_clusterQC_report.tsv --clusterVariants hifi_reads_variant_summary.tsv --clusterOnTargetVariants hifi_reads_variant_on_target_summary.tsv --prefix hifi_reads\n";
	exit;
}

my $unique_variant_summary=$prefix."_amplicon_workflow_variants_summary.tsv";
my $sequence_summary=$prefix."_amplicon_workflow_sequence_summary.tsv";


# #keys : {"file","format","type","num_seqs","sum_len","min_len","avg_len","max_len","Q1","Q2","Q3","sum_gap","N50","Q20%","Q30%","GC%"}
my %fastq_stats = parse_seqkit($fastq_seqkit_stats);
my %passed_cluster_stats = parse_seqkit($passed_cluster_seqkit_stats);
my %failed_cluster_stats = parse_seqkit($failed_cluster_seqkit_stats);

# #keys : {"target"=>"cluster"=>{"frequency","qualCount","qualMean","qualStd","qualMin","qual50%","qual90%","qualMax"}}
my %clusterQC_stats = parse_clusterQC($cluster_qc);

my ($count_snp,$count_indel,$count_ontarget_vars,$count_ontarget_snp,$count_ontarget_indel,@variants_tab) = parse_cluster_variants($cluster_variants_summ,$cluster_ontarget_variants_summ);

open (SEQ,">$sequence_summary") or die("Cannot write to - $sequence_summary");
print SEQ "file\tfastq_num_seqs\tfastq_sum_len\tfastq_min_len\tfastq_avg_len\tfastq_max_len\tfastq_Q1\tfastq_Q2\tfastq_Q3\tfastq_sum_gap\tfastq_N50\tfastq_Q20(%)\tfastq_Q30(%)\tfastq_GC(%)\tpbaa_passed_num_seq\tpbaa_failed_num_seq\ttotal_variants\ttotal_snps\ttotal_indels\ttotal_ontarget_variants\ttotal_ontarget_snps\ttotal_ontarget_indels\n";
print SEQ "$fastq_stats{'file'}\t$fastq_stats{'num_seqs'}\t$fastq_stats{'sum_len'}\t$fastq_stats{'min_len'}\t$fastq_stats{'avg_len'}\t$fastq_stats{'max_len'}\t$fastq_stats{'Q1'}\t$fastq_stats{'Q2'}\t$fastq_stats{'Q3'}\t$fastq_stats{'sum_gap'}\t$fastq_stats{'N50'}\t$fastq_stats{'Q20%'}\t$fastq_stats{'Q30%'}\t$fastq_stats{'GC%'}\t$passed_cluster_stats{'num_seqs'}\t$failed_cluster_stats{'num_seqs'}\t".scalar @variants_tab."\t$count_snp\t$count_indel\t$count_ontarget_vars\t$count_ontarget_snp\t$count_ontarget_indel\n";
close SEQ;

open (TMP,">temp.summary.tsv") or die("Cannot write to - temp.summary.tsv");
print TMP "Chr\tPos\tRef\tAlt\tis_on_target\tClinvar_consequence\tpbaa_guide_name\tcluster\n";
print TMP join("\n",@variants_tab)."\n";
close TMP;

my $sorted_variants = `head -n 1 temp.summary.tsv && tail -n +2 temp.summary.tsv | sort -k2,2n`;
open (VARS,">$unique_variant_summary") or die("Cannot write to - $unique_variant_summary");
print VARS $sorted_variants;
close VARS;

system("rm temp.summary.tsv");

sub parse_cluster_variants {
    my ($tsv_all,$tsv_ontarget) = (@_);
    my @variants;
    my %on_target;
    my %all_variants;
    my $count_ontarget_vars=0;
    my $count_snp=0;
    my $count_indel=0;
    my $count_ontarget_snp=0;
    my $count_ontarget_indel=0;

    open (OT,"$tsv_ontarget") or die("Cannot read table - $tsv_ontarget: $!");
    my $header_ontarget = <OT>;
    my @header_ontarget=split("\t",$header_ontarget);
    while(my $rec_ontarget=<OT>){
        chomp $rec_ontarget;
        # Sample	Target	Cluster	Numreads	CHROM	POS	ID	REF	ALT	BCSQ
        my @rec_ontarget_cols=split("\t",$rec_ontarget);
        my $variantKey="$rec_ontarget_cols[4]:$rec_ontarget_cols[5]:$rec_ontarget_cols[7]>$rec_ontarget_cols[8]";
        push @{$on_target{$variantKey}{'target'}}, $rec_ontarget_cols[1];
        push @{$on_target{$variantKey}{'cluster'}}, $rec_ontarget_cols[2];
        push @{$on_target{$variantKey}{'BCSQ'}}, $rec_ontarget_cols[-1];
    }

    open (ALL,"$tsv_all") or die("Cannot read table - $tsv_all: $!");
    my $header_all = <ALL>;
    my @header_all=split("\t",$header_all);
    while(my $rec_all=<ALL>){
        chomp $rec_all;
        # Sample	Target	Cluster	Numreads	CHROM	POS	ID	REF	ALT	BCSQ
        my @rec_all_cols=split("\t",$rec_all);
        my $variantKey="$rec_all_cols[4]:$rec_all_cols[5]:$rec_all_cols[7]>$rec_all_cols[8]";
        push @{$all_variants{$variantKey}{'target'}}, $rec_all_cols[1];
        push @{$all_variants{$variantKey}{'cluster'}}, $rec_all_cols[2];
        push @{$all_variants{$variantKey}{'BCSQ'}}, $rec_all_cols[-1];
        if (defined($on_target{$variantKey})) {
            $all_variants{$variantKey}{'is_ontarget'} = "Y";
        }
        else {
            $all_variants{$variantKey}{'is_ontarget'}= "N";
        }
    }
    close OT;
    close ALL;

    foreach my $var(keys %all_variants){
        my %seen_c;
        my %seen_t;
        my %seen_b;
        my @unique_cluster = grep { !$seen_c{$_}++ } @{$all_variants{$var}{'cluster'}};
        my @unique_target = grep { !$seen_t{$_}++ } @{$all_variants{$var}{'target'}};
        my @unique_BCSQ = grep { !$seen_b{$_}++ } @{$all_variants{$var}{'BCSQ'}};
        my ($chr,$pos,$ref,$alt) = split(/[:>]/, $var);
        my $var_string = "$chr\t$pos\t$ref\t$alt\t$all_variants{$var}{'is_ontarget'}\t".join("#",@unique_BCSQ)."\t".join(",",@unique_target)."\t".join(",",@unique_cluster);
        if (length($ref) == 1 && length($alt) == 1){
            $count_snp++;
        }
        else {
            $count_indel++;
        }
        if ($all_variants{$var}{'is_ontarget'} eq "Y"){
            $count_ontarget_vars++;
            if (length($ref) == 1 && length($alt) == 1){
                $count_ontarget_snp++;
            }
            else {
                $count_ontarget_indel++;
            }
        }
        push @variants,$var_string;
    }
    return($count_snp,$count_indel,$count_ontarget_vars,$count_ontarget_snp,$count_ontarget_indel,@variants);
}

sub parse_clusterQC {
    my ($tsv) = (@_);
    my %stats;
    open (TAB,"$tsv") or die("Cannot read table - $tsv: $!");
    my @headers=("target","cluster","frequency","qualCount","qualMean","qualStd","qualMin","qual50%","qual90%","qualMax");
    while(my $rec=<TAB>){
        chomp $rec;
        if ($rec =~ /^GBA/){
            my @rec_stats = split("\t",$rec);
            for(my $i=2;$i<=$#rec_stats;$i++){
                $stats{$rec_stats[0]}{abs($rec_stats[1])}{$headers[$i]}=$rec_stats[$i];
            }
        }
    }
    close TAB;
    return(%stats);
}

sub parse_seqkit {
    my ($tsv) = (@_);
    my %stats;
    open (TAB,"$tsv") or die("Cannot read table - $tsv: $!");
    my $header = <TAB>; chomp $header;
	my @seqkit_headers = split("\t",$header);
    
    my $stats = <TAB>; chomp $stats;
	my @seqkit_stats = split("\t",$stats);

    for(my $i=0; $i<=$#seqkit_headers;$i++){
        $seqkit_headers[$i] =~ s/[()]//g;
        $seqkit_stats[$i] =~ s/.hifi_reads.fastq|_pbaa_passed_cluster_sequences.fasta|_pbaa_failed_cluster_sequences.fasta//;

        $stats{$seqkit_headers[$i]}=$seqkit_stats[$i];
    }
    close TAB;
    return(%stats);
}