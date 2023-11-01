use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case);
use FindBin qw($Bin);
use Sys::Hostname;
use File::Basename;

my $fastq_seqkit_stats = '';
my $pbmm2_log = '';
my $dv_variants_summ = '';
my $dv_ontarget_variants_summ = '';
my $prefix = '';

my $help = '';

# GetOptions ('verbose' => \$verbose, 'all' => \$all);
GetOptions ('fastq=s' => \$fastq_seqkit_stats,
            'pbmmlog=s' => \$pbmm2_log,
            'allVariants=s' => \$dv_variants_summ,
            'onTargetVariants=s' => \$dv_ontarget_variants_summ,
            'prefix=s' => \$prefix,
			'help'  => \$help
			);

if (!$fastq_seqkit_stats || !$pbmm2_log || !$dv_variants_summ || !$dv_ontarget_variants_summ || !$prefix){
    $help = 1;
}

if ($help) {
    print "Usage:\nperl scripts/deepvariant_report.pl --fastq hifi_reads-3_raw_hifi_reads_fastq_stats.tab --pbmmlog hifi_reads-3_raw_hifi_to_reference_alignment.log --allVariants hifi_reads-3_raw_hifi_to_reference_alignment_pass_variants_annotated_summary.tsv --onTargetVariants hifi_reads-3_raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary.tsv --prefix test\n";
	exit;
}

my $unique_variant_summary=$prefix."_germline_workflow_variants_summary.tsv";
my $sequence_summary=$prefix."_germline_workflow_sequence_summary.tsv";


# #keys : {"file","format","type","num_seqs","sum_len","min_len","avg_len","max_len","Q1","Q2","Q3","sum_gap","N50","Q20%","Q30%","GC%"}
my %fastq_stats = parse_seqkit($fastq_seqkit_stats);

my $total_mapped_reads = parse_pbmm2log($pbmm2_log);

my $alignment_perc = 0;

my $total_unmapped_reads = $fastq_stats{'num_seqs'} - $total_mapped_reads;

my ($count_snp,$count_indel,$count_ontarget_vars,$count_ontarget_snp,$count_ontarget_indel,@variants_tab) = parse_dv_variants($dv_variants_summ,$dv_ontarget_variants_summ);

if ($total_mapped_reads != 0){
    $alignment_perc = sprintf("%.2f",($total_mapped_reads/$fastq_stats{'num_seqs'})*100);
}

open (SEQ,">$sequence_summary") or die("Cannot write to - $sequence_summary");
print SEQ "file\tfastq_num_seqs\tfastq_sum_len\tfastq_min_len\tfastq_avg_len\tfastq_max_len\tfastq_Q1\tfastq_Q2\tfastq_Q3\tfastq_sum_gap\tfastq_N50\tfastq_Q20(%)\tfastq_Q30(%)\tfastq_GC(%)\ttotal_mapped_reads\ttotal_unmapped_reads\ttotal_alignment%\ttotal_variants\ttotal_snps\ttotal_indels\ttotal_ontarget_variants\ttotal_ontarget_snps\ttotal_ontarget_indels\n";
print SEQ "$prefix\t$fastq_stats{'num_seqs'}\t$fastq_stats{'sum_len'}\t$fastq_stats{'min_len'}\t$fastq_stats{'avg_len'}\t$fastq_stats{'max_len'}\t$fastq_stats{'Q1'}\t$fastq_stats{'Q2'}\t$fastq_stats{'Q3'}\t$fastq_stats{'sum_gap'}\t$fastq_stats{'N50'}\t$fastq_stats{'Q20%'}\t$fastq_stats{'Q30%'}\t$fastq_stats{'GC%'}\t$total_mapped_reads\t$total_unmapped_reads\t$alignment_perc\t".scalar @variants_tab."\t$count_snp\t$count_indel\t$count_ontarget_vars\t$count_ontarget_snp\t$count_ontarget_indel\n";
close SEQ;

open (TMP,">temp.summary.tsv") or die("Cannot write to - temp.summary.tsv");
print TMP "Chr\tPos\tRef\tAlt\tis_on_target\tClinvar_consequence\tvcf_info\n";
print TMP join("\n",@variants_tab)."\n";
close TMP;

my $sorted_variants = `head -n 1 temp.summary.tsv && tail -n +2 temp.summary.tsv | sort -k2,2n`;
open (VARS,">$unique_variant_summary") or die("Cannot write to - $unique_variant_summary");
print VARS $sorted_variants;
close VARS;

system("rm temp.summary.tsv");

sub parse_dv_variants {
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
        # CHROM	POS	ID	REF	ALT	SAMPLE	BCSQ	GT:GQ:DP:AD:VAF:PL:BCSQ
        my @rec_ontarget_cols=split("\t",$rec_ontarget);
        my $variantKey="$rec_ontarget_cols[0]:$rec_ontarget_cols[1]:$rec_ontarget_cols[3]>$rec_ontarget_cols[4]";
        push @{$on_target{$variantKey}{'info'}}, $rec_ontarget_cols[-1];
        push @{$on_target{$variantKey}{'BCSQ'}}, $rec_ontarget_cols[-2];
    }

    open (ALL,"$tsv_all") or die("Cannot read table - $tsv_all: $!");
    my $header_all = <ALL>;
    my @header_all=split("\t",$header_all);
    while(my $rec_all=<ALL>){
        chomp $rec_all;
        # CHROM	POS	ID	REF	ALT	SAMPLE	BCSQ	GT:GQ:DP:AD:VAF:PL:BCSQ
        my @rec_all_cols=split("\t",$rec_all);
        my $variantKey="$rec_all_cols[0]:$rec_all_cols[1]:$rec_all_cols[3]>$rec_all_cols[4]";
        push @{$all_variants{$variantKey}{'info'}}, $rec_all_cols[-1];
        push @{$all_variants{$variantKey}{'BCSQ'}}, $rec_all_cols[-2];
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
        my %seen_i;
        my %seen_b;
        my @unique_info = grep { !$seen_i{$_}++ } @{$all_variants{$var}{'info'}};
        my @unique_BCSQ = grep { !$seen_b{$_}++ } @{$all_variants{$var}{'BCSQ'}};
        my ($chr,$pos,$ref,$alt) = split(/[:>]/, $var);
        my $var_string = "$chr\t$pos\t$ref\t$alt\t$all_variants{$var}{'is_ontarget'}\t".join(",",@unique_BCSQ)."\t".join(",",@unique_info);
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

sub parse_pbmm2log {
    my ($log) = (@_);
    my $total_mapped = 0;
    my $mapped_reads = `grep "Mapped Reads" $log`;
    chomp $mapped_reads;
    my @mapped_reads = split(" ",$mapped_reads);
    $total_mapped = $mapped_reads[-1];
    return $total_mapped;
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