# perl exon2GTF.pl $outfile $infile $source


use strict;
use warnings;
use Data::Dumper;



my $outfile = shift;
open(my $out, '>', $outfile) or die "Could not open file '$outfile' $!";

my $infile=shift;

my $source = shift;

#my $infile="test.txt"; 
my %data_store;
my %data_details;
open(my $in, '<', $infile ) or die "Can't open $infile: $!";
    while ( my $line = <$in> ) {
		my @pos = split /\t/, $line;
		
		my $chr = $pos[0];
		my $start = $pos[1];
		my $end = $pos[2];
		my $strand = $pos[3];

		my $gene = $pos[4];
		my $ORF = $pos[5];
		my $exon_id = $chr."_".$start."_".$end."_".$ORF;
		my $gene_name = $pos[7];
		my $orf_type = $pos[8];
		my $gene_type = $pos[9];	
		my $pept = $pos[10];

		$data_store{$gene}{$ORF}{$exon_id} = 1;
		
		if(not(exists($data_details{$gene}))){
		$data_details{$gene}{'chr'} = $chr;
		$data_details{$gene}{'start'} = $start;
		$data_details{$gene}{'end'} = $end;
		$data_details{$gene}{'strand'} = $strand;
		$data_details{$gene}{'info_strand'} = "";
		$data_details{$gene}{'gene_type'} =$gene_type;
		$data_details{$gene}{'gene_name'} = $gene_name;
		}else{
			if($data_details{$gene}{'start'} > $start){
				$data_details{$gene}{'start'} = $start;
			}
			if($data_details{$gene}{'end'} < $end){
                                $data_details{$gene}{'end'} = $end;
                        }
			if($data_details{$gene}{'strand'} ne $strand){
				$data_details{$gene}{'info_strand'} = 'mismatch_strand';
				$strand = $data_details{$gene}{'strand'};
			}
		}
		
		if(not(exists($data_details{$ORF}))){
                $data_details{$ORF}{'chr'} = $chr;
		$data_details{$ORF}{'start'} = $start;
                $data_details{$ORF}{'end'} = $end;
                $data_details{$ORF}{'gene_id'} = $gene;
		$data_details{$ORF}{'strand'} = $strand;
		$data_details{$ORF}{'info_chr'} = "";
		$data_details{$ORF}{'orf_type'} = $orf_type;
		$data_details{$ORF}{'pept'} = $pept;
		}else{
                        if($data_details{$ORF}{'start'} > $start){
                                $data_details{$ORF}{'start'} = $start;
                        }
                        if($data_details{$ORF}{'end'} < $end){
                                $data_details{$ORF}{'end'} = $end; 
                        }
                	if($data_details{$ORF}{'strand'} ne $strand){
                                $data_details{$ORF}{'info_strand'} = 'mismatch_strand';
                        }
			if($data_details{$ORF}{'pept'} ne $pept){
				$data_details{$ORF}{'info_strand'} = 'mismatch peptide';
			}
		}

		if(not(exists($data_details{$exon_id}))){
                	$data_details{$exon_id}{'chr'} = $chr;
			$data_details{$exon_id}{'start'} = $start;
                	$data_details{$exon_id}{'end'} = $end;
			$data_details{$exon_id}{'gene_id'} = $gene;
			$data_details{$exon_id}{'ORF_id'} = $ORF;
                	$data_details{$exon_id}{'strand'} = $strand;
			}
		else{
			print $exon_id;
		}
	}
    close $in;


print Dumper \%data_details;


foreach my $gene (keys %data_store){
	print $out $data_details{$gene}{'chr'}."\t".$source."\tgene\t".$data_details{$gene}{'start'}."\t".$data_details{$gene}{'end'}."\t.\t".$data_details{$gene}{'strand'}."\t.\tgene_id \"".$gene."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"; other info: ".$data_details{$gene}{'info_strand'}."\n";
	foreach my $ORF (keys %{$data_store{$gene}}){
		print $out $data_details{$ORF}{'chr'}."\t".$source."\tORF\t".$data_details{$ORF}{'start'}."\t".$data_details{$ORF}{'end'}."\t.\t".$data_details{$ORF}{'strand'}."\t.\tgene_id \"".$data_details{$ORF}{'gene_id'}."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"; ORF_id \"".$ORF."\"; ORF_type \"".$data_details{$ORF}{'orf_type'}."\"; ORF_pept \"".$data_details{$ORF}{'pept'}."\"; other info: ".$data_details{$gene}{'info_strand'}."\n";
	foreach my $exon (keys %{$data_store{$gene}{$ORF}}){
			print $out $data_details{$exon}{'chr'}."\t".$source."\texon\t".$data_details{$exon}{'start'}."\t".$data_details{$exon}{'end'}."\t.\t".$data_details{$exon}{'strand'}."\t.\tgene_id \"".$data_details{$exon}{'gene_id'}."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"; ORF_id \"".$data_details{$exon}{'ORF_id'}."\"; ORF_type \"".$data_details{$ORF}{'orf_type'}."\";\n";
		}
	}
}

close $out;
