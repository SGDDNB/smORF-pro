# perl collapseGTF.pl $outfile $infile

use strict;
use warnings;
use Data::Dumper;


my $outfile = shift;
open(my $out, '>', $outfile) or die "Could not open file '$outfile' $!";

my $infile=shift;
#my $infile="test.txt";
my %data_store;
my %data_details;
open(my $in, '<', $infile ) or die "Can't open $infile: $!";
    while ( my $line = <$in> ) {
	my @pos = split /\t/, $line;

	my $chr = $pos[0];
        my $start = $pos[3];
        my $end = $pos[4];
        my $strand = $pos[6];
	my $type = $pos[2];
	my $source = $pos[1];
	
	my @detail = split('"',$pos[8]);
	
	my $gene = $detail[1];
	my $gene_type = $detail[3];
	my $gene_name = $detail[5];
	my $ORFI_id;
	my $exon_id;
	my $ORFI_type;	
	if($type eq "gene"){
		$data_details{$gene}{'chr'} = $chr;
                $data_details{$gene}{'start'} = $start;
                $data_details{$gene}{'end'} = $end;
		$data_details{$gene}{'source'} = $source;
		$data_details{$gene}{'gene_type'} = $gene_type;
		$data_details{$gene}{'strand'} = $strand;
		$data_details{$gene}{'gene_name'} = $gene_name;
	}

	elsif($type eq "ORF"){
		$ORFI_id = $detail[7];
		$ORFI_type = $detail[9];
		$data_details{$ORFI_id}{'chr'} = $chr;
		$data_details{$ORFI_id}{'start'} = $start;
                $data_details{$ORFI_id}{'end'} = $end;
		$data_details{$ORFI_id}{'source'} = $source;
		$data_details{$ORFI_id}{'strand'} = $strand;	
		$data_details{$ORFI_id}{'ORFI_type'} = $ORFI_type;  
		$data_details{$ORFI_id}{'ORFI_pept'} = $detail[11];
	}
	
	elsif($type eq "exon"){
		$ORFI_id = $detail[7];
		$ORFI_type = $detail[9];
		$exon_id = $chr."_".$start."_".$end."_".$ORFI_id;
		$data_store{$gene}{$ORFI_id}{$exon_id} = 1;
		$data_details{$exon_id}{'chr'} = $chr;
                $data_details{$exon_id}{'start'} = $start;
                $data_details{$exon_id}{'end'} = $end;
		$data_details{$exon_id}{'source'} = $source;
		$data_details{$exon_id}{'strand'} = $strand;
		$data_details{$exon_id}{'ORFI_type'} = $ORFI_type;
	}
}
#print Dumper \%data_store;
#print Dumper \%data_details;

close $in;
my %data_store_groups;

foreach my $gene (keys %data_store){
	my @tx_3ps;
	foreach my $ORFI (keys %{$data_store{$gene}}){
		if($data_details{$ORFI}{'strand'} eq "-"){
			#print $data_details{$ORFI}{'start'}."\n";
			if($data_details{$ORFI}{'start'} ~~ @tx_3ps ){
				my $ORF_id = join("_",$gene,$data_details{$ORFI}{'start'});
				$data_store_groups{$gene}{$ORF_id}{$ORFI} = 1;
				if($data_details{$ORF_id}{'end'} < $data_details{$ORFI}{'end'}){
					$data_details{$ORF_id}{'end'}=$data_details{$ORFI}{'end'}
				}
			}
			else{
				push @tx_3ps, $data_details{$ORFI}{'start'};
				my $ORF_id = join("_",$gene,$data_details{$ORFI}{'start'});
				$data_store_groups{$gene}{$ORF_id}{$ORFI} = 1;
				$data_details{$ORF_id}{'chr'} = $data_details{$ORFI}{'chr'};
                		$data_details{$ORF_id}{'start'} = $data_details{$ORFI}{'start'};
                		$data_details{$ORF_id}{'end'} = $data_details{$ORFI}{'end'};
               		 	$data_details{$ORF_id}{'source'} = $data_details{$ORFI}{'source'};
                		$data_details{$ORF_id}{'strand'} = $data_details{$ORFI}{'strand'};
				$data_details{$ORF_id}{'ORF_type'} =  $data_details{$ORFI}{'ORFI_type'};
			}
		}
		elsif($data_details{$ORFI}{'strand'} eq "+"){
			if($data_details{$ORFI}{'end'} ~~ @tx_3ps ){
                                my $ORF_id = join("_",$gene,$data_details{$ORFI}{'end'});
                                $data_store_groups{$gene}{$ORF_id}{$ORFI} = 1;
                                if($data_details{$ORF_id}{'start'} > $data_details{$ORFI}{'start'}){
                                        $data_details{$ORF_id}{'start'}=$data_details{$ORFI}{'start'}
                                }
                        }
                        else{
                                push @tx_3ps, $data_details{$ORFI}{'end'};
                                my $ORF_id = join("_",$gene,$data_details{$ORFI}{'end'});
                                $data_store_groups{$gene}{$ORF_id}{$ORFI} = 1;
                                $data_details{$ORF_id}{'chr'} = $data_details{$ORFI}{'chr'};
                                $data_details{$ORF_id}{'start'} = $data_details{$ORFI}{'start'};
                                $data_details{$ORF_id}{'end'} = $data_details{$ORFI}{'end'};
                                $data_details{$ORF_id}{'source'} = $data_details{$ORFI}{'source'};
                                $data_details{$ORF_id}{'strand'} = $data_details{$ORFI}{'strand'};
				$data_details{$ORF_id}{'ORF_type'} =  $data_details{$ORFI}{'ORFI_type'};
                        }		
		}
	}
}

print Dumper \%data_store_groups;
print Dumper \%data_details;


#my $outfile2 = 'mismatch_ENS_Fantom.txt';
#open(my $mism, '>', $outfile2) or die "Could not open file '$outfile2' $!";
#my $outfile3 = 'duplicate_tx.txt';
#open(my $out3, '>', $outfile3) or die "Could not open file '$outfile3' $!";
#print Dumper \%data_store;
#print Dumper \%data_details;
#close $mism;

foreach my $gene (keys %data_store_groups){
	print $out $data_details{$gene}{'chr'}."\t".$data_details{$gene}{'source'}."\tgene\t".$data_details{$gene}{'start'}."\t".$data_details{$gene}{'end'}."\t.\t".$data_details{$gene}{'strand'}."\t.\tgene_id \"".$gene."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"\n"; 
      	foreach my $ORF (keys %{$data_store_groups{$gene}}){
		print $out $data_details{$ORF}{'chr'}."\t".$data_details{$ORF}{'source'}."\tORF\t".$data_details{$ORF}{'start'}."\t".$data_details{$ORF}{'end'}."\t.\t".$data_details{$ORF}{'strand'}."\t.\tgene_id \"".$gene."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"; ORF_id \"".$ORF."\"; ORF_type \"".$data_details{$ORF}{'ORF_type'}."\"\n";
		foreach my $ORFI (keys %{$data_store_groups{$gene}{$ORF}}){
			print $out $data_details{$ORFI}{'chr'}."\t".$data_details{$ORFI}{'source'}."\tiORF\t".$data_details{$ORFI}{'start'}."\t".$data_details{$ORFI}{'end'}."\t.\t".$data_details{$ORFI}{'strand'}."\t.\tgene_id \"".$gene."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"; ORF_id \"".$ORF."\"; iORF_id \"".$ORFI."\"; iORF_type \"".$data_details{$ORFI}{'ORFI_type'}."\"; iORF_pept \"".$data_details{$ORFI}{'ORFI_pept'}."\"\n";	
			foreach my $exon (keys %{$data_store{$gene}{$ORFI}}){
				 print $out $data_details{$exon}{'chr'}."\t".$data_details{$exon}{'source'}."\torfCDS\t".$data_details{$exon}{'start'}."\t".$data_details{$exon}{'end'}."\t.\t".$data_details{$exon}{'strand'}."\t.\tgene_id \"".$gene."\"; gene_biotype \"".$data_details{$gene}{'gene_type'}."\"; gene_name \"".$data_details{$gene}{'gene_name'}."\"; ORF_id \"".$ORF."\"; iORF_id \"".$ORFI."\"; iORF_type \"".$data_details{$exon}{'ORFI_type'}."\"\n";

			}
		}
        }
}

close $out;
