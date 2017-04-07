#!/usr/bin/perl

# information are stored in a hash table:
# 	key: gene name
# 	value: array where Y is present and "-" is absent
# 		array index:
# 			0: cancer_gene_census_06_26_2015.tsv
# 			1: cell_drug_FP_set.txt
# 			2: human_disease_textmining_filtered.csv
# 			3: NCG4_05-04-15_falsepositives_list.txt 
# 			4: NCG_CNV_gene.txt
# 			5: uniprot-keyword_oncogen.gene
# 			6: vogelstein_cnv_list.csv
# 			7: vogelstein_snv_list.csv


my ( $file_out, $file_in, %hash_matrix, $dir_root, @line, $geneName, $index_array );

$dir_root = "/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/GOLD_STANDARD";
$file_out = "$dir_root/gold_standard_matrix.xls";

%hash_matrix = ();	# initialize hash
open( OUT, ">$file_out" );	# initialize file out handler
open( UMB, ">$dir_root/ambiguous_genes_list.txt" );


# 0: cancer_gene_census_06_26_2015.tsv
$index_array = 0;
$file_in = "$dir_root/cancer_gene_census_06_26_2015.tsv";
open( IN, $file_in );
<IN>;	# skip header
while( <IN> ){
	chomp( @line = split( /\t/, $_ ) );
	$geneName = $line[0];
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 1: cell_drug_FP_set.txt
$index_array = 1;
$file_in = "$dir_root/cell_drug_FP_set.txt";
open( IN, $file_in );
<IN>;	# skip header
while( <IN> ){
	$geneName = $_;
	chomp( $geneName );
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 2: human_disease_textmining_filtered.csv
$index_array = 2;
$file_in = "$dir_root/human_disease_textmining_filtered.csv";
open( IN, $file_in );
<IN>;	# skip header
while( <IN> ){
	chomp( @line = split( /\t/, $_ ) );
	$geneName = $line[0];
	next unless( $line[3] eq "Cancer" && $line[5] >= 2.5 );
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 3: NCG4_05-04-15_falsepositives_list.txt
$index_array = 3;
$file_in = "$dir_root/NCG4_05-04-15_falsepositives_list.txt";
open( IN, $file_in );
<IN>;	# skip header
while( <IN> ){
	$geneName = $_;
	chomp( $geneName );
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 4: NCG_CNV_gene.txt
$index_array = 4;
$file_in = "$dir_root/NCG_CNV_gene.txt";
open( IN, $file_in );
<IN>;	# skip header
while( <IN> ){
	chomp( @line = split( /\t/, $_ ) );
	$geneName = $line[1];
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 5: uniprot-keyword_oncogen.gene
$index_array = 5;
$file_in = "$dir_root/uniprot-keyword_oncogen.gene";
open( IN, $file_in );
<IN>;	# skip header
while( <IN> ){
	$geneName = $_;
	chomp( $geneName );
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 6: vogelstein_cnv_list.csv
$index_array = 6;
$file_in = "$dir_root/vogelstein_cnv_list.csv";
open( IN, $file_in );
<IN>;	# skip header
<IN>;	# skip header
while( <IN> ){
	next if( /^*/ );
	chomp( @line = split( /\t/, $_ ) );
	$geneName = $line[0];
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );

# 7: vogelstein_snv_list.csv
$index_array = 7;
$file_in = "$dir_root/vogelstein_snv_list.csv";
open( IN, $file_in );
<IN>;	# skip header
<IN>;	# skip header
while( <IN> ){
	next if( /^*/ );
	chomp( @line = split( /\t/, $_ ) );
	$geneName = $line[0];
	unless( exists $hash_matrix{$geneName} ){
		my @array = ( "-" ) x 8;
		$hash_matrix{$geneName} = \@array;
	}
	$hash_matrix{$geneName}[$index_array] = "Y";
}
close( IN );


## generate matrix
print OUT "GeneName\tCGC_20150626\tcell_drug_FP_set\thuman_disease_textmining\tNCG4_20150504_false_positives\tNCG_CNV\tUniProt-oncogen\tvogelstein_cnv\tvogelstein_snv\n";
foreach $geneName ( sort keys %hash_matrix ){
	unless( ( grep( /Y/, @{$hash_matrix{$geneName}} ) >= 2 && $hash_matrix{$geneName}[1] eq "Y" && $hash_matrix{$geneName}[3] eq "-" ) # gene is in other list and [1] only
	|| ( grep( /Y/, @{$hash_matrix{$geneName}} ) >= 2 && $hash_matrix{$geneName}[1] eq "-" && $hash_matrix{$geneName}[3] eq "Y" ) # gene is in other list and [3] only
	|| ( grep( /Y/, @{$hash_matrix{$geneName}} ) >= 3 && $hash_matrix{$geneName}[1] eq "Y" && $hash_matrix{$geneName}[3] eq "Y" ) # gene is in other list and both [1] & [3] 
	){
		print OUT $geneName . "\t" . join( "\t", @{$hash_matrix{$geneName}} ) . "\n";
	} else{
		print UMB $geneName . "\n";
	}
}
close( OUT );
close( UMB );
