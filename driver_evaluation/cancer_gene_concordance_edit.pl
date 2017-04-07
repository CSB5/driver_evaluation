#!/usr/bin/perl
use warnings;

#require '/mnt/projects/bertrandd/oncoimpact/SCRIPT/MUTATION_BENCHMARK/ANALYSE_RESULT/common_functions_edit.pl';

my ($consolidated_result_dir, $data_gene_annotation_file, $gene_status_selection, $out_dir, $flag_plot, $additional_method_file, $use_known_method, $append_result, $script_dir, $file_ID) = @ARGV;

 require "$script_dir/common_functions_edit.pl";
our (@RANK_THRESHOLD);

#
#To update using that file:
#/mnt/pnsg10_projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/TEST_DATA/COAD/gene_mutation_frequency.txt
#

#Get the list of cancer genes
my %cancer_census;

cancer_annotation(\%cancer_census, $gene_status_selection, $script_dir);
#cancer_relaxed_annotation(\%cancer_census);

#Get the info from data_gene_annotation_file
my %sample_list = ();
my %data_gene_annotation = ();
gene_annotation(\%data_gene_annotation, \%sample_list, $data_gene_annotation_file, \%cancer_census);

#Construt the result data structures
my %method_ID = ();
my @ID_to_method = ();
my %method_result = ();
read_method_result(\%method_result, \%method_ID, \@ID_to_method, \%data_gene_annotation, $consolidated_result_dir, $additional_method_file, $use_known_method, $script_dir);

#print STDERR " *** THE ID TO METHOD @ID_to_method\n";<STDIN>;


#Read the prediction to compute the number of CC calls at different rank
my %method_concordance = ();
my @RANK_VAL = (@RANK_THRESHOLD, "ALL");

foreach $method (@ID_to_method){
    my @tab_concord = ();
    for($i = 0; $i < @RANK_VAL; $i++){
	$tab_concord[$i] = 0;
    }
    
    foreach $gene (keys %{$method_result{$method}}){
	$rank = $method_result{$method}->{$gene}->{"RANK"};
	if(exists $cancer_census{$gene}){
	    for($i = 0; $i < @RANK_VAL; $i++){
		$r = $RANK_VAL[$i];
		if($r eq "ALL" || $rank <= $r){
		    $tab_concord[$i]++;
		}
	    }
	}
    }
    
    $method_concordance{$method} = \@tab_concord;
}

#For the base line PPV values

my @base_line_concordance = ();for($i = 0; $i < @RANK_VAL; $i++){$base_line_concordance[$i] = 0;}

#Sort by mutation frequency
my $rank = 1;
foreach my $gene (sort {$data_gene_annotation{$b}->{"MUT_FREQ"} <=> $data_gene_annotation{$a}->{"MUT_FREQ"}} keys %data_gene_annotation){
    #$freq = $data_gene_annotation{$gene}->{"MUT_FREQ"};
    #print STDERR " *** ".$gene."\t".$freq."\n";<STDIN>;
    if(exists $cancer_census{$gene}){
	for($i = 0; $i < @RANK_VAL; $i++){
	    $r = $RANK_VAL[$i];
	    if($r eq "ALL" || $rank <= $r){
		$base_line_concordance[$i]++;
	    }
	}
    }
    $rank++;
}

my $nb_cancer_gene_mutated = $base_line_concordance[-1];
my $nb_gene_mutated = keys (%data_gene_annotation) + 0;
print STDERR " *** @base_line_concordance\n";

#For the sensitivity
#write_sensitivity_file(\%method_concordance, $nb_cancer_gene, "$out_dir/cancer_gene_pred_sens.dat");
#write_PPV_file(\%method_concordance, \@base_line_concordance, \%method_result, "$out_dir/cancer_gene_PPV.dat");

my $out_file = "$out_dir/cancer_gene_$gene_status_selection";
$out_file = "$out_dir/cancer_gene_$gene_status_selection\_$file_ID" if(defined($file_ID));
write_measure_file(\%method_concordance, \@base_line_concordance, \%method_result, $nb_cancer_gene_mutated, $nb_gene_mutated, $out_file);
write_rank_concordance_file(\%method_result, \%data_gene_annotation, 50);

sub write_rank_concordance_file{
    my ($method_result, $data_gene_annotation, $max_rank) = @_;
    my $nb_cc;my $rank;
    for(my $i = 0; $i < @ID_to_method; $i++){
	$method = $ID_to_method[$i];
	#print STDERR " *** write_rank_concordance_file $method\n";
	open(OUT, ">$out_dir/$method\_cc_precision_RANK_$max_rank.dat");
	#Sort by rank
	$nb_cc = 0;
	foreach $gene (sort {$method_result->{$method}->{$a}->{"RANK"} <=> $method_result->{$method}->{$b}->{"RANK"}} keys %{$method_result->{$method}}){
	    $rank = $method_result->{$method}->{$gene}->{"RANK"};
	    last if($rank > $max_rank);
	    $nb_cc++ if(exists $cancer_census{$gene});
	    print OUT $rank."\t".($nb_cc/$rank)."\n";
	}
	close(OUT);
    }

    #for the base line
    $method = "base_line";
    my $out_file = "$out_dir/$method\_cc_precision_RANK_$max_rank.dat";
    $out_file = "$out_dir/$method\_cc_precision_RANK_$max_rank\_$file_ID.dat" if(defined $file_ID);
    open(OUT, ">$out_dir/$method\_cc_precision_RANK_$max_rank.dat");
    $rank = 1;
    $nb_cc = 0;
    foreach my $gene (sort {$data_gene_annotation->{$b}->{"MUT_FREQ"} <=> $data_gene_annotation->{$a}->{"MUT_FREQ"}} keys %{$data_gene_annotation}){
	#$freq = $data_gene_annotation{$gene}->{"MUT_FREQ"};
	#print STDERR " *** ".$gene."\t".$freq."\n";<STDIN>;
	$nb_cc++ if(exists $cancer_census{$gene});
	print OUT $rank."\t".($nb_cc/$rank)."\n";
	$rank++;
    }
    close(OUT);
}



sub write_measure_file{
    my ($method_concordance, $base_line_concordance, $method_result, $nb_cancer_gene_in_data, $nb_gene_in_data, $out_file) = @_;

    my ($nb_driver_called, $nb_cancer_driver_called, $nb_cancer_gene);

    print STDERR " *** ".$out_file."\n";
    
    $write_type = ">"; $write_type = ">>" if($append_result eq "APPEND_RESULT");

    open(OUT_PPV, $write_type, "$out_file\_precision.dat");
    open(OUT_TPR, $write_type, "$out_file\_recall.dat");
    open(OUT_F1, $write_type, "$out_file\_F1.dat");

    #The header
    if(1 || $append_result ne "APPEND_RESULT"){
	print OUT_F1 "".(join("\t", @RANK_VAL))."\n";
	print OUT_TPR "".(join("\t", @RANK_VAL))."\n";
	print OUT_PPV "".(join("\t", @RANK_VAL))."\n";
    }
    
    for(my $i = -1; $i < @ID_to_method; $i++){
	#For the base line
	if($i == -1){
	    $method = "base_line";
	    $nb_method_driver = $nb_gene_in_data;
	}
	else{
	    $method = $ID_to_method[$i];
	    $nb_method_driver = keys (%{$method_result->{$method}});
	}

	print OUT_F1 $method;print OUT_PPV $method;	print OUT_TPR $method;
	for(my $j = 0; $j < @RANK_VAL; $j++){
	    
	    $nb_driver_called = $RANK_VAL[$j];
	    $nb_driver_called = $nb_method_driver if($nb_driver_called eq "ALL" || $nb_method_driver < $nb_driver_called);
	    
	    $nb_cancer_driver_called = $method_concordance->{$method}->[$j];
	    $nb_cancer_driver_called = $base_line_concordance->[$j] if($i == -1);#For the base line
	    
	    $nb_cancer_gene = $RANK_VAL[$j];
	    $nb_cancer_gene = $nb_cancer_gene_in_data if($nb_cancer_gene eq "ALL" || $nb_cancer_gene_in_data < $nb_cancer_gene);

	    $PPV= 0; $F1 = 0;	$TPR = 0;

	    $TPR = $nb_cancer_driver_called / $nb_cancer_gene if($nb_cancer_gene != 0);
	    $PPV = $nb_cancer_driver_called / $nb_driver_called if($nb_driver_called != 0);
	    $F1 = 2 * ($PPV * $TPR / ($PPV + $TPR)) if($PPV != 0 && $TPR != 0);

	    #To avoid the problem of missing data sets for the meta methods
	    if(($method eq "DriverNet" || $method eq "DawnRank" || $method eq "DriverDB" || $method eq "MutSig_meta"  || $method eq "IntOGen"
		|| $method eq "borda" || $method eq "borda_all" || $method eq "RRA" || $method eq "RRA_all"  || $method eq "random_forest" || $method eq "SVM" || $method eq "SVM_RBF" #This line should be removed after running the method on all samples
	       ) && $nb_driver_called == 0){
		$TPR = "NA";$PPV= "NA"; $F1 = "NA";
	    }

	    print OUT_F1 "\t".$F1;
	    print OUT_PPV "\t".$PPV;
	    print OUT_TPR "\t".$TPR;
	}
	print OUT_F1 "\n";print OUT_PPV "\n";	print OUT_TPR "\n";
    }
    close(OUT_PPV);close(OUT_F1);	close(OUT_TPR);
}

