#!/usr/bin/perl
use warnings;

#require '/mnt/projects/bertrandd/oncoimpact/SCRIPT/MUTATION_BENCHMARK/ANALYSE_RESULT/common_functions_edit.pl';

my ($consolidated_result_dir, $data_gene_annotation_file, $out_dir, $gene_status_selection, $additional_method_file, $use_known_method, $append_result, $flag_plot, $script_dir) = @ARGV;

 require "$script_dir/common_functions_edit.pl";

our (@RANK_THRESHOLD);

my (%pan_cancer, %cancer_census);
print STDERR $gene_status_selection."\n";
cancer_annotation(\%cancer_census, $gene_status_selection, $script_dir);

my %sample_list = ();
my %data_gene_annotation = ();
gene_annotation(\%data_gene_annotation, \%sample_list, $data_gene_annotation_file, \%cancer_census);


#Construt the result data structures
#Construt the result data structures
my %method_ID = ();
my @ID_to_method = ();
my %method_result = ();
read_method_result(\%method_result, \%method_ID, \@ID_to_method, \%data_gene_annotation, $consolidated_result_dir, $additional_method_file, $use_known_method, $script_dir);


#Init the method_driver_sample structure
my %method_sample_driver;
#foreach $method (("BASE_LINE", @ID_to_method)){
my %excluded_method = ();
foreach $method (@ID_to_method){
    if(! -e "$consolidated_result_dir/$method.result"){
	$excluded_method{$method} = 1;
    }
    print STDERR " *** Method used $method\n";
    my %map = ();
    $method_sample_driver{$method} = \%map;
    foreach $sample (keys %sample_list){
	my %map = ();
	$method_sample_driver{$method}->{$sample} = \%map;
    }
}

#Compute method_driver_sample structure based on gene

my $cc_rank;
#For the CC sample coverage of the base line
#$cc_rank = 1;
#foreach my $gene (sort {$data_gene_annotation{$b}->{"MUT_FREQ"} <=> $data_gene_annotation{$a}->{"MUT_FREQ"}} keys %data_gene_annotation){
#last if($rank_threshold ne "ALL" && $cc_rank > $rank_threshold);
#$freq = $data_gene_annotation{$gene}->{"MUT_FREQ"};
#print STDERR " *** ".$gene."\t".$freq."\n";<STDIN>;
#if(exists $cancer_census{$gene}){
#	$gene_sample_list = $data_gene_annotation{$gene}->{"SAMPLE"};

#print STDERR $gene."\t".$data_gene_annotation{$gene}->{"MUT_FREQ"}."\t".$cc_rank."\n" if($rank_threshold ne "ALL" && $rank_threshold == 10);

#	foreach $sample (@{$gene_sample_list}){
#	    $method_sample_driver{"BASE_LINE"}->{$sample}->{$gene} = 1;
#	}
#	$cc_rank++;
#    }
#}

#Writing 
if(0){
    my @RANK_VAL = (@RANK_THRESHOLD , "ALL");
    my $cc_rank;
    for(my $cmp_rank = 0; $cmp_rank < @RANK_VAL; $cmp_rank++){
	$rank_threshold = $RANK_VAL[$cmp_rank];
	#For the CC sample coverage of the base line
	$cc_rank = 1;
	foreach my $gene (sort {$data_gene_annotation{$b}->{"MUT_FREQ"} <=> $data_gene_annotation{$a}->{"MUT_FREQ"}} keys %data_gene_annotation){
	    last if($rank_threshold ne "ALL" && $cc_rank > $rank_threshold);
	    #$freq = $data_gene_annotation{$gene}->{"MUT_FREQ"};
	    #print STDERR " *** ".$gene."\t".$freq."\n";<STDIN>;
	    if(exists $cancer_census{$gene}){
		$gene_sample_list = $data_gene_annotation{$gene}->{"SAMPLE"};
		
		#print STDERR $gene."\t".$data_gene_annotation{$gene}->{"MUT_FREQ"}."\t".$cc_rank."\n" if($rank_threshold ne "ALL" && $rank_threshold == 10);
		
		#foreach $sample (@{$gene_sample_list}){
		#$method_sample_driver{"BASE_LINE"}->{$sample}->{$gene} = 1;
		#}
		$cc_rank++;
	    }
	}
	
	
	#For the method
	foreach $method (@ID_to_method){
	    foreach $gene (keys %{$method_result{$method}}){
		$rank = $method_result{$method}->{$gene}->{"RANK"};
		next if($rank_threshold ne "ALL" && $rank >  $rank_threshold);
		$gene_sample_list = $method_result{$method}->{$gene}->{"SAMPLE"};
		foreach $sample (@{$gene_sample_list}){
		    $method_sample_driver{$method}->{$sample}->{$gene} = 1;
		}
	    }
	}
	write_barplot_file(\%method_sample_driver, $rank_threshold, "ALL");
	write_PPV_file(\%method_sample_driver, \%data_gene_annotation, $rank_threshold, "ALL");
    }
}
	
#For the method
my $flag_use_sample_score = 1;
foreach $method (@ID_to_method){
    
    #print STDERR " **** $method\n";<STDIN>;

    foreach $gene (keys %{$method_result{$method}}){
	$gene_rank = $method_result{$method}->{$gene}->{"RANK"};
	$gene_sample_list = $method_result{$method}->{$gene}->{"SAMPLE"};
	$gene_sample_score_list = $method_result{$method}->{$gene}->{"SCORE"};
	$flag_use_sample_score = 1;
	if($gene_sample_score_list eq "-"){
	    $flag_use_sample_score = 0;
	}
	
	#print STDERR " SIZES : ".(@{$gene_sample_list})." = "."(@{$gene_sample_score_list})"."\n";

	for( my $i = 0; $i < @{$gene_sample_list}; $i++){
	    $sample = $gene_sample_list->[$i];
	    my @info_score = ("0", $gene_rank);
	    $info_score[0] = $gene_sample_score_list->[$i] if($flag_use_sample_score);
	    $info_score[0] = -1 * $info_score[0] if($method eq "fathmm");#Change the sign the score to have all the sample specific score in the higher the better type
	    $method_sample_driver{$method}->{$sample}->{$gene} = \@info_score;
	}
    }
}
write_barplot_file(\%method_sample_driver, "ALL", "ALL");
write_PPV_file(\%method_sample_driver, \%data_gene_annotation, 3, "ALL");
write_PPV_file(\%method_sample_driver, \%data_gene_annotation, 5, "ALL");
write_PPV_file(\%method_sample_driver, \%data_gene_annotation, 10, "ALL");


#For the sample results
if(0){    
    my $rank;
    my $method;
    my %sample_method_res = ();
    my @tmp = split(/\//, $consolidated_result_dir);
    my $cancer_type = $tmp[-1];

    foreach $sample (keys %sample_list){
	`mkdir -p "$consolidated_result_dir/../../SAMPLE_RESULT/$cancer_type/$sample"`;
	for(my $i = 0; $i < @ID_to_method; $i++){
	    $method = $ID_to_method[$i];
	    $sample_method_res = $method_sample_driver{$method}->{$sample};
	    $rank = 1;
	    
	    my $filename = "$consolidated_result_dir/../../SAMPLE_RESULT/$cancer_type/$sample/${method}.result";
	    open(my $fh, '>', $filename) or die "Could not open file $filename $!";

	    print $fh "Gene_name\tSample\tRank\tpValue\tInfo\n";
	    foreach $gene (sort by_score_and_rank (keys %{$sample_method_res})){
		print $fh "$gene\t-\t$rank\t-\t-\n";	
		$rank++;
	    }
	    
	}
    }
}


sub write_PPV_file{
    my ($method_sample_driver, $data_gene_annotation, $threshold_val, $gene_status_selection) = @_;
    my $PPV_file = "$out_dir/sample_precision_RANK_$threshold_val.dat";
    open(OUT, ">$PPV_file");open(OUT_F1, ">$PPV_file\_penalised");open(OUT_TPR, ">$out_dir/sample_recall_RANK_$threshold_val.dat");
    #The header
    print OUT "".join("\t", @ID_to_method)."\n";print OUT_TPR "".join("\t", @ID_to_method)."\n";print OUT_F1 "".join("\t", @ID_to_method)."\n";
    my ($PPV, $F1, $cmp_pred);
    my %sample_method_res = ();
    foreach $sample (keys %sample_list){
	print OUT $sample;print OUT_TPR $sample;print OUT_F1 $sample;
	for(my $i = 0; $i < @ID_to_method; $i++){
	    $method = $ID_to_method[$i];

	    #print STDERR " $sample $method\n";<STDIN>;
	    $nb_driver = 0;
	    #$nb_driver = keys(%{$method_sample_driver->{$method}->{$sample}});
	    #$nb_driver = $threshold_val if($nb_driver > $threshold_val); 
	    $nb_cc_driver = 0;
	    #Sort the prediction according to sample specifgic score
	    #In case of similar score or no sample specific score, the rank of the predicted gene in the agglomerative prediction is used
	    if(! exists $excluded_method{$method}){
		$sample_method_res = $method_sample_driver->{$method}->{$sample};
		$cmp_pred = 0;
		#print STDERR " *** ".join("\t", (keys %{$sample_method_res}));<STDIN>;
		foreach $gene (sort by_score_and_rank (keys %{$sample_method_res})){
		    if(! defined $data_gene_annotation->{$gene}->{"STATUS"}){
			#print STDERR " *** WARNING WEIRD gene in list:\t$method\t$sample\t $gene\n";
			next;
		    }
		    #print STDERR "\t $gene\n";
		    $nb_driver++;
		    $nb_cc_driver++ if(exists $cancer_census{$gene});
		    last if($nb_driver == $threshold_val);
		}
		$PPV = 0;#-0.1;
		$PPV = $nb_cc_driver/$nb_driver if($nb_driver != 0);
		#print OUT "\t".$PPV;
		
		#To penalise sample that did not output 5 driver per samples on their top50 driver list
		if($PPV < 0){
		    $F1 = $PPV;
		}
		else{
		    $F1 = 0;
		    $TPR = $nb_cc_driver/$threshold_val;#$TPR = 1 if($TPR > 1);
		    $F1 = 2 * ($PPV * $TPR / ($PPV + $TPR)) if($PPV != 0 && $TPR != 0);
		}
		
		if($method eq "DawnRank" && $nb_driver == 0){
		    $F1 = "NA";
		    $PPV = "NA";
		    $TPR = "NA";
		}
	    }
	    else{
		$F1 = "NA";
		$PPV = "NA";
		$TPR = "NA";
	    }
	    print OUT "\t".$PPV;
	    print OUT_TPR "\t".$TPR;
	    print OUT_F1 "\t".$F1;
	}
	print OUT "\n";	print OUT_TPR "\n";print OUT_F1 "\n";
    }
    close(OUT);close(OUT_F1);
}


sub by_score_and_rank {
    #print STDERR " ***by_score_and_rank\n----".$a."\n----".$b."\n";<STDIN>;
    #print STDERR " ***by_score_and_rank\n----".(join("\t", (@{$sample_method_res->{$a}})))."\n----".(join("\t", (@{$sample_method_res->{$b}})))."\n";<STDIN>;
    ($sample_method_res->{$b}->[0] <=> $sample_method_res->{$a}->[0] ||#sample spcefic gene score
     $sample_method_res->{$a}->[1] <=> $sample_method_res->{$b}->[1] #Gene rank 
    );
}

sub write_barplot_file{
     my ($method_sample_driver, $threshold_val, $gene_status_selection) = @_;

     @NUM_DRIVER_CAT = (0,1,3,8,15,25,100000);
     
     #Organise the data
     my %data_barplot;
     #for(my $i = -1; $i < @ID_to_method; $i++){
     for(my $i = 0; $i < @ID_to_method; $i++){
	 #Init the value 
	 $method = ($i == -1) ? "BASE_LINE" : $ID_to_method[$i];
	 
	 my %method_cat = ();foreach $v (@NUM_DRIVER_CAT){$method_cat{$v} = 0};
	 foreach $sample (keys %{$method_sample_driver->{$method}}){
	     
	     #Should not happen as the maf file should be filtered for the samples excluded
	     if(! defined $sample_list{$sample}){
		 #print STDERR " *** WARNING SAMPLE $sample NOT PRESENT IN THE SAMPLE LIST FOR METHOD $method\n";
		 next;
	     }

	     $nb_driver = keys(%{$method_sample_driver->{$method}->{$sample}});
	     for (my $i = 0; $i < @NUM_DRIVER_CAT; $i++){
		 $cat_val = $NUM_DRIVER_CAT[$i];
		 if($nb_driver <= $cat_val){
		     #print STDERR " *** $method $cat_val $sample\n" if($cat_val == 0);
		     $method_cat{$cat_val}++;
		     last;
		 }
	     }
	 }
	 $data_barplot{$method} = \%method_cat;
     }
     
     #Just to check if everything is fine
     my @sum = (); for(my $i = -1; $i < @ID_to_method; $i++){push(@sum, 0);}
     
     #Write the file
     #Write the header
     my $barplot_file = "$out_dir/sample_nb_driver_cat_RANK_$threshold_val.dat";
     open(OUT, ">$barplot_file");
     
     #print OUT "BASE_LINE\t".join("\t", @ID_to_method)."\n";
     print OUT join("\t", @ID_to_method)."\n";
     my $previous_val = "";
     my $current_val = "";
     for(my $j = 0; $j < @NUM_DRIVER_CAT; $j++){
	 $current_val = $NUM_DRIVER_CAT[$j];
	 if($j != 0 && $previous_val+1 != $current_val && $j != @NUM_DRIVER_CAT-1){
	     print OUT ($previous_val+1)."-".($current_val);
	 }
	 print OUT $current_val if($j == 0 || $previous_val+1 == $current_val);
	 print OUT ">".($previous_val+1) if($j == @NUM_DRIVER_CAT-1);

	 #for(my $i = -1; $i < @ID_to_method; $i++){
	 for(my $i = 0; $i < @ID_to_method; $i++){
	     $method = ($i == -1) ? "BASE_LINE" : $ID_to_method[$i];
	     $nb_sample_in_cat = $data_barplot{$method}->{$current_val};
	     if($excluded_method{$method}){
		 print OUT "\t"."NA";
	     }
	     else{
		 print OUT "\t".$nb_sample_in_cat;
	     }
	     $sum[$i+1] += $nb_sample_in_cat;
	 }
	 print OUT "\n";
	 
	 $previous_val = $current_val;
     }
     
     #print OUT "sum\t".(join("\t", @sum))."\n";
     
     close(OUT);
    
}

sub write_boxplot_file{
     my ($threshold_type, $threshold_val, $gene_status_selection) = @_;
     my $boxplot_file = "$out_dir/sample_cov_$threshold_type\_$threshold_val";
     if($gene_status_selection ne ""){
	$boxplot_file .= "\_$gene_status_selection";
    }
     
     open(OUT, ">$boxplot_file.dat");

     #The header
     print OUT "".(join("\t", @ID_to_method))."\n";
     
     foreach $sample (keys %sample_list){
	 print OUT $sample;
	 for(my $i = 0; $i < @ID_to_method; $i++){
	     $method = $ID_to_method[$i];
	     $res = 0;
	     if(exists $method_sample_driver{$method}->{$sample}){
		 $res = $method_sample_driver{$method}->{$sample};
	     }
	     print OUT "\t".$res;
	 }
	 print OUT "\n";
     }
     
     close(OUT);

     return $boxplot_file;

}

sub plot_boxplot{
    my ($matrix_file, $title) = @_;
    
    my $font_size = 3;
    my $note_font_size = 3;
    my $margin_size = 30; 

    open(OUT, ">$matrix_file.R");
    print OUT "pdf(file=\"$matrix_file.pdf\",
	paper=\"special\",
	width=25,
	height=25
	)\n";
    print OUT "profile <- read.table(\"$matrix_file.dat\", header = TRUE)\n";

    print OUT "palette <- rainbow(ncol(profile))\n";
    print OUT "boxplot.matrix(as.matrix(profile), col=palette, cex.axis=$font_size)\n";

    close(OUT);
    run_exe("R-3.0.0 --vanilla < $matrix_file.R");
}
