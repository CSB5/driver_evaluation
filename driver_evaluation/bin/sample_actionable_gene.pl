#!/usr/bin/perl
use warnings;

#require '/mnt/projects/bertrandd/oncoimpact/SCRIPT/MUTATION_BENCHMARK/ANALYSE_RESULT/common_functions.pl';

my ($consolidated_result_dir, $data_gene_annotation_file, $out_dir, $additional_method_file, $use_known_method, $append_result, $flag_plot, $script_dir) = @ARGV;

require "$script_dir/common_functions.pl";

our (@RANK_THRESHOLD);
#my (%pan_cancer, %cancer_census);
#print STDERR $gene_status_selection."\n";
my %cancer_census = ();
cancer_annotation(\%cancer_census, "CANCER_UNION", $script_dir);

my %actionable_gene = ();
open(FILE, "$script_dir/../ACTIONABLE_GENES/combine_target.dat");
while(<FILE>){
    @line = split(/\t/, $_);
    $gene = $line[0];
    $validation_type = $line[1];
    $actionable_gene{$gene} = $validation_type;
}
close(FILE);
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
foreach $method (@ID_to_method){
    my %map = ();
    $method_sample_driver{$method} = \%map;
    foreach $sample (keys %sample_list){
	my %map = ();
	$method_sample_driver{$method}->{$sample} = \%map;
    }
}

#Compute method_driver_sample structure based on gene

my $cc_rank;

	
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
#write_barplot_file(\%method_sample_driver, "ALL", "ALL");
#write_actionable_file(\%method_sample_driver, \%data_gene_annotation, 5, 1);
#write_actionable_file(\%method_sample_driver, \%data_gene_annotation, 5, 2);
#write_actionable_file(\%method_sample_driver, \%data_gene_annotation, 5, 3);
#
#write_actionable_file(\%method_sample_driver, \%data_gene_annotation, 10, 1);
#write_actionable_file(\%method_sample_driver, \%data_gene_annotation, 10, 2);
#write_actionable_file(\%method_sample_driver, \%data_gene_annotation, 10, 3);    
#
write_barplot_actionable_file(\%method_sample_driver, \%data_gene_annotation, 10);
write_barplot_actionable_file(\%method_sample_driver, \%data_gene_annotation, 5);


sub write_actionable_file{
    my ($method_sample_driver, $data_gene_annotation, $threshold_val, $actionable_type) = @_;
    my $res_file = "$out_dir/sample_actionable_gene_$actionable_type\_RANK_$threshold_val.dat";
    open(OUT, ">$res_file");
    #The header
    print OUT "".join("\t", @ID_to_method)."\n";
    my @actionable_gene = ();
    my %sample_method_res = ();
    foreach $sample (keys %sample_list){
	print OUT $sample;
	for(my $i = 0; $i < @ID_to_method; $i++){
	    $method = $ID_to_method[$i];

	    #print STDERR " $sample $method\n";<STDIN>;
	    $nb_driver = 0;
	    #$nb_driver = keys(%{$method_sample_driver->{$method}->{$sample}});
	    #$nb_driver = $threshold_val if($nb_driver > $threshold_val); 
	    $nb_ac_driver = 0;
	    #Sort the prediction according to sample specifgic score
	    #In case of similar score or no sample specific score, the rank of the predicted gene in the agglomerative prediction is used
	    $sample_method_res = $method_sample_driver->{$method}->{$sample};
	    #print STDERR " *** ".join("\t", (keys %{$sample_method_res}));<STDIN>;
	    foreach $gene (sort by_score_and_rank (keys %{$sample_method_res})){
		if(! defined $data_gene_annotation->{$gene}->{"STATUS"}){
		    #print STDERR " *** WARNING WEIRD gene in list:\t$method\t$sample\t $gene\n";
		    next;
		}
		#print STDERR "\t $gene\n";
		$nb_driver++;
		$nb_ac_driver++ if(exists $actionable_gene{$gene} && $actionable_gene{$gene} == $actionable_type);
		last if($nb_driver == $threshold_val);
	    }
	    print OUT "\t".$nb_ac_driver;
	    
	}
	print OUT "\n";
    }
    close(OUT);
}

sub write_barplot_actionable_file{
    my ($method_sample_driver, $data_gene_annotation, $threshold_val) = @_;
    my $profile_file = "$out_dir/sample_actionable_gene_profile_RANK_$threshold_val.dat";
    my $number_file = "$out_dir/sample_actionable_gene_number_RANK_$threshold_val.dat";
    my $number_file_cancer = "$out_dir/sample_cancer_gene_number_RANK_$threshold_val.dat";
    open(OUT_P, ">$profile_file");
    open(OUT_N, ">$number_file");
    open(OUT_N_C, ">$number_file_cancer");
    #
    my $gene_file = "$out_dir/sample_actionable_gene\_$threshold_val.dat";
    open(OUT_G, ">$gene_file");
    #The header
    #print OUT_P "".join("\t", @ID_to_method)."\n";
    #print OUT_N "".join("\t", @ID_to_method)."\n";
    my @NUM_CAT = (0,1,2,3,4,5);
    print OUT_P "1\t2\t3\t0\n";
    print OUT_N "5\t4\t3\t2\t1\t0\n";
    print OUT_N_C "5\t4\t3\t2\t1\t0\n";
    
    my @num_cat_res = ();my @num_cat_res_rev = ();
    my @num_cat_res_cancer = ();my @num_cat_res_cancer_rev = ();
    my @profile_res = ();
    
    #To compute the number of sample with at least an actionable gene mutated
    open(OUT_B, ">$out_dir/sample_actionable_gene.dat");
    my %sample_actionable = ();

    my %actionable_patient = ();
    my %actionable_patient_file = ("BRAF", ["$script_dir/../ACTIONABLE_GENES/BRAF_V600_patient.dat", "BRAFV600"],
				   "EGFR", ["$script_dir/../ACTIONABLE_GENES/EGFR_AMUT_patient.dat", "EGFR_MUT"],
				   "PIK3CA", ["$script_dir/../ACTIONABLE_GENES/PIK3CA_AMUT_patient.dat", "PIK3CA_MUT"]
	);

    foreach $g (keys %actionable_patient_file){
	open(FILE, $actionable_patient_file{$g}->[0]);
	$sample_actionable{$g} = [$actionable_patient_file{$g}->[1], {}];
	while(<FILE>){
	    chop $_;
	    @line = split(/\t/, $_);
	    $sample_actionable{$g}->[1]->{$line[1]} = 1;
	}
	close(FILE);
    }

    if($threshold_val == 5){
	foreach $g (keys %data_gene_annotation){
	    if($actionable_gene{$g}){
		foreach $s (@{$data_gene_annotation{$g}->{"SAMPLE"}}){
		    print OUT_B "ALL"."\t".$g."\t".$actionable_gene{$g}."\t".$s."\n";
		}
	    }
	}
    }
    close(OUT_B);

    for(my $i = 0; $i < @ID_to_method; $i++){
	$method = $ID_to_method[$i];
	print OUT_P $method."\t";
	print OUT_N $method."\t";
	print OUT_N_C $method."\t";
	#
	@profile_res = (0,0,0,0);
	#
	@num_cat_res = (0,0,0,0,0,0);
	@num_cat_res_cancer = (0,0,0,0,0,0);
	
	foreach $sample (keys %sample_list){
	    #print STDERR " $sample $method\n";<STDIN>;
	    $nb_driver = 0;
	    #$nb_driver = keys(%{$method_sample_driver->{$method}->{$sample}});
	    #$nb_driver = $threshold_val if($nb_driver > $threshold_val); 
	    $nb_ac_driver = 0;
	    $nb_cancer_driver = 0;
	    $best_ac_driver = 0;
	    #Sort the prediction according to sample specifgic score
	    #In case of similar score or no sample specific score, the rank of the predicted gene in the agglomerative prediction is used
	    $sample_method_res = $method_sample_driver->{$method}->{$sample};
	    #print STDERR " *** ".join("\t", (keys %{$sample_method_res}));<STDIN>;
	    foreach $gene (sort by_score_and_rank (keys %{$sample_method_res})){
		if(! defined $data_gene_annotation->{$gene}->{"STATUS"}){
		    #print STDERR " *** WARNING WEIRD gene in list:\t$method\t$sample\t $gene\n";
		    next;
		}
		#print STDERR "\t $gene\n";
		$nb_driver++;
		$nb_cancer_driver++ if(exists $cancer_census{$gene});
		if(exists $actionable_gene{$gene}){
		    $nb_ac_driver++;
		    print OUT_G $method."\t".$gene."\t".$actionable_gene{$gene}."\t".$sample."\n";
		    if(exists $sample_actionable{$gene} && exists $sample_actionable{$gene}->[1]->{$sample}){
			print OUT_G $method."\t".$sample_actionable{$gene}->[0]."\t".$actionable_gene{$gene}."\t".$sample."\n";
		    }
		}
		$best_ac_driver = $actionable_gene{$gene} if(exists $actionable_gene{$gene} && ($best_ac_driver == 0 || $best_ac_driver > $actionable_gene{$gene}));
		last if($nb_driver == $threshold_val);
	    }
	
	    #Profile
	    $profile_res[$best_ac_driver]++;
	    #
	    #Num
	    for (my $i = 0; $i < @NUM_CAT; $i++){
		if($nb_ac_driver <= $i){
		    #print STDERR " *** $method $cat_val $sample\n" if($cat_val == 0);
		    $num_cat_res[$i]++;
		    last;
		}
	    }
	    for (my $i = 0; $i < @NUM_CAT; $i++){
		if($nb_cancer_driver <= $i){
		    #print STDERR " *** $method $cat_val $sample\n" if($cat_val == 0);
		    $num_cat_res_cancer[$i]++;
		    last;
		}
	    }
	}
	push(@profile_res, $profile_res[0]);
	print OUT_P "".(join("\t", @profile_res[1..@profile_res-1]))."\n";
	#
	@num_cat_res_rev = ();
	@num_cat_res_cancer_rev = ();
	for (my $i = @num_cat_res-1; $i >= 0; $i--){
	    push(@num_cat_res_rev, $num_cat_res[$i]);
	    push(@num_cat_res_cancer_rev, $num_cat_res_cancer[$i]);
	}
	print OUT_N "".(join("\t", @num_cat_res_rev))."\n";
	print OUT_N_C "".(join("\t", @num_cat_res_cancer_rev))."\n";
    }
    close(OUT_N);
    close(OUT_N_C);
    close(OUT_P);
    close(OUT_G);
}


sub by_score_and_rank {
    #print STDERR " ***by_score_and_rank\n----".$a."\n----".$b."\n";<STDIN>;
    #print STDERR " ***by_score_and_rank\n----".(join("\t", (@{$sample_method_res->{$a}})))."\n----".(join("\t", (@{$sample_method_res->{$b}})))."\n";<STDIN>;
    ($sample_method_res->{$b}->[0] <=> $sample_method_res->{$a}->[0] ||#sample spcefic gene score
     $sample_method_res->{$a}->[1] <=> $sample_method_res->{$b}->[1] #Gene rank if the full list
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
	     print OUT "\t".$nb_sample_in_cat;
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
