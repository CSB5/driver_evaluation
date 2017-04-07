#!/usr/bin/perl
use warnings;


my @FIXED_METHOD_ORDER = (
    "MutationTaster", "SIFT", "PolyPhen2", "MutationAssessor", #Point mutation effect
    "transFIC", "fathmm", "CHASM",  #Point mutation effect Cancer related
    #"CHASM_average",  #Point mutation effect Cancer related
    "ActiveDriver_PTM", 
    #"ActiveDriver_phospho",
    "OncodriveFM", "OncodriveCLUST", "MutSigCV", #Point mutation bias
    "OncodriveCIS", "S2N", #E + C
    "NetBox", #Intergative method N + S + C 
    "HotNet2A",
    "DawnRank", 
    "DriverNet", #Intergative method N + E + S + C
    "oncoIMPACT-MEDIAN",
    "oncoIMPACT",#, "oncoIMPACT-v1"
    #Not if it should be there
    "oncoIMPACT-discovery"
    );

#Need to handle the ALL correctltly for the pairwise_comparision analysis
our (@RANK_THRESHOLD);
@RANK_THRESHOLD = (5, 10, 25, 50, 75, 100);#, "ALL");

#Construct the data base of results
sub read_method_result{
    my ($method_result, $method_ID, $ID_to_method, $data_gene_annotation, $consolidated_result_dir, $additional_method_file, $use_known_method, $script_dir) = @_;
    
    opendir(DIR, $consolidated_result_dir);
    my @all_result_file = readdir(DIR);
    close(DIR);
    
    #For the selected method for meta-analysis
    my %selected_method = ();
    my @selected_method_order = ();
    if($additional_method_file ne "NONE"){
#	@all_result_file = ();
	open(FILE, $additional_method_file);
	<FILE>;#to skip the header
	while(<FILE>){
	    chomp $_;
	    @line = split(/\s+/, $_);
	    #
	    $method = $line[0];
	    #
	    $method_file = "$method.result";
	    #
	    #For the file that are stored in the META directory
	    if(index($method,"/") != -1){
		@temp = split(/\//, $line[0]);
		$method = $temp[1];
	    }
	    
	    push(@all_result_file, $method_file);
	    #
	    $selected_method{$method_file} = $method;
	    push(@selected_method_order, $method);
	}
    }
    #print STDERR "\n@FIXED_METHOD_ORDER\n@selected_method_order\n";
    if($use_known_method){
	push(@FIXED_METHOD_ORDER, @selected_method_order);
    }
    else{
	#To only use the methods provided in the file
	@FIXED_METHOD_ORDER = @selected_method_order;
    }
    #print STDERR "\n@FIXED_METHOD_ORDER\n";<STDIN>;
    $cmp = 0;
    foreach $res_file (@all_result_file){
	if($res_file =~ m/(.*)\.result/){
	    print STDERR $cmp." -> ".$res_file."\n" if($cmp % 500 == 0);$cmp++;
	    #The method analysed
	    $method = $1;
	    #print STDERR " **** $method |$1|\n";
	    next if($method eq "oncoIMPACT-v2" || $method eq "oncoIMPACT-v1" || $method eq "oncoIMPACT-v2-p1000");

	    $method =  $selected_method{$res_file} if(defined($selected_method{$res_file}));		#if($additional_method_file ne "NONE");
	    
	    $method_ID->{$method} = 1;
	    
	    #print STDERR " *** $method\n";

	    #To get the info for the method
	    my %method_info = ();

	    $path_to_res_file = "$consolidated_result_dir/$res_file";
	    print STDERR " *** Analyse $path_to_res_file -> $method\n";#<STDIN>;
	    next if(! -e $path_to_res_file);
	    open(FILE, "$path_to_res_file");
	    <FILE>;#Skip the header
	    while(<FILE>){
		chop $_;
		@line = split(/\t/, $_);
		#@line = split(/\s+/, $_);
		$gene = $line[0];
		
		#Due to the problem with the annovar annotation
		next if($data_gene_annotation ne "NONE" && ! exists $data_gene_annotation->{$gene});
		
		$rank = $line[2];
		$gene_mutated_sample = "-";
		$gene_mutated_sample = $data_gene_annotation->{$gene}->{"SAMPLE"} if($data_gene_annotation ne "NONE");
		my %gene_info = ("RANK", $rank, "FREQ_CALL", 0, "USE", 0, "SAMPLE",  $gene_mutated_sample, "SCORE", "-");
		
		#For the methods that perform sample specific predictions
		$samples = $line[1];
		if($samples ne "ALL" && $samples ne ""  && $samples ne "NONE"){
		    my @sample_list = split(/\;/, $samples);
		    $scores = $line[5];
		    #print STDERR " ---- @line\n";<STDIN>;
		    my @score_list = split(/\;/, $scores);
		    ################
		    #To handle the problem with incoherant sample names
		    #for(my $i = 0; $i < @sample_list; $i++){
		    #print STDERR " ***  $sample_list[$i]\n";
		    #@tmp = split("-", $sample_list[$i]);
		    #$sample_list[$i] = join("-", @tmp[0..2]);
		    #print STDERR " ***  $sample_list[$i]\n";<STDIN>;
		    #}
		    ################
		    $gene_info{"SAMPLE"} = \@sample_list;
		    $gene_info{"SCORE"} = \@score_list;
		    
		}
		
		#print STDERR $gene_info{"SAMPLE"}."\t".$method."\t".$samples."\n";<STDIN>;

		$method_info{$gene} = \%gene_info;
	    }
	    close(FILE);
	    $method_result->{$method} = \%method_info;
	}
    }
    #update_method_ID_fixed($method_ID, $ID_to_method);
    
    if(@FIXED_METHOD_ORDER == 0){
	update_method_ID($method_ID, $ID_to_method);
    }
    else{
	update_method_ID_fixed($method_ID, $ID_to_method);
    }
}

sub update_method_ID_fixed{
    my ($method_ID, $ID_to_method) = @_;
    my $cmp_ID = 0;
    for(my $i = 0; $i < @FIXED_METHOD_ORDER; $i++){
	
	$method = $FIXED_METHOD_ORDER[$i];

	#print STDERR " *** method in fixed method order: |$method|\n";

	if(exists $method_ID->{$method}){
	    #print STDERR " *** method with a valid ID: $method\n";
	    $method_ID->{$method} = $cmp_ID;
	    $ID_to_method->[$cmp_ID] =  $method;
	    $cmp_ID++;
	}
    }
}

sub update_method_ID{
    my ($method_ID, $ID_to_method) = @_;
    my $cmp_ID = 0;
    foreach $method  (keys %{$method_ID}){
	$method_ID->{$method} = $cmp_ID;
	$ID_to_method->[$cmp_ID] =  $method;
	$cmp_ID++;
    }
}


sub gene_annotation{

    my ($data_gene_annotation, $sample_list, $data_gene_annotation_file, $cancer_annotation) = @_;
    

    print STDERR " *** Read gene annotation file $data_gene_annotation_file\n";

    open(FILE, $data_gene_annotation_file);

    <FILE>;#To skip the header
    my %filter_map = ();
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$gene =  $line[0];
	
	#Gene annotation
	$mut_freq = $line[1];
	
	$gene_status = "-";
	$gene_status = "POSITIVE" if(exists $cancer_annotation->{$gene});
	
	#Sample annotation
	@sample_info_tmp = split(/\;/, $line[2]);
	my @sample_info = ();
	%filter_map = ();
	my %mutation_type = ("CNA", 0, "MISSENSE", 0, "OTHER_SNV", 0);
	foreach $s_m (@sample_info_tmp){
	    @tmp = split(/\:/, $s_m);
	    #
	    $mut_type = $tmp[1];#Only CNA and SNV probably need to use the bam file to get the missense vs other SNV type
	    $mut_type = "MISSENSE" if($mut_type eq "MUT");
	    $mut_type = "CNA" if($mut_type eq "AMPL" || $mut_type eq "DEL");
	    $mutation_type{$mut_type}++;
	    #
	    $sample = $tmp[0];
	    if(! defined $filter_map{$sample}){
		$sample_list->{$sample} = 1;
		push(@sample_info, $sample);
		$filter_map{$sample} = 1;
	    }
	}
	
	#Update the geneannotation
	my %gene_annot = ("MUT_FREQ", $mut_freq, "STATUS", $gene_status, "SAMPLE", \@sample_info, "MUTATION_TYPE", \%mutation_type);
	$data_gene_annotation->{$gene} = \%gene_annot;
    }
    close(FILE);
    
    print STDERR " *** Completed\n";

}

sub cancer_relaxed_annotation{
    my ($map_gene, $script_dir) = @_;
    my $NGC4_false_positive_file = "$script_dir/ANALYSE_RESULT/NCG4_05-04-15_falsepositives_list.txt";
    my $NGC4_file = "$script_dir/ANALYSE_RESULT/ngc_4_cancergenes.txt";
    
    my %FP_gene = ();
    open(FILE, $NGC4_false_positive_file);
    <FILE>;
    while(<FILE>){
	chop $_;
	$FP_gene{$_} = 1
    }
    close(FILE);

    open(FILE, $NGC4_file);
    <FILE>;
    while(<FILE>){
	@line = split(/\t/, $_);
	$map_gene->{$line[1]} = 1;
    }
    close(FILE)
}


sub cancer_annotation{
    
    my ($cancer_census, $gene_status_selection, $script_dir) = @_;
    
    my $cancer_gene_census_file;
    my $gold_standard_dir = "$script_dir/GOLD_STANDARD/";
    
    #For the gene list to exclude
    my %exclude_gene_list = ();
    #CHASM list
    #if(index($gene_status_selection, "CHASM") != -1 || index($gene_status_selection, "fathmm") != -1 ){
    if(index($gene_status_selection, "CHASM") != -1 ){
	$cancer_gene_census_file = "$gold_standard_dir/CHASM_training_genes.txt";
	open(FILE, $cancer_gene_census_file);
	#<FILE>;
	while(<FILE>){
	    chop $_;
	    $gene = $_;
	    #chomp($line[0]);
	    $exclude_gene_list{$gene} = 1;
	}
	close(FILE);
	print STDERR $gene_status_selection." + CHASM list size = ".(keys %exclude_gene_list)."\n";
    }
    
    #if(index($gene_status_selection, "fathmm") != -1 || index($gene_status_selection, "CHASM") != -1){
    if(index($gene_status_selection, "fathmm") != -1){
	$cancer_gene_census_file = "$gold_standard_dir/fathmm_training_genes.txt";
	open(FILE, $cancer_gene_census_file);
	#<FILE>;
	while(<FILE>){
	    chop $_;
	    $gene = $_;
	    #chomp($line[0]);
	    $exclude_gene_list{$gene} = 1;
	}
	close(FILE);
	print STDERR $gene_status_selection." + fathmm list size = ".(keys %exclude_gene_list)."\n";
    }

    #We went to remove all gene for wchich we have an annotation from the training set
    if(index($gene_status_selection, "NO_ANNOTATION") != -1){
	foreach $gene (keys %exclude_gene_list){
	    $cancer_census->{$gene} = $exclude_gene_list{$gene};
	}
    }
    
    #cancer gene census
    #my $cancer_gene_census_file = "/mnt/projects/bertrandd/oncoimpact/SCRIPT/oncoIMPACT/cancer_gene_census.csv";
    if(index($gene_status_selection, "CANCER") != -1 && index($gene_status_selection, "FP") == -1){
	$cancer_gene_census_file = "$gold_standard_dir/cancer_gene_census_06_26_2015.tsv";
	open(FILE, $cancer_gene_census_file);
	#<FILE>;
	while(<FILE>){
	    chop $_;
	    @line = split(/\s+/, $_);    
	    $gene = $line[0];
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    $cancer_census->{$gene} = 1;
	}
	close(FILE);
	print STDERR $gene_status_selection." + CGS list size = ".(keys %{$cancer_census})."\n";
    }
        
    if(index($gene_status_selection, "CGC_CNA") != -1 && index($gene_status_selection, "FP") == -1){
	#$cancer_gene_census_file = "$gold_standard_dir/cancer_gene_census_06_26_2015.tsv";
	open(FILE, "cat $gold_standard_dir/Census_amp_13_07_2015.csv $gold_standard_dir/Census_del_13_07_2015.csv |");
	#open(FILE, "cat $gold_standard_dir/Census_del_13_07_2015.csv |");
	#<FILE>;
	while(<FILE>){
	    chop $_;
	    @line = split(/\t/, $_);    
	    $gene = $line[0];
	    next if($gene eq "Gene Symbol");
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    $cancer_census->{$gene} = 1;
	}
	close(FILE);
	print STDERR $gene_status_selection." + CGS_CN list size = ".(keys %{$cancer_census})."\n";
	
	$cancer_gene_census_file = "$gold_standard_dir/NCG_CNV_gene.txt";
	open(FILE, $cancer_gene_census_file);
	<FILE>;
	while(<FILE>){
	    chop $_;
	    @line = split(/\s+/, $_);    
	    $gene = $line[1];
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    $cancer_census->{$gene} = 1 if(!exists $only_cs{$gene});
	}
	close(FILE);
	print STDERR $gene_status_selection."+ CGC_CN list size = ".(keys %{$cancer_census})."\n";
	
    }


    #print STDERR " *** Nb CC gene:".(keys %{$cancer_census})."\n";
    
    if(index($gene_status_selection, "CANCER_NO_F1") != -1 || index($gene_status_selection, "CANCER_UNION") != -1){

	my %only_cs = %{$cancer_census};
	%{$cancer_census} = () if($gene_status_selection ne "CANCER_UNION");
	#uniprot oncogene protein key word
	$cancer_gene_census_file = "$gold_standard_dir/uniprot-keyword_oncogen.gene";
	open(FILE, $cancer_gene_census_file);
	<FILE>;
	while(<FILE>){
	    chop $_;
	    @line = split(/\s+/, $_);    
	    #chomp($line[0]);
	    $gene = $line[0];
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    $cancer_census->{$gene} = 1 if(!exists $only_cs{$gene} || $gene_status_selection eq "CANCER_UNION");
	}
	close(FILE);
	print STDERR $gene_status_selection." + uniprot list size = ".(keys %{$cancer_census})."\n";
	
	#print STDERR " *** Nb CC + PROTgene:".(keys %{$cancer_census})."\n";<STDIN>;
	
	#Need to find the paper
	$cancer_gene_census_file = "$gold_standard_dir/NCG_CNV_gene.txt";
	open(FILE, $cancer_gene_census_file);
	<FILE>;
	while(<FILE>){
	    chop $_;
	    @line = split(/\s+/, $_);    
	    $gene = $line[1];
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    $cancer_census->{$gene} = 1 if(!exists $only_cs{$gene} || $gene_status_selection eq "CANCER_UNION");
	}
	close(FILE);
	print STDERR $gene_status_selection."+ NCG_CNV list size = ".(keys %{$cancer_census})."\n";
	
	#cut -f1 ~/vogelstein_* | grep -v Table - | grep -v Gene - | grep -v number > ~/vogelstein_list.dat
	$cancer_gene_census_file = "$gold_standard_dir/vogelstein_list.dat";
	open(FILE, $cancer_gene_census_file);
	<FILE>;
	while(<FILE>){
	    chop $_;
	    @line = split(/\s+/, $_);    
	    #chomp($line[0]);
	    $gene = $line[0];
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    $cancer_census->{$gene} = 1 if(!exists $only_cs{$gene} || $gene_status_selection eq "CANCER_UNION");
	}
	close(FILE);
	print STDERR $gene_status_selection." + vogelstein list size = ".(keys %{$cancer_census})."\n";
	
	#data from DISEASES: Test mining and data integration of disease-gene association used in xxx paper
	$cancer_gene_census_file = "$gold_standard_dir/human_disease_textmining_filtered.csv";
	open(FILE, $cancer_gene_census_file);
	while(<FILE>){
	    chop $_;
	    @line = split(/\t/, $_);    
	    #chomp($line[0]);
	    $gene = $line[1];
	    #
	    if(exists $exclude_gene_list{$gene} && $cancer_census->{$gene}){
		delete  $cancer_census->{$gene};
		next;
	    }
	    next if(exists $exclude_gene_list{$gene});
	    #
	    
	    $disease_type = $line[3];
	    $star_score = $line[5];
	    #print STDERR " $gene $disease_type $star_score\n";<STDIN>; 
	    if($disease_type eq "Cancer" && $star_score >= 2.5 && (!exists $only_cs{$gene} || $gene_status_selection eq "CANCER_UNION")){
		$cancer_census->{$gene} = 1;
		#print STDERR " *** $gene\n";
	    }
	}
	close(FILE);
	print STDERR $gene_status_selection." + human_disease list size = ".(keys %{$cancer_census})."\n";#<STDIN>;

	#cut -f2,4,6 GOLD_STANDARD/human_disease_textmining_filtered.csv | grep -e Cancer  | awk '{if($3 > 2.5) print $1}' | wc -l

    }

    if(index($gene_status_selection, "CANCER_FP") != -1){
	$cancer_gene_census_file = "$gold_standard_dir/NCG4_05-04-15_falsepositives_list.txt";
	open(FILE, $cancer_gene_census_file);
	while(<FILE>){
	    chop $_;
	    @line = split(/\s+/, $_);    
	    $gene = $line[0];
	    $cancer_census->{$gene} = 1;
	}
	close(FILE);
	print STDERR $gene_status_selection." + NGC4 = ".(keys %{$cancer_census})."\n";#<STDIN>;

	$cancer_gene_census_file = "$gold_standard_dir/cell_drug_FP_set.txt";
	open(FILE, $cancer_gene_census_file);
	while(<FILE>){
	    chop $_;
	    $gene = $_;
	    $cancer_census->{$gene} = 1;
	}
	close(FILE);
	print STDERR $gene_status_selection." + cancer cell = ".(keys %{$cancer_census})."\n";#<STDIN>;
    }

    print STDERR $gene_status_selection." FINAL list size = ".(keys %{$cancer_census})."\n";

}

#cut -f7,11-12,46 DATA/OV/annovar.hg18_multianno.txt|head -n1;cut -f2,7,11-12,46 DATA/OV/annovar.hg18_multianno.txt | grep 7518263
