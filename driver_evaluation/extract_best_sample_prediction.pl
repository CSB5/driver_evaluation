#!/usr/bin/perl
use warnings;

#require '/mnt/projects/bertrandd/oncoimpact/SCRIPT/MUTATION_BENCHMARK/ANALYSE_RESULT/common_functions_edit.pl';

my ($consolidated_result_dir, $data_gene_annotation_file, $max_rank, $additional_method_file, $script_dir) = @ARGV;

 require "$script_dir/common_functions_edit.pl";

my $out_dir = "$consolidated_result_dir\_$max_rank";
run_exe("mkdir $out_dir") if(! -d $out_dir);#<STDIN>;

my (%pan_cancer, %cancer_census);
cancer_annotation(\%cancer_census, "CANCER", $script_dir);

my %sample_list = ();
my %data_gene_annotation = ();
gene_annotation(\%data_gene_annotation, \%sample_list, $data_gene_annotation_file, \%cancer_census);


#Construt the result data structures
#Construt the result data structures
my %method_ID = ();
my @ID_to_method = ();
my %method_result = ();
read_method_result(\%method_result, \%method_ID, \@ID_to_method, \%data_gene_annotation, $consolidated_result_dir, $additional_method_file, 1, $script_dir);


#Init the method_driver_sample structure
my %method_sample_driver;
foreach $method (("BASE_LINE", @ID_to_method)){
    my %map = ();
    $method_sample_driver{$method} = \%map;
    foreach $sample (keys %sample_list){
	my %map = ();
	$method_sample_driver{$method}->{$sample} = \%map;
    }
}

#Compute method_driver_sample structure based on gene

#For the method
my $flag_use_sample_score = 1;
my $nb_driver;
my %gene_result = ();
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

    #Write the output file
    %gene_result = ();
    open(OUT, ">$out_dir/$method.result");
    foreach $sample (keys %sample_list){
	$sample_method_res = $method_sample_driver{$method}->{$sample};
	$nb_driver = 0;
	#print STDERR " *** ".join("\t", (keys %{$sample_method_res}));<STDIN>;
	foreach $gene (sort by_score_and_rank (keys %{$sample_method_res})){
	    #if(! defined $data_gene_annotation->{$gene}->{"STATUS"}){
	    #print STDERR " *** WARNING WEIRD gene in list:\t$method\t$sample\t $gene\n";
	    #next;
	    #}
	    if(! exists $gene_result{$gene}){
		my @sample_info = ("","");
		$gene_result{$gene} = \@sample_info;
	    }
	    $gene_result{$gene}->[0] .= $sample.";";
	    $gene_result{$gene}->[1] .= "NA".";";
	    $nb_driver++;
	    last if($nb_driver == $max_rank);
	}
    }
    open(OUT, ">$out_dir/$method.result");
    print OUT "Gene_name\tSample\tRank\tScore\tInfo\tSample-specific_score\n";
    foreach $gene (keys %gene_result){
	print OUT $gene."\t".$gene_result{$gene}->[0]."\t".$method_result{$method}->{$gene}->{"RANK"}."\t".$method_result{$method}->{$gene}->{"RANK"}."\t"."-"."\n";
    }
    close(OUT);
}

sub by_score_and_rank {
    #print STDERR " ***by_score_and_rank\n----".$a."\n----".$b."\n";<STDIN>;
    #print STDERR " ***by_score_and_rank\n----".(join("\t", (@{$sample_method_res->{$a}})))."\n----".(join("\t", (@{$sample_method_res->{$b}})))."\n";<STDIN>;
    ($sample_method_res->{$b}->[0] <=> $sample_method_res->{$a}->[0] ||#sample spcefic gene score
     $sample_method_res->{$a}->[1] <=> $sample_method_res->{$b}->[1] #Gene rank 
    );
}

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";#<STDIN>;
    print STDERR `$exe` if($run);
}
