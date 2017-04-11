#!/usr/bin/perl
use warnings;

#require '/mnt/projects/bertrandd/oncoimpact/SCRIPT/MUTATION_BENCHMARK/ANALYSE_RESULT/common_functions_edit.pl';

my ($consolidated_result_dir, $data_gene_annotation_file, $out_dir, $gene_status_selection, $flag_plot, $additional_method_file, $use_known_method, $append_result, $flag_only_driver_number, $script_dir) = @ARGV;

 require "$script_dir/common_functions.pl";

#$flag_meta = 0 if( ! defined $flag_meta );

#`mkdir $out_dir`;
our (@RANK_THRESHOLD);
my @CAT_MUT_SHARE = (0, 25, 50, 75);

#*** TO DO ***
#Remove the MUT_FREQ threshold type 
#Handle the ALL in the RANK_THRESHOLD list
my %threshold_list = 
    ("MUT_FREQ", [0], #[0.10, 0.05, 0.02, 0.01, 0],
     "RANK", \@RANK_THRESHOLD);
     #"MUT_FREQ",  [0]);
     #"RANK", [5]);

$gene_status_selection = "" if(! defined $gene_status_selection);


#
#To update using that file:
#/mnt/pnsg10_projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/TEST_DATA/COAD/gene_mutation_frequency.txt
#


my %cancer_census;
cancer_annotation(\%cancer_census, $gene_status_selection, $script_dir);

my %sample_list = ();
my %data_gene_annotation = ();
gene_annotation(\%data_gene_annotation, \%sample_list, $data_gene_annotation_file, \%cancer_census);


#Construt the result data structures
my %method_ID = ();
my @ID_to_method = ();
my %method_result = ();
read_method_result(\%method_result, \%method_ID, \@ID_to_method, \%data_gene_annotation, $consolidated_result_dir, $additional_method_file, $use_known_method, $script_dir);
#read_method_result(\%method_result, \%method_ID, \@ID_to_method, \%data_gene_annotation, $consolidated_result_dir);

#Write the file that contain the number of calls of each methods
#if($flag_meta){
#    write_driver_number_file(\%method_result, "$out_dir/driver_number_meta.dat");
#}
#else{
write_driver_number_file(\%method_result, "$out_dir/driver_number.dat");
#}
exit(0) if($flag_only_driver_number);

#Init the pairwise matrise method distance
#print STDERR " *** @ID_to_method\n";<STDIN>;
my @matrix = ();
for(my $i = 0; $i < @ID_to_method; $i++){
    my @tab = ();
    for(my $j = 0; $j < @ID_to_method; $j++){
	$val = 0;
	push(@tab, $val);
    }
    push(@matrix, \@tab);
}

#For the gene name
my $gene_diretory = $out_dir."/../GENE";
#run_exe("mkdir -p $gene_diretory") if(! -d $gene_diretory);
#<STDIN>;

$gene_diretory .= "/$gene_status_selection/";
run_exe("rm -r $gene_diretory") if(-d $gene_diretory);
run_exe("mkdir -p $gene_diretory") if(! -d $gene_diretory);

#Write the pairwise comparison matrises
my $inter;
my $matrix_file;
foreach $threshold_type (keys %threshold_list){
    foreach $threshold_val (@{$threshold_list{$threshold_type}}){
	#Only process the 50 rank !!!
	next if(!($threshold_type eq "RANK" && $threshold_val == 50));
	#
	clear_matrix(\@matrix);
	foreach $method_1 (@ID_to_method){
	    foreach $method_2 (@ID_to_method){
		#next if($method_1 eq $method_2);
		$inter = 0;
		$list_size = 0;
		#print STDERR $threshold_type."\t".$threshold_val."\t".$method_1."\t".$method_2."\t".$gene."\n";#<STDIN>;
		foreach $gene (keys %{$method_result{$method_1}}){
		    
		    #To test the gene status selection
		    next if(! exists $data_gene_annotation{$gene} || #should be replaced
			    ($gene_status_selection ne "ALL" && $data_gene_annotation{$gene}->{"STATUS"} ne "POSITIVE"));

		    #The must pass the threshold for method_1
		    #print $threshold_type."\t".$threshold_val."\t".$method_1."\t".$method_2."\t".$gene."\n";#<STDIN>;
		    if(test_threshold($threshold_type, $threshold_val, $gene, $method_1)){
			
			#Single method analysis
			next if($method_1 eq $method_2);
			    
			$list_size++;
			$method_result{$method_1}->{$gene}->{"USE"} = 1;
			
			#print STDERR $threshold_type."\t".$threshold_val."\t".$method_1."\t".$method_2."\t".$gene."\n";#<STDIN>;
			if(exists $method_result{$method_2}->{$gene} && test_threshold($threshold_type, $threshold_val, $gene, $method_2)){
			    $inter++;
			    $method_result{$method_1}->{$gene}->{"FREQ_CALL"}++;
			}
		    }
		}
		#$matrix[$method_ID{$method_1}]->[$method_ID{$method_2}] = sprintf("%.2f", $inter/$list_size) if($list_size != 0);
		$matrix[$method_ID{$method_1}]->[$method_ID{$method_2}] = $inter if($list_size != 0);
	    }
	}
	
	if($gene_status_selection eq "ALL" && $threshold_val == 50){
	    write_gene_heatmap_file($threshold_val);
	    #exit(0);
	}
	
	#For the violinplots of the driver mutation frequencies
	$matrix_file = write_violin_plot_file($threshold_type, $threshold_val, $gene_status_selection);
	print STDERR " *** violin plot file $matrix_file\n";
	
	#For the barplot pairwise comparision files
	$gene_file = "NONE";
	$gene_file = "$gene_diretory/$threshold_type\_$threshold_val" if($threshold_val == 100);
	$matrix_file = write_barplot_file($threshold_type, $threshold_val, $gene_status_selection, $gene_file);
	
	plot_barplot($matrix_file, "$threshold_type\_$threshold_val", @ID_to_method+0, 1) if($flag_plot eq "PLOT");
	plot_barplot($matrix_file, "$threshold_type\_$threshold_val", @ID_to_method+0, 0) if($flag_plot eq "PLOT");
	
	####For the HEAT map pairwise comparision files
	#Write the file
	$matrix_file = write_heat_map_file(\@matrix, $threshold_type, $threshold_val, $gene_status_selection);
	plot_heat_map($matrix_file, "$threshold_type\_$threshold_val") if($flag_plot eq "PLOT");
	
	#Delete the gene value for the next test
	for(my $i = 0; $i < @ID_to_method; $i++){
	    $method = $ID_to_method[$i];
	    #print STDERR " *** $method\n";<STDIN>;
	    foreach $gene (keys %{$method_result{$method}}){
		#print STDERR " *** $gene\n";<STDIN>;
		$method_result{$method}->{$gene}->{"FREQ_CALL"} = 0;
		$method_result{$method}->{$gene}->{"USE"} = 0;
	    }
	}

    }
}

sub write_gene_heatmap_file{
    my ($threshold_val) = @_;
    my $method_file = "$out_dir/gene_method_pred_RANK_$threshold_val.dat";
    my $gene_file = "$out_dir/gene_mut_type_file.dat";
    open(OUT_M, ">$method_file");
    open(OUT_G, ">$gene_file");
    
    my $flag_gene_with_pred = 0;

    #Get the gene list used: union of the top $threshold_val calls
    my %gene_list = ();
    for(my $i = 0; $i < @ID_to_method; $i++){
	$method = $ID_to_method[$i];
	foreach $gene (keys %{$method_result{$method}}){
	    $gene_list{$gene} = 1 if($method_result{$method}->{$gene}->{"USE"} == 1);
	}
    }
    
    #Out the file for each gene
    #The headers
    print OUT_G "MISSENSE"."\t"."OTHER-SNV"."\t"."CNA"."\n";
    print OUT_M "".join("\t", @ID_to_method)."\n";
    foreach $gene (keys %gene_list){

	$flag_gene_with_pred = 0;

	$cna_freq = $data_gene_annotation{$gene}->{"MUTATION_TYPE"}->{"CNA"};
	$miss_freq = $data_gene_annotation{$gene}->{"MUTATION_TYPE"}->{"MISSENSE"};
	$other_freq = $data_gene_annotation{$gene}->{"MUTATION_TYPE"}->{"OTHER_SNV"};

	#The gene file
	print OUT_G 
	    $gene.
	    "\t".$data_gene_annotation{$gene}->{"MUT_FREQ"}.
	    "\t".$miss_freq.
	    "\t".$other_freq.
	    "\t".$cna_freq.
	    "\n";
	
	#The method file
	print OUT_M $gene;
	for(my $i = 0; $i < @ID_to_method; $i++){
	    $method = $ID_to_method[$i];
	    $res = -0.1;
	    if(exists $method_result{$method}->{$gene} && $method_result{$method}->{$gene}->{"USE"} == 1){
		#$res = 55;
		$mutation_freq = $cna_freq + $miss_freq + $other_freq;
		$res = ($cna_freq) / $mutation_freq if($mutation_freq != 0);
		$flag_gene_with_pred = 1;
	    }
	    print OUT_M "\t".$res;
	}
	if(!$flag_gene_with_pred){
	    print STDERR " *** WARNING There is no prediction of $gene\n";<STDIN>;
	}
	print OUT_M "\n";
    }

    close(OUT_G);close(OUT_M);

}

sub write_driver_number_file{
    my ($method_result, $out_file) = @_;
    print STDERR " *** ".$out_file."\n";
    $write_type = ">"; $write_type = ">>" if($append_result eq "APPEND_RESULT");
    open(OUT, $write_type, "$out_file");
    for(my $i = 0; $i < @ID_to_method; $i++){
	$method = $ID_to_method[$i];
	$nb_driver = keys (%{$method_result->{$method}});
	print OUT $method."\t".$nb_driver."\n";
    }
    close(OUT);
}

sub write_violin_plot_file{
    my ($threshold_type, $threshold_val, $gene_status_selection) = @_;
    my %method_driver_freq = ();
    $max_nb_driver = 0;
    for(my $i = 0; $i < @ID_to_method; $i++){
	$method = $ID_to_method[$i];
	my @tab = ();
	$method_driver_freq{$method} = \@tab;
	foreach $gene (keys %{$method_result{$method}}){
	    #Update the matrix
	    if($method_result{$method}->{$gene}->{"USE"}){
		#$driver_mut_freq = @{$data_gene_annotation{$gene}->{"SAMPLE"}};
		#patient speciffic driver predictions !!!
		$driver_mut_freq = @{$method_result{$method}->{$gene}->{"SAMPLE"}}+0;
		#print STDERR " *** gene ".$gene."\t".$driver_mut_freq."\n";<STDIN>;# if($method eq "NetBox" && $driver_mut_freq < 0.02);
		push(@{$method_driver_freq{$method}}, $driver_mut_freq);
	    }
	}
	$nb_driver = @{$method_driver_freq{$method}}+0;
	$max_nb_driver = $nb_driver if($nb_driver > $max_nb_driver);
    }

    #Write the file
    my $matrix_file = "$out_dir/driver_freq_distribution_$threshold_type\_$threshold_val";
    if($gene_status_selection ne ""){
	$matrix_file .= "\_$gene_status_selection";
    }
    $write_type = ">"; $write_type = ">>" if($append_result eq "APPEND_RESULT");
    open(OUT, $write_type, "$matrix_file.dat");

    for(my $i = 0; $i < @ID_to_method; $i++){
	$method = $ID_to_method[$i];
	$nb_driver = @{$method_driver_freq{$method}}+0;
	#print STDERR "$method -> $nb_driver\n";<STDIN>;
	print OUT $method."\t".(join("\t", @{$method_driver_freq{$method}}));
	for(my $j = $nb_driver; $j < $max_nb_driver; $j++){
	    print OUT "\t"."NA";
	}
	print OUT "\n";
    }
    close(OUT);

    return $matrix_file;
    #profile <- read.table("../AAA_RES/driver_freq_distribution_RANK_50_ALL.dat", header = FALSE, row.names=1)
    #t_p = t(profile)
    #boxplot(t_p)
}


#Plot the barplot
sub write_barplot_file{
    my ($threshold_type, $threshold_val, $gene_status_selection, $gene_list_pref) = @_;

    $gene_list_pref = "NONE" if(! defined $gene_list_pref);
    #print STDERR " *** Write $gene_list_pref\n";<STDIN>;

    #Contruct the frequency call matrix
    #Init the call frequence matrice method
    my @matrix_freq_call = ();
    for(my $i = 0; $i < @ID_to_method; $i++){
	#Init the value 
	$method = $ID_to_method[$i];
	open(OUT, ">$gene_list_pref\_$method.gene") if($gene_list_pref ne "NONE");
	my @tab = ();
	for(my $j = 0; $j < @ID_to_method; $j++){
	    $val = 0;
	    push(@tab, $val);
	}

	foreach $gene (keys %{$method_result{$method}}){
	    #Update the matrix
	    if($method_result{$method}->{$gene}->{"USE"}){
		$tab[$method_result{$method}->{$gene}->{"FREQ_CALL"}]++;
		print OUT $gene."\t".$method_result{$method}->{$gene}->{"RANK"}."\n" if($gene_list_pref ne "NONE");
	    }
	}
	close(OUT) if($gene_list_pref ne "NONE");

	push(@matrix_freq_call, \@tab);
	
    }
    close(OUT);

    my $matrix_file = "$out_dir/freq_call_$threshold_type\_$threshold_val";
    if($gene_status_selection ne ""){
	$matrix_file .= "\_$gene_status_selection";
    }
    #if($flag_meta){
    #$matrix_file .= "_meta";
    #}

    $write_type = ">"; $write_type = ">>" if($append_result eq "APPEND_RESULT");
    open(OUT, $write_type, "$matrix_file.dat");
    #header
    print OUT "".join("\t", @ID_to_method)."\n";
    for(my $j = 0; $j < @{$matrix_freq_call[0]}; $j++){
	print OUT $j;
	for(my $i = 0; $i < @ID_to_method; $i++){
	    print OUT "\t".$matrix_freq_call[$i]->[$j];
	}
	print OUT "\n";
    }
    close(OUT);

    return $matrix_file;

}

sub plot_barplot{
    my ($matrix_file, $title, $nb_method, $raw_val) = @_;
    
    my $font_size = 2.5;
    my $note_font_size = 2.5;
    my $margin_size = 30; 
    

    $out_file = $matrix_file;
    $out_file .= "_RAW" if($raw_val);

    open(OUT, ">$out_file.R");
    
    print OUT "pdf(file=\"$out_file.pdf\",
	paper=\"special\",
	width=25,
	height=25
	)\n";
    print OUT "profile <- read.table(\"$matrix_file.dat\", header = TRUE)\n";

    #Compute the freqency
    if(! $raw_val){
	foreach $method (@ID_to_method){
	    print OUT "profile\$$method = profile\$$method/sum(profile\$$method)\n";
	}
    }

    print OUT "profile_mat = as.matrix(profile)\n";
    print OUT "palette <- colorRampPalette(c('#f0f3ff','#0033BB'))($nb_method)\n";
    print OUT "barplot(profile_mat, main=\"$title\", col=palette,  cex.lab=$font_size, cex.names=2, cex.axis = $font_size, cex.main= $font_size)\n";
    
    if(index($matrix_file, "RANK") == -1 && $raw_val){
	print OUT "legend(\"topright\", legend = rownames(profile), fill = palette, cex=2);"
    }

    close(OUT);
    
    run_exe("R-3.0.0 --vanilla < $out_file.R");
}

#Plot the heatmap
sub write_heat_map_file{
    my ($mat, $threshold_type, $threshold_val, $gene_status_selection) = @_;

    #need update for different threshold values
    #Work only for MUT_FREQ_0
    my $nb_mutated_genes = keys(%data_gene_annotation);

    my $matrix_file = "$out_dir/pairwise_comparision_$threshold_type\_$threshold_val";
    if($gene_status_selection ne ""){
	$matrix_file .= "\_$gene_status_selection";
    }

    $write_type = ">"; $write_type = ">>" if($append_result eq "APPEND_RESULT");
    open(OUT, $write_type, "$matrix_file.dat");
    #header
    print OUT "".join("\t", @ID_to_method)."\n";

    my ($union_size, $nb_common_driver, $min_driver) = 0;

    for(my $i = 0; $i < @ID_to_method; $i++){
	$method_i = $ID_to_method[$i];
	#
	$nb_driver_i = keys (%{$method_result{$method_i}});
	$nb_driver_i = $threshold_val if($threshold_type eq "RANK" && $threshold_val < $nb_driver_i);
	#
	print OUT $method_i;
	for(my $j = 0; $j < @ID_to_method; $j++){

	    if($j == $i){
		print OUT "\t1";
		next;
	    }

	    $method_j = $ID_to_method[$j];
	    #
	    $nb_driver_j = keys (%{$method_result{$method_j}});
	    $nb_driver_j = $threshold_val if($threshold_type eq "RANK" && $threshold_val < $nb_driver_j);
	    #
	    $min_driver = $nb_driver_j;$min_driver = $nb_driver_i if($nb_driver_i < $nb_driver_j);
	    #
	    $nb_common_driver = $mat->[$i]->[$j];
	    #
	    $union_size = $nb_driver_i  + $nb_driver_j - $nb_common_driver;
	    $comparison_value = 0;
	    
	    #Jaccard index: Fraction of the union shared
	    $comparison_value = $nb_common_driver / $union_size if($union_size != 0);
	    
	    #Fractin of the prediction shared
	    #$comparison_value = $nb_common_driver / $nb_driver_i if($nb_driver_i != 0);
	    
	    print OUT "\t".$comparison_value;
	    #print OUT "\t".((2 * $nb_common_driver) / ($nb_driver_i  + $nb_driver_j));
	    #print OUT "\t".((($nb_common_driver / $nb_driver_i) + ($nb_common_driver / $nb_driver_j))/2);
	    #print OUT "\t".($nb_common_driver / $min_driver);
	    #
	    #This method did not work due to the high number of passenger genes compared to the number of drivers
	    #
	    #I: common driver
	    #A: nb_driver_i
	    #B: nb_driver_j
	    #N: nb_mutated_genes
	    #P: nb_common_passenger
	    #P = |N| - (|A| + |B|) + |I| = | N - (A U B) |
	    #$nb_common_passenger = $nb_mutated_genes - ($nb_driver_i + $nb_driver_j) + $nb_common_driver;
	    #print STDERR $method_i."\t".$method_j."\n";
	    #print STDERR $nb_mutated_genes."\t".$nb_driver_i."\t".$nb_common_driver."\t".$nb_driver_j."\t".$nb_common_driver."\t = ".$nb_common_passenger."\n";<STDIN>;
	    #Aggrement
	    #print OUT "\t".(sprintf("%.2f",($nb_common_driver + $nb_common_passenger) / $nb_mutated_genes));
	}
	print OUT "\n";
    }
    close(OUT);
    
    return $matrix_file;

}

sub plot_heat_map{
    my ($matrix_file, $title) = @_;

    my $font_size = 4;
    my $note_font_size = 4;
    my $margin_size = 30; 

    open(OUT, ">$matrix_file.R");
    
    print OUT "pdf(file=\"$matrix_file.pdf\",
	paper=\"special\",
	width=35,
	height=35
	)\n";

    print OUT "library(\"gplots\")\n";
    print OUT "palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(10)\n";
    print OUT "profile <- read.table(\"$matrix_file.dat\", header = TRUE)\n";
    print OUT "profile_mat = as.matrix(profile)\n";
#print OUT "heatmap.2(profile_mat, col=redgreen, margin=c(5, 20), key=TRUE, scale=\"row\", density.info=\"none\", trace=\"none\")\n";
    
    
    my $str_key = "key = FALSE,
lmat=rbind(c(2),c(3),c(1),c(3)), 
    lhei=c(1,1,9,0), 
    lwid=c(1)";
    
    if($matrix_file eq "$out_dir/pairwise_comparision_MUT_FREQ_0" || $matrix_file eq "$out_dir/pairwise_comparision_MUT_FREQ_0_CANCER"){
	$str_key = "key=TRUE, keysize=1.3";
	$note_font_size = 4;
	$font_size = 6;
	$margin_size = 30; 
    }
    
    print OUT 
	"heatmap.2(profile_mat,
main = \"$title\", 
scale=\"none\",
Rowv=FALSE,Colv=FALSE,
breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), #symbreaks = TRUE,
#
#key.par=list(mgp=c(1.5, 0.5, 0),
#mar=c(2.5, 2.5, 1, 0)),
#
dendrogram = \"none\", 
#$str_key,
cellnote=as.matrix(profile_mat),notecol=\"black\",notecex=$note_font_size,
               #hclustfun = function(x) hclust(x,method = 'complete'),
               #distfun = function(x) dist(x,method = 'euclidean'),
               margin=c($margin_size, $margin_size), 
col=palette, cexRow=$font_size, cexCol=$font_size, 
               density.info=\"none\", trace=\"none\"";
    #print OUT ",ColSideColors = $color_subtype" if($subtype_file ne "NONE");
    print OUT ")\n";
    close(OUT);

    run_exe("R-3.0.0 --vanilla < $matrix_file.R");
    run_exe("pdfcrop  $matrix_file.pdf $matrix_file\_temp.pdf");
    run_exe("mv  $matrix_file\_temp.pdf $matrix_file.pdf");

}

sub clear_matrix{
    my ($mat) = @_;
    for(my $i = 0; $i < @{$mat}; $i++){
	for(my $j = 0; $j < @{$mat->[$i]}; $j++){
	    $val = 0;
	    $val = 1 if($i == $j); 
	    $mat->[$i]->[$j] = $val;
	}
    }
}

sub test_threshold{
    my ($t_type, $t_val, $gene, $meth) = @_;
    if(#Gene annotation base threshold
       ($t_type =~ m/MUT_FREQ/ && $data_gene_annotation{$gene}->{$t_type} >= $t_val) ||
       #Method output based threshold
       ($t_type =~ m/RANK/ && $method_result{$meth}->{$gene}->{$t_type} <= $t_val)){
	return 1;
    }
    return 0;
}


sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";
    print STDERR `$exe` if($run);
}

#heatmap.2(profile_mat, scale="none", col=palette, density.info="none", trace="none", hclustfun = function(x) hclust(x,method = 'complete'), distfun = function(x) dist(x,method = 'euclidean'), dendrogram = "col", key = FALSE, labRow = NULL, ylab = NULL)


