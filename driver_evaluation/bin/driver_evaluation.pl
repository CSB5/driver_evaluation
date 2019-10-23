#!/usr/bin/perl
use warnings;
use Config::Simple;
#use Getopt::Long;
use POSIX 'strftime';
use 5.010;


#########
#my ($analysis_dir, $final_out_dir, $script_dir, $selected_method_file) = @ARGV;
#
my $FLAG_ANANLIZE_DATA_SET = 1;
my $FLAG_PLOT_SAMPLE_BASED_ANALYSIS = 1;
my $FLAG_PLOT_DRIVER_MUTATION_FREQUENCY = 1;
my $FLAG_PLOT_METHOD_COMPARISON_HEAT_MAP = 0;
my $FLAG_COMPUTE_DRIVER_NUMBER = 1;
my $FLAG_COMPARISON_WITH_CANCER_GENE = 1;
my $FLAG_PLOT_CONCORDANCE_ANALYSIS = 1;
my $FLAG_PLOT_SHARED_METHOD_CALL = 0;
my $FLAG_ONLY_DRIVER_NUMBER = 0;
my $FLAG_ANALIZE_DOWNSAMPLING_DATA_SET   =0;
my $FLAG_ANALIZE_FP_DATA_SET   =0;
my $FLAG_PLOT_FP_RESULT   =0;
my $FLAG_PLOT_SAMPLE_ACTIONABLE_ANALYSIS =1;



my $help_message = "

Usage:
        result_evaluation.pl <out_dir> <method_dir>

Options:
    --out_dir: path to the evaluation result directory
    --method_dir: path to the method prediction
\n";

if ( @ARGV == 0 ) {
        print $help_message;
        exit 0;
}

my $configFile; my $flag_help;

my $index = 0;
foreach my $argument (@ARGV){

	if($argument eq "--out_dir"){
		$final_out_dir = $ARGV[$index+1];
	}
	if($argument eq "--method_dir"){
                $new_method_dir = $ARGV[$index+1];
        }
	$index++;
}
if(!defined($final_out_dir)){
	print "Path to the evaluation result directory, not defined!\n";
	exit 0;
}
if(!defined($new_method_dir)){
        print "Path to the method prediction, not defined!\n";
        exit 0;
}

#print STDERR "result_evaluation.pl --out_dir $final_out_dir --method_dir $new_method_dir \n";

#VARIABLES THAT THE PATH TO THE DATA AND SCRIPT
#THOSE VARIABLE ARE INITIALAZED DURING THE INSTALLATION

#my $analysis_dir = "/mnt/projects/bertrandd/opera_lg/kohjy_tmp/driver_evaluation_2/driver_evaluation/EVALUATION_DATA_SET/RESULTS";
#my $script_dir   = "/mnt/projects/bertrandd/opera_lg/kohjy_tmp/driver_evaluation_2/driver_evaluation/bin";

my $analysis_dir = "XX_ANALYSIS_DIR";
my $script_dir   = "YY_SCRIPT_DIR";


if(index($analysis_dir, "XX_ANALYSIS") == 0 || index($script_dir, "YY_SCRIPT") == 0){
	print "Please run the installation script, install.pl !!\n";
	exit 0;
}

#print STDERR " *** Pass\n";exit(0);



run_exe("mkdir $final_out_dir") unless (-d $final_out_dir);

#print STDERR " *** $new_method_dir\n";
opendir(DIR, $new_method_dir);
my @all_cancer_type_prediction = readdir(DIR);
close(DIR);

my @tmp = split(/\//, $new_method_dir);
my $additional_method_name = $tmp[@tmp-1];
$additional_method_name = $new_method_dir if(@tmp == 0);

#$selected_method_file    =config{'general.selected_method_file'}; 
$selected_method_file    = "$final_out_dir/additional_method_name.txt";

#run_exe("rm $selected_method_file;touch $selected_method_file");
run_exe("echo -e \"Method\" > $selected_method_file");
run_exe("echo -e \"$additional_method_name\" >> $selected_method_file");

foreach $file(@all_cancer_type_prediction){
    if($file =~ /(.+)\.result/){
	$cancer_type = $1;
	run_exe("cp $new_method_dir/$file $analysis_dir/$cancer_type/$additional_method_name.result");
    }
}

#exit(0);


if (!defined($analysis_dir)) {
        print "analysis dir option need to be specified\n";
        exit 0;
}
if (!defined($final_out_dir)) {
        print "final out dir option need to be specified\n";
        exit 0;
}
if (!defined($script_dir)) {
        print "script dir option need to be specified\n";
        exit 0;
}

require "$script_dir/common_functions.pl";



#Order used for all the plots rpresenting the data sets
my @DATA_SET_ORDER_PLOT = ("READ", "PAAD", "BLCA", "KIRP", "STAD", 
			   "LUAD", "LUSC", "LIHC", "COAD", "GBM", 
			   "PRAD", "OV", "THCA", "KIRC", "BRCA");

my %CAT_COLOR = (
    "BASE_LINE", "\"black\"", #"\"blue\"",
    "FUNC_SNV", "\"purple\"", #"\"blue\"",
    "FUNC_SNV_C", "\"darkorange2\"", #\"hotpink4\"", #"\"blue\"",
    "BIAS_SNV", "\"darkblue\"",
    "INTEG_E_C", "\"darkgreen\"",
    "INTEG_N_E_S_C", "\"darkred\"",
    #
    "META_METHOD", "\"purple4\"",
    "ADHOC", "\"orangered\""#,"\"slategrey\"",
    );
    

my $BASE_LINE_FULL_NAME = "base_line";
#[Category, short_name
my %METHOD_INFO = (
    #Base line
    $BASE_LINE_FULL_NAME, ["BASE_LINE", "BL"],
    
    #SNV functional prediction
    "MutationTaster", ["FUNC_SNV", "MT"],
    "SIFT", ["FUNC_SNV", "SIFT"],
    "PolyPhen2", ["FUNC_SNV", "PP2"],
    "MutationAssessor", ["FUNC_SNV", "MA"],
    #
    "transFIC", ["FUNC_SNV_C", "TF"],
    "fathmm", ["FUNC_SNV_C", "FH"],
    "CHASM", ["FUNC_SNV_C", "CM"],
    #
    #SNV mutation bias
    "ActiveDriver", ["BIAS_SNV", "AD"],
    "OncodriveFM", ["BIAS_SNV", "OFM"],
    "OncodriveCLUST", ["BIAS_SNV", "OCL"], 
    "MutSigCV", ["BIAS_SNV", "MCV"],
    #
    #CNV + Expression
    "OncodriveCIS", ["INTEG_E_C", "OCI"],
    "S2N", ["INTEG_E_C", "S2N"],
    #Intergative method N + E + S + C
    #N + S + C 
    "HotNet2A", ["INTEG_N_E_S_C", "HN2"],
    "NetBox", ["INTEG_N_E_S_C", "NB"],
    #N + E + S + C
    "DriverNet", ["INTEG_N_E_S_C", "DN"],
    "DawnRank", ["INTEG_N_E_S_C", "DR"],
    "oncoIMPACT", ["INTEG_N_E_S_C", "OI"],
    #
    #The CNA only analysis
    #
    #
    #Name for the meta method
    "DriverDB", ["META_METHOD", "DDB"],
    "MutSig", ["META_METHOD", "MutSig"],
    
    );


############################################################################################
###################### Run the analysis of the real data sets ###############################
#############################################################################################
#Some global variables
my @tested_rank = @RANK_THRESHOLD;
my @tested_rank_all = (@RANK_THRESHOLD, "ALL");
my %measure_tested = ("F1", "Combined score",
		      "recall", "Sensitivity",
		      "precision", "Concordance");

opendir(DIR, $analysis_dir) || die "can't opendir $analysis_dir: $!";
@all_data_set = readdir(DIR);
#print join(':', @all_data_set);
close(DIR);

#Run the analysis for data set
@DATA_SET_ORDER = ();

$selected_method_file = "NONE" if(!defined($selected_method_file));

#my $selected_method_file = $ARGV[0];	#"/home/bertrandd/PROJECT_LINK/ONCOIMPACT/MUTATION_BENCHMARK/PAPER_RELEASE/extra_methods.txt";		#"NONE";
#my $selected_method_file = "/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/SOFTWARE_TESTBED/Borda/Top3_medianF1_consensus_methods.txt";
#my $selected_method_file = "/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/SOFTWARE_TESTBED/Borda/Top3_no_CHASM_medianF1_consensus_methods.txt";

my @new_method;
#print STDERR "$selected_method_file";
if($selected_method_file ne "NONE"){
	open(my $fh, "<", $selected_method_file)
		or die "Can't open < input.txt: $!";


<$fh>;
while(<$fh>)  {   
    chomp $_;
    my @line = split(/\s+/, $_);
    #print STDERR "!!!!! $line[0]\n";
    push(@new_method, $line[0]);
}
close($fh);
}


print STDERR " *** Analyse method results. Please wait ...\n";
foreach $data_set (@all_data_set){
	
    #next if($data_set eq "BRCA" || $data_set eq "UCEC");# || $data_set eq "LUSC");#Data set to exclude
    next if($data_set eq "UCEC");# || $data_set eq "LUSC");#Data set to exclude
     next if($data_set eq "." || $data_set eq ".." || $data_set =~ /^\./ );

    print STDERR "\n";
    $data_set_result_dir = "$analysis_dir/$data_set";
    #run_exe("mkdir -p $data_set_result_dir");
    if(-e "$data_set_result_dir"){
#####	print STDERR " *** ANALYSE $data_set\n";
	push(@DATA_SET_ORDER, $data_set);
	#$data_set_data_dir = "$analysis_dir/../DATA_TCGA/$data_set/gene_mutation_frequency.txt";
	$data_set_data_dir = "$analysis_dir/../DATA/$data_set/gene_mutation_frequency.txt";

	$out_dir = "$data_set_result_dir/PLOT";
	if( $FLAG_ANANLIZE_DATA_SET || ! -d $out_dir){
	    #if( ! -d $out_dir){
	    run_exe("mkdir -p $out_dir") if(! -d $out_dir);;
	    
	    #Run the pairwise comparision analysis
	    if($FLAG_COMPUTE_DRIVER_NUMBER || $FLAG_PLOT_METHOD_COMPARISON_HEAT_MAP || $FLAG_PLOT_SHARED_METHOD_CALL || $FLAG_PLOT_DRIVER_MUTATION_FREQUENCY){
 		run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir ALL NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
		if($FLAG_COMPARISON_WITH_CANCER_GENE){
		    #The Cancer plots
		    print STDERR " *** Analyze cancer type, $data_set ...\n";
		    run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_UNION NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
		    run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
		    run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CGC_CNA NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
		    run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_NO_CS NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
		    run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_FP NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
		    #
		    if($FLAG_PLOT_SHARED_METHOD_CALL){
			#The data exluding the CHASM test data set
			print STDERR " *** The data exluding the CHASM test data set ...\n";
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_CHASM NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_NO_CS_CHASM NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_FP_CHASM NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir NO_ANNOTATION_CHASM NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
			#
			#fathmm
			print STDERR " *** fathmm ...\n";
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_fathmm NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_NO_CS_fathmm NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_FP_fathmm NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
			run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir NO_ANNOTATION_fathmm NO_PLOT NONE 1 NEW_RESULT $FLAG_ONLY_DRIVER_NUMBER $script_dir");#<STDIN>;
		    }
		}
		 run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir ALL NO_PLOT $selected_method_file 1 APPEND_RESULT 1 $script_dir") if($selected_method_file ne "NONE");
		#run_exe("$script_dir/pair_wise_comparision.pl $meta_data_res_dir $data_set_data_dir $out_dir ALL NO_PLOT $selected_method_file 1 APPEND_RESULT");
	    }
	    
	    #Run the sample based analysis
	    print STDERR " *** Run the sample based analysis ...\n";
	    run_exe("$script_dir/sample_evaluation.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_UNION NONE 1 NEW_RESULT NO_PLOT $script_dir") if($FLAG_PLOT_SAMPLE_BASED_ANALYSIS && ($selected_method_file eq "NONE")) ;#<STDIN>;
	    run_exe("$script_dir/sample_evaluation.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER_UNION $selected_method_file 1 APPEND_RESULT NO_PLOT $script_dir") if($FLAG_PLOT_SAMPLE_BASED_ANALYSIS && ($selected_method_file ne "NONE"));	   
 	
	    #
	    #print STDERR "$FLAG_PLOT_SAMPLE_ACTIONABLE_ANALYSIS $selected_method_file\n";
	    print STDERR " *** Sample actionable gene ...\n";
	    run_exe("$script_dir/sample_actionable_gene.pl $data_set_result_dir $data_set_data_dir $out_dir NONE 1 NEW_RESULT NO_PLOT $script_dir") if($FLAG_PLOT_SAMPLE_ACTIONABLE_ANALYSIS  && ($selected_method_file eq "NONE")) ;#<STDIN>;
            run_exe("$script_dir/sample_actionable_gene.pl $data_set_result_dir $data_set_data_dir $out_dir $selected_method_file 1 APPEND_RESULT NO_PLOT $script_dir") if($FLAG_PLOT_SAMPLE_ACTIONABLE_ANALYSIS  && ($selected_method_file ne "NONE"));

	    #Run the concordance comparision analysis
	    #run_exe("$script_dir/pair_wise_comparision.pl $data_set_result_dir $data_set_data_dir $out_dir CANCER NO_PLOT");
	    if($FLAG_PLOT_CONCORDANCE_ANALYSIS){
		print STDERR " *** Run the concordance comparision analysis ...\n";
		run_exe("$script_dir/cancer_gene_concordance.pl $data_set_result_dir $data_set_data_dir CANCER_UNION $out_dir NO_PLOT NONE 1 NEW_RESULT $script_dir") if($selected_method_file eq "NONE");#<STDIN>;
		run_exe("$script_dir/cancer_gene_concordance.pl $data_set_result_dir $data_set_data_dir CANCER_UNION $out_dir NO_PLOT $selected_method_file 1 APPEND_RESULT $script_dir") if($selected_method_file ne "NONE");
	    }
	}
	#last;
    }
}
#exit(0);

#############################################################################################################################################################################################

if(0){
#False positive analysis
my @noise_rate_FP = (0.02, 0.05, 0.1);
#my @rank_FP = (3, 8, 15);
#my @rank_FP = (1, 5, 15);
my @rank_FP = (5, 10, 5000);
my $nb_replicate_FP = 3;

my $FP_analysis_script = "$script_dir/FPanalysis.py";
my $best_sample_prediction_extraction_script = "$script_dir/extract_best_sample_prediction.pl";
my $FP_analysis_dir = "$analysis_dir/../FP_ANALYSIS/";

my $count=0;
foreach my $nm (@new_method) {
	$new_method[$count]="$nm.result";
	$count++;
}

my $newmtd = join( "__", @new_method );




print STDERR " *** Analysis of $FP_analysis_dir\n";

opendir(DIR, $FP_analysis_dir) || die "can't opendir $FP_analysis_dir: $!";
@all_data_set = readdir(DIR);
close(DIR);

my @DATA_SET_ORDER_FP;

foreach $data_set (@all_data_set){

    next if($data_set eq "." || $data_set eq ".." ||
            $data_set eq "READ" ||
            #$data_set eq "BRCA" ||
            #$data_set eq "LUAD" ||
            #$data_set eq "KIRP" ||
            #$data_set eq "LUSC" ||
            #$data_set eq "STAD" ||
            #$data_set eq "THCA" ||
            $data_set eq "OV" ||
            $data_set eq "COAD"
            #$data_set eq "LIHC" ||
            #$data_set eq "LUSC"
        );#Data set to exclude

    #print STDERR " *** ANALYSE $data_set\n";<STDIN>;

    push(@DATA_SET_ORDER_FP, $data_set);
    #
    $real_data_set_dir = "$analysis_dir/../DATA/$data_set/";
    $real_data_set_result_dir = "$analysis_dir/../RESULTS/$data_set/";
    #
    if($FLAG_ANALIZE_FP_DATA_SET){
        for(my $cmp_param = 0; $cmp_param < @noise_rate_FP; $cmp_param++){
            $n_rate = $noise_rate_FP[$cmp_param];
            for(my $rep = 0; $rep < $nb_replicate_FP; $rep++){
                #
                $sim_data_set_dir = "$FP_analysis_dir/$data_set/noise_rate_$n_rate/REP_$rep/DATA";
                #$sim_data_set_result_dir = "$sim_data_set_dir/../ANALYSIS/CONSOLIDATED_RESULTS/$CONSOLIDATED_DATA_SET_ID/";
                $sim_data_set_result_dir = "$sim_data_set_dir/../RESULTS/";
		print STDERR "\n *** *** ANALYSE DIRECTORIES:\n\t$sim_data_set_dir\n\t$sim_data_set_result_dir\n";
                #
                run_exe("mkdir -p $sim_data_set_result_dir");
                if(-e "$sim_data_set_result_dir"){

                    #Analysis of the directory
                    $out_dir = "$sim_data_set_result_dir/PLOT";
                    if( $FLAG_ANALIZE_FP_DATA_SET){# || ! -d $out_dir){
                        #if( ! -d $out_dir){
                        run_exe("mkdir -p $out_dir") if(! -d $out_dir);

                        #Problem with actideDriver run
                        if(! -e "$sim_data_set_result_dir/ActiveDriver_PTM.result"){
                            print STDERR " *** Need to construt an empty file: $sim_data_set_result_dir/ActiveDriver_PTM.result\n";
                            run_exe("touch $sim_data_set_result_dir/ActiveDriver_PTM.result");#<STDIN>;
                        }

                        #run_exe("$FP_analysis_script $sim_data_set_dir/GDAC_somatic_mutations.filtered.maf $real_data_set_result_dir $sim_data_set_result_dir $out_dir/FP_rate.dat $out_dir/FDR_rate.dat $newmtd $out_dir/FP_num.dat > /dev/null");
                        #For the rank selected
                        foreach $rank (@rank_FP){
                            run_exe("$best_sample_prediction_extraction_script ".(substr($sim_data_set_result_dir, 0, -1))." $sim_data_set_dir/gene_mutation_frequency.txt $rank NONE $script_dir") if($selected_method_file eq "NONE");
			    run_exe("$best_sample_prediction_extraction_script ".(substr($sim_data_set_result_dir, 0, -1))." $sim_data_set_dir/gene_mutation_frequency.txt $rank $selected_method_file $script_dir") if($selected_method_file ne "NONE");
                            run_exe("$FP_analysis_script $sim_data_set_dir/GDAC_somatic_mutations.filtered.maf $real_data_set_result_dir ".(substr($sim_data_set_result_dir, 0, -1))." $out_dir/FP_rate_$rank.dat $out_dir/FDR_rate_$rank.dat $newmtd $out_dir/FP_num_$rank.dat> /dev/null");
                            #run_exe("rm -r ".(substr($sim_data_set_result_dir, 0, -1))."\_$rank");
                        }
                        #<STDIN>;
                    }
                }
            }
        }
    }
}
#exit(0);


}
#############################################################################################################################################################################################

if(0){
#Downsampling analysis
@noise_rate = (20, 50, 100);
#@noise_rate = (20, 100);
#
@nb_replicate = (10, 3, 3);
#
my $FREQ_COMMON_DRIVER = 0.1;
#
my $downsampling_analysis_script = "$script_dir/ds_analysis.py";
my $downsampling_analysis_dir = "$analysis_dir/../DOWNSAMPLING_ANALYSIS/";

print STDERR " *** Analysis of $downsampling_analysis_dir\n";

opendir(DIR, $downsampling_analysis_dir) || die "can't opendir $analysis_dir: $!";
@all_data_set = readdir(DIR);
close(DIR);

my @DATA_SET_ORDER_DOWNSAMPLING;

foreach $data_set (@all_data_set){

    next if($data_set eq "PRAD" || $data_set eq "." || $data_set eq "..");#Data set to exclude
    #print STDERR " *** ANALYSE $data_set\n";<STDIN>;
    push(@DATA_SET_ORDER_DOWNSAMPLING, $data_set);
    if( $FLAG_ANALIZE_DOWNSAMPLING_DATA_SET){# || ! -d $out_dir){
        #
        $real_data_set_dir = "$analysis_dir/../DATA/$data_set/";
        $real_data_set_result_dir = "$analysis_dir/../RESULTS/$data_set/";
        #

        for(my $cmp_param = 0; $cmp_param < @noise_rate; $cmp_param++){
            $n_rate = $noise_rate[$cmp_param];
            for(my $rep = 0; $rep < $nb_replicate[$cmp_param]; $rep++){
                #
                $sim_data_set_dir = "$downsampling_analysis_dir/$data_set/noise_rate_$n_rate/REP_$rep/DATA";
                #$sim_data_set_result_dir = "$sim_data_set_dir/../ANALYSIS/CONSOLIDATED_RESULTS/$CONSOLIDATED_DATA_SET_ID/";
                $sim_data_set_result_dir = "$sim_data_set_dir/../RESULTS/";
		run_exe("mkdir -p $sim_data_set_result_dir");
		print STDERR "\n *** *** ANALYSE DIRECTORIES:\n\t$sim_data_set_dir\n\t$sim_data_set_result_dir\n";
                #
                if(-e "$sim_data_set_result_dir"){
                    #Analysis of the directory
                    $out_dir = "$sim_data_set_result_dir/PLOT";
                    if( $FLAG_ANALIZE_DOWNSAMPLING_DATA_SET){# || ! -d $out_dir){
                        #if( ! -d $out_dir){
                        run_exe("mkdir -p $out_dir") if(! -d $out_dir);

                        #Problem with actideDriver run
                        if(! -e "$sim_data_set_result_dir/ActiveDriver_PTM.result"){
                            print STDERR " *** Need to construt an empty file: $sim_data_set_result_dir/ActiveDriver_PTM.result\n";
                            run_exe("touch $sim_data_set_result_dir/ActiveDriver_PTM.result");#<STDIN>;
                        }

                        #For the downsampling we are using 2 types of threshold on the gene list
                        #Whole gene list
                        run_exe("$downsampling_analysis_script $sim_data_set_dir/ $real_data_set_dir $sim_data_set_result_dir $real_data_set_result_dir $out_dir/stab_rec.dat 0 $FREQ_COMMON_DRIVER $newmtd> /dev/null");
                        #Recovery: Top 20 genes of the simulated list are present in the full data list
                        #Stability: Top 20 genes of the full data are present in the simulated list
                        run_exe("$downsampling_analysis_script $sim_data_set_dir/ $real_data_set_dir $sim_data_set_result_dir $real_data_set_result_dir $out_dir/stab_rec_20.dat 20 $FREQ_COMMON_DRIVER $newmtd> /dev/null");
                        #<STDIN>;
                    }
                }
            }
        }
    }
}
}
#############################################################################################


my @METHOD_COLOR = ();
my @METHOD_ORDER = ();
my $R_PATH = "R";

#Construst the METHOD_ORDER and METHOD_COLOR global variables !!!
##Plot the boxplot of the number of predictions
plot_nb_call_box_plot("$final_out_dir/driver_number");#<STDIN>;#To change to add the meta

#For the plot that need the baseline data
my @METHOD_ORDER_BL = ($BASE_LINE_FULL_NAME, @METHOD_ORDER);
my @METHOD_COLOR_BL = ($CAT_COLOR{$METHOD_INFO{$BASE_LINE_FULL_NAME}->[0]}, @METHOD_COLOR);
my @METHOD_NAME_BL = ($METHOD_INFO{$BASE_LINE_FULL_NAME}->[1], @METHOD_NAME);

if($FLAG_PLOT_FP_RESULT){
    #
    $ana_file = "FP_rate";
    $out_file = "$final_out_dir/FP_rate";
    combine_FP($FP_analysis_dir, $ana_file, \@noise_rate_FP, $nb_replicate_FP, $out_file);
#    plot_FP($out_file, "False Positive Rate");
    #
    #plot_stab_rc($out_file, "Downsampling analysis\\n$n_rate samples");
    #
    $ana_file = "FDR_rate";
    $out_file = "$final_out_dir/FDR_rate";
    combine_FP($FP_analysis_dir, $ana_file, \@noise_rate_FP, $nb_replicate_FP, $out_file);
#    plot_FP($out_file, "Fraction of False Positives");
    #
    $ana_file = "FP_num";
    $out_file = "$final_out_dir/FP_num";
    combine_FP($FP_analysis_dir, $ana_file, \@noise_rate_FP, $nb_replicate_FP, $out_file);
#    plot_FP($out_file, "");
#    plot_FP("$out_file\_top_x", "");

    #
    if(1){
        $out_file_legend = "$final_out_dir/FP_legend";
#        plot_FP($out_file, $out_file_legend);
#        run_exe("pdfcrop  $out_file_legend.pdf $out_file_legend\_temp.pdf");
#        run_exe("mv  $out_file_legend\_temp.pdf $out_file_legend.pdf");
        #
        $out_file_legend = "$final_out_dir/FP_top_x_legend";
#        plot_FP("$out_file\_top_x", $out_file_legend);
#        run_exe("pdfcrop  $out_file_legend.pdf $out_file_legend\_temp.pdf");
#        run_exe("mv  $out_file_legend\_temp.pdf $out_file_legend.pdf");
    }
}



#plot the measure boxplot
if($FLAG_PLOT_CONCORDANCE_ANALYSIS){
    #
    #$out_file = "$final_out_dir/reverse_CANCER_F1_RANK_50";
    #run_exe("cp $final_out_dir/cancer_gene\_CANCER_F1\_RANK_50.dat $out_file.dat");
    #plot_measure_box_plot_reverse("$out_file", $measure_tested{"F1"}, "Top 50 Predictions");

    #exit(0);
    #For the rank concordance analysis
    combine_rank_concordance_file(50);
    #@method_in_rank_concordance = ("PolyPhen2", "CHASM", "fathmm", "MutSigCV", "DriverNet", "oncoIMPACT");
    #@method_in_rank_concordance = ("PolyPhen2", "fathmm", "MutSigCV", "S2N", "oncoIMPACT");
    #plot_rank_concordance_file("$final_out_dir/cc_precision_RANK_50", \@method_in_rank_concordance, 50);
    #
    #
    #exit(0);
    foreach $m (keys %measure_tested){
        $out_file = "$final_out_dir/cancer_gene";
	#if($m ne "TPR"){
		#if($m eq "F1") {$m = "CS"}
        	combine_sample_measure_file("CANCER_UNION_$m", $out_file);
	#}
	#foreach $rank (@tested_rank){
	    #plot_measure_box_plot("$out_file\_CANCER_UNION_$m\_RANK_".$rank, $measure_tested{$m}, "Top $rank Predictions");
	#}
    }
}

#For the sample based analysis
#For the plot that need the baseline data
my $title = "";
if($FLAG_PLOT_SAMPLE_BASED_ANALYSIS){
    #print STDERR " ***  ***|"."@tested_rank"."|\n";
    foreach $rank (@tested_rank_all){
	#
	$rank = 3 if($rank eq 50);
	#print STDERR " *** $rank ***\n";

	$out_file = "$final_out_dir/sample_nb_driver_cat_RANK_$rank";
	combined_sample_based_category("sample_nb_driver_cat_RANK_$rank.dat", $out_file.".dat") if($rank eq "ALL");

	$title = "Top $rank Predictions";$title = "All Predictions" if($rank eq "ALL");
	

#	plot_bar_plot_sample_file("$out_file", $title, 1);
	
	#For the legend
	if($rank eq "ALL"){
#	    plot_bar_plot_sample_file("$out_file", "$final_out_dir/sample_nb_driver_cat_legend", 1) ;
#	    run_exe("pdfcrop  $out_file.pdf $out_file\_temp.pdf");
#	    run_exe("mv  $out_file\_temp.pdf $out_file.pdf");
	}
	#
	if($rank eq 3 || $rank eq 5 || $rank eq 10){
	    $out_file = "$final_out_dir/sample_precision_RANK_$rank";
	    combined_sample_based_PPV("sample_precision_RANK_$rank.dat", $out_file.".dat");
#	    plot_sample_based_PPV($out_file, $title);
	    #
	    $out_file = "$final_out_dir/sample_F1_RANK_$rank";
	    combined_sample_based_PPV("sample_precision_RANK_$rank.dat_penalised", $out_file.".dat");
	    #
	    $out_file = "$final_out_dir/sample_recall_RANK_$rank";
	    combined_sample_based_PPV("sample_recall_RANK_$rank.dat", $out_file.".dat");
	    #plot_measure_box_plot("sample_precision_RANK_$rank.dat_penalised", $measure_tested{"F1"}, "Top $rank Predictions");
	}
    }
}

if($FLAG_PLOT_SAMPLE_ACTIONABLE_ANALYSIS){
    #print STDERR " ***  ***|"."@tested_rank"."|\n";
    my $dir_actionable_gene = "$final_out_dir/ACTIONABLE_GENE";
    my $out_file;
    run_exe("mkdir -p $dir_actionable_gene") if(! -d $dir_actionable_gene);
    

    foreach $rank (@tested_rank_all){
	if($rank eq 5 || $rank eq 10){
	#if($rank eq 10){
	    $out_file = "$final_out_dir/sample_actionable_profile\_$rank";
	    
	    combined_sample_actionable_profile("sample_actionable_gene_profile_RANK_$rank.dat", $out_file);#Construct the average value for bar plot and matrix for cancer type view
	    #
#	    plot_bar_plot_sample_actionable_gene_file("$out_file\_all", "toto", 1);
	    #WEIRD
	    #$out_file = "$final_out_dir/sample_actionable_profile\_$rank";
#####	    heat_map_actionable_gene_by_cancer_type("$out_file\_cancer_type", "toto", "HIGH");
	    
	    #For the gene heat map
	    $out_file = "$dir_actionable_gene/actionable\_$rank";
	    combined_sample_actionable_gene("sample_actionable_gene_$rank.dat", $out_file);#Construct the average value for bar plot and matrix for cancer type view
	    foreach $m (@METHOD_ORDER){
#		#plot_sample_actionable_gene("$out_file\_$m", $m, 1);
	    }
	    
#	    foreach $ag (@ACTIONABLE_GENE_STUDIED){
#		#heat_map_actionable_gene_by_cancer_type("$out_file\_$ag", $ag, "LOW");
#		plot_bar_plot_sample_actionable_gene_file("$out_file\_$ag\_bar", $ag, 0);
#	    }

	}
    }
}

#Clean the result directory
run_exe("rm $analysis_dir/*/$additional_method_name.result");
run_exe("rm -r $analysis_dir/*/PLOT/");
run_exe("rm -r $analysis_dir/*/GENE/");
run_exe("rm $selected_method_file");

sub combined_sample_actionable_gene{
    my ($file_type, $out_file) = @_;
    #
    #Init the data site
    my %method_best_gene = ();
    foreach $m (@METHOD_ORDER){
        $best_gene_file = "$out_file\_$m\_best_list.dat";
        #run_exe("grep -h $m $analysis_dir/*/CONSOLIDATED_RESULTS/$CONSOLIDATED_DATA_SET_ID/PLOT/$file_type | cut -f1-3  | sort | uniq -c| awk '{print \$2\"\t\"\$3\"\t\"\$4\"\t\"\$1}' | sort -k4,4 -nr | head -n20 | cut -f2,3 > $best_gene_file");
        run_exe("grep -h $m $analysis_dir/*/PLOT/$file_type | cut -f1-3  | sort | uniq -c| awk '{print \$2\"\t\"\$3\"\t\"\$4\"\t\"\$1}' | sort -k4,4 -nr | head -n20 | cut -f2,3 > $best_gene_file");
	$method_best_gene{$m} = {};
        $method_best_gene{$m}->{"ORDER"} = ();
        open(FILE, $best_gene_file);
        while(<FILE>){
            chop $_;
            @line = split(/\t/, $_);
            $gene = $line[0];
            $method_best_gene{$m}->{$gene} = {};
            push(@{$method_best_gene{$m}->{"ORDER"}}, $gene);
        }
        close(FILE);
    }

    #Read all the gene file
    my %data_set_nb_sample;
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
#        $data_set_plot_dir = "$analysis_dir/$data_set/CONSOLIDATED_RESULTS/$CONSOLIDATED_DATA_SET_ID/PLOT";
	$data_set_plot_dir = "$analysis_dir/$data_set/PLOT";

        if(-e $data_set_plot_dir){
            #get the number of sample
            $line_with_sample = `head -n 2 $data_set_plot_dir/sample_actionable_gene_profile_RANK_5.dat | tail -n 1`;
            @tmp = split(/\t/, $line_with_sample);
            $data_set_nb_sample{$data_set} = $tmp[1] + $tmp[2] + $tmp[3] + $tmp[4];

            open(FILE, "$data_set_plot_dir/$file_type");
            #print STDERR " *** Read file $data_set_plot_dir/$file_type\n";<STDIN>;
            <FILE>;#Skip the header
            while(<FILE>){
                chop $_;
                @line = split(/\t/, $_);
                $m = $line[0];
                $gene = $line[1];
                if(exists $method_best_gene{$m}->{$gene}){
                    if(! exists $method_best_gene{$m}->{$gene}->{$data_set}){
                        $method_best_gene{$m}->{$gene}->{$data_set} = 0;
                    }
                    $method_best_gene{$m}->{$gene}->{$data_set}++;
                }
            }
            close(FILE);
        }
    }

    #Get the method heat map file
    my $nb_sample;
    foreach $m (@METHOD_ORDER){
        $best_gene_file = "$out_file\_$m.dat";
        open(OUT, ">$best_gene_file");
        print OUT "".join("\t", @DATA_SET_ORDER)."\n";
        foreach $g (@{$method_best_gene{$m}->{"ORDER"}}){
            print OUT $g;
            for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
                $data_set = $DATA_SET_ORDER[$cmp_data];
                $nb_sample = 0;
                $nb_sample = $method_best_gene{$m}->{$g}->{$data_set} if(exists $method_best_gene{$m}->{$g}->{$data_set});
                print OUT "\t".$nb_sample/$data_set_nb_sample{$data_set};
            }
            print OUT "\n";
        }
        close(OUT);
    }
}



sub heat_map_actionable_gene_by_cancer_type{
    ($matrix_file, $title, $plot_type) = @_;
    my $color_list = join(",", @METHOD_COLOR);
    my $font_size = 1.5;
    my $title_font_size = 3;
    my $note_font_size = 4;
    my $margin_size = 30;

    open(OUT, ">$matrix_file.R");

    print OUT "pdf(file=\"$matrix_file.pdf\",
        paper=\"special\",
        width=7,
        height=7
        )\n";

    print OUT "library(\"gplots\")\n";
    print OUT "pairs.breaks <- seq(0, 1, by=0.01)\n" if($plot_type eq "HIGH");
    print OUT "pairs.breaks <- seq(0, 0.8, by=0.008)\n" if($plot_type eq "LOW");

    #Choice of the color
    #print OUT "palette = c( \"black\", colorRampPalette(c(\"darkred\", \"indianred\"))(5), colorRampPalette(c(\"indianred\", \"lightblue\"))(15), colorRampPalette(c(\"lightblue\", \"darkblue\"))(30), colorRampPalette(c(\"darkblue\", \"darkgreen\"))(49))\n";# if($plot_type eq "LOW");
    print OUT "palette = c( \"black\", colorRampPalette(c(\"darkred\", \"indianred\"))(19), colorRampPalette(c(\"indianred\", \"lightblue\"))(30), colorRampPalette(c(\"lightblue\", \"darkblue\"))(50))\n";# if($plot_type eq "HIGH");

    print OUT "profile_mat = read.table(\"$matrix_file.dat\")\n";

    print OUT "
        heatmap.2(as.matrix(profile_mat),
                  main = \"$title\",
                  cex.main = $title_font_size,
                  scale=\"none\",

                  #For the colors
                  ColSideColors = c($color_list),
                  #RowSideColors = c($color_list),
                  #Rowv=FALSE,
                  Colv=FALSE,
                  #breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), #symbreaks = TRUE,
#
#key.par=list(mgp=c(1.5, 0.5, 0),
#mar=c(2.5, 2.5, 1, 0)),

                  #For the tree
                  dendrogram = \"row\",
                  hclustfun = function(x) hclust(x,method = 'complete'),
                  distfun = function(x) dist(x,method = 'euclidean'),
                                  #cellnote=as.matrix(profile_mat),notecol=\"black\",notecex=$note_font_size,
                  #margins =c(5,5),
                  #margin=c($margin_size, $margin_size),
                  col=palette, breaks = pairs.breaks, cexRow=$font_size, cexCol=$font_size,
                  density.info=\"none\", trace=\"none\"";
    #print OUT ",ColSideColors = $color_subtype" if($subtype_file ne "NONE");
    print OUT ")\n";
    close(OUT);

    run_exe("$R_PATH --vanilla < $matrix_file.R");
}



#Construct the average value for bar plot and matrix for cancer type view
sub combined_sample_actionable_profile{
    ($file_type, $out_file) = @_;
    #
    #Init the data site
    my %method_profile = ();
    foreach $m (@METHOD_ORDER){
        $method_profile{$m} = {};
        $method_profile{$m}->{"ALL"} =  [0, 0, 0, 0, 0];
    }
    #Read the file
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
#        $data_set_plot_dir = "$analysis_dir/$data_set/CONSOLIDATED_RESULTS/$CONSOLIDATED_DATA_SET_ID/PLOT";

	 $data_set_plot_dir = "$analysis_dir/$data_set/PLOT";

        if(-e $data_set_plot_dir){

            open(FILE, "$data_set_plot_dir/$file_type");
            <FILE>;#Skip the header
            while(<FILE>){
                chop $_;
                @line = split(/\t/, $_);
                $method = $line[0];
                $method_profile{$method}->{"$data_set"} =  [0, 0];
                for(my $i = 1; $i < @line-1; $i++){
                    $method_profile{$method}->{"ALL"}->[$i-1] += $line[$i];$method_profile{$method}->{"ALL"}->[@line-1] += $line[$i];
                    $method_profile{$method}->{$data_set}->[0] += $line[$i];
                }
                if(!($data_set eq "BRCA" && $method eq "DawnRank")){
                    $method_profile{$method}->{"ALL"}->[@line-1-1] += $line[@line-1];
                    $method_profile{$method}->{"ALL"}->[@line-1] += $line[@line-1];
                }
                $method_profile{$method}->{$data_set}->[1] += $line[@line-1];
            }
            close(FILE);
        }
    }

    #Write the average file
    open(OUT_ALL, ">$out_file\_all.dat");
    #print OUT_ALL "1\t2\t3\t0\n";
    print OUT_ALL "".(join("\t", @METHOD_NAME))."\n";
    for (my $cmp_cat = 0; $cmp_cat < 4; $cmp_cat++){
        print OUT_ALL $cmp_cat;
        for (my $cmp_method = 0; $cmp_method < @METHOD_ORDER; $cmp_method++){
            $method = $METHOD_ORDER[$cmp_method];
            #
            #print STDERR $method."\n";<STDIN>;
            #
            print OUT_ALL "\t".$method_profile{$method}->{"ALL"}->[$cmp_cat]/$method_profile{$method}->{"ALL"}->[@line-1];
        }
        print OUT_ALL "\n";
    }

    #Write the cancer type file
    open(OUT_TYPE, ">$out_file\_cancer_type.dat");
    print OUT_TYPE "".(join("\t", @METHOD_NAME))."\n";
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
        print OUT_TYPE $data_set;
        for (my $cmp_method = 0; $cmp_method < @METHOD_ORDER; $cmp_method++){
            $method = $METHOD_ORDER[$cmp_method];
            print OUT_TYPE "\t".($method_profile{$method}->{$data_set}->[0]/($method_profile{$method}->{$data_set}->[1]+$method_profile{$method}->{$data_set}->[0]));
        }
        print OUT_TYPE "\n";
    }
}



sub combine_FP{
    ($dir, $ana_file, $noise_rate_list, $nb_replicate, $out_file) = @_;

    #VERY BAD, NEED TO CHANGE THAT
    my %data_set_size = ();
    my $total_sample_size = 0;
    #wc -l DATA/*/complete_samples_list.txt | tr '/' '\t' | awk '{if($2 != "total") print $3"\t"$1}' >
    open(FILE, "DATA/data_set_size.dat");
    while(<FILE>){
        chop $_;
        @line = split(/\t/, $_);
        $nb_sample = $line[1];
        $data_set_size{$line[0]} = $nb_sample;
        $total_sample_size += $nb_sample;
    }
    close(FILE);

    print STDERR " *** combine_FP $ana_file ".(join(",", @{$noise_rate_list}))." -> $out_file\n";#<STDIN>;
    my $noise_rate;
    my %method_res = ();
    foreach $m (@METHOD_ORDER){
        my %map;
        $method_res{$m} = \%map;
    }

    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER_FP; $cmp_data++){
        $data_set = $DATA_SET_ORDER_FP[$cmp_data];
        foreach $m (keys %method_res){
            my %map = ();
            $method_res{$m}->{$data_set} = \%map;
        }
        for(my $cmp_noise_rate = 0; $cmp_noise_rate < @{$noise_rate_list}; $cmp_noise_rate++){
            $noise_rate = $noise_rate_list->[$cmp_noise_rate];
            foreach my $m (keys %method_res){
                my @tab = ();
                for($i = -1; $i < @rank_FP; $i++){push(@tab, 0);}
                $method_res{$m}->{$data_set}->{$noise_rate} = \@tab;
            }

            for(my $cmp_rank = -1; $cmp_rank < @rank_FP; $cmp_rank++){
                for(my $rep = 0; $rep < $nb_replicate; $rep++){
                    #
                    #$rep_result_file = "$dir/$data_set/noise_rate_$noise_rate/REP_$rep/ANALYSIS/CONSOLIDATED_RESULTS/$CONSOLIDATED_DATA_SET_ID/PLOT/$ana_file";
                    $rep_result_file = "$dir/$data_set/noise_rate_$noise_rate/REP_$rep/RESULTS/PLOT/$ana_file";
		    $rep_result_file .= "_$rank_FP[$cmp_rank]" if($cmp_rank != -1);
                    $rep_result_file .= ".dat";
                    #print STDERR " *** Read result file: $rep_result_file\n";<STDIN>;

                    open(FILE, $rep_result_file);
                    <FILE>;#To get rid of the header
                    while(<FILE>){
                        chop $_;
                        @line = split(/\t/, $_);
                        $method = $line[0];
                        #
                        next if(! exists $METHOD_ANALYSED_FP{$method});
                        #
                        $val = $line[1];
                        #$method_res{$method}->{$data_set}->{$noise_rate}->[$cmp_rank+1] += $val;
                        $method_res{$method}->{$data_set}->{$noise_rate}->[$cmp_rank+1] += $val * $data_set_size{$data_set};
                        #print STDERR $method." ".$data_set_size{$data_set}." ".$val." ".$method_res{$method}->{$data_set}->{$noise_rate}->[$cmp_rank+1]."\n";<STDIN>;
                        }
                    close(FILE);
                    #
                }
            }
        }
    }

    #Write the average file
    open(OUT, ">$out_file.dat");
    print OUT "".(join("\t", @{$noise_rate_list}))."\n";

    open(OUT_TOP, ">$out_file\_top_x.dat");
    print OUT_TOP "".(join("\t", @{$noise_rate_list}))."\n";

    #print OUT "stability\trecovery\trecovery_rare\trecovery_common\n";
    my $nb_data_set = @DATA_SET_ORDER_FP;
    my $nb_noise_rate = @{$noise_rate_list} + 0;
    my $nb_rank = @rank_FP;
    print STDERR "nb_data_set:$nb_data_set\tnb_replicate:$nb_replicate_FP\n";


    my @avg_res;
    for(my $cmp_method = 0; $cmp_method < @METHOD_ORDER; $cmp_method++){
        $method = $METHOD_ORDER[$cmp_method];
        #
        next if(! exists $METHOD_ANALYSED_FP{$method});
        #
#####        print STDERR " *** method analysised $method\n";
        @avg_res = ();
        for($i = -$nb_noise_rate; $i < $nb_noise_rate*@rank_FP; $i++){push(@avg_res, 0);}

        for($cmp_rank = 0; $cmp_rank <= $nb_rank; $cmp_rank++){
            for(my $cmp_noise_rate = 0; $cmp_noise_rate < $nb_noise_rate; $cmp_noise_rate++){
                $noise_rate = $noise_rate_list->[$cmp_noise_rate];
                for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER_FP; $cmp_data++){
                    $data_set = $DATA_SET_ORDER_FP[$cmp_data];
                    #$val = $method_res{$method}->{$data_set}->{$noise_rate}->[$cmp_rank]/$nb_replicate;
                    $val = $method_res{$method}->{$data_set}->{$noise_rate}->[$cmp_rank];#/$nb_replicate;
                    #print STDERR "".($cmp_noise_rate + ($cmp_rank*$nb_rank))." -> ".$val." ".$avg_res[$cmp_noise_rate + $cmp_rank]."\n";
                    $avg_res[$cmp_noise_rate + ($cmp_rank*$nb_rank)] += $val;
                }
                #<STDIN>;
                #$avg_res[$cmp_noise_rate] += ($method_res{$method}->{$data_set}->{$noise_rate}->[0]/$nb_replicate);
                #$avg_res[$i+3] += ($method_res{$method}->{$data_set}->{$noise_rate}->[1]/$nb_replicate);
            }
        }

        for(my $i = 0; $i < @avg_res; $i++){
            #$avg_res[$i] = $avg_res[$i]/$nb_data_set;
            $avg_res[$i] = $avg_res[$i]/($nb_replicate * $total_sample_size);
        }
        #print OUT $METHOD_INFO{$method}->[1]."\t".(join("\t", @avg_res))."\n";
        #print OUT $METHOD_INFO{$method}->[1]."\t".(join("\t", @avg_res[0..2]))."\n";
        print OUT $METHOD_INFO{$method}->[1]."\t".(join("\t", @avg_res[9..11]))."\n";
        #print OUT $METHOD_INFO{$method}->[1]."\t".(join("\t", @avg_res[16..18]))."\n";
        #print OUT $METHOD_INFO{$method}->[1]."\t".(join("\t", @avg_res))."\n";
        #print OUT_TOP $METHOD_INFO{$method}->[1]."\t".$avg_res[4]."\t".$avg_res[7]."\t".$avg_res[10]."\t".$avg_res[1]."\n";
        #print OUT_TOP $METHOD_INFO{$method}->[1]."\t".$avg_res[5]."\t".$avg_res[8]."\t".$avg_res[11]."\t".$avg_res[14]."\n";
        print OUT_TOP $METHOD_INFO{$method}->[1]."\t".$avg_res[5]."\t".$avg_res[8]."\t".$avg_res[11]."\n";
    }
    close(OUT);
    close(OUT_TOP);
}


sub combined_sample_based_PPV{
    ($file_type, $out_file) = @_;

    #print STDERR " *** combined_sample_based_precision\n";

    open(OUT, ">$out_file");
    #Write the header
    print OUT join("\t", @METHOD_NAME)."\n";;
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
        $data_set_plot_dir = "$analysis_dir/$data_set/PLOT";

        if(-e $data_set_plot_dir){
            open(FILE, "$data_set_plot_dir/$file_type");
            $first = 1;
            <FILE>;#Skip the header
            print OUT  <FILE>;
            close(FILE);
        }
    }
    close(OUT);
}


sub combined_sample_based_category{
    ($file_type, $out_file) = @_;

    print STDERR " *** combined_sample_based_category\n";

    my $first_dataset =1;
    %method_res = ();
    my @category_order = ();
    my @nb_sample_analysed = ();for(my $cmp_method = 0; $cmp_method < @METHOD_ORDER; $cmp_method++){$nb_sample_analysed[$cmp_method] = 0};
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
        $data_set_plot_dir = "$analysis_dir/$data_set/PLOT";
        if(-e $data_set_plot_dir){
            open(FILE, "$data_set_plot_dir/$file_type");
            $first = 1;
            <FILE>;#Skip the header
            while(<FILE>){
                chop $_;
                @line = split(/\t/, $_);
                $category = $line[0];
                push(@category_order, $category) if($first_dataset);

                #Initialize the method result table
                for(my $cmp_method = 0; $cmp_method < @METHOD_ORDER; $cmp_method++){
                    $method = $METHOD_ORDER[$cmp_method];
                    $nb_sample_in_catagory = $line[$cmp_method+1];
                    if(! exists $method_res{$method}){
                        my %map = ();$method_res{$method} = \%map;
                    }
                    if(! exists $method_res{$method}->{$category}){
                        $method_res{$method}->{$category} = 0;
                    }
		    if($nb_sample_in_catagory ne "NA"){
			$method_res{$method}->{$category} += $nb_sample_in_catagory;
			$nb_sample_analysed[$cmp_method] += $nb_sample_in_catagory;# if($cmp_method == @METHOD_ORDER-1);
		    }
                }
            }
            $first_dataset = 0;
        }
        else{
            print STDERR " *** WARNING $data_set WITHOUT PLOTS\n";
        }
    }

    #Construt the combine  file and the categ file
    open(OUT, ">$out_file");
    #Write the header
    print OUT join("\t", @METHOD_NAME)."\n";;
    for(my $cmp_cat = 0; $cmp_cat < @category_order; $cmp_cat++){
        $cat = $category_order[$cmp_cat];
        print OUT $cat;
        for(my $cmp_method = 0; $cmp_method < @METHOD_ORDER; $cmp_method++){
            $method = $METHOD_ORDER[$cmp_method];
	    #print STDERR $method." ->".$nb_sample_analysed[$cmp_method]."\n";
            print OUT "\t".(sprintf("%.4f",($method_res{$method}->{$cat}/$nb_sample_analysed[$cmp_method])));
        }
        print OUT "\n";
    }
    close(OUT);


}

sub combine_rank_concordance_file{
    my ($rank, $specific_method_order) = @_;

    $specific_method_order = \@METHOD_ORDER if(! defined $specific_method_order);

    my %method_res = ();
    my %method_res_nb_data_set = ();
    foreach $m (@{$specific_method_order}){
        my @tab = ();for (my $i = 1; $i <= $rank; $i++){$tab[$i] = 0;}
        $method_res{$m} = \@tab;
        my @tab1 = ();for (my $i = 1; $i <= $rank; $i++){$tab1[$i] = 0;}
        $method_res_nb_data_set{$m} = \@tab1;
    }
    my $nb_data_set;
    my $curr_rank;
    my $curr_concord;
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
        #
        foreach $m (@{$specific_method_order}){
            $result_file = "$analysis_dir/$data_set/PLOT/$m\_cc_precision_RANK_$rank.dat";
            #print STDERR " *** Read result file: $result_file\n";<STDIN>;
            open(FILE, $result_file);
            while(<FILE>){
                chop $_;
                @line = split(/\t/, $_);
                $curr_rank = $line[0];
                $curr_concord = $line[1];
                $method_res{$m}->[$curr_rank] += $curr_concord;
                $method_res_nb_data_set{$m}->[$curr_rank]++;
            }
            #To complte value in case of method having less the 50 calls
            #The final concordance is given to all the remaining ranks
            #for(my $i = $curr_rank+1; $i <= $rank; $i++){
            #$method_res{$m}->[$i] += $curr_concord;
            #}
        }
        close(FILE);
    }

    foreach $m (@{$specific_method_order}){
        open(OUT, ">$final_out_dir/$m\_cc_precision_RANK_$rank.dat");
        for (my $i = 1; $i <= $rank; $i++){
            $nb_data_set = $method_res_nb_data_set{$m}->[$i];
            last if($nb_data_set == 0 ||
                    #80% of the data set they analyzed
                    ($nb_data_set < 12 && $m eq "MutSigCV") ||
                    ($nb_data_set < 6 && $m eq "MutSig_meta")
                );
            print OUT $i."\t".($method_res{$m}->[$i]/$nb_data_set)."\n";
        }
        close(OUT);
    }

}



sub combine_sample_measure_file{

    my ($file_type, $out_file, $specific_method_order, $specific_method_name) = @_;

    $specific_method_order = \@METHOD_ORDER_BL if(!defined $specific_method_order);
    $specific_method_name = \@METHOD_NAME_BL if(!defined $specific_method_name);

    #print STDERR " *** combine_sample_measure_file $file_type $out_file\n";<STDIN>;

    my %method_res = ();
    #Read the result file
    my @rank_tested = ();
    for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
        $data_set = $DATA_SET_ORDER[$cmp_data];
        $data_set_plot_dir = "$analysis_dir/$data_set/PLOT";
        #print STDERR " *** $data_set_plot_dir\n";<STDIN>;
        if(-e $data_set_plot_dir){
            #push(@DATA_SET_ORDER, $data_set);
            open(FILE, "$data_set_plot_dir/cancer_gene_$file_type.dat");
            $first = 1;
            while(<FILE>){
                chop $_;
                @line = split(/\t/, $_);
                #To get the rank tested at the first read
                if($first == 1){
                    @rank_tested = @line;
                    #print STDERR " *** RANK_TESTED: @rank_tested\n";<STDIN>;
                    $first = 0;
                }

                #Initialize the method result table
                $method = $line[0];
                if(! exists $method_res{$method}){
                    my @tab_1 = ();for(my $k = 0; $k < @rank_tested; $k++){my @tab_2 = (); $tab_1[$k] = \@tab_2;}
                    $method_res{$method} = \@tab_1;
                }
                for(my $cmp_rank = 0; $cmp_rank < @rank_tested; $cmp_rank++){
                    #Update the value
                    push(@{$method_res{$method}->[$cmp_rank]}, $line[$cmp_rank+1]);
                }
            }
        }
        else{
            print STDERR " *** WARNING $data_set WITHOUT PLOTS\n";
        }
    }

    #Construt the combine concordance evalutation file for each rank tested
    for(my $cmp_rank = 0; $cmp_rank < @rank_tested; $cmp_rank++){
	if($rank_tested[$cmp_rank] eq "10" || $rank_tested[$cmp_rank] eq "50"){
	    
        open(OUT, ">$out_file\_$file_type\_RANK_".$rank_tested[$cmp_rank].".dat");
        print OUT "".join("\t", @{$specific_method_name})."\n";

        for(my $cmp_data = 0; $cmp_data < @DATA_SET_ORDER; $cmp_data++){
            print OUT $DATA_SET_ORDER[$cmp_data];

            for(my $cmp_method = 0; $cmp_method < @{$specific_method_order}; $cmp_method++){
                $method = $specific_method_order->[$cmp_method];
                $res = $method_res{$method}->[$cmp_rank]->[$cmp_data];
                $res = "NA" if(! defined $res);
                print OUT "\t".$res;
            }
            print OUT "\n";
        }
        close(OUT);
	
	}
    }
}



sub plot_nb_call_box_plot{
    my ($out_file) = @_;
    my %method_res = ();
    my $cmp_meta = 1;
    #Read the result file
    for(my $i = 0; $i < @DATA_SET_ORDER; $i++){
        $data_set = $DATA_SET_ORDER[$i];
        $data_set_plot_dir = "$analysis_dir/$data_set/PLOT";
        if(-e $data_set_plot_dir){
            #push(@DATA_SET_ORDER, $data_set);
            open(FILE, "$data_set_plot_dir/driver_number.dat");
            while(<FILE>){
                chop $_;
                @line = split(/\t/, $_);
                $method = $line[0];
                $call = $line[1];
                if(! exists $method_res{$method}){
                    my @tab = ();
                    $method_res{$method} = \@tab;
                    push(@METHOD_ORDER, $method);

                    #update the method info for the meta method
                    if(! exists $METHOD_INFO{$method}){
        ###                my @tab_info = ("META_METHOD", "MT_$cmp_meta");
       			my @tab_info = ("META_METHOD", "$new_method[$cmp_meta-1]");
	                 $METHOD_INFO{$method} = \@tab_info;
                        $cmp_meta++;
                    }

                    push(@METHOD_NAME, $METHOD_INFO{$method}->[1]);
                    push(@METHOD_COLOR, $CAT_COLOR{$METHOD_INFO{$method}->[0]});
                }
                push(@{$method_res{$method}}, $call);
            }
        }
        else{
            print STDERR " *** WARNING $data_set WITHOUT PLOTS\n";
        }
    }

    #Construct the file to plot
    open(OUT, ">$out_file.dat");
    print OUT "".(join("\t", @METHOD_NAME))."\n";
    for(my $j = 0; $j < @DATA_SET_ORDER; $j++){
        print OUT $DATA_SET_ORDER[$j];
        for(my $i = 0; $i < @METHOD_ORDER; $i++){
            print OUT "\t".$method_res{$METHOD_ORDER[$i]}->[$j];
        }
        print OUT "\n";
    }
    close(OUT);

#    #plot the data
#    open(OUT, ">$out_file.R");
#    my $font_size = 1.5;
#    my $note_font_size = 3;
#    my $margin_size = 30;
#
#    open(OUT, ">$out_file.R");
#    print OUT "pdf(file=\"$out_file.pdf\",
#        paper=\"special\",
#        width=18,
#        height=7
#        )\n";
#    print OUT "profile <- read.table(\"$out_file.dat\", header = TRUE)\n";
#    print OUT "l10_profile = log10(profile+1)\n";
#
#    print OUT "par(mar = c(5,10,4,2)+0.1)\n";
#    print OUT "par(mgp=c(3, 2.5, 0))\n";
#
#
#    print OUT "palette <- c(".join(",", @METHOD_COLOR).")\n";
#    print OUT "
#        boxplot(
#            as.matrix(l10_profile),
#            col = \"white\",
#            ylim = c(0,max(l10_profile)),
#            border = palette,
#            boxwex = 0.5,
#            medlwd = 5,
#            boxlwd = 5,
#            outpch = NA,
#            axes = FALSE,
#            )\n";
#
#    $str_list = ();
#    for(my $i = 0; $i < @METHOD_ORDER; $i++){
#        $str_list .= "l10_profile\$$METHOD_NAME[$i],"
#    }
#    chop $str_list;
#
#    print OUT "stripchart(
#        list($str_list),
#        vertical = TRUE,
#        method = \"jitter\",
#        pch = 21,
#        col = \"black\",
#        bg = \"black\",
#        cex = 1.5,
#        add = TRUE)\n";
#
#    ## add axes
#    my $str = "tick_pos = c(";for(my $i= 0; $i < @METHOD_ORDER; $i++){$str .= ($i+1).",";}chop $str;$str .= ")";
#    print OUT $str."\n";
#    #print OUT "".(plot_color_axis(\@METHOD_COLOR, \@METHOD_NAME, $font_size))."\n";
#    #print OUT "axis(1, at = 1:".(@METHOD_ORDER+0).", labels = c(\"".(join("\",\"",@METHOD_NAME))."\"), cex.axis = $font_size)\n";
#    print OUT "axis(2, at = 0:4, labels = c(0,10, 100, 1000, 10000), cex.axis = $font_size)\n";
#    ## now draw the y-axis annotation on a different line out from the plot
#    ## using the extra margin space:
#    print OUT "title(ylab = \"Number of Driver Predicted\", cex.lab = $font_size, line = 7)\n";
#
#
#    close(OUT);
#    run_exe("$R_PATH --vanilla < $out_file.R");

}

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    #print STDERR $exe."\n";#<STDIN>;
    print STDERR `$exe` if($run);
}

