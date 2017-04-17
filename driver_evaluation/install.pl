#!/usr/bin/perl
use warnings;
use Cwd;

my $script_dir = getcwd;

print STDERR " *** Install driver_prediction\n";
`sed -i 's|YYY_SCRIPT_DIR|'$script_dir/bin'|' bin/driver_evaluation.pl`;
#run_exe("sed -i 's|YYY_SCRIPT_DIR|'$script_dir'|' bin/driver_evaluation.pl");
`sed -i 's|XXX_ANALYSIS_DIR|'$script_dir/EVALUATION_DATA_SET/RESULTS'|' bin/driver_evaluation.pl`;

print STDERR " *** Download the genomic and prediction data ... Please wait\n";
run_exe("wget ftp://ftp2.gis.a-star.edu.sg/ConsensusDriver/evaluation_data_set.tar.gz");
print STDERR " *** Extract the genomic and prediction data ... Please wait\n";
run_exe("tar -xvzf evaluation_data_set.tar.gz");
run_exe("rm evaluation_data_set.tar.gz");

sub run_exe{
    my ($exe) = @_;
    $run = 0;
    print STDERR $exe."\n";#<STDIN>;
    print STDERR `$exe` if($run);
}



