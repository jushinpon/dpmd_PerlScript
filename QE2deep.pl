use warnings;
use strict;
use Data::Dumper;
my $output = './output';
`rm -rf $output`;
`mkdir $output`;

open my $fileID ,"> $output/00file_summary.txt";#keep the corresponding fileid and filename
print $fileID "***fileID  => filename***\n\n";
my $ry2eV = 13.605693009;
my $bohr2ang = 0.52917721067;
my $kbar2bar = 1000.;
my $kbar2evperang3 = 1e3 / 1.602176621e6;
my $force_convert = $ry2eV / $bohr2ang;

my @out = <./*.out>;# all slurm or qe output files
@out = sort @out;
for my $id (0..$#out){
	open my $all ,"< $out[$id]";
	my @all = <$all>;
	close($all);
my $natom = `cat $out[$id]|sed -n '/number of atoms\\/cell/p'|awk '{print \$5}'`;	
chomp $natom;
if(!$natom){die "You don't get the Atom Number!!!\n";}

#open required output files	
	my $temp = sprintf("%04d", $id + 1);
	print $fileID "$temp $out[$id]\n";	
	open my $eraw ,"> $output/energy_$temp.raw";
	open my $vraw ,"> $output/virial__$temp.raw";
	open my $fraw ,"> $output/force_$temp.raw";
	open my $craw ,"> $output/coord_$temp.raw";
	open my $braw ,"> $output/box_$temp.raw";
	
################# energy ############
#!    total energy              =    (-158.01049803) Ry
	my @totalenergy = grep {if(m/^\s*!\s*total energy\s*=\s*([-+]?\d*\.?\d*)/){$_ = $1*$ry2eV;}} @all;
    for (0..$#totalenergy){
		print $eraw "$totalenergy[$_]";
		print $eraw "\n" if ($_ ne $#totalenergy);
	}
###virial (kbar)
#   0.00000058  -0.00000001  -0.00000003            (0.09)       (-0.00)       (-0.00)
	my @virial = grep {if(m/^\s*[-+]?\d+\.?\d+\s+[-+]?\d+\.?\d+\s+[-+]?\d+\.?\d+
		\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)/x){
			$_ = [$1*$kbar2bar,$2*$kbar2bar,$3*$kbar2bar];}} @all;
	for (1..@virial){
		chomp @{$virial[$_ -1]}[0..2];
		print $vraw "@{$virial[$_ -1]}[0..2] ";
		print $vraw "\n" if ($_% 3 == 0 and $_ ne scalar @virial);
	}
############## force ############   #Ry/au
##     atom    1 type  1   force =     0.00000466    0.00000837    0.00000332
	my @force = grep {if(m/^.+force =\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)/){
			$_ = [$1*$force_convert,$2*$force_convert,$3*$force_convert];}} @all;
	for (1..@force){
		#print "$_ @{$force[$_ -1]}[0..2]\n";
		chomp @{$force[$_ -1]}[0..2];
		print $fraw "@{$force[$_ -1]}[0..2] ";
		print $fraw "\n" if ($_% $natom == 0 and $_ ne scalar @force);
	}
############## coord ############
##ATOMIC_POSITIONS (angstrom)        
##Al           -0.0000004209       -0.0000004098       -0.0000002490
my @coord = grep {if(m/^\w+\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)/){
			$_ = [$1,$2,$3];}} @all;
	for (1..@coord){
		#print "$_ @{$coord[$_ -1]}[0..2]\n";
		my $mod = $_ % ($natom*3);
		if(! grep /$mod/,(1..$natom*2)){
			chomp @{$coord[$_ -1]}[0..2];
			print $craw "@{$coord[$_ -1]}[0..2] ";
			print $craw "\n" if ($mod == 0 and $_ ne scalar @coord);
		}
	}
############## box ############
##CELL_PARAMETERS (angstrom)
##4.031848986   0.000000009   0.000000208
my @box = grep {if(m/^\s{1,3}([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)\s+([-+]?\d+\.?\d+)$/){
			$_ = [$1,$2,$3];}} @all;
	for (1..@box){
		#print "$_ @{$box[$_ -1]}[0..2]\n";
		my $mod = $_ % 6;
		if(! grep /$mod/,(1..3)){
			chomp @{$box[$_ -1]}[0..2];
			print $braw "@{$box[$_ -1]}[0..2] ";
			print $braw "\n" if ($mod == 0 and $_ ne scalar @box);
		}
	}
# end of qe2dpmd convertion
	close ($eraw);
	close ($vraw);
	close ($fraw);
	close ($craw);
	close ($braw);
} # main loop
