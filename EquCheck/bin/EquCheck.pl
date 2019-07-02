#############################################################
# 
#        EquCheck 
#        1.0.0
#
#        Shamik Bhat, Sahitya Talachutla, Hai Lin 
#
#        University of Colorado Denver 
#
#        Dec. 19, 2018
#
#    This program is used to check whether or not a recorded 
#  time series from an molecular dynamics simulation has 
#  reached equilibrium. 
#
#    This program is provided as it is and free of charge.
#  No liability is accepted under any condition.
#
#    Questions and Comments please send to 
#         Prof. Hai Lin (hai.lin@ucdenver.edu)
#
##############################################################

###### use #######     
use strict;
use warnings;
use diagnostics;
use Cwd qw(abs_path);
use Cwd;
use Scalar::Util 'looks_like_number';

############################################################## 
#
#              Interactive input 
#
##############################################################

######  Welcome information  ######
my @ucdenver = 'Foo';
print "Welcome to EquCheck 1.0.0. \n";
print "This program examines if a time series recorded in MD simulation reaches equilibrium. \n";
print "Written by Shamik Bhat, Sahitya Talachutla, Hai Lin. \n";
print "Last code change at 12/19/2018.\n"; 
print 'Please send any comments to Prof. Hai Lin (hai.lin@ucdenver.edu).', "\n";

######  Set the working directory  #######
my $dirfil;
my $direc;
my $cwd;
$cwd= getcwd();
print "\nThe current working directory of EquCheck.pl is\n$cwd\n";	
while ((my $t=0) == 0){
	print "\nWould you like to change the directory? Answer Y or N (default).\n";
	$dirfil=<STDIN>;
	$dirfil=uc($dirfil);
	chomp $dirfil;
	if ($dirfil eq ""){
		$direc= getcwd();
		$dirfil="N";
		chomp $dirfil;
		$t=1;
		last;
	}
	if ($dirfil eq "Y"){
		while ((my $tr=0) == 0){
			print "What is the new working directory? (Do not include forward slash at end)\n";
			$direc=<STDIN>;
			chomp $direc;
			if (-d $direc and (substr $direc, -1) ne "/"){
				chdir $direc;
				my $newd=getcwd();
				print "The new working directory is \n$newd\n";
				$tr=1;
				last;
			}
			else{
				print "Sorry that is not a valid directry. Please give a valid directory.\n";
				next;
			}
		}
		$t=1;
		last;	
	}
	elsif ($dirfil eq "N"){
		$direc= getcwd();
		$t=1;
		last;
	}
	else{
		print "Sorry that is not a valid value. Please try again. \n";
		next;
	}
}

######  Set the data file  ######
my $tempdir;
my $fnh;
my @files;
$tempdir=$direc;
my $filename;
while ((my $tr=0) == 0){
	print "\nType the file name within directory!\n";
	opendir DIR, $direc or die "Error in opnening directory $direc";
	while( ($fnh = readdir(DIR))) {
  	 print("$fnh\n");
  	 chomp $fnh;
   	push @files, $fnh;
	}
	closedir DIR;
	$filename=<STDIN>;
	my $pass=0;
	chomp $filename;
	for (my $l=0; $l<scalar(@files); $l++){
		if ($filename eq $files[$l]){
			$pass=1;
			last;
		}
		else{
			next;
		}
	}
	if ($pass==1){
		$tr=1;
		last;
	}
	else{
		print "Sorry that is not a valid file name. Please choose a filename within the directory.\n";
	}
}	


######  define parameters  ######
my @finaldata;  # defies the final data array to be used for the rest of the program
my $Zscore;
my $Tscore;
my $yes="Y";
my $no="N";
my $samplevar=0;
my $fn= "EquCheck.out";
my $bas= "1";
chomp $bas;
my $inter= "2";
chomp $inter;
my $adv= "3";
chomp $adv;
open (FH, ">", $fn) or die "Can't open file\n";	

######  Global defines z score  ######

######  Read the time series from the data file and save to an array ######
sub DATAMOD () {
	# Provides user input for the file name
	chomp $filename;
	# Removes excess lines from file
	open (my $DATA, '<', $filename) 
		or die "Can't open $filename: $!";
	# Opens the user speficied data file; show an error if file is not found.
	my @rmsd;
	# defines the array for the raw data.
	while (<$DATA>) {
		chomp;
		push @rmsd, $_;
	}
	# this loop removes excess \n and adds each line of data as one array variable to @rmsd.
	my @commasplitter;
	# defines the array for which we add commas to the data, instead of spaces.
	for my $line (@rmsd) {
		push (@commasplitter, join (",", split (/[:,\s\/]+/, $line)));
	}
	# The loop removes extra spaces from @rmsd, and inserts commas into a new array called @commasplitter.
	my @dataset;
	# This defines the data set in data point form.
	for my $dataline (@commasplitter) {
		my $reversedata= reverse $dataline; # reverses each number so that we can remove the excess commas before each line
		my $tust = substr $reversedata, -1;
		if ($tust eq ","){
			chop $reversedata; # removes the excess comma at the end	
		}
		my $datapoints= reverse $reversedata; # gives the final data points without extra commas.
		push (@dataset, $datapoints); # all of the data is in a new array.
	}
	for my $finalvalue (@dataset) {
		my @valueset = split (",", $finalvalue); # removes any commas from the data, and stores the variables into a temporary array @valueset.
		push (@finaldata, [$valueset[0], $valueset[1]]); # each variable is then inputted into a final array, so that each data point is its own individual array of valeus.
	}
	close $DATA;
}


######  Input parameters  ######
my $fail=0;
# Adds a failsafe in case a start time that is not compatible is chosen.
my $n;
# Asks for user input on sample size n
my $sl;
# Asks for user input for segment length of each sample length.
my $alevel; #<STDIN>;
my $ulevel;
# Asks for specified confidence level to reject null hypothesis.
my $stime; #<STDIN>;
my $timestep; #<STDIN>;
# This gives us the number of time units between each data point taken.
my $ststep;
# Defines the start time without units ($ststep).
my $tor;
# Defines the number of points per segment.
my @sampleval;
# Defines new array of sample values.
my $startcode;
# Defines the first data point index of the first segment of the sample.
my $endcode;
# Defines the last data point index of the first segment of the sample.
my $sum=0;
# defines sum variable.
my $var=0;
# defines variance as a variable.
my $AZPass=0;
my $VZPass=0;
my $otor;
my $oststep;
my $response;
my $n2;
my $nor;
my $lvlr;

######  Interactive parameters input  ######
sub DATACALL(){

#####  Ask for sample size  ######
	while ((my $t=0) == 0){
		print "What is the sample size? For the sake of statistical power, this must be greater than or equal to 24. \n";
		print "The default sample size is 24. Press enter for default sample size. \n";
		$n=<STDIN>;
		chomp $n;
		if ($n eq ""){
			$n=24;
			$t=1;
			last;
		}
		if (looks_like_number($n)) {
			if ($n>=24){
				$t=1;
				last;
			}
		elsif ($n<24){
				print "The sample size is less than 24. Please choose a sample size greater than or equal to 24. \n";
				next;
			}
		}
		else{
			print "Sorry that is not a valid value. Please try again\n";
		}
	}
	$nor=$n;

######  Ask for sample size cap  ######
	while ((my $tr=0) == 0){
		print "Would you like to impose a limit on sample size? (Y or N)\n";
		print "The default answer is N. Press Enter for default.\n";
		$response=<STDIN>;
		$response=uc($response); # lowercase becomes uppercase.
		chomp $response;
		if ($response eq ""){
			$response = "N";
			chomp $response;
		}
		if ($response eq $yes){
			print "What is the sample size cap?\n";
			$n2=<STDIN>;
			chomp $n2;
		if (looks_like_number($n2)) {
			if ($n<=$n2){
				$tr=1;
				last;
			}
		}
			else {
				print "Sorry. The cap you chose is lower than the original sample size or does not exist. Try choosing a value greater than your sample size.\n";
				next;
			}
		}
		elsif ($response eq $no){
			$tr=1;
			last;
		}
		else{
			print "Sorry that is not a valid response. Please try again.\n";
			next;
		}
	}

######   Gives time interval  ######

		my $intsp = $finaldata[2][0]-$finaldata[1][0];
		print "The time interval (in time steps) between adjacent data points is $intsp. \n\n";
		$timestep= $intsp;
		chomp $timestep;
		# This gives us the number of time units between each data point taken.


######  Ask for segment length  ######
	while ((my $tr=0) == 0){
		print "What is your segment length (fluctuation time) in time steps? For example, if the time series is recorded every 10 time steps, and if the fluctuation time is 200 time steps, each fragment will have 200/10 = 20 data points. \n";
		print "The default segment length is 2*(time interval).\n";
		$sl=<STDIN>;
		chomp $sl;
		if ($sl eq ""){
			$sl = 2*$timestep;
			chomp $sl;
			$tr=1;
			last;
		}
		#Illegal division 
		if ( looks_like_number($sl) ){
		if ($sl <= 0){
			$sl = 2*$timestep;
			print "WARNING: Negative value or 0 was entered. A multiple of time interval must be entered. The segment length is set to 2*(time interval).\n";
			print FH "WARNING: Negative value or 0 was entered. A multiple of time interval must be entered. The segment length is set to 2*(time interval).\n";
			chomp $sl;
			$tr=1;
			last;
		}
			$tr=1;
			last;
		}
		else{
			print "Sorry that is not a valid input for segment length. Please choose again. \n";
			next;
		}
	}
while ( (my $tr=0)==0){
		my $seglength=$sl/$timestep;
		if (not $seglength =~ /^\d+$/ or $seglength == 1){
			print "ERROR IN USER INPUT\n";
			print "That is not a valid value for segment length. Ensure your segment length is a mutliple of time interval.\n";
			print "What is your segment length (fluctuation time) in time steps? For example, if the time series is recorded every 10 time steps, and if the fluctuation time is 200 time steps, each fragment will have 200/10 = 20 data points. \n";
			print "The default segment length is 2*(time interval).\n";
			$sl=<STDIN>;
			chomp $sl;
		if ($sl eq ""){
			$sl = 2*$timestep;
			$tr=1;
			last;
			}
		if ( looks_like_number($sl)){
		if ($sl < 0){
			$sl = 2*$timestep;
			print "WARNING: Negative value was entered. The segment length is set to 2*(time interval).\n";
			print FH "WARNING: Negative value was entered. The segment length is set to 2*(time interval).\n";
			$tr=1;
			last;
				}
			}
		else{
			print "Sorry that is not a valid input for segment length. Please choose again. \n";
			next;
			}
		}
		else{
			$tr=1;
			last;
		}
	}
	

######  Ask for the alpha level  ######
	while ((my $tr=0) == 0){
		print "What is the alpha level? This is probability of making the wrong decision when the null hypothesis is true. Allowed values are 0.01, 0.02, 0.05 (default), 0.10, and 0.50. \n";
		$alevel=<STDIN>;
		chomp $alevel;
		if ($alevel eq ""){
			$alevel = "0.05";
			chomp $alevel;
			$tr=1;
			last;
		}
		if ( looks_like_number($alevel)){
			chomp $alevel;
			if ($alevel < 0){
				$alevel=0.05;
				print "WARNING: Negative value was entered. The alpha value is set to 0.05.\n";
				print FH "WARNING: Negative value was entered. The alpha value is set to 0.05.\n";
				$tr=1;
				last;
			}
			$ulevel= 1-($alevel/2);
			my $uzp=($ulevel-0.50)*100;
  			chomp $uzp;
  			my $utp=$uzp;
  			chomp $utp;
			if ($utp == 0 or $utp == 25 or $utp == 45 or $utp == 47.5 or $utp == 49 or $utp == 49.5 or $alevel == 0.10){
				$tr=1;
				last;
			}
			else{
				print "Sorry that is not a valid input for alpha level. Please choose again. \n";
				next;
			}
		}
		else{
			print "Sorry that is not a valid input for alpha level. Please choose again. \n";
			next;
		}
	}




######   Ask for the start time  ######
	while ((my $tr=0) == 0){
		print "What is the estimated time (in time steps) when the system reaches equilibrium? \n";
        print "Default value is 0 (i.e. the system is already in equilibrium). \n";
		$stime=<STDIN>;
		chomp $stime;
		if ($stime eq ""){
			$stime = "0";
			chomp $stime;
			$tr=1;
			last;
		}
		if ( looks_like_number($stime)){
		if ($stime < 0){
			$stime=0;
			print "WARNING: Negative value was entered. The start time is set to 0. \n";
			print FH "WARNING: Negative value was entered. The start time is set to 0. \n";
			$tr=1;
			last;
		}
			$tr=1;
			last;
		}
		else{
			print "Sorry that is not a valid input for start time. Please choose again. \n";
			next;
		}
	}

while ( (my $tr=0)==0){
		my $starttest=$stime/$timestep;
		if ($starttest =~ /\./){
			print "ERROR IN USER INPUT\n";
			print "That is not a valid value for start time. Ensure your start time is a mutliple of time interval.\n";
			print "What is the estimated time (in time steps) when the system reaches equilibrium? \n";
        		print "Default value is 0 (i.e. the system is already in equilibrium). \n";
			$stime=<STDIN>;
			chomp $stime;
		if ($stime eq ""){
			$stime = "0";
			chomp $stime;
			$tr=1;
			last;
		}

		if ( looks_like_number($stime)){
		if ($stime < 0){
			$stime = 0;
			print "WARNING: Negative value was entered. The start time is set to 0.\n";
			print FH "WARNING: Negative value was entered. The start time is set to 0. \n";
			$tr=1;
			last;
					}
			next;
		}
		else{
			print "Sorry that is not a valid input for start time. Please choose again. \n";
			next;
		}
		}
		else {
			$tr=1;
			last;
		}
}
	$ststep=$stime/$timestep;
	$oststep=$ststep;


######  Ask for output verbose level  ######
	while ((my $tr=0)==0){
                print "Please select the level (1, 2, or 3) of output details? Level 1 is brief, level 2 (default) is modest, and level 3 is very detailed.\n";
		# print "What level of information would you like in the output? \n(1= Final test information, 2= Includes every time test passes and increments n, 3= includes all diagnostic informaiton if test fails.\n ";
		# print "The default level is level 2. Press Enter for default level.\n";
		$lvlr=<STDIN>;
		chomp $lvlr;
		if ($lvlr eq ""){
			$lvlr="2";
			chomp $lvlr;
			$tr=1;
			last;
		}
		if ($lvlr eq $bas or $lvlr eq $inter or $lvlr eq $adv){
			$tr++;
			last;
		}
		else {
			print "Sorry. That is not a valid level. Please Choose 1, 2, or 3.\n";
		}
	}


#################################################
#
#       Prepare the parameters
#
#################################################

######  Convert time variables into the unit of time interval  ###### HL10072018
	# Dividing by the time interval($timeststep) here, gives us the start time without units ($ststep).
	$tor = $sl / $timestep;
	$otor = $tor;
	# Dividing segment length ($sl) by the time interval($timeststep) here, gives us the number of points per segment.
	$startcode = $ststep - 1;
	# Defines the first data point index of the first segment of the sample.
	$endcode = $ststep + $tor - 2;
	# Defines the last data point index of the first segment of the sample.
	$sum = 0;
	# defines sum variable.
	$var = 0;
	# defines variance as a variable.
}

######  Define Z values ###### HL10072018
my $ztrue;
my $ttrue;
# Defines another failsafe for Z conversion.
print "SAMPLE INFORMATION, SHAPE TEST BELOW \n";
sub Z_and_T_CONVERSION(){
	######  Converts Alpha Level to Z and T input Parameters ###### HL10072018
	$ulevel= 1-($alevel/2);
	my @_Pctile_Z_map =
    ([ 0,     0 ],
    [ 0.5,   0.01253347 ],
    [ 1,     0.02506891 ],
    [ 1.5,   0.03760829 ],
    [ 2,     0.05015358 ],
    [ 2.5,   0.06270678 ],
    [ 3,     0.07526986 ],
    [ 3.5,   0.08784484 ],
    [ 4,     0.1004337 ],
    [ 4.5,   0.1130385 ],
    [ 5,     0.1256613 ],
    [ 5.5,   0.1383042 ],
    [ 6,     0.1509692 ],
    [ 6.5,   0.1636585 ],
    [ 7,     0.1763742 ],
    [ 7.5,   0.1891184 ],
    [ 8,     0.2018935 ],
    [ 8.5,   0.2147016 ],
    [ 9,     0.227545 ],
    [ 9.5,   0.240426 ],
    [ 10,    0.2533471 ],
    [ 10.5,  0.2663106 ],
    [ 11,    0.279319 ],
    [ 11.5,  0.2923749 ],
    [ 12,    0.3054808 ],
    [ 12.5,  0.3186394 ],
    [ 13,    0.3318533 ],
    [ 13.5,  0.3451255 ],
    [ 14,    0.3584588 ],
    [ 14.5,  0.3718561 ],
    [ 15,    0.3853205 ],
    [ 15.5,  0.3988551 ],
    [ 16,    0.4124631 ],
    [ 16.5,  0.426148 ],
    [ 17,    0.4399132 ],
    [ 17.5,  0.4537622 ],
    [ 18,    0.4676988 ],
    [ 18.5,  0.4817268 ],
    [ 19,    0.4958503 ],
    [ 19.5,  0.5100735 ],
    [ 20,    0.5244005 ],
    [ 20.5,  0.538836 ],
    [ 21,    0.5533847 ],
    [ 21.5,  0.5680515 ],
    [ 22,    0.5828415 ],
    [ 22.5,  0.5977601 ],
    [ 23,    0.612813 ],
    [ 23.5,  0.628006 ],
    [ 24,    0.6433454 ],
    [ 24.5,  0.6588377 ],
    [ 25,    0.6744898 ],
    [ 25.5,  0.6903088 ],
    [ 26,    0.7063026 ],
    [ 26.5,  0.7224791 ],
    [ 27,    0.7388468 ],
    [ 27.5,  0.755415 ],
    [ 28,    0.7721932 ],
    [ 28.5,  0.7891917 ],
    [ 29,    0.8064212 ],
    [ 29.5,  0.8238936 ],
    [ 30,    0.8416212 ],
    [ 30.5,  0.8596174 ],
    [ 31,    0.8778963 ],
    [ 31.5,  0.8964734 ],
    [ 32,    0.9153651 ],
    [ 32.5,  0.9345893 ],
    [ 33,    0.9541653 ],
    [ 33.5,  0.9741139 ],
    [ 34,    0.9944579 ],
    [ 34.5,  1.015222 ],
    [ 35,    1.036433 ],
    [ 35.5,  1.058122 ],
    [ 36,    1.080319 ],
    [ 36.5,  1.103063 ],
    [ 37,    1.126391 ],
    [ 37.5,  1.150349 ],
    [ 38,    1.174987 ],
    [ 38.5,  1.200359 ],
    [ 39,    1.226528 ],
    [ 39.5,  1.253565 ],
    [ 40,    1.281552 ],
    [ 40.5,  1.310579 ],
    [ 41,    1.340755 ],
    [ 41.5,  1.372204 ],
    [ 42,    1.405072 ],
    [ 42.5,  1.439531 ],
    [ 43,    1.475791 ],
    [ 43.5,  1.514102 ],
    [ 44,    1.554774 ],
    [ 44.5,  1.598193 ],
    [ 45,    1.644854 ],
    [ 45.5,  1.695398 ],
    [ 46,    1.750686 ],
    [ 46.5,  1.811911 ],
    [ 47,    1.880794 ],
    [ 47.5,  1.959964 ],
    [ 48,    2.053749 ],
    [ 48.5,  2.17009 ],
    [ 49,    2.326348 ],
    [ 49.5,  2.575829 ],
    [ 49.9,  3.090232 ],
    [ 49.95, 3.290527 ],
    [ 49.99, 3.719016 ]);
   # This is the table of Z statistics

my @_Pctile_T_map =
    ([ 0,     0 ],
    [ 25,    0.858 ],
    [ 30,    0.8416212 ],
    [ 35,    1.060 ],
    [ 40,    1.319 ],
    [ 45,    1.714 ],
    [ 47.5,  2.069 ],
    [ 49,    2.500 ],
    [ 49.5,  2.807 ],
    [ 49.9,  3.485 ],
    [ 49.95, 3.768 ]);
  	my $pt;
  	my $uzp=($ulevel-0.50)*100;
  	chomp $uzp;
  	my $utp=$uzp;
  	chomp $utp;
  	# This converts the alphalevel given to that the tables accept.
  	######  Inputs Alpha Level and Finds Corresponding Z and T Values ###### HL10072018
  	for (my $i=0; $i<scalar(@_Pctile_Z_map); $i++) {
  		if ($uzp == $_Pctile_Z_map[$i][0]){
  			$ztrue=1;
  			$Zscore = $_Pctile_Z_map[$i][1];
  			last;
  		}
  		else{
  			$ztrue=0;
  			next;
  		}
  		}
  		
 # This calls for the z statistic by cycling though the the z table and comparing to the alpha level to find the right value. This then gets stored.
	for (my $i=0; $i<scalar(@_Pctile_T_map); $i++) {
  		if ($utp == $_Pctile_T_map[$i][0]){
  			$ttrue=1;
  			$Tscore= $_Pctile_T_map[$i][1];
  			last;

  		}
  		else{
  			$ttrue=0;
  			next;
  # This calls for the t statistic by cycling through the t table and comparing to the alpha level to find the right value. This then gets stored.
  		}
  	  		}
  	if ($alevel == 0.10){
  		$uzp=45.0;
  		$utp=45.0;
  		$Zscore=1.644854;
  		$Tscore=1.714;
  		$ztrue=1;
  		$ttrue=1;
  		}
}	


######  Define Global Variables for Thermodynamic Averages ###### HL10072018
my $sampleavgsum=0;
my $sampleaverage=0;
my $samvarsum=0;
######  Obtains Statistical Averages and Standard Deviations ###### HL10072018
sub SAMPLE(){
	# sample averages, std deviaitons, n, k, tor, a need to be defined in a new subroutine.
	@sampleval= [0];
	$sampleavgsum=0;
	$sampleaverage=0;
	$samvarsum=0;
	$samplevar=0;
	pop @sampleval;
	$fail=0;
	# Redefines failure parameters and sample parameters for each new run to prevent cumulative summation accross samples.
	$startcode = $ststep-1;
	$endcode = $ststep+$tor-2;
	# Defines the first and last data points for the segment.
	 for (my $j=0; $j<=($n-1); $j++){
		if($endcode < scalar(@finaldata) or $startcode < 0){ # CHANGE TO ACCOUNT FOR GOING OVER.)
			for (my $i=$startcode; $i<=$endcode; $i++) {
				$sum += $finaldata[$i][1];
			}
			# Cycles through the points within the segment ands all of the averages for that point within $sum.
			my $segav= $sum / $tor;
			# Here we get the segment average for the $jth segment.
			for (my $k=$startcode; $k<=$endcode; $k++) {
				my $var1 = ($finaldata[$k][1] - $segav)**2;
				$var +=$var1;
			}
			# Segment Variances are obtaines and added to the $var value in a similar manner.
			my $stdev= ($var/($tor-1))**(1/2);
			# We get the segement standard deviation of the $jth segment.
			push (@sampleval, [$j+1, $segav, $stdev]);
			# This adds the segment's id ($j+1), the segment's average ($segav), and the stand deviaiton of the segment ($stdev).
			$var=0;
			$sum=0;
			# These reset the sum and variance of the data for the next segment $j+1.
			$startcode += $tor;
			$endcode += $tor;		
			# These reset the start and end incides to recall from the dat for the next segment $j+1.
			}
		else {
			last;
		}
	}
	for (my $i=0; $i<scalar(@sampleval); $i++) {
		$sampleavgsum = $sampleavgsum + $sampleval[$i][1]; 
		
	}
	$sampleaverage=$sampleavgsum/$n; 
	# The sample average was obtained by cyling through the segment averages in @sampleval and adding thing to $sampleavgsum. Then the number of segments was divided.
	for (my $i=0; $i<scalar(@sampleval); $i++){
		$samvarsum+=($sampleval[$i][1] - $sampleaverage)**2; 
	}
	$samplevar = $samvarsum / ($n - 1); 
	# The sample variance was obtained by cycling though the segment variances and adding the difference in means to the $samvarsum. Then sample variance was calcualted.
	if ($endcode>scalar(@finaldata)){
		$fail = 1;
	}
	else{
		$fail = 0;
	}
	#This is a failsafe in case the sampled data goes beond what is available in the set. This leads to an immediate stop in sampling.
}

######  Conducts the Shapiro-Wilk Test for $n<=50 ###### HL10072018
my $W=30;
my $WPass=0;
my $SWfail=0;
my $Wcomparison;
# Defines parameters in case test does not pass initial analysis.
sub SW(){
	$SWfail=13;
	my @SWtable =
   ([ 24, 1, 0.4493],
    [ 24, 2, 0.3098],
    [ 24, 3, 0.2554],
    [ 24, 4, 0.2145],
    [ 24, 5, 0.1807],
    [ 24, 6, 0.1512],
    [ 24, 7, 0.1245],
    [ 24, 8, 0.0878],
    [ 24, 9, 0.0764],
    [ 24, 10, 0.0459],
    [ 24, 11, 0.0321],
    [ 24, 12, 0.0107],
    [ 25, 1,  0.4450],
    [ 25, 2,  0.3069],
    [ 25, 3, 0.2533],
    [ 25, 4, 0.2148],
    [ 25, 5, 0.1822],
    [ 25, 6, 0.1539],
    [ 25, 7, 0.1283],
    [ 25, 8, 0.1046],
    [ 25, 9, 0.0823],
    [ 25, 10, 0.0610],
    [ 25, 11, 0.0403],
    [ 25, 12, 0.02],  
    [ 25, 13, 0.0],
    [ 26, 1, 0.4407],
    [ 26, 2, 0.3043],
    [ 26, 3, 0.2533],
    [ 26, 4, 0.2151],
    [ 26, 5, 0.1822],
    [ 26, 6, 0.1563],
    [ 26, 7, 0.1316],
    [ 26, 8, 0.1089],
    [ 26, 9, 0.0823], 
    [ 26, 10, 0.0672],
    [ 26, 11, 0.0476],
    [ 26, 12, 0.0284],
    [ 26, 13, 0.0094],
    [ 27, 1, 0.4366],
    [ 27, 2, 0.3018],
    [ 27, 3, 0.2522],
    [ 27, 4, 0.2152],
    [ 27, 5, 0.1848],
    [ 27, 6, 0.1584],
    [ 27, 7, 0.1346],
    [ 27, 8, 0.1128],
    [ 27, 9, 0.0923],
    [ 27, 10, 0.0728],
    [ 27, 11, 0.0540],
    [ 27, 12, 0.0358],
    [ 27, 13, 0.0178],
    [ 27, 14, 0.0],
    [ 28, 1, 0.4328],
    [ 28, 2, 0.2992],
    [ 28, 3, 0.2510],
    [ 28, 4, 0.2152],
    [ 28, 5, 0.1857],
    [ 28, 6, 0.1601],
    [ 28, 7, 0.1372],
    [ 28, 8, 0.1162],
    [ 28, 9, 0.0965],
    [ 28, 10, 0.0778],
    [ 28, 11, 0.0598],
    [ 28, 12, 0.0424],
    [ 28, 13, 0.0253],
    [ 28, 14, 0.0084],
    [ 29, 1, 0.4291],
    [ 29, 2, 0.2968],
    [ 29, 3, 0.2499],
    [ 29, 4, 0.2150],
    [ 29, 5, 0.1864],
    [ 29, 6, 0.1616],
    [ 29, 7, 0.1395],
    [ 29, 8, 0.1192],
    [ 29, 9, 0.1002], 
    [ 29, 10, 0.0822],
    [ 29, 11, 0.0650],
    [ 29, 12, 0.0483],
    [ 29, 13, 0.0320],
    [ 29, 14, 0.0159],
    [ 29, 15, 0.00000],
    [ 30, 1, 0.4254],
    [ 30, 2, 0.2944],
    [ 30, 3, 0.2487],
    [ 30, 4, 0.2148],
    [ 30, 5, 0.1870],
    [ 30, 6, 0.1630],
    [ 30, 7, 0.1415],
    [ 30, 8, 0.1219],
    [ 30, 9, 0.1036],
    [ 30, 10, 0.0862],
    [ 30, 11, 0.0697],
    [ 30, 12, 0.0537],
    [ 30, 13, 0.0381],
    [ 30, 14, 0.0227],
    [ 30, 15, 0.0076],
    [ 31, 1, 0.4220],
    [ 31, 2, 0.2921],
    [ 31, 3, 0.2475],
    [ 31, 4, 0.2145],
    [ 31, 5, 0.1874],
    [ 31, 6, 0.1641],
    [ 31, 7, 0.1433],
    [ 31, 8, 0.1243],
    [ 31, 9, 0.1066],
    [ 31, 10, 0.0899],
    [ 31, 11, 0.0739],
    [ 31, 12, 0.0585],
    [ 31, 13, 0.0435],
    [ 31, 14, 0.0289],
    [ 31, 15, 0.0144],
    [ 31, 16, 0.00000],
    [ 32, 1, 0.4188],
    [ 32, 2, .2898],
    [ 32, 3, .2463],
    [ 32, 4, .2141],
    [ 32, 5, .1878],
    [ 32, 6, .1651],
    [ 32, 7, .1449],
    [ 32, 8, .1265],
    [ 32, 9, .1093],
    [ 32, 10, .0931],
    [ 32, 12, 0.0629],
    [ 32, 13, 0.4857],
    [ 32, 14, 0.0344],
    [ 32, 15, 0.0206],
    [ 32, 16, 0.0068],
    [ 33, 1, 0.4156],
    [ 33, 2, 0.2898],
    [ 33, 3, 0.2463],
    [ 33, 4, 0.2137],
    [ 33, 5, 0.1880],
    [ 33, 6, 0.1660],
    [ 33, 7, 0.1463],
    [ 33, 8, 0.1284],
    [ 33, 9, 0.1118],
    [ 33, 10, 0.0961],
    [ 33, 11, 0.0812],
    [ 33, 12, 0.0669],
    [ 33, 13, 0.0530],
    [ 33, 14, 0.0395],
    [ 33, 15, 0.0262],
    [ 33, 16, 0.0131],
    [ 33, 17, 0.0000],
    [ 34, 1, 0.4127],
    [ 34, 2, 0.2854],
    [ 34, 3, 0.2439],
    [ 34, 4, 0.2132],
    [ 34, 5, 0.1882],
    [ 34, 6, 0.1667],
    [ 34, 7, 0.1475],
    [ 34, 8, 0.1301],
    [ 34, 9, 0.1140],
    [ 34, 10, 0.0988],
    [ 34, 11, 0.0844],
    [ 34, 12, 0.0706],
    [ 34, 13, 0.0572],
    [ 34, 14, 0.0441],
    [ 34, 15, 0.0314],
    [ 34, 16, 0.0187],
    [ 34, 17, 0.0062],
    [ 35, 1, 0.4096],
    [ 35, 2, 0.2834],
    [ 35, 3, 0.2427],
    [ 35, 4, 0.2127],
    [ 35, 5, 0.1883],
    [ 35, 6, 0.1673],
    [ 35, 7, 0.1487],
    [ 35, 8, 0.1317],
    [ 35, 9, 0.1160],
    [ 35, 10, 0.1013],
    [ 35, 11, 0.0873],
    [ 35, 12, 0.0739],
    [ 35, 13, 0.0610],
    [ 35, 14, 0.0484],
    [ 35, 15, 0.0361],
    [ 35, 16, 0.0239],
    [ 35, 17, 0.0119],
    [ 35, 18, 0.0000],
    [ 36, 1, 0.4068],
    [ 36, 2, 0.2813],
    [ 36, 3, 0.2415],
    [ 36, 4, 0.2121],
    [ 36, 5, 0.1883],
    [ 36, 6, 0.1678],
    [ 36, 7, 0.1496],
    [ 36, 8, 0.1331],
    [ 36, 9, 0.1179],
    [ 36, 10, 0.1036],
    [ 36, 11, 0.0900],
    [ 36, 12, 0.0770],
    [ 36, 13, 0.0645],
    [ 36, 14, 0.0523],
    [ 36, 15, 0.0404],
    [ 36, 16, 0.0287],
    [ 36, 17, 0.0172],
    [ 36, 18, 0.0057],
    [ 37, 1, 0.4040],
    [ 37, 2, 0.2794],
    [ 37, 3, 0.2403],
    [ 37, 4, 0.2116],
    [ 37, 5, 0.1883],
    [ 37, 6, 0.1683],
    [ 37, 7, 0.1505],
    [ 37, 8, 0.1344],
    [ 37, 9, 0.1196],
    [ 37, 10, 0.1056],
    [ 37, 11, 0.0924],
    [ 37, 12, 0.0798],
    [ 37, 13, 0.0677],
    [ 37, 14, 0.0559],
    [ 37, 15, 0.0444],
    [ 37, 16, 0.0331],
    [ 37, 17, 0.0220],
    [ 37, 18, 0.0110],
    [ 37, 19, 0.0000],
    [ 38, 1, 0.4015],
    [ 38, 2, 0.2774],
    [ 38, 3, 0.2391],
    [ 38, 4, 0.2110],
    [ 38, 5, 0.1881],
    [ 38, 6, 0.1686],
    [ 38, 7, 0.1513],
    [ 38, 8, 0.1356],
    [ 38, 9, 0.1211],
    [ 38, 10, 0.1075],
    [ 38, 11, 0.0947],
    [ 38, 12, 0.0824],
    [ 38, 13, 0.0706],
    [ 38, 14, 0.0592],
    [ 38, 15, 0.0481],
    [ 38, 16, 0.0372],
    [ 38, 17, 0.0264],
    [ 38, 18, 0.0158],
    [ 38, 19, 0.0053],
    [ 39, 1, 0.3989],
    [ 39, 2, 0.2755],
    [ 39, 3, 0.2380],
    [ 39, 4, 0.2104],
    [ 39, 5, 0.1880],
    [ 39, 6, 0.1689],
    [ 39, 7, 0.1520],
    [ 39, 8, 0.1366],
    [ 39, 9, 0.1225],
    [ 39, 10, 0.1092],
    [ 39, 11, 0.0967],
    [ 39, 12, 0.0848],
    [ 39, 13, 0.0733],
    [ 39, 14, 0.0622],
    [ 39, 15, 0.0515],
    [ 39, 16, 0.0409],
    [ 39, 17, 0.0305],
    [ 39, 18, 0.0203],
    [ 39, 19, 0.0101],
    [ 39, 20, 0.000], 
    [ 40, 1, 0.3964],
    [ 40, 2, 0.2737],
    [ 40, 3, 0.2368],
    [ 40, 4, 0.2098],
    [ 40, 5, 0.1878],
    [ 40, 6, 0.1691],
    [ 40, 7, 0.1526],
    [ 40, 8, 0.1376],
    [ 40, 9, 0.1237],
    [ 40, 10, 0.1108],
    [ 40, 11, 0.0986],
    [ 40, 12, 0.0870],
    [ 40, 13, 0.0759],
    [ 40, 14, 0.0651],
    [ 40, 15, 0.0546],
    [ 40, 16, 0.0444],
    [ 40, 17, 0.0343],
    [ 40, 18, 0.0244],
    [ 40, 19, 0.0146],
    [ 40, 20, 0.0049],
    [ 41, 1, 0.3940],
    [ 41, 2, 0.2719],
    [ 41, 3, 0.2357],
    [ 41, 4, 0.2091],
    [ 41, 5, 0.1876],
    [ 41, 6, 0.1693],
    [ 41, 7, 0.1531],
    [ 41, 8, 0.1384],
    [ 41, 9, 0.1249],
    [ 41, 10, 0.1123],
    [ 41, 11, 0.1004],
    [ 41, 12, 0.0891],
    [ 41, 13, 0.0782],
    [ 41, 14, 0.0677], 
    [ 41, 15, 0.0575],
    [ 41, 16, 0.0476],
    [ 41, 17, 0.0379],
    [ 41, 18, 0.0283],
    [ 41, 19, 0.0188],
    [ 41, 20, 0.0094],
    [ 41, 21, 0.000],
    [ 42, 1, 0.3917],
    [ 42, 2, 0.2701],
    [ 42, 3, 0.2345],
    [ 42, 4, 0.2085],
    [ 42, 5, 0.1874],
    [ 42, 6, 0.1694],
    [ 42, 7, 0.1535],
    [ 42, 8, 0.1392],
    [ 42, 9, 0.1259],
    [ 42, 10, 0.1136],
    [ 42, 11, .1004],
    [ 42, 12, .0891],
    [ 42, 13, .0782],
    [ 42, 14, .0677],
    [ 42, 15, .0575],
    [ 42, 16, .0476],
    [ 42, 17, .0379],
    [ 42, 18, .0283],
    [ 42, 19, .0188],
    [ 42, 20, .0094],
    [ 42, 21, .0000],
    [ 43, 1, 0.3894],
    [ 43, 2, 0.2684],
    [ 43, 3, 0.2334],
    [ 43, 4, 0.2078],
    [ 43, 5, 0.1871],
    [ 43, 6, 0.1695],
    [ 43, 7, 0.1539],
    [ 43, 8, 0.1398],
    [ 43, 9, 0.1269],
    [ 43, 10, .1149],
    [ 43, 11, .1035],
    [ 43, 12, .0927],
    [ 43, 13, .0824],
    [ 43, 14, .0724],
    [ 43, 15, .0628],
    [ 43, 16, .0534],
    [ 43, 17, .0442],
    [ 43, 18, .0352],
    [ 43, 19, .0263],
    [ 43, 20, .0175],
    [ 43, 21, .0087],
    [ 43, 22, .0000],
    [ 44, 1, 0.3872],
    [ 44, 2, 0.2667],
    [ 44, 3, 0.2323],
    [ 44, 4, 0.2072],
    [ 44, 5, 0.1868],
    [ 44, 6, 0.1695],
    [ 44, 7, 0.1542],
    [ 44, 8, 0.1405],
    [ 44, 9, 0.1278],
    [ 44, 10, .1160],
    [ 44, 11, .1049],
    [ 44, 12, .0943],
    [ 44, 13, .0842],
    [ 44, 14, .0745],
    [ 44, 15, .0651],
    [ 44, 16, 0.0560],
    [ 44, 17, .0471],
    [ 44, 18, .0383],
    [ 44, 19, .0296],
    [ 44, 20, .0211],
    [ 44, 21, .0126],
    [ 44, 22, .0042],
    [ 45, 1, 0.3850],
    [ 45, 2, 0.2651],
    [ 45, 3, 0.2313],
    [ 45, 4, 0.2065],
    [ 45, 5, 0.1865],
    [ 45, 6, 0.1695],
    [ 45, 7, 0.1545],
    [ 45, 8, 0.1410],
    [ 45, 9, 0.1286],
    [ 45, 10, .1170],
    [ 45, 11, .1062],
    [ 45, 12, .0959],
    [ 45, 13, .0860],
    [ 45, 14, .0765],
    [ 45, 15, .0673],
    [ 45, 16, .0584],
    [ 45, 17, .0497],
    [ 45, 18, .0412],
    [ 45, 19, .0328],
    [ 45, 20, .0245],
    [ 45, 21, .0163],
    [ 45, 22, .0081],
    [ 45, 23, .0000],
    [ 46, 1, 0.3830],
    [ 46, 2, 0.2635],
    [ 46, 3, 0.2302],
    [ 46, 4, 0.2058],
    [ 46, 5, 0.1862],
    [ 46, 6, 0.1695],
    [ 46, 7, 0.1548],
    [ 46, 8, 0.1415],
    [ 46, 9, 0.1293],
    [ 46, 10, .1180],
    [ 46, 11, .1073],
    [ 46, 12, .0972],
    [ 46, 13, .0876],
    [ 46, 14, .0783],
    [ 46, 15, .0694],
    [ 46, 16, .0607],
    [ 46, 17, .0522],
    [ 46, 18, .0439],
    [ 46, 19, .0357],
    [ 46, 20, .0277],
    [ 46, 21, .0197],
    [ 46, 22, .0118],
    [ 46, 23, 0.0039],
    [ 47, 1, 0.3808],
    [ 47, 2, 0.2620],
    [ 47, 3, 0.2291],
    [ 47, 4, 0.2052],
    [ 47, 5, 0.1859],
    [ 47, 6, 0.1695],
    [ 47, 7, 0.1550],
    [ 47, 8, 0.1420],
    [ 47, 9, 0.1300],
    [ 47, 10, .1189],
    [ 47, 11, .1085],
    [ 47, 12, .0986],
    [ 47, 13, .0892],
    [ 47, 14, .0801],
    [ 47, 15, .0713],
    [ 47, 16, .0628],
    [ 47, 17, .0546],
    [ 47, 18, .0465],
    [ 47, 19, .0385],
    [ 47, 20, .0307],
    [ 47, 21, .0229],
    [ 47, 22, .0153],
    [ 47, 23, .0076],
    [ 47, 24, .0000],
    [ 48, 1, 0.3789],
    [ 48, 2, 0.2604],
    [ 48, 3, 0.2281],
    [ 48, 4, 0.2045],
    [ 48, 5, 0.1855],
    [ 48, 6, 0.1693],
    [ 48, 7, 0.1551],
    [ 48, 8, 0.1423],
    [ 48, 9, 0.1306],
    [ 48, 10, .1197],
    [ 49, 11, .1095],
    [ 48, 12, .0998],
    [ 48, 13, .0906],
    [ 48, 14, .0817],
    [ 48, 15, .0731],
    [ 48, 16, .0648],
    [ 48, 17, .0568],
    [ 48, 18, .0489],
    [ 48, 19, .0411],
    [ 48, 20, .0335],
    [ 48, 21, .0259],
    [ 48, 22, .0185],
    [ 48, 23, .0111],
    [ 48, 24, .0037],
    [ 49, 1, 0.3770],
    [ 49, 2, 0.2589],
    [ 49, 3, 0.2271],
    [ 49, 4, 0.2038],
    [ 49, 5, 0.1851],
    [ 49, 6, 0.1692],
    [ 49, 7, 0.1553],
    [ 49, 8, 0.1427],
    [ 49, 9, 0.1312],
    [ 49, 10, .1205],
    [ 49, 11, .1105],
    [ 49, 12, .1010],
    [ 49, 13, .0919],
    [ 49, 14, .0832],
    [ 49, 15, .0748],
    [ 49, 16, .0667],
    [ 49, 17, .0588],
    [ 49, 18, .0511],
    [ 49, 19, .0436],
    [ 49, 20, .0288],
    [ 49, 21, .0215],
    [ 49, 22, .0143],
    [ 49, 23, .0071],
    [ 49, 24, .0000],
    [ 50, 1, 0.3751],
    [ 50, 2, .2574],
    [ 50, 3, 0.2260],
    [ 50, 4, 0.2032],
    [ 50, 5, 0.1847],
    [ 50, 6, 0.1691],
    [ 50, 7, 0.1554],
    [ 50, 8, 0.1430],
    [ 50, 9, 0.1317],
    [ 50, 10, .1212],
    [ 50, 11, .1113],
    [ 50, 12, .1020],
    [ 50, 13, .0932],
    [ 50, 14, .0846],
    [ 50, 15, .0764],
    [ 50, 16, .0685],
    [ 50, 17, .0608],
    [ 50, 18, .0532],
    [ 50, 19, .0459],
    [ 50, 20, .0386],
    [ 50, 21, .0314],
    [ 50, 22, .0244],
    [ 50, 23, .0174],
    [ 50, 24, .0104],
    [ 50, 25, .0035]);
    # Makes a table with A-values, based on n, and k. Based on the orignal 1965 paper, check User Guide for more information.
    my @samplemeans; 
    my $A;
    my $A_argument;
    my $B=0;
    for (my $i=0; $i<scalar(@sampleval); $i++) {
    	push (@samplemeans, $sampleval[$i][1]); 
    }
	my @sortedmeans = sort { $a <=> $b } @samplemeans;
	# Added all segment means to @samplemeans into an array by cycling through @sampleval, and the means were sorted by number.
	my $l=int($n/2); #defines l as in the Shapiro-Wilk paper.
	for (my $i=0; $i<$l; $i++) { 
		for (my $j=0; $j<scalar(@SWtable); $j++) { #a loop that matches both "n" and "k" in the above array.
			if ($n == $SWtable[$j][0] and $i == $SWtable[$j][1]) { #if a match is found, stores A-value in a variable.
				$A=$SWtable[$j][2]; 
				$A_argument = $A*($samplemeans[$n-$i]-$samplemeans[$i-1]); #follows the summation equation in the paper based on A-value.
			}
			else {
				next; 
			}
			$B+=$A_argument; #sums the argument based on all the A-values and segment mean data.
		}
	}
######  Compares Outputed Variable to Test Statistics based on $alevel ###### HL10072018
		$W=((($B)**2)/(($samplevar)*($n-1))); #Uses B and sample variance to calculate the Shapiro-Wilk test statistic.
		my @Wtable =  #defines a table based on alpha level and W critical values for n=24, which can be used for all n-values up to 50.
		([0.01, .884],
		[0.02, .898],
		[0.05, .916],
		[0.10, .930],
		[0.50, .963],
		[0.90, .981],
		[0.95, .984],
		[0.98, .987],
		[0.99, .989]);
		for (my $i=0; $i<scalar(@Wtable); $i++) { #matches user-given alpha level with one in the table.
			if ($alevel == $Wtable[$i][0]) {
				$Wcomparison = $Wtable[$i][1]; #if alpha level match is found, stores the respective critical value.
			}
			else {
				next;
			}
		}
		# $Wtable is cycled through to find the appropriate alpha level $alevel, then  the corresponding critical value is stored.
		if ($W<=$Wcomparison) { 
			$WPass = 1;
		}
		else {
			$WPass = 0;
		}	
		# Tests whether the null hypothesis that there is normality in segment means is not rejected. If so, Shapiro Wilk Test Pasts	
}	

######  Conducts the Mann-Kendal Test for Averages ###### HL10072018
my $I=0;
my $semI=0;
my $MKSTD=0;
my $S=0;
my $semS=0;
my $deviate=30;
my $semdeviate=30;
# Defines some of the global variables needed in the upcoming subroutine that define testing variables
my $AMK=0;
sub MKA(){
	$AMK=12;
	$AZPass=0;
	$S=0;
	$I=0;
	$MKSTD=0;
	$semI=0;
	# Resets the global variables if the subroutine is changed again.
	my $factorial=1; #defines one of the factorials needed for the calculation of a mathematical combination
	my $factorial2=1; #defines the other factorial needed for the calculation of a mathematical combination
	for (my $j=0; $j<=($n-1); $j++){ #this for loop runs through each segment
		for (my $i=0; $i<$j; $i++){ #this for loop picks a segment mean and compares it to array term the j for loop is on
			if ($sampleval[$i][1]>$sampleval[$j][1]){
				$semI=$I+1;
				next; #this skips the points where there is a downward change
			}
			elsif ($sampleval[$i][1]<$sampleval[$j][1]){
				$I=$I+1; #counts the number of comparisons where there is an upward change
			}
			else {
				next; #skips any points that are equal to each other 
			}	
			}
		}		
	$MKSTD=(($n*($n-1)*(2*$n+5))/18)**.5; #this is the standard deviation of the test statistic S
	for (my $i=1; $i<=$n; $i++){ #this calculates the factorial of the sample size, n
		$factorial = $factorial*$i;
	}
	for (my $i=1; $i<=($n-2); $i++){ #this calculates the factorial of n-2, required to calculate the combination value from 2 to n
		$factorial2=$factorial2*$i;
	}
	my $combination=($factorial)/(2*$factorial2); #this finalizes the combination value from 2 to n
	$S=(2*$I)-$combination; #calculates the S test statistic
	$semS=$I-$semI;
	$deviate=$S/$MKSTD;
	$semdeviate=$semS/(($n*($n-1))/2); #calculates the standard deviate test statistic needed for analysis #this will be compared to test statistic values for percentiles (e.g. 1.96, 1.5, etc.)

	if (abs($deviate)<$Zscore){
		$AZPass=1;		
	}
	# Determines if Test Passes. The null hypothesis reveals no trend in segment means. The test passes if the null hypothesis is not rejected.
if ($lvlr eq $adv){
	print "\nDIAGNOSTIC INFORMATION\n";
	print "The start time step is now $ststep\n";
	print "The sample size is now $n\n";
	print "The segment length is now $tor\n";
	print FH "\nDIAGNOSTIC INFORMATION\n";
	print FH "The start time step is now $ststep\n";
	print FH "The sample size is now $n\n";
	print FH "The segment length is now $tor\n";
}
# Prints Diagnostic Test Information if Level 3 is called.
}
######  Conducts the Mann-Kenall Test for Variances ###### HL10072018
	my $vS=0;
	my $vI=0;
	my $semvI=0;
	my $vMKSTD=0;
	my $vfactorial=1;
	my $vfactorial2=1;
	my $vdeviate=30;
	my $VMK=0;
# These are variables used for the subroutine that determine the test statistic for test pass or fail that is required for later subroutines.
sub MKV(){ #refer to the above comments for the below subroutine; it is repeated for the Mann Kendall variance test
	$VMK=14;
	$vS=0;
	$vI=0;
	$vMKSTD=0;
	$vfactorial=1;
	$vfactorial2=1;
	$VZPass=0;
	for (my $j=0; $j<=($n-1); $j++){
		for (my $i=0; $i<$j; $i++){
			if ($sampleval[$i][2]>$sampleval[$j][2]){
				#$I=$I-1;
				next;
			}
			elsif ($sampleval[$i][2]<$sampleval[$j][2]){
				$vI=$vI+1;
			}
			else {
				next;
			}
		}
	}	
	$vMKSTD=(($n*($n-1)*(2*$n+5))/18)**.5;
	for (my $i=1; $i<=$n; $i++){
		$vfactorial=$vfactorial*$i;
	}
	for (my $i=1; $i<=($n-2); $i++){
		$vfactorial2=$vfactorial2*$i;
	}
	my $vcombination=($vfactorial)/(2*$vfactorial2);
	$vS=(2*$vI)-$vcombination;
	$vdeviate=$vS/$vMKSTD;
	if (abs($vdeviate)<$Zscore){
		$VZPass=1;
	}
# Determines if Test Passes. The null hypothesis reveals no trend in segment means. The test passes if the null hypothesis is not rejected.


}
#The subroutine MKV() is the same as MKA(), except it works for variances.
######  Conducts the Shape Test for Normality at $n>=50 ###### HL10072018
my $sfail=0;
my $kfail=0;
my $SkewTEST=30;
my $KurTEST=30;
my $ST=0;
# Globally defines variables that determine pass or fail for Shape Test.
sub SHAPE(){
	$ST=15;
	$sfail=0;
	$kfail=0;
	#Gives variables for failure, and resets failure variables.
	#M2
	my $M2shapeargument=0;
	my $M2=0;
	for (my $i=0; $i<scalar(@sampleval); $i++){
		$M2shapeargument = $M2shapeargument+(($sampleval[$i][1]-$sampleaverage)**2); #finds the 2nd moment according to the Schiferl and Wallace paper equation for the Shape Test
		#this is done by finding difference between each segment mean and sample average, squared, and then finding the summation of those values
	}
	$M2=$M2shapeargument/$n;
	#M3
	my $M3shapeargument=0; #finds the 3rd moment according to the Schiferl and Wallace paper equation for the Shape Test
	my $M3=0;
	for (my $i=0; $i<scalar(@sampleval); $i++){
		$M3shapeargument = $M3shapeargument+(($sampleval[$i][1]-$sampleaverage)**3);
	}
	$M3=$M3shapeargument/$n;
	#M4 
	my $M4shapeargument=0; #finds the 4th moment according to the Schiferl and Wallace paper equation for the Shape Test
	my $M4=0;
	for (my $i=0; $i<scalar(@sampleval); $i++){
		$M4shapeargument = $M4shapeargument+(($sampleval[$i][1]-$sampleaverage)**4);
	}
	$M4=$M4shapeargument/$n;
	#Test Statistic calculations
	my $rootB1=$M3/($M2**1.5); #calculates the skewness as defined by Schiferl and Wallace
	my $B2=$M4/($M2**2); #calculates the kurtosis as defined by Schiferl and Wallace
	my $excessB2=$B2-3;
	#Z Statistics for Skewness and Kurtosis
	my $SES=0;
	my $sampleskew=0;
	my $G1=0;
	$SkewTEST=0;
	#Defines standard error of skewness and the sample skewness constant.
	$sampleskew=((($n)*($n-1))**0.5)/($n-2);
	$SES=((6*$n*($n-1))/(($n-2)*($n+1)*($n+3)))**0.5;
	$G1=$sampleskew*$rootB1;
	$SkewTEST= $G1/$SES;
	my $G2=(($n-1)*(($n+1)*$excessB2+6))/(($n-2)*($n-3));
	my $SEK=2*$SES*(($n**2-1)/(($n-3)*($n+5)))**0.5;
	$KurTEST= $G2/$SEK;
	# Calculates the standard normal deviates for Kurtosis Test and Skeness Test according to User's Guide.
	if (abs($SkewTEST)<$Zscore){
		$sfail=1;
	}
	if (abs($KurTEST)<$Zscore){
		$kfail=1;
	}
	# Null hypothesis for both Skewness and Kurtosis Test is that there is no skewness or excess kurtosis that is not normal. If the test passes, there is no evidence of non-normality when considering skewness or excess kurtosis.
	# https://brownmath.com/stat/shape.htm#Skewness
}

######  Conducts the Von-Neumann Test for lack of Correlation ###### HL10072018
my $VNfail=0;
my $VNnormaldeviate=30;
my $VNF=0;
# Global Variables used to give information on tests variables that indicate success or failure.
sub VN(){
	$VNF=16;
	$VNfail=0;
	my $VNargument;
	my $VNq2;
	my $VNr;
	my $VNstdev;
	$VNstdev=((1+(1/($n-1)))*(1/($n+1)))**0.5;
	#This is the standard deviation of the von neumann test
	for (my $i=0; $i<(scalar(@sampleval)-1); $i++) {
		$VNargument+=($sampleval[$i+1][1]-$sampleval[$i][1])**2;
	}
	#This cycles through the get the VN argument
	$VNq2=($VNargument/((2*$n)-2));
	#This gives us the first test statistic.
	$VNr=$VNq2/($samplevar);
	#This converts the test statistic to r that is used.
	$VNnormaldeviate=($VNr-1)/$VNstdev;
	#This is the standard normal deviate to determine if it passes the VN test. Basically the test statistics.
	if (abs($VNnormaldeviate)<$Zscore){
		$VNfail++;
	}	
	#The null hypothesis is that there is no correlation between the segment averages. If the test passes, there is no evidence that the null hypothesis is rejected.
	
}
DATAMOD();
DATACALL();
my $fail2=0;
######  Conducts Steps A and B ###### HL10072018
sub AUTOTEST(){
	$ststep=$oststep;
	Z_and_T_CONVERSION();
	while ($ttrue==0 or $ztrue==0){
		print "Sorry. That alpha level cannot be used with The program. Please choose a new alpha level. \n";
		print "What is the alpha level? This is probability of making the wrong decision when the null hypothesis is true. Allowed values are 0.01, 0.02, 0.05 (default), 0.10, and 0.50. \n";
		$alevel=<STDIN>;
		if ($alevel < 0){
			$alevel=0.05;
			print "Negative value was entered. The alpha value is set to 0.05.\n";
			last;
		}
		if ($alevel eq ""){
			$alevel = "0.05";
			chomp $alevel;
			last;
		}
		Z_and_T_CONVERSION();
	}
	#CHANGE	
		SAMPLE();
		MKA();
		MKV();
	while ((my $i=0) ==0){
		if ($n<24) {
			$fail2=1;
			last;
		}
		if ($fail>0){
			last;
		}
		if ($AZPass ==1 and $VZPass == 1){
			$i++;
			last;
		}
		elsif ($AZPass !=1 or $VZPass !=1){
			$ststep+=1;
			SAMPLE();
			MKA();
			MKV();
		}
		else {
			next;
		}	
		}
	if ($fail>0){
	print "\nError! Your Start Time was too late, please choose an earlier start time.\n";
		$fail2++;
	}
	elsif ($AZPass !=1 or $VZPass !=1){
	print "\nError! Equilibration either did not occur, or the program did not detect equilibration. \n";
	}
	elsif ($AZPass == 1 and $VZPass ==1){
	#print "Hooray! Your data is equilibrated. There is some diagnostic information above, and relevant data in the EquCheck.out in the folder where this program ran from. \n";
	}
	#print "The Start Time is NOW $ststep time intervals\n";
}
#This subroutine is an automated method to get start time. It covers tests a and b and determines whether or not both have passed as in the powerpoint.
#If the man-kendall tests for averages and means failed. It indicates that the start time determined was too early. This delays it.
######  Conducts Tests C and D ###### HL10072018
my $lop;
my $ntor=0;
sub TEST(){
	AUTOTEST();
	$lop=0;
	while ($lop==0){
	if ($n<24){
		$fail=1;
		$lop=1;
		last;
	}
	elsif ($n>=24){
	if ($AZPass == 1 and $VZPass ==1){
		if ($n<=50){
			SW();
			VN();
		}
		if ($n>50){
			SHAPE();
			VN();
		}
		if ($n>50){
			while (($kfail!=1 or $sfail!=1 or $VNfail==0) and $fail==0){
				$tor++;
				AUTOTEST();
				SHAPE();
				VN();
			}
		$lop=1;
		last;
		}
		if ($n<=50){
			while (($WPass!=1 or $VNfail==0) and $fail==0){
				$tor++;
				AUTOTEST();
				SW();
				VN();
			}
		$lop=1;
		last;
		}
	}
	else{
		$lop=1;
		last;
	}
	}
	}
	$ntor=$tor;
}
#This subroutine automates the process for steps c and d. It covers von neumann tests and the normality tests. If any fail, it indicates the length of each segment is


#too small and we need to increase it, automated here. It appends the AUTOTEST() subroutene with new informaiton.
my $Error;
# Globally defined value for the Error within the sample average.
######  Automates the EquCheck Program###### HL10072018
sub OVERALL(){
	TEST();
	my $negZ=$Zscore*-1;
	$Error= $Tscore*($samplevar/$n)**0.5;
	if ($n<24){
		$fail=1;
		last;
	}
	elsif ($n>=24){
	if ($response eq $no){
		while ((my $stop=0) ==0){
			if ($fail==0){
		if ($lvlr eq $inter or $lvlr eq $adv){
				print FH "\n\nNEW RUN \n";
				chomp $n;
				print FH "Sample Size = $n\n";
				print FH "Segment Length (divided by Time Interval) = $tor data points\n";
				print FH "The Start Time Interval = $ststep time intervals\n";
				print FH "Alpha Level = $alevel\n";
				print FH "The Z Score used is $Zscore\n";
				print FH "The T Score used is $Tscore\n";
		if ($deviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print FH "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print FH "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print FH "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print FH "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print FH "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print FH "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print FH "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print FH "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print FH "The Test Statistic for Shapiro Wilk is $W";
				print FH "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print FH "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print FH "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}		
	}	
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
		print "\n\nNEW RUN \n";
				chomp $n;
				print "Sample Size = $n\n";
				print "Segment Length (divided by Time Interval) = $tor data points\n";
				print "The Start Time Interval = $ststep time intervals\n";
				print "Alpha Level = $alevel\n";
				print "The Z Score used is $Zscore\n";
				print "The T Score used is $Tscore\n";
		if ($deviate==30){
			print "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print "The Test Statistic for Shapiro Wilk is $W";
				print "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
}	
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#	
			$deviate=30;
			$vdeviate=30;
			if ($n<=50){
				$W=30;
			}
			if ($n>50){
				$KurTEST=30;
				$SkewTEST=30;
			}
			$VNnormaldeviate=30;			
				$n++;
				$tor=$otor;
				TEST();
			}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
		
			if ($fail!=0){
			if ($lvlr eq $inter or $lvlr eq $adv){	
				print FH "\n\nNEW RUN \n";
				chomp $n;
				print FH "Sample Size = $n\n";
				print FH "Segment Length (divided by Time Interval) = $tor data points\n";
				print FH "The Start Time Interval = $ststep time intervals\n";
				print FH "Alpha Level = $alevel\n";
				print FH "The Z Score used is $Zscore\n";
				print FH "The T Score used is $Tscore\n";
		if ($deviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print FH "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH	 "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH	 "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print FH "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH	 "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH	 "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print FH "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print FH "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH	 "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH	 "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print FH "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print FH "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH	 "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH	 "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print FH "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print FH "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print FH "The Test Statistic for Shapiro Wilk is $W";
				print FH	 "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print FH "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print FH "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH	 "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH	 "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#	
		print "\n\nNEW RUN \n";
				chomp $n;
				print "Sample Size = $n\n";
				print "Segment Length (divided by Time Interval) = $tor data points\n";
				print "The Start Time Interval = $ststep time intervals\n";
				print "Alpha Level = $alevel\n";
				print "The Z Score used is $Zscore\n";
				print "The T Score used is $Tscore\n";
		if ($deviate==30){
			print "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print "The Test Statistic for Shapiro Wilk is $W";
				print "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
}	
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
			$deviate=30;
			$vdeviate=30;
			if ($n<=50){
				$W=30;
			}
			if ($n>50){
				$KurTEST=30;
				$SkewTEST=30;
			}
			$VNnormaldeviate=30;
				$n--;
				$tor=$otor;
				TEST();
				$stop++;
				last;
			}
		}
	}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
	elsif ($response eq $yes){
		while (((my $stop=0) ==0) and ($n<$n2)){
		if ($fail==0){
		if ($lvlr eq $inter or $lvlr eq $adv){
				print FH "\n\nNEW RUN \n";
				chomp $n;
				print FH "Sample Size = $n\n";
				print FH "Segment Length (divided by Time Interval) = $tor data points\n";
				print FH "The Start Time Interval = $ststep time intervals\n";
				print FH "Alpha Level = $alevel\n";
				print FH "The Z Score used is $Zscore\n";
				print FH "The T Score used is $Tscore\n";
		if ($deviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print FH "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print FH "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print FH "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print FH "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print FH "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print FH "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print FH "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print FH "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print FH "The Test Statistic for Shapiro Wilk is $W";
				print FH "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print FH "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print FH "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#
		print "\n\nNEW RUN \n";
				chomp $n;
				print "Sample Size = $n\n";
				print "Segment Length (divided by Time Interval) = $tor data points\n";
				print "The Start Time Interval = $ststep time intervals\n";
				print "Alpha Level = $alevel\n";
				print "The Z Score used is $Zscore\n";
				print "The T Score used is $Tscore\n";
		if ($deviate==30){
			print "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print "The Test Statistic for Shapiro Wilk is $W";
				print "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
	}	
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#	
			$deviate=30;
			$vdeviate=30;
			if ($n<=50){
				$W=30;
			}
			if ($n>50){
				$KurTEST=30;
				$SkewTEST=30;
			}
			$VNnormaldeviate=30;			
				$n++;
				$tor=$otor;
				TEST();
		}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#

			if ($fail!=0){
			if ($lvlr eq $inter or $lvlr eq $adv){
			print FH "\n\nNEW RUN \n";
				chomp $n;
				print FH "Sample Size = $n\n";
				print FH "Segment Length (divided by Time Interval) = $tor data points\n";
				print FH "The Start Time Interval = $ststep time intervals\n";
				print FH "Alpha Level = $alevel\n";
				print FH "The Z Score used is $Zscore\n";
				print FH "The T Score used is $Tscore\n";
		if ($deviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print FH "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print FH "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print FH "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print FH "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print FH "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print FH "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print FH "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print FH "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print FH "The Test Statistic for Shapiro Wilk is $W";
				print FH "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print FH "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print FH "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
		print "\n\nNEW RUN \n";
				chomp $n;
				print "Sample Size = $n\n";
				print "Segment Length (divided by Time Interval) = $tor data points\n";
				print "The Start Time Interval = $ststep time intervals\n";
				print "Alpha Level = $alevel\n";
				print "The Z Score used is $Zscore\n";
				print "The T Score used is $Tscore\n";
		if ($deviate==30){
			print "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print "The Test Statistic for Shapiro Wilk is $W";
				print "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
}	
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#	
			$deviate=30;
			$vdeviate=30;
			if ($n<=50){
				$W=30;
			}
			if ($n>50){
				$KurTEST=30;
				$SkewTEST=30;
			}
			$VNnormaldeviate=30;			
				$n--;
				$tor=$otor;
				TEST();
				$stop++;
				last;
			}
		}
	}
	}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
my $startline=$ststep+1;
my $starttime=($ststep+1)*$timestep;
if ($n>=24){
if ($lvlr eq "1" or $lvlr eq $inter or $lvlr eq $adv){	
	if ($fail==0){
		print FH "\n\n\nSUMMARY OF FINAL RESULTS\n";
		print FH "USER INPUT\n";
		print FH "Original Sample Size = $nor\n";
		if ($response eq $yes){
			print FH "Sample Size Cap = $n2\n";
		}
		print FH "Time Interval = $timestep time units\n";
		print FH "Segment Length = $sl time units\n";
		print FH "Alpha Level = $alevel\n";
		print FH "The User Estimated Start Time = $stime time units\n"; 
		print FH "\n\nYour data is equilibrated! RELEVANT STATISTICAL INFORMATION BELOW \n\n";
		print FH "The start time is $startline (line number) or $starttime (actual time step)\n";
		print FH "The sample average is $sampleaverage with an error of $Error \n";
		print FH "The sample variance is $samplevar\n";
		print FH "The total number of segments is $n\n";
		print FH "The length of each segment is $tor data points";
		print FH "\n\nTHE TEST STATISTICS BELOW \n\n";
		print FH "The alpha level you chose is $alevel\n";
		print FH "The Z Score used is $Zscore\n";
		print FH "The T Score used is $Tscore\n";
		my $negZ=$Zscore*-1;
		if ($deviate<0){
			print FH "\nMann-Kendall Test for Averages\n";
			print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
			print FH	 "$deviate > $negZ, so the test passes.\n";
		}
		elsif ($deviate>=0){
			print FH "\nMann-Kendall Test for Averages\n";
			print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
			print FH	 "$deviate < $Zscore, so the test passes.\n";
		}
		if ($vdeviate<0){
			print FH "\nMann-Kendall Test for Variance\n";
			print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
			print FH	 "$vdeviate > $negZ, so the test passes.\n";
		}
		elsif ($vdeviate>=0){
			print FH "\nMann-Kendall Test for Variance\n";
			print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
			print FH	 "$vdeviate < $Zscore, so the test passes.\n";
		}
		if ($n<=50){
			print FH "\nShapiro-Wilk Test for Normality\n";
			print FH "The Test Statistic for Shapiro Wilk is $W\n";
			print FH	 "$W < $Wcomparison, so the test passes.\n";
		}
		if ($n>50){
			if ($SkewTEST<0){
				print FH "\nShape Test for Normality\n";
				print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
				print FH	 "$SkewTEST > $negZ, so the test passes.\n";
			}
			elsif ($SkewTEST>=0){
				print FH "\nShape Test for Normality\n";
				print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
				print FH	 "$SkewTEST < $Zscore, so the test passes.\n";
			}
			if ($KurTEST<0){
				print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
				print FH	 "$KurTEST > $negZ, so the test passes.\n";
			}
			elsif ($KurTEST>=0){
				print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
				print FH	 "$KurTEST < $Zscore, so the test passes.\n";
			}
		}
			if ($VNnormaldeviate<0){
				print FH "\nThe Von Neumann Test for Trend in Segment Averages\n";
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH	 "$VNnormaldeviate > $negZ, so the test passes.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print FH "\nThe Von Neumann Test for Trend in Segment Averages\n";
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH	 "$VNnormaldeviate < $Zscore, so the test passes.\n";
			}
		close FH;

	
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
print "\n\n\nSUMMARY OF FINAL RESULTS\n";
print "USER INPUT\n";
		print "Original Sample Size = $nor\n";
		if ($response eq $yes){
			print "Sample Size Cap = $n2\n";
		}
		print "Time Interval = $timestep time units\n";
		print "Segment Length = $sl time units\n";
		print "Alpha Level = $alevel\n";
		print "The User Estimated Start Time = $stime time units\n"; 
		print "\n\nYour data is equilibrated! RELEVANT STATISTICAL INFORMATION BELOW \n\n";
		print "The start time is $startline (line number) or $starttime (actual time step)\n";
		print "The sample average is $sampleaverage with an error of $Error \n";
		print "The sample variance is $samplevar\n";
		print "The total number of segments is $n\n";
		print "The length of each segment is $tor data points";
		print "\n\nTHE TEST STATISTICS BELOW \n\n";
		print "The alpha level you chose is $alevel\n";
		print "The Z Score used is $Zscore\n";
		print "The T Score used is $Tscore\n";
		$negZ=$Zscore*-1;
		if ($deviate<0){
			print "\nMann-Kendall Test for Averages\n";
			print "The Mann-Kendall Statistic for Averages is $deviate\n";
			print "$deviate > $negZ, so the test passes.\n";
		}
		elsif ($deviate>=0){
			print "\nMann-Kendall Test for Averages\n";
			print "The Mann-Kendall Statistic for Averages is $deviate\n";
			print "$deviate < $Zscore, so the test passes.\n";
		}
		if ($vdeviate<0){
			print "\nMann-Kendall Test for Variance\n";
			print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
			print "$vdeviate > $negZ, so the test passes.\n";
		}
		elsif ($vdeviate>=0){
			print "\nMann-Kendall Test for Variance\n";
			print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
			print "$vdeviate < $Zscore, so the test passes.\n";
		}
		if ($n<=50){
			print "\nShapiro-Wilk Test for Normality\n";
			print "The Test Statistic for Shapiro Wilk is $W\n";
			print "$W < $Wcomparison, so the test passes.\n";
		}
		if ($n>50){
			if ($SkewTEST<0){
				print "\nShape Test for Normality\n";
				print "The Shape Test Statistic for Skewness is $SkewTEST\n";
				print "$SkewTEST > $negZ, so the test passes.\n";
			}
			elsif ($SkewTEST>=0){
				print "\nShape Test for Normality\n";
				print "The Shape Test Statistic for Skewness is $SkewTEST\n";
				print "$SkewTEST < $Zscore, so the test passes.\n";
			}
			if ($KurTEST<0){
				print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
				print "$KurTEST > $negZ, so the test passes.\n";
			}
			elsif ($KurTEST>=0){
				print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
				print "$KurTEST < $Zscore, so the test passes.\n";
			}
		}
			if ($VNnormaldeviate<0){
				print "\nThe Von Neumann Test for Trend in Segment Averages\n";
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate > $negZ, so the test passes.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print "\nThe Von Neumann Test for Trend in Segment Averages\n";
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate < $Zscore, so the test passes.\n";
			}
}
}
}
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
if ($lvlr eq "1" or $lvlr eq $inter or $lvlr eq $adv){	
	if ($fail!=0){
		my $negZ=$Zscore*-1;
		print FH "\n\n\nSUMMARY OF FINAL RESULTS\n";
		print FH "USER INPUT\n";
		print FH "Sample Size = $nor\n";
		if ($response eq $yes){
			print FH "Sample Size Cap = $n2\n";
		}
		print FH "Time Interval = $timestep time units\n";
		print FH "Segment Length = $sl time units\n";
		print FH "Alpha Level = $alevel\n";
		print FH "The User Estimated Start Time = $stime time units\n"; 
		print FH "\nError. It seems equilibration may not have occurred. RELEVANT DATA BELOW.\n\n";
		print FH "The alpha level you chose is $alevel\n";
		print FH "The Z Score used is $Zscore\n";
		print FH "The T Score used is $Tscore\n";
		if ($deviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print FH "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH	 "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print FH "The Mann-Kendall Statistic for Averages is $deviate\n";
				print FH	 "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print FH "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print FH "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH	 "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print FH "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print FH	 "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print FH "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print FH "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH	 "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print FH "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print FH	 "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print FH "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print FH "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH	 "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print FH "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print FH	 "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print FH "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print FH "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print FH "The Test Statistic for Shapiro Wilk is $W";
				print FH	 "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print FH "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print FH "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH	 "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print FH "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print FH	 "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}

		close FH;
	}

#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#

	if ($fail!=0){
		my $negZ=$Zscore*-1;
	print "\n\n\nSUMMARY OF FINAL RESULTS\n";
		print "USER INPUT\n";
		print "Sample Size = $nor\n";
		if ($response eq $yes){
			print "Sample Size Cap = $n2\n";
		}
		print "Time Interval = $timestep time units\n";
		print "Segment Length = $sl time units\n";
		print "Alpha Level = $alevel\n";
		print "The User Estimated Start Time = $stime time units\n"; 
		print "\nError. It seems equilibration may not have occurred. RELEVANT DATA BELOW.\n\n";
		print "The alpha level you chose is $alevel\n";
		print "The Z Score used is $Zscore\n";
		print "The T Score used is $Tscore\n";
		if ($deviate==30){
			print "\nThe Mann-Kendall Test for Segment Means did not Occur, as a previous test may have failed.\n";
		}
		elsif ($AZPass!=1 and $AMK==12){
			print "\nThe Mann-Kendall Test for Segment Means has failed. Refer to the User Guide for more information. \n";
			if ($deviate<0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate < $negZ, so the test fails.\n";
			}
			elsif ($deviate>=0){
				print "The Mann-Kendall Statistic for Averages is $deviate\n";
				print "$deviate > $Zscore, so the test fails.\n";
			}
		}
		if ($vdeviate==30){
			print "\nThe Mann-Kendall Test for Segment Variances did not occur, as a previous test may have failed.\n";
		}
		elsif ($VZPass!=1 and $VMK==14){
			print "\nThe Mann-Kendall Test for Segment Variances has failed. Refer to the User Guide for more information.  \n";
			if ($vdeviate<0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate < $negZ, so the test fails.\n";
			}
			elsif ($vdeviate>=0){
				print "The Mann-Kendall Statistic for Variance is $vdeviate\n";
				print "$vdeviate > $Zscore, so the test fails.\n";
			}
		}
		if ($n>50){
			if ($KurTEST==30){
				print "\nThe Kurtosis Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($kfail!=1 and $ST=15){
				print "\nThe Kurtosis Test has failed. Refer to the User Guide for more information.  \n";
				if ($KurTEST<0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST < $negZ, so the test fails.\n";
				}
				elsif ($KurTEST>=0){
					print "The Shape Test Statistic for Kurtosis is $KurTEST\n";
					print "$KurTEST > $Zscore, so the test fails.\n";
				}
			}
			if ($SkewTEST==30){
				print "\nThe Skewness Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($sfail!=1 and $ST=15){
				print "\nThe Skewness Test has failed. Refer to the User Guide for more information.  \n";
				if ($SkewTEST<0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST < $negZ, so the test fails.\n";
				}
				elsif ($SkewTEST>=0){
					print "The Shape Test Statistic for Skewness is $SkewTEST\n";
					print "$SkewTEST > $Zscore, so the test fails.\n";
				}
			}
		}	
		if ($n<=50){
			if ($W==30){
				print "\nThe Shapiro-Wilk Test did not occur, as a previous test may have failed.\n";
			}
			elsif ($WPass!=1 and $SWfail==13){
				print "\nThe Shapiro-Wilk Test has failed. Refer to the User Guide for more information.  \n";
				print "The Test Statistic for Shapiro Wilk is $W";
				print "$W > $Wcomparison, so the test fails.\n";
			}
		}
		if ($VNnormaldeviate==30){
				print "\nThe Von-Neumann Test did not occur, as a previous test may have failed.\n";
			}
		elsif ($VNfail==0 and $VNF==16){
			print "\nThe Von-Neumann Test has failed. Refer to the User Guide for more information. \n";
			if ($VNnormaldeviate<0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate < $negZ, so the test fails.\n";
			}
			elsif ($VNnormaldeviate>=0){
				print "The Von Neumann Test Statistic is $VNnormaldeviate\n";
				print "$VNnormaldeviate > $Zscore, so the test fails.\n";
			}
		}
	}
}
}
#This subroutine overall takes each parts of the test, and conducts it automatically given the data that was called earlier.
#OVERALL(), takes tests a,b,c, and d and increases the number of segments accordingly to determine whether or not the data is complete.
OVERALL();
