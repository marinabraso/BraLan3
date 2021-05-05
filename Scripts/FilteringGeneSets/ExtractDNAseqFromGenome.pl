#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];
        $id =~ s/\.[0-9]+//g;
	    # print "saw chromosome ".$id."\n";
	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

#########################################################################################

sub reverseComplement{
    my $sequence=$_[0];
    
    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
}

###############################################################################################

sub readCoordinates{
    my $pathin=$_[0];
    my $blocks=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $gene=$s[0];
	my $chr=$s[1];
	my $start=$s[2];
	my $end=$s[3];
	my $strand=$s[4];
    if($strand==1){$strand="+";}
    if($strand==-1){$strand="-";}

	if(!(exists $blocks->{$gene})){
	    $blocks->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[], "end"=>[]};
	}

	push(@{$blocks->{$gene}{"start"}}, $start);
	push(@{$blocks->{$gene}{"end"}}, $end);
	
	$line=<$input>;
    }

    close($input);
} 

#############################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts cDNA fasta sequences. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################################
##########################################################################################

my %parameters;

$parameters{"pathCoordinates"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathCoordinates", "pathGenomeSequence", "pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

#####################################################################################
#####################################################################################

print "Reading coordinates...\n";

my %Coordinates;
readCoordinates($parameters{"pathCoordinates"}, \%Coordinates);
my $nbg=keys %Coordinates;

print "Found ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Reading genome sequence...\n";
my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);
print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

foreach my $geneid (keys %Coordinates){
    my $chr=$Coordinates{$geneid}{"chr"};
    my $strand=$Coordinates{$geneid}{"strand"};
    
    if(exists $genome{$chr}){
	my $nbexons=@{$Coordinates{$geneid}{"start"}};

	my $sequence="";
     	
	for(my $i=0; $i<$nbexons; $i++){
	    my $start=${$Coordinates{$geneid}{"start"}}[$i];
	    my $end=${$Coordinates{$geneid}{"end"}}[$i];
	    
	    my $thisseq=substr $genome{$chr}, ($start-1), ($end-$start+1);
	    $sequence.=$thisseq;
	}
	
	if($strand eq "-"){
	    $sequence=reverseComplement($sequence);
	}elsif($strand ne "+"){
        print "Error: strand or wrong annotation of strand for gene ".$geneid.".\n";
        exit;
    }
	
	print $output ">".$geneid."\n";
    print $output $sequence."\n";

    } else{
	print "Error: cannot find ".$chr." in genome.\n";
    }
}

close($output);

print "Done.\n";

##############################################################
