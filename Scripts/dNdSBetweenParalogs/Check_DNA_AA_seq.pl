#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub BuildingGeneticCode{
    my $hashDNA=$_[0];
    
	my   @AAs=qw/F F L L S S S S Y Y 1 1 C C 1 W L L L L P P P P H H Q Q R R R R I I I M T T T T N N K K S S R R V V V V A A A A D D E E G G G G/;
	my @Base1=qw/T T T T T T T T T T T T T T T T C C C C C C C C C C C C C C C C A A A A A A A A A A A A A A A A G G G G G G G G G G G G G G G G/;
	my @Base2=qw/T T T T C C C C A A A A G G G G T T T T C C C C A A A A G G G G T T T T C C C C A A A A G G G G T T T T C C C C A A A A G G G G/;
	my @Base3=qw/T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G/;
	
	for(my$c=0;$c<=scalar(@AAs);$c++){
		#print $Base1[$c].$Base2[$c].$Base3[$c]."\t".$AAs[$c]."\n";
		$hashDNA->{$Base1[$c].$Base2[$c].$Base3[$c]}=$AAs[$c];
	}
}

sub ReadFasta{
    my $fasta=$_[0];
    my $pathin=$_[1];

    my $input;

    my @ext=split("\\.",$pathin);
    if($ext[-1] eq "gz"){
		open($input,"zcat $pathin |");
    }else{
		open($input, $pathin);
    }
    
    my $line=<$input>;
    while($line){
    	my $Sequence = "";
    	my $ID = "";
		chomp $line;
	    my $fstchar=substr $line, 0,1;
	    if($fstchar eq ">"){
			$ID = substr $line, 1;
			$line=<$input>;
			$fstchar=substr $line, 0,1;
	    }
	    while($fstchar ne ">"){
	    	chomp $line;
	    	$Sequence = $Sequence.$line;
			$line=<$input>;
			$fstchar=substr $line, 0,1;	    	
	    }
	    $fasta->{$ID}=$Sequence;
	}
    
    close($input);
}

sub ExtractSequenceFromFasta{
	my $id=$_[0];
	my $seq=$_[1];
	my $pathfasta=$_[2];

	print $id."\n".$seq."\n".$pathfasta."\n";
#	my $cmd = "awk -v id=\"".$id."\" '{if(\$1 == \">\"id){print \$0; valid=1;next;} if(valid==1){if(\$1 ~ />/){valid=0}else{print \$0}}}' ".$pathfasta;    
#	$seq = system_bash($cmd);
#	chomp $seq;
#	print $seq."\n";
#	return($seq);
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

sub system_bash {
  my @args = ( "bash", "-c", shift );
  system(@args);
}

##########################################################################
##########################################################################

## parameters 

my %parameters;
$parameters{"DNAfasta"}="NA";
$parameters{"AAfasta"}="NA";
$parameters{"output"}="NA";

my @defaultpars=("DNAfasta", "AAfasta", "output");


my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
	
    if(exists $parameters{$parname}){
		$parameters{$parname}=$parval;
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

##############################################################
##############################################################

print "Building genetic code hash...\n";

my %GeneticCode;

BuildingGeneticCode(\%GeneticCode);

print "Genetic code build.\n\n";

##############################################################

print "Reading DNA fasta file...\n";

my %DNAseq;

#ReadFasta(\%DNAseq, $parameters{"DNAfasta"});

print "DNA fasta charged\n\n";

#print "BLAG02000806\n".$DNAseq{"BLAG02000806"}."\n";

##############################################################

print "For each Protein...\n";

my $cmd = "grep '>' ".$parameters{"DNAfasta"}." | sed 's/>//g' | head"; 
my @Protlist = system_bash($cmd);
chomp @Protlist;

foreach my $protID (@Protlist){
	#print $protID."\n";
	my $DNAseq = "buu";
	#ExtractSequenceFromFasta($protID, $DNAseq, $parameters{"DNAfasta"});
	#print $DNAseq."\n";
}

##############################################################

print "\nDone.\n";

##############################################################

