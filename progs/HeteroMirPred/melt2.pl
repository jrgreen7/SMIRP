#! /usr/bin/perl -w

use strict;
use warnings;


use Getopt::Long;
Getopt::Long::Configure 'gnu_getopt', 'no_auto_abbrev', 'no_ignore_case';

use constant R => 0.0019872;

sub systemError ($) {
    die ($? == -1 ? "Error: $! from $_[0]\n" : 'Exit status ' . ($? >> 8) . " from $_[0]\n")
}

sub safelyOpen (*$$) {
    open $_[0], $_[1], $_[2] or die "Can't open $_[2]: $!\n";
}


my ($temp, $Ct, %options) = (37);
GetOptions \%options, 'version|V' => sub { version('hybrid2.pl') }, 'help|h' => sub { usage() }, 'NA|n=s', 'temperature|t=f' => \$temp, 'maxloop=i', 'allpairs', 'sodium|N=f', 'magnesium|M=f', 'Ct|C=f' => \$Ct, 'polymer|p', 'nodangle', 'simple', 'prefilter=s', 'noisolate|I', 'm|maxbp=i', 'circular' or die $!;

my @args;
foreach (keys %options) {
    if ($_ eq 'allpairs' or $_ eq 'circular' or $_ eq 'nodangle' or $_ eq 'noisolate' or $_ eq 'polymer' or $_ eq 'simple' or $_ eq 'zip') {
	push @args, "--$_";
    } else {
	push @args, "--$_";
	push @args, $options{$_} if defined $options{$_};
    }
}

my $file1 = shift;
unless ($file1) {
    print STDERR "Error: file not specified\nRun 'melt.pl -h' for help\n";
    exit 1;
}
my $file2 = shift;
my $prefix = $file1;
$prefix =~ s/\.seq$//;

if ($file2) {
    unless (defined $Ct) {
	die "Error: two files found but no strand concentration given\n";
    }

    $prefix .= "-$file2";
    $prefix =~ s/\.seq$//;
    system('hybrid-min', -t => $temp, -T => $temp, @args, $file1, $file2) == 0 or systemError('hybrid-min');
    print "dG\tdH\tdS\tTm\n";
    safelyOpen *IN, '<', "$prefix.ct";
    while (<IN>) {
	my ($len, $dG, $dH) = /(\d+)\sdG = ([^\s]+)\sdH = ([^\s]+)/ or next;
	my $dS = ($dH - $dG) / (273.15 + $temp);
	my $homo = not ($len % 2);
	my @bases;
	for (my $i = 1; $i <= $len; ++$i) {
	    my $line = scalar <IN>;
	    my (undef, $base, $prev, $next) = split /\s+/, $line;
	    if (($i == $len / 2 and $next) or ($i == $len / 2 + 1 and $prev)) {
		$homo = 0;
	    }
	    if ($i <= $len / 2) {
		push @bases, $base;
	    } else {
		$homo = 0 unless @bases and $base eq shift @bases;
	    }
	}
	my $factor = $homo ? 1 : 4;
	my $Tm = $dH / ($dS + R * log($Ct / $factor));
	printf "%.1f\t%.1f\t%.1f\t%.1f\n", $dG, $dH, 1000 * $dS, $Tm - 273.15;
    }
    close IN or die $!;
} else {
    my $suffix = 'DH';
    $suffix = 'DHD' if exists $options{NA} and $options{NA} eq 'DNA';
    system('./hybrid-ss-min', -t => $temp, -T => $temp, @args, $file1) == 0 or systemError('hybrid-ss-min');
    print "dG\tdH\tdS\tTm\n";
    safelyOpen *DG, '<', "$prefix.dG";
    open DH, '-|', './ct-energy', -s => $suffix, "$prefix.ct" or die "Can't execute ct-energy: $!";
    scalar <DG>;
    while (defined(my $dG = <DG>) and defined (my $dH = <DH>)) {
	$dG = (split /\s+/, $dG)[1];
	chomp $dH;
	my $dS = ($dH - $dG) / (273.15 + $temp);
	my $Tm = $dH / $dS;
	printf "\t%.1f\t%.1f\t%.1f\t%.1f\n", $dG, $dH, 1000 * $dS, $Tm - 273.15;
	#######
	open(FILE3, "> tempmelt.txt");
	printf (FILE3 "%.1f\t%.1f\t%.1f\t%.1f", $dG, $dH, 1000 * $dS, $Tm - 273.15);
	######
   }
    close DH or die $!;
    close DG or die $!;
}
