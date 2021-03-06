# $Id: BPpsilite.pm,v 1.22 2002/10/22 07:38:45 lapp Exp $
# Bioperl module Bio::Tools::BPpsilite
############################################################
#	based closely on the Bio::Tools::BPlite modules
#	Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
#	Lorenz Pollak (lorenz@ist.org, bioperl port)
#
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# October 20, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::BPpsilite - Lightweight BLAST parser for (iterated) psiblast reports

=head1 SYNOPSIS

  use Bio::Tools::BPpsilite;
  open FH, "t/psiblastreport.out";
  $report = Bio::Tools::BPpsilite->new(-fh=>\*FH);

  # determine number of iterations executed by psiblast
  $total_iterations = $report->number_of_iterations;
  $last_iteration = $report->round($total_iterations);

  # Process only hits found in last iteration ...
   $oldhitarray_ref = $last_iteration->oldhits;
   HIT: while($sbjct = $last_iteration->nextSbjct) {
  	  $id = $sbjct->name;
  	  $is_old =  grep  /\Q$id\E/, @$oldhitarray_ref;
  	  if ($is_old ){next HIT;}
  #  do something with new hit...
  }


=head1 DESCRIPTION

BPpsilite is a package for parsing multiple iteration PSIBLAST
reports.  It is based closely on Ian Korf's BPlite.pm module for
parsing single iteration BLAST reports (as modified by Lorenz Pollak).

Two of the four basic objects of BPpsilite.pm are identical to the
corresponding objects in BPlite - the "HSP.pm" and "Sbjct.pm" objects.
This DESCRIPTION documents only the one new object, the "iteration",
as well as the additional methods that are implemented in BPpsilite
that are not in BPlite.  See the BPlite documentation for information
on the BPlite, SBJCT and HSP objects.

The essential difference between PSIBLAST and the other BLAST programs
(in terms of report parsing) is that PSIBLAST performs multiple
iterations of the BLASTing of the database and the results of all of
these iterations are stored in a single PSIBLAST report.  (For general
information on PSIBLAST see the README.bla file in the standalone
BLAST distribution and references therein). PSIBLAST's use of multiple
iterations imposes additional demands on the report parser: * There
are several iterations of hits.  Many of those hits will be repeated
in more than one iteration.  Often only the last iteration will be of
interest.  * Each iteration will list two different kinds of hits -
repeated hits that were used in the model and newly identified hits -
which may need to be processed in different manners * The total number
of iterations performed is not displayed in the report until (almost)
the very end of the report.  (The user can specify a maximum number of
iterations for the PSIBLAST search, but the program may perform fewer
iterations if convergence is reached)

BPpsilite addresses these issues by offering the following methods:

* The total number of iteration used is given by the method
   number_of_iterations as in:

	$total_iterations = $report->number_of_iterations;

* Results from an arbitrary iteration round can be accessed by using
  the 'round' method:

	$iteration3_report = $report->round(3);

* The ids of the sequences which passed the significance threshold for
  the first time in the "nth" iteration can be identified by using the
  newhits method.  Previously identified hits are identified by using
  the oldhits method, as in:

 	$oldhitarray_ref = $iteration3_report->oldhits;
 	$newhitarray_ref = $iteration3_report->newhits;

BPpsilite.pm should work equally well on reports generated by the
StandAloneBlast.pm local BLAST module as with reports generated by
remote psiblast searches. For examples of usage of BPpsilite.pm, the
user is referred to the BPpsilite.t script in the "t" directory.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Peter Schattner

Email: schattner@alum.mit.edu

=head1 CONTRIBUTORS

Jason Stajich, jason@cgt.mc.duke.edu

=head1 ACKNOWLEDGEMENTS

Based on work of:
Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
Lorenz Pollak (lorenz@ist.org, bioperl port)

=head1 COPYRIGHT

BPlite.pm is copyright (C) 1999 by Ian Korf. 

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

package Bio::Tools::BPpsilite;

use strict;
use vars qw(@ISA);
use Bio::Tools::BPlite::Iteration; #
use Bio::Tools::BPlite::Sbjct; #   Debug code
use Bio::Root::Root; # root interface to inherit from
use Bio::Root::IO;
use Bio::Tools::BPlite; 

@ISA = qw(Bio::Root::Root Bio::Root::IO);

sub new {
  my ($class, @args) = @_; 
  my $self = $class->SUPER::new(@args);
  
  # initialize IO
  $self->_initialize_io(@args);
  $self->{'_tempdir'} = $self->tempdir('CLEANUP' => 1);
  $self->{'QPATLOCATION'} = [];  # Anonymous array of query pattern locations for PHIBLAST
  $self->{'NEXT_ITERATION_NUMBER'} = 1;
  $self->{'TOTAL_ITERATION_NUMBER'} = -1;  # -1 indicates preprocessing not yet done

  if ($self->_parseHeader) {$self->{'REPORT_DONE'} = 0} # there are alignments
  else                     {$self->{'REPORT_DONE'} = 1} # empty report
  
  return $self; # success - we hope!
}

=head2 query

 Title    : query
 Usage    : $query = $obj->query();
 Function : returns the query object
 Returns  : query object
 Args     :

=cut

sub query    {shift->{'QUERY'}}

=head2 qlength

 Title    : qlength
 Usage    : $len = $obj->qlength();
 Function : returns the length of the query 
 Returns  : length of query
 Args     :

=cut

sub qlength  {shift->{'LENGTH'}}

=head2 database

 Title    : database
 Usage    : $db = $obj->database();
 Function : returns the database used in this search
 Returns  : database used for search
 Args     :

=cut

sub database {shift->{'DATABASE'}}

=head2 number_of_iterations

 Title    : number_of_iterations
 Usage    : $total_iterations = $obj-> number_of_iterations();
 Function : returns the total number of iterations used in this search
 Returns  : total number of iterations used for search
 Args     : none

=cut


=head2 pattern

 Title    : database
 Usage    : $pattern = $obj->pattern();
 Function : returns the pattern used in a PHIBLAST search

=cut

sub pattern {shift->{'PATTERN'}}

=head2 query_pattern_location

 Title    : query_pattern_location
 Usage    : $qpl = $obj->query_pattern_location();
 Function : returns reference to array of locations in the query sequence
            of pattern used in a PHIBLAST search

=cut

sub query_pattern_location {shift->{'QPATLOCATION'}}




sub number_of_iterations {
	my $self = shift;
  	if ($self->{'TOTAL_ITERATION_NUMBER'} == -1){&_preprocess($self);}
	$self->{'TOTAL_ITERATION_NUMBER'};
}

=head2 round

 Title    : round
 Usage    : $Iteration3 = $report->round(3);
 Function : Method of retrieving data from a specific iteration 
 Example  :  
 Returns  : reference to requested Iteration object or null if argument
		is greater than total number of iterations
 Args     : number of the requested iteration

=cut

sub round {
  my $self = shift;
  my $iter_num = shift;
  $self->_initialize_io(-file => Bio::Root::IO->catfile
			($self->{'_tempdir'},"iteration".$iter_num.".tmp"));
  if( ! $self->_fh ) {
      $self->throw("unable to re-open iteration file for round ".$iter_num);
  }
  return Bio::Tools::BPlite::Iteration->new(-round=>$iter_num,
					    -parent=>$self);
}

# begin private routines

sub _parseHeader {
  my ($self) = @_;

  
  while(defined ($_ = $self->_readline) ) {
    if ($_ =~ /^Query=\s+([^\(]+)/)    {
      my $query = $1;
      while(defined ($_ = $self->_readline)) {
        last if $_ !~ /\S/;
	$query .= $_;
      }
      $query =~ s/\s+/ /g;
      $query =~ s/^>//;
      $query =~ /\((\d+)\s+\S+\)\s*$/;
      my $length = $1;
      $self->{'QUERY'} = $query;
      $self->{'LENGTH'} = $length;
    }
    elsif ($_ =~ /^Database:\s+(.+)/) {$self->{'DATABASE'} = $1}
    elsif ($_ =~ /^\s*pattern\s+(\S+).*position\s+(\d+)\D/) 
    {   # For PHIBLAST reports
	$self->{'PATTERN'} = $1;
	push (@{$self->{'QPATLOCATION'}}, $2);
    } elsif ($_ =~ /^>|^Results from round 1/)    {
	$self->_pushback($_); 
	return 1;
    } elsif ($_ =~ /^Parameters|^\s+Database:/) {
	$self->_pushback($_); 
	return 0; # there's nothing in the report
    }
  }
}

=head2 _preprocess

 Title    : _preprocess
 Usage    : internal routine, not called directly
 Function : determines number of iterations in report and prepares
	    data so individual iterations canbe parsed in non-sequential 
            order 
 Example  :  
 Returns  : nothing. Sets TOTAL_ITERATION_NUMBER in object's hash
 Args     : reference to calling object

=cut

#'
sub _preprocess {
    my $self = shift;
#	$self->throw(" PSIBLAST report preprocessing not implemented yet!");

    my  $oldround = 0;
    my ($currentline, $currentfile, $round);

# open output file for data from iteration round #1
    $round = 1;
    $currentfile = Bio::Root::IO->catfile($self->{'_tempdir'}, 
					  "iteration$round.tmp");
    open (FILEHANDLE, ">$currentfile") || 
	$self->throw("cannot open filehandle to write to file $currentfile");

    while(defined ($currentline = $self->_readline()) ) {
	if ($currentline =~ /^Results from round\s+(\d+)/) {
	    if ($oldround) { close (FILEHANDLE) ;}
	    $round = $1;
	    $currentfile = Bio::Root::IO->catfile($self->{'_tempdir'}, 
						  "iteration$round.tmp");

	    close FILEHANDLE;
	    open (FILEHANDLE, ">$currentfile") || 
		$self->throw("cannot open filehandle to write to file $currentfile");
	    $oldround = $round;
	}elsif ($currentline =~ /CONVERGED/){ # This is a fix for psiblast parsing with -m 6 /AE
	    $round--;
	}
	print FILEHANDLE $currentline ;
	
    }
    $self->{'TOTAL_ITERATION_NUMBER'}= $round;
# It is necessary to close filehandle otherwise the whole
# file will not be read later !!
    close FILEHANDLE;
}

1;

__END__
