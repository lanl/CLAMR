package TestHarness::Grammar;
use Parse::RecDescent;

{ my $ERRORS;


package Parse::RecDescent::namespace000001;
use strict;
use vars qw($skip $AUTOLOAD  );
$skip = '\s*';

    use Carp;
    use Data::Dumper;
    use vars qw($suite $debug);

    local $| = 1;
    $debug = 0;   # Set to 1 for debugging.
  ;


{
local $SIG{__WARN__} = sub {0};
# PRETEND TO BE IN Parse::RecDescent NAMESPACE
*Parse::RecDescent::namespace000001::AUTOLOAD	= sub
{
	no strict 'refs';
	$AUTOLOAD =~ s/^Parse::RecDescent::namespace000001/Parse::RecDescent/;
	goto &{$AUTOLOAD};
}
}

push @Parse::RecDescent::namespace000001::ISA, 'Parse::RecDescent';
# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::set_param_line
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"set_param_line"};
	
	Parse::RecDescent::_trace(q{Trying rule: [set_param_line]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{set_param_line})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/[ \\t]* ([\\w_]+) [ \\t]* [=:]* [ \\t*]* ([-._\\w\\d]+) [ \\t]* \\n/x]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{set_param_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{set_param_line});
		%item = (__RULE__ => q{set_param_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/[ \\t]* ([\\w_]+) [ \\t]* [=:]* [ \\t*]* ([-._\\w\\d]+) [ \\t]* \\n/x]}, Parse::RecDescent::_tracefirst($text),
					  q{set_param_line})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[ \t]* ([\w_]+) [ \t]* [=:]* [ \t*]* ([-._\w\d]+) [ \t]* \n)//x)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{set_param_line})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
      my $func = lc $1;
      my $value = $2;
      $suite->$func($value);
   };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/[ \\t]* ([\\w_]+) [ \\t]* [=:]* [ \\t*]* ([-._\\w\\d]+) [ \\t]* \\n/x]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{set_param_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{set_param_line})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{set_param_line})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{set_param_line});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{set_param_line})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::_alternation_1_of_production_2_of_rule_include
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"_alternation_1_of_production_2_of_rule_include"};
	
	Parse::RecDescent::_trace(q{Trying rule: [_alternation_1_of_production_2_of_rule_include]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{_alternation_1_of_production_2_of_rule_include})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/\\+ | include/ix]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_2_of_rule_include})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_2_of_rule_include});
		%item = (__RULE__ => q{_alternation_1_of_production_2_of_rule_include});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/\\+ | include/ix]}, Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_2_of_rule_include})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\+ | include)//ix)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/\\+ | include/ix]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_2_of_rule_include})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{_alternation_1_of_production_2_of_rule_include})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{_alternation_1_of_production_2_of_rule_include})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{_alternation_1_of_production_2_of_rule_include});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{_alternation_1_of_production_2_of_rule_include})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::include
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"include"};
	
	Parse::RecDescent::_trace(q{Trying rule: [include]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{include})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/\\+ | include/ix SUITENAME /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{include})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{include});
		%item = (__RULE__ => q{include});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying repeated subrule: [/\\+ | include/ix]},
				  Parse::RecDescent::_tracefirst($text),
				  q{include})
					if defined $::RD_TRACE;
		$expectation->is(q{})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::_alternation_1_of_production_1_of_rule_include, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [/\\+ | include/ix]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{include})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [_alternation_1_of_production_1_of_rule_include]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$item{q{/\\+ | include/ix}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying subrule: [SUITENAME]},
				  Parse::RecDescent::_tracefirst($text),
				  q{include})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{SUITENAME})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::SUITENAME($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [SUITENAME]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{include})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [SUITENAME]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$item{q{SUITENAME}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    # We are going to parse another suite file. $suite is a package variable
    # which will be changed to the new suite. We must use local to protect the
    # value of $suite.
    my $new_suite;
    {
      my $old_suite = $suite;
      local $suite;
      local $_;
      $new_suite = $old_suite->add_testsuite($item{SUITENAME});
    }
    $return = [@item];
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/\\+ | include/ix SUITENAME /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/\\+ | include/ix TESTNAME CPUS_RUNTIME /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{include})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[1];
		$text = $_[1];
		my $_savetext;
		@item = (q{include});
		%item = (__RULE__ => q{include});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying repeated subrule: [/\\+ | include/ix]},
				  Parse::RecDescent::_tracefirst($text),
				  q{include})
					if defined $::RD_TRACE;
		$expectation->is(q{})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::_alternation_1_of_production_2_of_rule_include, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [/\\+ | include/ix]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{include})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [_alternation_1_of_production_2_of_rule_include]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$item{q{/\\+ | include/ix}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying subrule: [TESTNAME]},
				  Parse::RecDescent::_tracefirst($text),
				  q{include})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{TESTNAME})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::TESTNAME($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [TESTNAME]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{include})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [TESTNAME]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$item{q{TESTNAME}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying repeated subrule: [CPUS_RUNTIME]},
				  Parse::RecDescent::_tracefirst($text),
				  q{include})
					if defined $::RD_TRACE;
		$expectation->is(q{CPUS_RUNTIME})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::CPUS_RUNTIME, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [CPUS_RUNTIME]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{include})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [CPUS_RUNTIME]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$item{q{CPUS_RUNTIME}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    my $test = $suite->add_test($item{TESTNAME});
    #print "CPUS_RUNTIME: ", Dumper ($item{"CPUS_RUNTIME(?)"});
    my $cpus_runtime = $item{"CPUS_RUNTIME(?)"}->[0];
    while (my ($func, $value) = each %$cpus_runtime) {
      #print "func:$func:  value:$value\n";
      $func = lc $func; # functions are always lower case
      $test->$func($value);
    }

    $return = [@item];
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/\\+ | include/ix TESTNAME CPUS_RUNTIME /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{include})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{include})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{include})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{include});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{include})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::SUITENAME
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"SUITENAME"};
	
	Parse::RecDescent::_trace(q{Trying rule: [SUITENAME]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{SUITENAME})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/[-\\/\\w0-9_]+\\.suite/i]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{SUITENAME})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{SUITENAME});
		%item = (__RULE__ => q{SUITENAME});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/[-\\/\\w0-9_]+\\.suite/i]}, Parse::RecDescent::_tracefirst($text),
					  q{SUITENAME})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[-\/\w0-9_]+\.suite)//i)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{SUITENAME})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $item[1];
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/[-\\/\\w0-9_]+\\.suite/i]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{SUITENAME})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{SUITENAME})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{SUITENAME})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{SUITENAME});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{SUITENAME})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::COMMENT
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"COMMENT"};
	
	Parse::RecDescent::_trace(q{Trying rule: [COMMENT]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{COMMENT})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/[ \\t]* ([\\#;\\$] | \\/\\/) [^\\n]*/x]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{COMMENT})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{COMMENT});
		%item = (__RULE__ => q{COMMENT});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/[ \\t]* ([\\#;\\$] | \\/\\/) [^\\n]*/x]}, Parse::RecDescent::_tracefirst($text),
					  q{COMMENT})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[ \t]* ([\#;\$] | \/\/) [^\n]*)//x)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/[ \\t]* ([\\#;\\$] | \\/\\/) [^\\n]*/x]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{COMMENT})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{COMMENT})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{COMMENT})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{COMMENT});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{COMMENT})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::suitefile
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"suitefile"};
	
	Parse::RecDescent::_trace(q{Trying rule: [suitefile]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{suitefile})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [<reject>]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{suitefile})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{suitefile});
		%item = (__RULE__ => q{suitefile});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {$suite = $arg[0]};
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		

		Parse::RecDescent::_trace(q{>>Rejecting production<< (found <reject>)},
					 Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		

		$_tok = undef;
		
		last unless defined $_tok;


		Parse::RecDescent::_trace(q{>>Matched production: [<reject>]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		local $skip = defined($skip) ? $skip : $Parse::RecDescent::skip;
		Parse::RecDescent::_trace(q{Trying production: [<skip: ""> if_else, or line /\\Z/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{suitefile})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[1];
		$text = $_[1];
		my $_savetext;
		@item = (q{suitefile});
		%item = (__RULE__ => q{suitefile});
		my $repcount = 0;


		

		Parse::RecDescent::_trace(q{Trying directive: [<skip: "">]},
					Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE; 
		$_tok = do { my $oldskip = $skip; $skip= ""; $oldskip };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying repeated subrule: [if_else, or line]},
				  Parse::RecDescent::_tracefirst($text),
				  q{suitefile})
					if defined $::RD_TRACE;
		$expectation->is(q{if_else, or line})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::_alternation_1_of_production_2_of_rule_suitefile, 1, 100000000, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [if_else, or line]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{suitefile})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [_alternation_1_of_production_2_of_rule_suitefile]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		$item{q{if_else, or line}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying terminal: [/\\Z/]}, Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\Z/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\Z)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $suite;
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [<skip: ""> if_else, or line /\\Z/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{suitefile})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{suitefile})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{suitefile})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{suitefile});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{suitefile})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::STRING
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"STRING"};
	
	Parse::RecDescent::_trace(q{Trying rule: [STRING]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{STRING})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/["''"]? ([^"''"\\n]*) ["''"]?/x]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{STRING})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{STRING});
		%item = (__RULE__ => q{STRING});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/["''"]? ([^"''"\\n]*) ["''"]?/x]}, Parse::RecDescent::_tracefirst($text),
					  q{STRING})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:["''"]? ([^"''"\n]*) ["''"]?)//x)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{STRING})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $1;
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/["''"]? ([^"''"\\n]*) ["''"]?/x]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{STRING})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{STRING})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{STRING})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{STRING});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{STRING})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::NUM_CPUS
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"NUM_CPUS"};
	
	Parse::RecDescent::_trace(q{Trying rule: [NUM_CPUS]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{NUM_CPUS})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/[0-9]+/ /:/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{NUM_CPUS})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{NUM_CPUS});
		%item = (__RULE__ => q{NUM_CPUS});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/[0-9]+/]}, Parse::RecDescent::_tracefirst($text),
					  q{NUM_CPUS})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[0-9]+)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying terminal: [/:/]}, Parse::RecDescent::_tracefirst($text),
					  q{NUM_CPUS})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/:/})->at($text);
		$_savetext = $text;

		if ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?::)//)
		{
			$text = $_savetext;
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN2__}=$&;
		$text = $_savetext;

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{NUM_CPUS})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
   #print "NUM_CPUS: ", Dumper(@item);
   $return = {$item[0] => $item[1]};
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/[0-9]+/ /:/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{NUM_CPUS})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{NUM_CPUS})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{NUM_CPUS})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{NUM_CPUS});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{NUM_CPUS})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::if_block
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"if_block"};
	
	Parse::RecDescent::_trace(q{Trying rule: [if_block]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{if_block})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/if/i <commit> paren_enclosed bracket_enclosed]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{if_block})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{if_block});
		%item = (__RULE__ => q{if_block});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/if/i]}, Parse::RecDescent::_tracefirst($text),
					  q{if_block})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:if)//i)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		

		Parse::RecDescent::_trace(q{Trying directive: [<commit>]},
					Parse::RecDescent::_tracefirst($text),
					  q{if_block})
						if defined $::RD_TRACE; 
		$_tok = do { $commit = 1 };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying subrule: [paren_enclosed]},
				  Parse::RecDescent::_tracefirst($text),
				  q{if_block})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{paren_enclosed})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::paren_enclosed($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [paren_enclosed]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{if_block})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [paren_enclosed]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{if_block})
						if defined $::RD_TRACE;
		$item{q{paren_enclosed}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying subrule: [bracket_enclosed]},
				  Parse::RecDescent::_tracefirst($text),
				  q{if_block})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{bracket_enclosed})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::bracket_enclosed($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [bracket_enclosed]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{if_block})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [bracket_enclosed]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{if_block})
						if defined $::RD_TRACE;
		$item{q{bracket_enclosed}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{if_block})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    print "\nIF_BLOCK: ", Dumper (%item), "\n\n" if $debug;
    # Foreach bracket_enclosed assign the Judge created in the paren_enclosed.
    my $judge = $item{paren_enclosed};
    foreach my $test (@{$item{bracket_enclosed}}) {
      $test->add_judge($judge);
    }

    # These tests are all taken care of so I don't need to return them, but
    # I need to return the judge.
    $return = $judge;
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/if/i <commit> paren_enclosed bracket_enclosed]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{if_block})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{if_block})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{if_block})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{if_block});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{if_block})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::CPUS_RUNTIME
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"CPUS_RUNTIME"};
	
	Parse::RecDescent::_trace(q{Trying rule: [CPUS_RUNTIME]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{CPUS_RUNTIME})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [NUM_CPUS TIME_LIMIT]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{CPUS_RUNTIME});
		%item = (__RULE__ => q{CPUS_RUNTIME});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [NUM_CPUS]},
				  Parse::RecDescent::_tracefirst($text),
				  q{CPUS_RUNTIME})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::NUM_CPUS($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [NUM_CPUS]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{CPUS_RUNTIME})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [NUM_CPUS]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$item{q{NUM_CPUS}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying repeated subrule: [TIME_LIMIT]},
				  Parse::RecDescent::_tracefirst($text),
				  q{CPUS_RUNTIME})
					if defined $::RD_TRACE;
		$expectation->is(q{TIME_LIMIT})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::TIME_LIMIT, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [TIME_LIMIT]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{CPUS_RUNTIME})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [TIME_LIMIT]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$item{q{TIME_LIMIT}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
       my $time_limit = $item{"TIME_LIMIT(?)"}->[0] || {};
       $return = {%{$item{NUM_CPUS}}, %$time_limit};
     };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [NUM_CPUS TIME_LIMIT]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [TIME_LIMIT NUM_CPUS]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[1];
		$text = $_[1];
		my $_savetext;
		@item = (q{CPUS_RUNTIME});
		%item = (__RULE__ => q{CPUS_RUNTIME});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [TIME_LIMIT]},
				  Parse::RecDescent::_tracefirst($text),
				  q{CPUS_RUNTIME})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::TIME_LIMIT($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [TIME_LIMIT]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{CPUS_RUNTIME})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [TIME_LIMIT]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$item{q{TIME_LIMIT}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying repeated subrule: [NUM_CPUS]},
				  Parse::RecDescent::_tracefirst($text),
				  q{CPUS_RUNTIME})
					if defined $::RD_TRACE;
		$expectation->is(q{NUM_CPUS})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::NUM_CPUS, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [NUM_CPUS]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{CPUS_RUNTIME})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [NUM_CPUS]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$item{q{NUM_CPUS}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
       my $cpus =  $item{"NUM_CPUS(?)"}->[0] || {};
       $return = {%$cpus, %{$item{TIME_LIMIT}}};
     };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [TIME_LIMIT NUM_CPUS]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{CPUS_RUNTIME})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{CPUS_RUNTIME})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{CPUS_RUNTIME});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{CPUS_RUNTIME})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::_alternation_1_of_production_1_of_rule_include
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"_alternation_1_of_production_1_of_rule_include"};
	
	Parse::RecDescent::_trace(q{Trying rule: [_alternation_1_of_production_1_of_rule_include]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{_alternation_1_of_production_1_of_rule_include})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/\\+ | include/ix]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_include})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_include});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_include});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/\\+ | include/ix]}, Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_include})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\+ | include)//ix)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/\\+ | include/ix]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_include})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{_alternation_1_of_production_1_of_rule_include})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{_alternation_1_of_production_1_of_rule_include})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{_alternation_1_of_production_1_of_rule_include});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{_alternation_1_of_production_1_of_rule_include})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::line
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"line"};
	
	Parse::RecDescent::_trace(q{Trying rule: [line]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{line})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		local $skip = defined($skip) ? $skip : $Parse::RecDescent::skip;
		Parse::RecDescent::_trace(q{Trying production: [<skip: "[ \t]*"> print_, or include, or exclude, or BLANK_LINE, or comment_line, or set_param_line]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{line});
		%item = (__RULE__ => q{line});
		my $repcount = 0;


		

		Parse::RecDescent::_trace(q{Trying directive: [<skip: "[ \t]*">]},
					Parse::RecDescent::_tracefirst($text),
					  q{line})
						if defined $::RD_TRACE; 
		$_tok = do { my $oldskip = $skip; $skip= "[ \t]*"; $oldskip };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying subrule: [_alternation_1_of_production_1_of_rule_line]},
				  Parse::RecDescent::_tracefirst($text),
				  q{line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{print_, or include, or exclude, or BLANK_LINE, or comment_line, or set_param_line})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::_alternation_1_of_production_1_of_rule_line($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [_alternation_1_of_production_1_of_rule_line]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [_alternation_1_of_production_1_of_rule_line]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{line})
						if defined $::RD_TRACE;
		$item{q{_alternation_1_of_production_1_of_rule_line}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{line})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {1};
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [<skip: "[ \t]*"> print_, or include, or exclude, or BLANK_LINE, or comment_line, or set_param_line]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{line})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{line})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{line});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{line})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::else_block
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"else_block"};
	
	Parse::RecDescent::_trace(q{Trying rule: [else_block]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{else_block})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/else/i bracket_enclosed]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{else_block})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{else_block});
		%item = (__RULE__ => q{else_block});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/else/i]}, Parse::RecDescent::_tracefirst($text),
					  q{else_block})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:else)//i)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying subrule: [bracket_enclosed]},
				  Parse::RecDescent::_tracefirst($text),
				  q{else_block})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{bracket_enclosed})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::bracket_enclosed($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [bracket_enclosed]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{else_block})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [bracket_enclosed]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{else_block})
						if defined $::RD_TRACE;
		$item{q{bracket_enclosed}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{else_block})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    print "\nELSE_BLOCK: ", Dumper (%item), "\n\n" if $debug;
    # return the tests enclosed in brackets
    $return = $item{bracket_enclosed};
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/else/i bracket_enclosed]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{else_block})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{else_block})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{else_block})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{else_block});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{else_block})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::_alternation_1_of_production_1_of_rule_line
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"_alternation_1_of_production_1_of_rule_line"};
	
	Parse::RecDescent::_trace(q{Trying rule: [_alternation_1_of_production_1_of_rule_line]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [print_]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_line});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [print_]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::print_($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [print_]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_1_of_rule_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [print_]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$item{q{print_}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [print_]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [include]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[1];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_line});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [include]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::include($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [include]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_1_of_rule_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [include]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$item{q{include}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [include]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [exclude]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[2];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_line});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [exclude]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::exclude($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [exclude]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_1_of_rule_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [exclude]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$item{q{exclude}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [exclude]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [BLANK_LINE]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[3];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_line});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [BLANK_LINE]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::BLANK_LINE($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [BLANK_LINE]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_1_of_rule_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [BLANK_LINE]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$item{q{BLANK_LINE}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [BLANK_LINE]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [comment_line]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[4];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_line});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [comment_line]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::comment_line($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [comment_line]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_1_of_rule_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [comment_line]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$item{q{comment_line}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [comment_line]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [set_param_line]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[5];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_1_of_rule_line});
		%item = (__RULE__ => q{_alternation_1_of_production_1_of_rule_line});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [set_param_line]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::set_param_line($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [set_param_line]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_1_of_rule_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [set_param_line]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$item{q{set_param_line}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [set_param_line]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{_alternation_1_of_production_1_of_rule_line})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{_alternation_1_of_production_1_of_rule_line})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{_alternation_1_of_production_1_of_rule_line});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{_alternation_1_of_production_1_of_rule_line})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::bracket_enclosed
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"bracket_enclosed"};
	
	Parse::RecDescent::_trace(q{Trying rule: [bracket_enclosed]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{bracket_enclosed})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		local $skip = defined($skip) ? $skip : $Parse::RecDescent::skip;
		Parse::RecDescent::_trace(q{Trying production: [<skip: ""> /\{/ COMMENT /[ \\t]*\\n/ line /\}/ COMMENT /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{bracket_enclosed});
		%item = (__RULE__ => q{bracket_enclosed});
		my $repcount = 0;


		

		Parse::RecDescent::_trace(q{Trying directive: [<skip: "">]},
					Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE; 
		$_tok = do { my $oldskip = $skip; $skip= ""; $oldskip };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying terminal: [/\{/]}, Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\{/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:{)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying repeated subrule: [COMMENT]},
				  Parse::RecDescent::_tracefirst($text),
				  q{bracket_enclosed})
					if defined $::RD_TRACE;
		$expectation->is(q{COMMENT})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::COMMENT, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [COMMENT]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{bracket_enclosed})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [COMMENT]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$item{q{COMMENT}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying terminal: [/[ \\t]*\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/[ \\t]*\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[ \t]*\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN2__}=$&;
		

		Parse::RecDescent::_trace(q{Trying repeated subrule: [line]},
				  Parse::RecDescent::_tracefirst($text),
				  q{bracket_enclosed})
					if defined $::RD_TRACE;
		$expectation->is(q{line})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::line, 0, 100000000, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [line]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{bracket_enclosed})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [line]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$item{q{line}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying terminal: [/\}/]}, Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\}/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:})//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN3__}=$&;
		

		Parse::RecDescent::_trace(q{Trying repeated subrule: [COMMENT]},
				  Parse::RecDescent::_tracefirst($text),
				  q{bracket_enclosed})
					if defined $::RD_TRACE;
		$expectation->is(q{COMMENT})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::COMMENT, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [COMMENT]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{bracket_enclosed})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [COMMENT]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$item{q{COMMENT}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN4__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    # Just return the list of problems/suites. These are objects.
    $return = $item{'line(s?)'};
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [<skip: ""> /\{/ COMMENT /[ \\t]*\\n/ line /\}/ COMMENT /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{bracket_enclosed})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{bracket_enclosed})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{bracket_enclosed});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{bracket_enclosed})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::TIME_LIMIT
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"TIME_LIMIT"};
	
	Parse::RecDescent::_trace(q{Trying rule: [TIME_LIMIT]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{TIME_LIMIT})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/(?:(\\d+) :)? (\\d+)/x]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{TIME_LIMIT})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{TIME_LIMIT});
		%item = (__RULE__ => q{TIME_LIMIT});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/(?:(\\d+) :)? (\\d+)/x]}, Parse::RecDescent::_tracefirst($text),
					  q{TIME_LIMIT})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:(?:(\d+) :)? (\d+))//x)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{TIME_LIMIT})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
   #print "TIME_LIMIT: ", Dumper(@item);
   $return = {$item[0] => $item[1]};
 };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/(?:(\\d+) :)? (\\d+)/x]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{TIME_LIMIT})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{TIME_LIMIT})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{TIME_LIMIT})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{TIME_LIMIT});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{TIME_LIMIT})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::print_
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"print_"};
	
	Parse::RecDescent::_trace(q{Trying rule: [print_]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{print_})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/print|error|warn(ing)?/i STRING /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{print_})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{print_});
		%item = (__RULE__ => q{print_});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/print|error|warn(ing)?/i]}, Parse::RecDescent::_tracefirst($text),
					  q{print_})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:print|error|warn(ing)?)//i)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying subrule: [STRING]},
				  Parse::RecDescent::_tracefirst($text),
				  q{print_})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{STRING})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::STRING($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [STRING]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{print_})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [STRING]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{print_})
						if defined $::RD_TRACE;
		$item{q{STRING}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{print_})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN2__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{print_})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $suite->add_test('PRINT', $item{STRING});
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/print|error|warn(ing)?/i STRING /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{print_})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{print_})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{print_})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{print_});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{print_})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::paren_enclosed
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"paren_enclosed"};
	
	Parse::RecDescent::_trace(q{Trying rule: [paren_enclosed]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{paren_enclosed})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		local $skip = defined($skip) ? $skip : $Parse::RecDescent::skip;
		Parse::RecDescent::_trace(q{Trying production: [<skip: "">]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{paren_enclosed})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{paren_enclosed});
		%item = (__RULE__ => q{paren_enclosed});
		my $repcount = 0;


		

		Parse::RecDescent::_trace(q{Trying directive: [<skip: "">]},
					Parse::RecDescent::_tracefirst($text),
					  q{paren_enclosed})
						if defined $::RD_TRACE; 
		$_tok = do { my $oldskip = $skip; $skip= ""; $oldskip };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{paren_enclosed})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {extract_bracketed($text,'()')};
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{paren_enclosed})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    # create a new Judge based on the condition/testproblem set contained in the
    # parens
    my $judge = TestHarness::Judge->new();
    my $string = $item[1];
    # $string should look lik this
    # (!noh, regression.suite, simple)

    # Remove surrounding parens.
    $string =~ s/\s*[()]\s*//g;
 
    # Split comma delimited list.
    my @tests = split /\s*,\s*/;

    # process each test/testsuite
    foreach my $test(@tests) {
      my $not = "true";
      my $name = "";

      # parse the test to look for condition/testname
      ($not, $name) = ($test =~ /^\s* (\-|!|not|fail|pass)? \s* ([\w.]+) \s*$/ix);
      $judge->add_criterion($name, $not);
    }

    $return = $judge;
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION2__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [<skip: "">]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{paren_enclosed})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{paren_enclosed})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{paren_enclosed})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{paren_enclosed});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{paren_enclosed})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::if_else
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"if_else"};
	
	Parse::RecDescent::_trace(q{Trying rule: [if_else]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{if_else})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		local $skip = defined($skip) ? $skip : $Parse::RecDescent::skip;
		Parse::RecDescent::_trace(q{Trying production: [<skip: ""> if_block <commit> elsif_block else_block]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{if_else})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{if_else});
		%item = (__RULE__ => q{if_else});
		my $repcount = 0;


		

		Parse::RecDescent::_trace(q{Trying directive: [<skip: "">]},
					Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE; 
		$_tok = do { my $oldskip = $skip; $skip= ""; $oldskip };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying subrule: [if_block]},
				  Parse::RecDescent::_tracefirst($text),
				  q{if_else})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{if_block})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::if_block($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [if_block]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{if_else})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [if_block]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE;
		$item{q{if_block}} = $_tok;
		push @item, $_tok;
		
		}

		

		Parse::RecDescent::_trace(q{Trying directive: [<commit>]},
					Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE; 
		$_tok = do { $commit = 1 };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE2__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying repeated subrule: [elsif_block]},
				  Parse::RecDescent::_tracefirst($text),
				  q{if_else})
					if defined $::RD_TRACE;
		$expectation->is(q{elsif_block})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::elsif_block, 0, 100000000, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [elsif_block]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{if_else})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [elsif_block]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE;
		$item{q{elsif_block}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying repeated subrule: [else_block]},
				  Parse::RecDescent::_tracefirst($text),
				  q{if_else})
					if defined $::RD_TRACE;
		$expectation->is(q{else_block})->at($text);
		
		unless (defined ($_tok = $thisparser->_parserepeat($text, \&Parse::RecDescent::namespace000001::else_block, 0, 1, $_noactions,$expectation,undef))) 
		{
			Parse::RecDescent::_trace(q{<<Didn't match repeated subrule: [else_block]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{if_else})
							if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched repeated subrule: [else_block]<< (}
					. @$_tok . q{ times)},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE;
		$item{q{else_block}} = $_tok;
		push @item, $_tok;
		


		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    print "\nIF_ELSE: ", Dumper (%item), "\n\n" if $debug;
    # Foreach testproblem/suite from the elsif_blocks, assign the anti-judge
    # created in the if_block. 
    # 
    # Foreach testproblem/suite from the else_block, assign the anti-judge from
    # the if_block and from each else_block.
    #
    # Get a list of judges and files.
    #
    # If_block returns a judge.
    my @judges = ($item{if_block});

    # Elsif_block returns an array (judge, (lists of tests/testsuites)). We
    # must process each elsif_block. Each block has to worry about the judges
    # created in all previous elsif_blocks as well as the one created in the
    # if_block.
    foreach (@{$item{'elsif_block(s?)'}}) {
      my @tests = @{$_[1]};
      
      # apply the current set of anti-judges
      foreach my $test(@tests) {
	foreach my $judge(@judges) {
	  $test->add_anti_judge($judge);
	}
      }
      
      # Add the locally produced judge
      push @judges, $_[0];
    }

    # Else_block returns a list of tests/testsuites. There can only be one such
    # block but the parser will return it in an array if it exists.
    foreach (@{$item{'else_block(?)'}}) {
      # apply the current set of anti-judges
      foreach my $test(@$_) {
	foreach my $judge(@judges) {
	  $test->add_anti_judge($judge);
	}
      }
    }
    #
    $return = 1;
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [<skip: ""> if_block <commit> elsif_block else_block]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{if_else})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{if_else})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{if_else})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{if_else});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{if_else})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::_alternation_1_of_production_2_of_rule_suitefile
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"_alternation_1_of_production_2_of_rule_suitefile"};
	
	Parse::RecDescent::_trace(q{Trying rule: [_alternation_1_of_production_2_of_rule_suitefile]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{_alternation_1_of_production_2_of_rule_suitefile})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [if_else]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_2_of_rule_suitefile});
		%item = (__RULE__ => q{_alternation_1_of_production_2_of_rule_suitefile});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [if_else]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_2_of_rule_suitefile})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::if_else($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [if_else]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_2_of_rule_suitefile})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [if_else]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		$item{q{if_else}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [if_else]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [line]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[1];
		$text = $_[1];
		my $_savetext;
		@item = (q{_alternation_1_of_production_2_of_rule_suitefile});
		%item = (__RULE__ => q{_alternation_1_of_production_2_of_rule_suitefile});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying subrule: [line]},
				  Parse::RecDescent::_tracefirst($text),
				  q{_alternation_1_of_production_2_of_rule_suitefile})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::line($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [line]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{_alternation_1_of_production_2_of_rule_suitefile})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [line]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		$item{q{line}} = $_tok;
		push @item, $_tok;
		
		}


		Parse::RecDescent::_trace(q{>>Matched production: [line]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{_alternation_1_of_production_2_of_rule_suitefile})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{_alternation_1_of_production_2_of_rule_suitefile})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{_alternation_1_of_production_2_of_rule_suitefile});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{_alternation_1_of_production_2_of_rule_suitefile})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::BLANK_LINE
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"BLANK_LINE"};
	
	Parse::RecDescent::_trace(q{Trying rule: [BLANK_LINE]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{BLANK_LINE})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/[ \\t]*\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{BLANK_LINE})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{BLANK_LINE});
		%item = (__RULE__ => q{BLANK_LINE});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/[ \\t]*\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{BLANK_LINE})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[ \t]*\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/[ \\t]*\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{BLANK_LINE})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{BLANK_LINE})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{BLANK_LINE})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{BLANK_LINE});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{BLANK_LINE})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::elsif_block
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"elsif_block"};
	
	Parse::RecDescent::_trace(q{Trying rule: [elsif_block]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{elsif_block})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/else?\\s*if/i paren_enclosed bracket_enclosed]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{elsif_block})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{elsif_block});
		%item = (__RULE__ => q{elsif_block});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/else?\\s*if/i]}, Parse::RecDescent::_tracefirst($text),
					  q{elsif_block})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:else?\s*if)//i)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying subrule: [paren_enclosed]},
				  Parse::RecDescent::_tracefirst($text),
				  q{elsif_block})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{paren_enclosed})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::paren_enclosed($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [paren_enclosed]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{elsif_block})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [paren_enclosed]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{elsif_block})
						if defined $::RD_TRACE;
		$item{q{paren_enclosed}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying subrule: [bracket_enclosed]},
				  Parse::RecDescent::_tracefirst($text),
				  q{elsif_block})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{bracket_enclosed})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::bracket_enclosed($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [bracket_enclosed]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{elsif_block})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [bracket_enclosed]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{elsif_block})
						if defined $::RD_TRACE;
		$item{q{bracket_enclosed}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{elsif_block})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    print "\nELSIF_BLOCK: ", Dumper (%item), "\\nn" if $debug;
    # Foreach bracket_enclosed assign the Judge created in the
    # paren_enclosed.
    my $judge = $item{paren_enclosed};
    foreach my $test (@{$item{bracket_enclosed}}) {
      $test->add_judge($judge);
    }

    # I need to return the judge and the included tests.
    my $return = [$judge, $item{bracket_enclosed}];
 };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/else?\\s*if/i paren_enclosed bracket_enclosed]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{elsif_block})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{elsif_block})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{elsif_block})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{elsif_block});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{elsif_block})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::TESTNAME
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"TESTNAME"};
	
	Parse::RecDescent::_trace(q{Trying rule: [TESTNAME]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{TESTNAME})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/[-\\/\\w0-9_]+(\\.test)?/i]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{TESTNAME})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{TESTNAME});
		%item = (__RULE__ => q{TESTNAME});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/[-\\/\\w0-9_]+(\\.test)?/i]}, Parse::RecDescent::_tracefirst($text),
					  q{TESTNAME})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:[-\/\w0-9_]+(\.test)?)//i)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{TESTNAME})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $item[1];
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/[-\\/\\w0-9_]+(\\.test)?/i]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{TESTNAME})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{TESTNAME})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{TESTNAME})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{TESTNAME});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{TESTNAME})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::comment_line
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"comment_line"};
	
	Parse::RecDescent::_trace(q{Trying rule: [comment_line]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{comment_line})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		local $skip = defined($skip) ? $skip : $Parse::RecDescent::skip;
		Parse::RecDescent::_trace(q{Trying production: [<skip: ""> COMMENT /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{comment_line})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{comment_line});
		%item = (__RULE__ => q{comment_line});
		my $repcount = 0;


		

		Parse::RecDescent::_trace(q{Trying directive: [<skip: "">]},
					Parse::RecDescent::_tracefirst($text),
					  q{comment_line})
						if defined $::RD_TRACE; 
		$_tok = do { my $oldskip = $skip; $skip= ""; $oldskip };
		if (defined($_tok))
		{
			Parse::RecDescent::_trace(q{>>Matched directive<< (return value: [}
						. $_tok . q{])},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		else
		{
			Parse::RecDescent::_trace(q{<<Didn't match directive>>},
						Parse::RecDescent::_tracefirst($text))
							if defined $::RD_TRACE;
		}
		
		last unless defined $_tok;
		push @item, $item{__DIRECTIVE1__}=$_tok;
		

		Parse::RecDescent::_trace(q{Trying subrule: [COMMENT]},
				  Parse::RecDescent::_tracefirst($text),
				  q{comment_line})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{COMMENT})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::COMMENT($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [COMMENT]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{comment_line})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [COMMENT]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{comment_line})
						if defined $::RD_TRACE;
		$item{q{COMMENT}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{comment_line})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		


		Parse::RecDescent::_trace(q{>>Matched production: [<skip: ""> COMMENT /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{comment_line})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{comment_line})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{comment_line})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{comment_line});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{comment_line})
	}
	$_[1] = $text;
	return $return;
}

# ARGS ARE: ($parser, $text; $repeating, $_noactions, \@args)
sub Parse::RecDescent::namespace000001::exclude
{
	my $thisparser = $_[0];
	$ERRORS = 0;
	my $thisrule = $thisparser->{"rules"}{"exclude"};
	
	Parse::RecDescent::_trace(q{Trying rule: [exclude]},
				  Parse::RecDescent::_tracefirst($_[1]),
				  q{exclude})
					if defined $::RD_TRACE;

	
	my $err_at = @{$thisparser->{errors}};

	my $score;
	my $score_return;
	my $_tok;
	my $return = undef;
	my $_matched=0;
	my $commit=0;
	my @item = ();
	my %item = ();
	my $repeating =  defined($_[2]) && $_[2];
	my $_noactions = defined($_[3]) && $_[3];
 	my @arg =        defined $_[4] ? @{ &{$_[4]} } : ();
	my %arg =        ($#arg & 01) ? @arg : (@arg, undef);
	my $text;
	my $lastsep="";
	my $expectation = new Parse::RecDescent::Expectation($thisrule->expected());
	$expectation->at($_[1]);
	
	my $thisline;
	tie $thisline, q{Parse::RecDescent::LineCounter}, \$text, $thisparser;

	

	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/ \\- |  ! |  not |  exclude/ix SUITENAME /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{exclude})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[0];
		$text = $_[1];
		my $_savetext;
		@item = (q{exclude});
		%item = (__RULE__ => q{exclude});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/ \\- |  ! |  not |  exclude/ix]}, Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?: \- |  ! |  not |  exclude)//ix)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying subrule: [SUITENAME]},
				  Parse::RecDescent::_tracefirst($text),
				  q{exclude})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{SUITENAME})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::SUITENAME($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [SUITENAME]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{exclude})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [SUITENAME]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$item{q{SUITENAME}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN2__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $suite->exclude_suite($item{SUITENAME});
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/ \\- |  ! |  not |  exclude/ix SUITENAME /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


	while (!$_matched && !$commit)
	{
		
		Parse::RecDescent::_trace(q{Trying production: [/\\- | ! | not | exclude/ix TESTNAME /\\n/]},
					  Parse::RecDescent::_tracefirst($_[1]),
					  q{exclude})
						if defined $::RD_TRACE;
		my $thisprod = $thisrule->{"prods"}[1];
		$text = $_[1];
		my $_savetext;
		@item = (q{exclude});
		%item = (__RULE__ => q{exclude});
		my $repcount = 0;


		Parse::RecDescent::_trace(q{Trying terminal: [/\\- | ! | not | exclude/ix]}, Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\- | ! | not | exclude)//ix)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN1__}=$&;
		

		Parse::RecDescent::_trace(q{Trying subrule: [TESTNAME]},
				  Parse::RecDescent::_tracefirst($text),
				  q{exclude})
					if defined $::RD_TRACE;
		if (1) { no strict qw{refs};
		$expectation->is(q{TESTNAME})->at($text);
		unless (defined ($_tok = Parse::RecDescent::namespace000001::TESTNAME($thisparser,$text,$repeating,$_noactions,undef)))
		{
			
			Parse::RecDescent::_trace(q{<<Didn't match subrule: [TESTNAME]>>},
						  Parse::RecDescent::_tracefirst($text),
						  q{exclude})
							if defined $::RD_TRACE;
			$expectation->failed();
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched subrule: [TESTNAME]<< (return value: [}
					. $_tok . q{]},
					  
					  Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$item{q{TESTNAME}} = $_tok;
		push @item, $_tok;
		
		}

		Parse::RecDescent::_trace(q{Trying terminal: [/\\n/]}, Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$lastsep = "";
		$expectation->is(q{/\\n/})->at($text);
		

		unless ($text =~ s/\A($skip)/$lastsep=$1 and ""/e and   $text =~ s/\A(?:\n)//)
		{
			
			$expectation->failed();
			Parse::RecDescent::_trace(q{<<Didn't match terminal>>},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;

			last;
		}
		Parse::RecDescent::_trace(q{>>Matched terminal<< (return value: [}
						. $& . q{])},
						  Parse::RecDescent::_tracefirst($text))
					if defined $::RD_TRACE;
		push @item, $item{__PATTERN2__}=$&;
		

		Parse::RecDescent::_trace(q{Trying action},
					  Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		

		$_tok = ($_noactions) ? 0 : do {
    $return = $suite->exclude_test($item{TESTNAME});
  };
		unless (defined $_tok)
		{
			Parse::RecDescent::_trace(q{<<Didn't match action>> (return value: [undef])})
					if defined $::RD_TRACE;
			last;
		}
		Parse::RecDescent::_trace(q{>>Matched action<< (return value: [}
					  . $_tok . q{])}, $text)
						if defined $::RD_TRACE;
		push @item, $_tok;
		$item{__ACTION1__}=$_tok;
		


		Parse::RecDescent::_trace(q{>>Matched production: [/\\- | ! | not | exclude/ix TESTNAME /\\n/]<<},
					  Parse::RecDescent::_tracefirst($text),
					  q{exclude})
						if defined $::RD_TRACE;
		$_matched = 1;
		last;
	}


        unless ( $_matched || defined($return) || defined($score) )
	{
		

		$_[1] = $text;	# NOT SURE THIS IS NEEDED
		Parse::RecDescent::_trace(q{<<Didn't match rule>>},
					 Parse::RecDescent::_tracefirst($_[1]),
					 q{exclude})
					if defined $::RD_TRACE;
		return undef;
	}
	if (!defined($return) && defined($score))
	{
		Parse::RecDescent::_trace(q{>>Accepted scored production<<}, "",
					  q{exclude})
						if defined $::RD_TRACE;
		$return = $score_return;
	}
	splice @{$thisparser->{errors}}, $err_at;
	$return = $item[$#item] unless defined $return;
	if (defined $::RD_TRACE)
	{
		Parse::RecDescent::_trace(q{>>Matched rule<< (return value: [} .
					  $return . q{])}, "",
					  q{exclude});
		Parse::RecDescent::_trace(q{(consumed: [} .
					  Parse::RecDescent::_tracemax(substr($_[1],0,-length($text))) . q{])}, 
					  Parse::RecDescent::_tracefirst($text),
					  , q{exclude})
	}
	$_[1] = $text;
	return $return;
}
}
package TestHarness::Grammar; sub new { my $self = bless( {
                 '_AUTOTREE' => undef,
                 'localvars' => '',
                 'startcode' => '',
                 '_check' => {
                               'thisoffset' => '',
                               'itempos' => '',
                               'prevoffset' => '',
                               'prevline' => '',
                               'prevcolumn' => '',
                               'thiscolumn' => ''
                             },
                 'namespace' => 'Parse::RecDescent::namespace000001',
                 '_AUTOACTION' => undef,
                 'rules' => {
                              'set_param_line' => bless( {
                                                           'impcount' => 0,
                                                           'calls' => [],
                                                           'opcount' => 0,
                                                           'prods' => [
                                                                        bless( {
                                                                                 'number' => '0',
                                                                                 'strcount' => 0,
                                                                                 'dircount' => 0,
                                                                                 'uncommit' => undef,
                                                                                 'error' => undef,
                                                                                 'patcount' => 1,
                                                                                 'actcount' => 1,
                                                                                 'items' => [
                                                                                              bless( {
                                                                                                       'description' => '/[ \\\\t]* ([\\\\w_]+) [ \\\\t]* [=:]* [ \\\\t*]* ([-._\\\\w\\\\d]+) [ \\\\t]* \\\\n/x',
                                                                                                       'rdelim' => '/',
                                                                                                       'pattern' => '[ \\t]* ([\\w_]+) [ \\t]* [=:]* [ \\t*]* ([-._\\w\\d]+) [ \\t]* \\n',
                                                                                                       'hashname' => '__PATTERN1__',
                                                                                                       'lookahead' => 0,
                                                                                                       'ldelim' => '/',
                                                                                                       'mod' => 'x',
                                                                                                       'line' => 257
                                                                                                     }, 'Parse::RecDescent::Token' ),
                                                                                              bless( {
                                                                                                       'hashname' => '__ACTION1__',
                                                                                                       'lookahead' => 0,
                                                                                                       'line' => 258,
                                                                                                       'code' => '{
      my $func = lc $1;
      my $value = $2;
      $suite->$func($value);
   }'
                                                                                                     }, 'Parse::RecDescent::Action' )
                                                                                            ],
                                                                                 'line' => undef
                                                                               }, 'Parse::RecDescent::Production' )
                                                                      ],
                                                           'name' => 'set_param_line',
                                                           'vars' => '',
                                                           'changed' => 0,
                                                           'line' => 257
                                                         }, 'Parse::RecDescent::Rule' ),
                              '_alternation_1_of_production_2_of_rule_include' => bless( {
                                                                                           'impcount' => 0,
                                                                                           'calls' => [],
                                                                                           'opcount' => 0,
                                                                                           'prods' => [
                                                                                                        bless( {
                                                                                                                 'number' => '0',
                                                                                                                 'strcount' => 0,
                                                                                                                 'dircount' => 0,
                                                                                                                 'uncommit' => undef,
                                                                                                                 'error' => undef,
                                                                                                                 'patcount' => 1,
                                                                                                                 'actcount' => 0,
                                                                                                                 'items' => [
                                                                                                                              bless( {
                                                                                                                                       'description' => '/\\\\+ | include/ix',
                                                                                                                                       'rdelim' => '/',
                                                                                                                                       'pattern' => '\\+ | include',
                                                                                                                                       'hashname' => '__PATTERN1__',
                                                                                                                                       'lookahead' => 0,
                                                                                                                                       'ldelim' => '/',
                                                                                                                                       'mod' => 'ix',
                                                                                                                                       'line' => 263
                                                                                                                                     }, 'Parse::RecDescent::Token' )
                                                                                                                            ],
                                                                                                                 'line' => undef
                                                                                                               }, 'Parse::RecDescent::Production' )
                                                                                                      ],
                                                                                           'name' => '_alternation_1_of_production_2_of_rule_include',
                                                                                           'vars' => '',
                                                                                           'changed' => 0,
                                                                                           'line' => 263
                                                                                         }, 'Parse::RecDescent::Rule' ),
                              'include' => bless( {
                                                    'impcount' => 1,
                                                    'calls' => [
                                                                 '_alternation_1_of_production_1_of_rule_include',
                                                                 'SUITENAME',
                                                                 '_alternation_1_of_production_2_of_rule_include',
                                                                 'TESTNAME',
                                                                 'CPUS_RUNTIME'
                                                               ],
                                                    'opcount' => 0,
                                                    'prods' => [
                                                                 bless( {
                                                                          'number' => '0',
                                                                          'strcount' => 0,
                                                                          'dircount' => 0,
                                                                          'uncommit' => undef,
                                                                          'error' => undef,
                                                                          'patcount' => 1,
                                                                          'actcount' => 1,
                                                                          'items' => [
                                                                                       bless( {
                                                                                                'subrule' => '_alternation_1_of_production_1_of_rule_include',
                                                                                                'expected' => '/\\\\+ | include/ix',
                                                                                                'min' => 0,
                                                                                                'argcode' => undef,
                                                                                                'max' => 1,
                                                                                                'matchrule' => 0,
                                                                                                'repspec' => '?',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 184
                                                                                              }, 'Parse::RecDescent::Repetition' ),
                                                                                       bless( {
                                                                                                'subrule' => 'SUITENAME',
                                                                                                'matchrule' => 0,
                                                                                                'implicit' => undef,
                                                                                                'argcode' => undef,
                                                                                                'lookahead' => 0,
                                                                                                'line' => 184
                                                                                              }, 'Parse::RecDescent::Subrule' ),
                                                                                       bless( {
                                                                                                'description' => '/\\\\n/',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => '\\n',
                                                                                                'hashname' => '__PATTERN1__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => '',
                                                                                                'line' => 184
                                                                                              }, 'Parse::RecDescent::Token' ),
                                                                                       bless( {
                                                                                                'hashname' => '__ACTION1__',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 185,
                                                                                                'code' => '{
    # We are going to parse another suite file. $suite is a package variable
    # which will be changed to the new suite. We must use local to protect the
    # value of $suite.
    my $new_suite;
    {
      my $old_suite = $suite;
      local $suite;
      local $_;
      $new_suite = $old_suite->add_testsuite($item{SUITENAME});
    }
    $return = [@item];
  }'
                                                                                              }, 'Parse::RecDescent::Action' )
                                                                                     ],
                                                                          'line' => undef
                                                                        }, 'Parse::RecDescent::Production' ),
                                                                 bless( {
                                                                          'number' => '1',
                                                                          'strcount' => 0,
                                                                          'dircount' => 0,
                                                                          'uncommit' => undef,
                                                                          'error' => undef,
                                                                          'patcount' => 1,
                                                                          'actcount' => 1,
                                                                          'items' => [
                                                                                       bless( {
                                                                                                'subrule' => '_alternation_1_of_production_2_of_rule_include',
                                                                                                'expected' => '/\\\\+ | include/ix',
                                                                                                'min' => 0,
                                                                                                'argcode' => undef,
                                                                                                'max' => 1,
                                                                                                'matchrule' => 0,
                                                                                                'repspec' => '?',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 198
                                                                                              }, 'Parse::RecDescent::Repetition' ),
                                                                                       bless( {
                                                                                                'subrule' => 'TESTNAME',
                                                                                                'matchrule' => 0,
                                                                                                'implicit' => undef,
                                                                                                'argcode' => undef,
                                                                                                'lookahead' => 0,
                                                                                                'line' => 198
                                                                                              }, 'Parse::RecDescent::Subrule' ),
                                                                                       bless( {
                                                                                                'subrule' => 'CPUS_RUNTIME',
                                                                                                'expected' => undef,
                                                                                                'min' => 0,
                                                                                                'argcode' => undef,
                                                                                                'max' => 1,
                                                                                                'matchrule' => 0,
                                                                                                'repspec' => '?',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 198
                                                                                              }, 'Parse::RecDescent::Repetition' ),
                                                                                       bless( {
                                                                                                'description' => '/\\\\n/',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => '\\n',
                                                                                                'hashname' => '__PATTERN1__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => '',
                                                                                                'line' => 198
                                                                                              }, 'Parse::RecDescent::Token' ),
                                                                                       bless( {
                                                                                                'hashname' => '__ACTION1__',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 199,
                                                                                                'code' => '{
    my $test = $suite->add_test($item{TESTNAME});
    #print "CPUS_RUNTIME: ", Dumper ($item{"CPUS_RUNTIME(?)"});
    my $cpus_runtime = $item{"CPUS_RUNTIME(?)"}->[0];
    while (my ($func, $value) = each %$cpus_runtime) {
      #print "func:$func:  value:$value\\n";
      $func = lc $func; # functions are always lower case
      $test->$func($value);
    }

    $return = [@item];
  }'
                                                                                              }, 'Parse::RecDescent::Action' )
                                                                                     ],
                                                                          'line' => 198
                                                                        }, 'Parse::RecDescent::Production' )
                                                               ],
                                                    'name' => 'include',
                                                    'vars' => '',
                                                    'changed' => 0,
                                                    'line' => 183
                                                  }, 'Parse::RecDescent::Rule' ),
                              'SUITENAME' => bless( {
                                                      'impcount' => 0,
                                                      'calls' => [],
                                                      'opcount' => 0,
                                                      'prods' => [
                                                                   bless( {
                                                                            'number' => '0',
                                                                            'strcount' => 0,
                                                                            'dircount' => 0,
                                                                            'uncommit' => undef,
                                                                            'error' => undef,
                                                                            'patcount' => 1,
                                                                            'actcount' => 1,
                                                                            'items' => [
                                                                                         bless( {
                                                                                                  'description' => '/[-\\\\/\\\\w0-9_]+\\\\.suite/i',
                                                                                                  'rdelim' => '/',
                                                                                                  'pattern' => '[-\\/\\w0-9_]+\\.suite',
                                                                                                  'hashname' => '__PATTERN1__',
                                                                                                  'lookahead' => 0,
                                                                                                  'ldelim' => '/',
                                                                                                  'mod' => 'i',
                                                                                                  'line' => 162
                                                                                                }, 'Parse::RecDescent::Token' ),
                                                                                         bless( {
                                                                                                  'hashname' => '__ACTION1__',
                                                                                                  'lookahead' => 0,
                                                                                                  'line' => 163,
                                                                                                  'code' => '{
    $return = $item[1];
  }'
                                                                                                }, 'Parse::RecDescent::Action' )
                                                                                       ],
                                                                            'line' => undef
                                                                          }, 'Parse::RecDescent::Production' )
                                                                 ],
                                                      'name' => 'SUITENAME',
                                                      'vars' => '',
                                                      'changed' => 0,
                                                      'line' => 159
                                                    }, 'Parse::RecDescent::Rule' ),
                              'COMMENT' => bless( {
                                                    'impcount' => 0,
                                                    'calls' => [],
                                                    'opcount' => 0,
                                                    'prods' => [
                                                                 bless( {
                                                                          'number' => '0',
                                                                          'strcount' => 0,
                                                                          'dircount' => 0,
                                                                          'uncommit' => undef,
                                                                          'error' => undef,
                                                                          'patcount' => 1,
                                                                          'actcount' => 0,
                                                                          'items' => [
                                                                                       bless( {
                                                                                                'description' => '/[ \\\\t]* ([\\\\#;\\\\$] | \\\\/\\\\/) [^\\\\n]*/x',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => '[ \\t]* ([\\#;\\$] | \\/\\/) [^\\n]*',
                                                                                                'hashname' => '__PATTERN1__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => 'x',
                                                                                                'line' => 253
                                                                                              }, 'Parse::RecDescent::Token' )
                                                                                     ],
                                                                          'line' => undef
                                                                        }, 'Parse::RecDescent::Production' )
                                                               ],
                                                    'name' => 'COMMENT',
                                                    'vars' => '',
                                                    'changed' => 0,
                                                    'line' => 252
                                                  }, 'Parse::RecDescent::Rule' ),
                              'suitefile' => bless( {
                                                      'impcount' => 1,
                                                      'calls' => [
                                                                   '_alternation_1_of_production_2_of_rule_suitefile'
                                                                 ],
                                                      'opcount' => 0,
                                                      'prods' => [
                                                                   bless( {
                                                                            'number' => '0',
                                                                            'strcount' => 0,
                                                                            'dircount' => 1,
                                                                            'uncommit' => undef,
                                                                            'error' => undef,
                                                                            'patcount' => 0,
                                                                            'actcount' => 1,
                                                                            'items' => [
                                                                                         bless( {
                                                                                                  'hashname' => '__ACTION1__',
                                                                                                  'lookahead' => 0,
                                                                                                  'line' => 30,
                                                                                                  'code' => '{$suite = $arg[0]}'
                                                                                                }, 'Parse::RecDescent::Action' ),
                                                                                         bless( {
                                                                                                  'hashname' => '__DIRECTIVE1__',
                                                                                                  'name' => '<reject>',
                                                                                                  'lookahead' => 0,
                                                                                                  'line' => 30
                                                                                                }, 'Parse::RecDescent::UncondReject' )
                                                                                       ],
                                                                            'line' => undef
                                                                          }, 'Parse::RecDescent::Production' ),
                                                                   bless( {
                                                                            'number' => '1',
                                                                            'strcount' => 0,
                                                                            'dircount' => 1,
                                                                            'uncommit' => undef,
                                                                            'error' => undef,
                                                                            'patcount' => 1,
                                                                            'actcount' => 1,
                                                                            'items' => [
                                                                                         bless( {
                                                                                                  'hashname' => '__DIRECTIVE1__',
                                                                                                  'name' => '<skip: "">',
                                                                                                  'lookahead' => 0,
                                                                                                  'line' => 31,
                                                                                                  'code' => 'my $oldskip = $skip; $skip= ""; $oldskip'
                                                                                                }, 'Parse::RecDescent::Directive' ),
                                                                                         bless( {
                                                                                                  'subrule' => '_alternation_1_of_production_2_of_rule_suitefile',
                                                                                                  'expected' => 'if_else, or line',
                                                                                                  'min' => 1,
                                                                                                  'argcode' => undef,
                                                                                                  'max' => 100000000,
                                                                                                  'matchrule' => 0,
                                                                                                  'repspec' => 's',
                                                                                                  'lookahead' => 0,
                                                                                                  'line' => 31
                                                                                                }, 'Parse::RecDescent::Repetition' ),
                                                                                         bless( {
                                                                                                  'description' => '/\\\\Z/',
                                                                                                  'rdelim' => '/',
                                                                                                  'pattern' => '\\Z',
                                                                                                  'hashname' => '__PATTERN1__',
                                                                                                  'lookahead' => 0,
                                                                                                  'ldelim' => '/',
                                                                                                  'mod' => '',
                                                                                                  'line' => 31
                                                                                                }, 'Parse::RecDescent::Token' ),
                                                                                         bless( {
                                                                                                  'hashname' => '__ACTION1__',
                                                                                                  'lookahead' => 0,
                                                                                                  'line' => 32,
                                                                                                  'code' => '{
    $return = $suite;
  }'
                                                                                                }, 'Parse::RecDescent::Action' )
                                                                                       ],
                                                                            'line' => 31
                                                                          }, 'Parse::RecDescent::Production' )
                                                                 ],
                                                      'name' => 'suitefile',
                                                      'vars' => '',
                                                      'changed' => 0,
                                                      'line' => 29
                                                    }, 'Parse::RecDescent::Rule' ),
                              'STRING' => bless( {
                                                   'impcount' => 0,
                                                   'calls' => [],
                                                   'opcount' => 0,
                                                   'prods' => [
                                                                bless( {
                                                                         'number' => '0',
                                                                         'strcount' => 0,
                                                                         'dircount' => 0,
                                                                         'uncommit' => undef,
                                                                         'error' => undef,
                                                                         'patcount' => 1,
                                                                         'actcount' => 1,
                                                                         'items' => [
                                                                                      bless( {
                                                                                               'description' => '/["\'\'"]? ([^"\'\'"\\\\n]*) ["\'\'"]?/x',
                                                                                               'rdelim' => '/',
                                                                                               'pattern' => '["\'\'"]? ([^"\'\'"\\n]*) ["\'\'"]?',
                                                                                               'hashname' => '__PATTERN1__',
                                                                                               'lookahead' => 0,
                                                                                               'ldelim' => '/',
                                                                                               'mod' => 'x',
                                                                                               'line' => 245
                                                                                             }, 'Parse::RecDescent::Token' ),
                                                                                      bless( {
                                                                                               'hashname' => '__ACTION1__',
                                                                                               'lookahead' => 0,
                                                                                               'line' => 246,
                                                                                               'code' => '{
    $return = $1;
  }'
                                                                                             }, 'Parse::RecDescent::Action' )
                                                                                    ],
                                                                         'line' => undef
                                                                       }, 'Parse::RecDescent::Production' )
                                                              ],
                                                   'name' => 'STRING',
                                                   'vars' => '',
                                                   'changed' => 0,
                                                   'line' => 243
                                                 }, 'Parse::RecDescent::Rule' ),
                              'NUM_CPUS' => bless( {
                                                     'impcount' => 0,
                                                     'calls' => [],
                                                     'opcount' => 0,
                                                     'prods' => [
                                                                  bless( {
                                                                           'number' => '0',
                                                                           'strcount' => 0,
                                                                           'dircount' => 0,
                                                                           'uncommit' => undef,
                                                                           'error' => undef,
                                                                           'patcount' => 2,
                                                                           'actcount' => 1,
                                                                           'items' => [
                                                                                        bless( {
                                                                                                 'description' => '/[0-9]+/',
                                                                                                 'rdelim' => '/',
                                                                                                 'pattern' => '[0-9]+',
                                                                                                 'hashname' => '__PATTERN1__',
                                                                                                 'lookahead' => 0,
                                                                                                 'ldelim' => '/',
                                                                                                 'mod' => '',
                                                                                                 'line' => 225
                                                                                               }, 'Parse::RecDescent::Token' ),
                                                                                        bless( {
                                                                                                 'description' => '/:/',
                                                                                                 'rdelim' => '/',
                                                                                                 'pattern' => ':',
                                                                                                 'hashname' => '__PATTERN2__',
                                                                                                 'lookahead' => -1,
                                                                                                 'ldelim' => '/',
                                                                                                 'mod' => '',
                                                                                                 'line' => 225
                                                                                               }, 'Parse::RecDescent::Token' ),
                                                                                        bless( {
                                                                                                 'hashname' => '__ACTION1__',
                                                                                                 'lookahead' => 0,
                                                                                                 'line' => 226,
                                                                                                 'code' => '{
   #print "NUM_CPUS: ", Dumper(@item);
   $return = {$item[0] => $item[1]};
  }'
                                                                                               }, 'Parse::RecDescent::Action' )
                                                                                      ],
                                                                           'line' => undef
                                                                         }, 'Parse::RecDescent::Production' )
                                                                ],
                                                     'name' => 'NUM_CPUS',
                                                     'vars' => '',
                                                     'changed' => 0,
                                                     'line' => 225
                                                   }, 'Parse::RecDescent::Rule' ),
                              'if_block' => bless( {
                                                     'impcount' => 0,
                                                     'calls' => [
                                                                  'paren_enclosed',
                                                                  'bracket_enclosed'
                                                                ],
                                                     'opcount' => 0,
                                                     'prods' => [
                                                                  bless( {
                                                                           'number' => '0',
                                                                           'strcount' => 0,
                                                                           'dircount' => 1,
                                                                           'uncommit' => undef,
                                                                           'error' => undef,
                                                                           'patcount' => 1,
                                                                           'actcount' => 1,
                                                                           'items' => [
                                                                                        bless( {
                                                                                                 'description' => '/if/i',
                                                                                                 'rdelim' => '/',
                                                                                                 'pattern' => 'if',
                                                                                                 'hashname' => '__PATTERN1__',
                                                                                                 'lookahead' => 0,
                                                                                                 'ldelim' => '/',
                                                                                                 'mod' => 'i',
                                                                                                 'line' => 83
                                                                                               }, 'Parse::RecDescent::Token' ),
                                                                                        bless( {
                                                                                                 'hashname' => '__DIRECTIVE1__',
                                                                                                 'name' => '<commit>',
                                                                                                 'lookahead' => 0,
                                                                                                 'line' => 83,
                                                                                                 'code' => '$commit = 1'
                                                                                               }, 'Parse::RecDescent::Directive' ),
                                                                                        bless( {
                                                                                                 'subrule' => 'paren_enclosed',
                                                                                                 'matchrule' => 0,
                                                                                                 'implicit' => undef,
                                                                                                 'argcode' => undef,
                                                                                                 'lookahead' => 0,
                                                                                                 'line' => 83
                                                                                               }, 'Parse::RecDescent::Subrule' ),
                                                                                        bless( {
                                                                                                 'subrule' => 'bracket_enclosed',
                                                                                                 'matchrule' => 0,
                                                                                                 'implicit' => undef,
                                                                                                 'argcode' => undef,
                                                                                                 'lookahead' => 0,
                                                                                                 'line' => 83
                                                                                               }, 'Parse::RecDescent::Subrule' ),
                                                                                        bless( {
                                                                                                 'hashname' => '__ACTION1__',
                                                                                                 'lookahead' => 0,
                                                                                                 'line' => 84,
                                                                                                 'code' => '{
    print "\\nIF_BLOCK: ", Dumper (%item), "\\n\\n" if $debug;
    # Foreach bracket_enclosed assign the Judge created in the paren_enclosed.
    my $judge = $item{paren_enclosed};
    foreach my $test (@{$item{bracket_enclosed}}) {
      $test->add_judge($judge);
    }

    # These tests are all taken care of so I don\'t need to return them, but
    # I need to return the judge.
    $return = $judge;
  }'
                                                                                               }, 'Parse::RecDescent::Action' )
                                                                                      ],
                                                                           'line' => undef
                                                                         }, 'Parse::RecDescent::Production' )
                                                                ],
                                                     'name' => 'if_block',
                                                     'vars' => '',
                                                     'changed' => 0,
                                                     'line' => 83
                                                   }, 'Parse::RecDescent::Rule' ),
                              'CPUS_RUNTIME' => bless( {
                                                         'impcount' => 0,
                                                         'calls' => [
                                                                      'NUM_CPUS',
                                                                      'TIME_LIMIT'
                                                                    ],
                                                         'opcount' => 0,
                                                         'prods' => [
                                                                      bless( {
                                                                               'number' => '0',
                                                                               'strcount' => 0,
                                                                               'dircount' => 0,
                                                                               'uncommit' => undef,
                                                                               'error' => undef,
                                                                               'patcount' => 0,
                                                                               'actcount' => 1,
                                                                               'items' => [
                                                                                            bless( {
                                                                                                     'subrule' => 'NUM_CPUS',
                                                                                                     'matchrule' => 0,
                                                                                                     'implicit' => undef,
                                                                                                     'argcode' => undef,
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 213
                                                                                                   }, 'Parse::RecDescent::Subrule' ),
                                                                                            bless( {
                                                                                                     'subrule' => 'TIME_LIMIT',
                                                                                                     'expected' => undef,
                                                                                                     'min' => 0,
                                                                                                     'argcode' => undef,
                                                                                                     'max' => 1,
                                                                                                     'matchrule' => 0,
                                                                                                     'repspec' => '?',
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 213
                                                                                                   }, 'Parse::RecDescent::Repetition' ),
                                                                                            bless( {
                                                                                                     'hashname' => '__ACTION1__',
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 214,
                                                                                                     'code' => '{
       my $time_limit = $item{"TIME_LIMIT(?)"}->[0] || {};
       $return = {%{$item{NUM_CPUS}}, %$time_limit};
     }'
                                                                                                   }, 'Parse::RecDescent::Action' )
                                                                                          ],
                                                                               'line' => undef
                                                                             }, 'Parse::RecDescent::Production' ),
                                                                      bless( {
                                                                               'number' => '1',
                                                                               'strcount' => 0,
                                                                               'dircount' => 0,
                                                                               'uncommit' => undef,
                                                                               'error' => undef,
                                                                               'patcount' => 0,
                                                                               'actcount' => 1,
                                                                               'items' => [
                                                                                            bless( {
                                                                                                     'subrule' => 'TIME_LIMIT',
                                                                                                     'matchrule' => 0,
                                                                                                     'implicit' => undef,
                                                                                                     'argcode' => undef,
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 218
                                                                                                   }, 'Parse::RecDescent::Subrule' ),
                                                                                            bless( {
                                                                                                     'subrule' => 'NUM_CPUS',
                                                                                                     'expected' => undef,
                                                                                                     'min' => 0,
                                                                                                     'argcode' => undef,
                                                                                                     'max' => 1,
                                                                                                     'matchrule' => 0,
                                                                                                     'repspec' => '?',
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 218
                                                                                                   }, 'Parse::RecDescent::Repetition' ),
                                                                                            bless( {
                                                                                                     'hashname' => '__ACTION1__',
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 219,
                                                                                                     'code' => '{
       my $cpus =  $item{"NUM_CPUS(?)"}->[0] || {};
       $return = {%$cpus, %{$item{TIME_LIMIT}}};
     }'
                                                                                                   }, 'Parse::RecDescent::Action' )
                                                                                          ],
                                                                               'line' => 218
                                                                             }, 'Parse::RecDescent::Production' )
                                                                    ],
                                                         'name' => 'CPUS_RUNTIME',
                                                         'vars' => '',
                                                         'changed' => 0,
                                                         'line' => 211
                                                       }, 'Parse::RecDescent::Rule' ),
                              '_alternation_1_of_production_1_of_rule_include' => bless( {
                                                                                           'impcount' => 0,
                                                                                           'calls' => [],
                                                                                           'opcount' => 0,
                                                                                           'prods' => [
                                                                                                        bless( {
                                                                                                                 'number' => '0',
                                                                                                                 'strcount' => 0,
                                                                                                                 'dircount' => 0,
                                                                                                                 'uncommit' => undef,
                                                                                                                 'error' => undef,
                                                                                                                 'patcount' => 1,
                                                                                                                 'actcount' => 0,
                                                                                                                 'items' => [
                                                                                                                              bless( {
                                                                                                                                       'description' => '/\\\\+ | include/ix',
                                                                                                                                       'rdelim' => '/',
                                                                                                                                       'pattern' => '\\+ | include',
                                                                                                                                       'hashname' => '__PATTERN1__',
                                                                                                                                       'lookahead' => 0,
                                                                                                                                       'ldelim' => '/',
                                                                                                                                       'mod' => 'ix',
                                                                                                                                       'line' => 263
                                                                                                                                     }, 'Parse::RecDescent::Token' )
                                                                                                                            ],
                                                                                                                 'line' => undef
                                                                                                               }, 'Parse::RecDescent::Production' )
                                                                                                      ],
                                                                                           'name' => '_alternation_1_of_production_1_of_rule_include',
                                                                                           'vars' => '',
                                                                                           'changed' => 0,
                                                                                           'line' => 263
                                                                                         }, 'Parse::RecDescent::Rule' ),
                              'line' => bless( {
                                                 'impcount' => 1,
                                                 'calls' => [
                                                              '_alternation_1_of_production_1_of_rule_line'
                                                            ],
                                                 'opcount' => 0,
                                                 'prods' => [
                                                              bless( {
                                                                       'number' => '0',
                                                                       'strcount' => 0,
                                                                       'dircount' => 1,
                                                                       'uncommit' => undef,
                                                                       'error' => undef,
                                                                       'patcount' => 0,
                                                                       'actcount' => 1,
                                                                       'items' => [
                                                                                    bless( {
                                                                                             'hashname' => '__DIRECTIVE1__',
                                                                                             'name' => '<skip: "[ \\t]*">',
                                                                                             'lookahead' => 0,
                                                                                             'line' => 156,
                                                                                             'code' => 'my $oldskip = $skip; $skip= "[ \\t]*"; $oldskip'
                                                                                           }, 'Parse::RecDescent::Directive' ),
                                                                                    bless( {
                                                                                             'subrule' => '_alternation_1_of_production_1_of_rule_line',
                                                                                             'matchrule' => 0,
                                                                                             'implicit' => 'print_, or include, or exclude, or BLANK_LINE, or comment_line, or set_param_line',
                                                                                             'argcode' => undef,
                                                                                             'lookahead' => 0,
                                                                                             'line' => 157
                                                                                           }, 'Parse::RecDescent::Subrule' ),
                                                                                    bless( {
                                                                                             'hashname' => '__ACTION1__',
                                                                                             'lookahead' => 0,
                                                                                             'line' => 158,
                                                                                             'code' => '{1}'
                                                                                           }, 'Parse::RecDescent::Action' )
                                                                                  ],
                                                                       'line' => undef
                                                                     }, 'Parse::RecDescent::Production' )
                                                            ],
                                                 'name' => 'line',
                                                 'vars' => '',
                                                 'changed' => 0,
                                                 'line' => 154
                                               }, 'Parse::RecDescent::Rule' ),
                              'else_block' => bless( {
                                                       'impcount' => 0,
                                                       'calls' => [
                                                                    'bracket_enclosed'
                                                                  ],
                                                       'opcount' => 0,
                                                       'prods' => [
                                                                    bless( {
                                                                             'number' => '0',
                                                                             'strcount' => 0,
                                                                             'dircount' => 0,
                                                                             'uncommit' => undef,
                                                                             'error' => undef,
                                                                             'patcount' => 1,
                                                                             'actcount' => 1,
                                                                             'items' => [
                                                                                          bless( {
                                                                                                   'description' => '/else/i',
                                                                                                   'rdelim' => '/',
                                                                                                   'pattern' => 'else',
                                                                                                   'hashname' => '__PATTERN1__',
                                                                                                   'lookahead' => 0,
                                                                                                   'ldelim' => '/',
                                                                                                   'mod' => 'i',
                                                                                                   'line' => 111
                                                                                                 }, 'Parse::RecDescent::Token' ),
                                                                                          bless( {
                                                                                                   'subrule' => 'bracket_enclosed',
                                                                                                   'matchrule' => 0,
                                                                                                   'implicit' => undef,
                                                                                                   'argcode' => undef,
                                                                                                   'lookahead' => 0,
                                                                                                   'line' => 111
                                                                                                 }, 'Parse::RecDescent::Subrule' ),
                                                                                          bless( {
                                                                                                   'hashname' => '__ACTION1__',
                                                                                                   'lookahead' => 0,
                                                                                                   'line' => 112,
                                                                                                   'code' => '{
    print "\\nELSE_BLOCK: ", Dumper (%item), "\\n\\n" if $debug;
    # return the tests enclosed in brackets
    $return = $item{bracket_enclosed};
  }'
                                                                                                 }, 'Parse::RecDescent::Action' )
                                                                                        ],
                                                                             'line' => undef
                                                                           }, 'Parse::RecDescent::Production' )
                                                                  ],
                                                       'name' => 'else_block',
                                                       'vars' => '',
                                                       'changed' => 0,
                                                       'line' => 111
                                                     }, 'Parse::RecDescent::Rule' ),
                              '_alternation_1_of_production_1_of_rule_line' => bless( {
                                                                                        'impcount' => 0,
                                                                                        'calls' => [
                                                                                                     'print_',
                                                                                                     'include',
                                                                                                     'exclude',
                                                                                                     'BLANK_LINE',
                                                                                                     'comment_line',
                                                                                                     'set_param_line'
                                                                                                   ],
                                                                                        'opcount' => 0,
                                                                                        'prods' => [
                                                                                                     bless( {
                                                                                                              'number' => '0',
                                                                                                              'strcount' => 0,
                                                                                                              'dircount' => 0,
                                                                                                              'uncommit' => undef,
                                                                                                              'error' => undef,
                                                                                                              'patcount' => 0,
                                                                                                              'actcount' => 0,
                                                                                                              'items' => [
                                                                                                                           bless( {
                                                                                                                                    'subrule' => 'print_',
                                                                                                                                    'matchrule' => 0,
                                                                                                                                    'implicit' => undef,
                                                                                                                                    'argcode' => undef,
                                                                                                                                    'lookahead' => 0,
                                                                                                                                    'line' => 262
                                                                                                                                  }, 'Parse::RecDescent::Subrule' )
                                                                                                                         ],
                                                                                                              'line' => undef
                                                                                                            }, 'Parse::RecDescent::Production' ),
                                                                                                     bless( {
                                                                                                              'number' => '1',
                                                                                                              'strcount' => 0,
                                                                                                              'dircount' => 0,
                                                                                                              'uncommit' => undef,
                                                                                                              'error' => undef,
                                                                                                              'patcount' => 0,
                                                                                                              'actcount' => 0,
                                                                                                              'items' => [
                                                                                                                           bless( {
                                                                                                                                    'subrule' => 'include',
                                                                                                                                    'matchrule' => 0,
                                                                                                                                    'implicit' => undef,
                                                                                                                                    'argcode' => undef,
                                                                                                                                    'lookahead' => 0,
                                                                                                                                    'line' => 262
                                                                                                                                  }, 'Parse::RecDescent::Subrule' )
                                                                                                                         ],
                                                                                                              'line' => 262
                                                                                                            }, 'Parse::RecDescent::Production' ),
                                                                                                     bless( {
                                                                                                              'number' => '2',
                                                                                                              'strcount' => 0,
                                                                                                              'dircount' => 0,
                                                                                                              'uncommit' => undef,
                                                                                                              'error' => undef,
                                                                                                              'patcount' => 0,
                                                                                                              'actcount' => 0,
                                                                                                              'items' => [
                                                                                                                           bless( {
                                                                                                                                    'subrule' => 'exclude',
                                                                                                                                    'matchrule' => 0,
                                                                                                                                    'implicit' => undef,
                                                                                                                                    'argcode' => undef,
                                                                                                                                    'lookahead' => 0,
                                                                                                                                    'line' => 262
                                                                                                                                  }, 'Parse::RecDescent::Subrule' )
                                                                                                                         ],
                                                                                                              'line' => 262
                                                                                                            }, 'Parse::RecDescent::Production' ),
                                                                                                     bless( {
                                                                                                              'number' => '3',
                                                                                                              'strcount' => 0,
                                                                                                              'dircount' => 0,
                                                                                                              'uncommit' => undef,
                                                                                                              'error' => undef,
                                                                                                              'patcount' => 0,
                                                                                                              'actcount' => 0,
                                                                                                              'items' => [
                                                                                                                           bless( {
                                                                                                                                    'subrule' => 'BLANK_LINE',
                                                                                                                                    'matchrule' => 0,
                                                                                                                                    'implicit' => undef,
                                                                                                                                    'argcode' => undef,
                                                                                                                                    'lookahead' => 0,
                                                                                                                                    'line' => 262
                                                                                                                                  }, 'Parse::RecDescent::Subrule' )
                                                                                                                         ],
                                                                                                              'line' => 262
                                                                                                            }, 'Parse::RecDescent::Production' ),
                                                                                                     bless( {
                                                                                                              'number' => '4',
                                                                                                              'strcount' => 0,
                                                                                                              'dircount' => 0,
                                                                                                              'uncommit' => undef,
                                                                                                              'error' => undef,
                                                                                                              'patcount' => 0,
                                                                                                              'actcount' => 0,
                                                                                                              'items' => [
                                                                                                                           bless( {
                                                                                                                                    'subrule' => 'comment_line',
                                                                                                                                    'matchrule' => 0,
                                                                                                                                    'implicit' => undef,
                                                                                                                                    'argcode' => undef,
                                                                                                                                    'lookahead' => 0,
                                                                                                                                    'line' => 263
                                                                                                                                  }, 'Parse::RecDescent::Subrule' )
                                                                                                                         ],
                                                                                                              'line' => 262
                                                                                                            }, 'Parse::RecDescent::Production' ),
                                                                                                     bless( {
                                                                                                              'number' => '5',
                                                                                                              'strcount' => 0,
                                                                                                              'dircount' => 0,
                                                                                                              'uncommit' => undef,
                                                                                                              'error' => undef,
                                                                                                              'patcount' => 0,
                                                                                                              'actcount' => 0,
                                                                                                              'items' => [
                                                                                                                           bless( {
                                                                                                                                    'subrule' => 'set_param_line',
                                                                                                                                    'matchrule' => 0,
                                                                                                                                    'implicit' => undef,
                                                                                                                                    'argcode' => undef,
                                                                                                                                    'lookahead' => 0,
                                                                                                                                    'line' => 263
                                                                                                                                  }, 'Parse::RecDescent::Subrule' )
                                                                                                                         ],
                                                                                                              'line' => 263
                                                                                                            }, 'Parse::RecDescent::Production' )
                                                                                                   ],
                                                                                        'name' => '_alternation_1_of_production_1_of_rule_line',
                                                                                        'vars' => '',
                                                                                        'changed' => 0,
                                                                                        'line' => 262
                                                                                      }, 'Parse::RecDescent::Rule' ),
                              'bracket_enclosed' => bless( {
                                                             'impcount' => 0,
                                                             'calls' => [
                                                                          'COMMENT',
                                                                          'line'
                                                                        ],
                                                             'opcount' => 0,
                                                             'prods' => [
                                                                          bless( {
                                                                                   'number' => '0',
                                                                                   'strcount' => 0,
                                                                                   'dircount' => 1,
                                                                                   'uncommit' => undef,
                                                                                   'error' => undef,
                                                                                   'patcount' => 4,
                                                                                   'actcount' => 1,
                                                                                   'items' => [
                                                                                                bless( {
                                                                                                         'hashname' => '__DIRECTIVE1__',
                                                                                                         'name' => '<skip: "">',
                                                                                                         'lookahead' => 0,
                                                                                                         'line' => 149,
                                                                                                         'code' => 'my $oldskip = $skip; $skip= ""; $oldskip'
                                                                                                       }, 'Parse::RecDescent::Directive' ),
                                                                                                bless( {
                                                                                                         'description' => '/\\{/',
                                                                                                         'rdelim' => '/',
                                                                                                         'pattern' => '{',
                                                                                                         'hashname' => '__PATTERN1__',
                                                                                                         'lookahead' => 0,
                                                                                                         'ldelim' => '/',
                                                                                                         'mod' => '',
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Token' ),
                                                                                                bless( {
                                                                                                         'subrule' => 'COMMENT',
                                                                                                         'expected' => undef,
                                                                                                         'min' => 0,
                                                                                                         'argcode' => undef,
                                                                                                         'max' => 1,
                                                                                                         'matchrule' => 0,
                                                                                                         'repspec' => '?',
                                                                                                         'lookahead' => 0,
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Repetition' ),
                                                                                                bless( {
                                                                                                         'description' => '/[ \\\\t]*\\\\n/',
                                                                                                         'rdelim' => '/',
                                                                                                         'pattern' => '[ \\t]*\\n',
                                                                                                         'hashname' => '__PATTERN2__',
                                                                                                         'lookahead' => 0,
                                                                                                         'ldelim' => '/',
                                                                                                         'mod' => '',
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Token' ),
                                                                                                bless( {
                                                                                                         'subrule' => 'line',
                                                                                                         'expected' => undef,
                                                                                                         'min' => 0,
                                                                                                         'argcode' => undef,
                                                                                                         'max' => 100000000,
                                                                                                         'matchrule' => 0,
                                                                                                         'repspec' => 's?',
                                                                                                         'lookahead' => 0,
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Repetition' ),
                                                                                                bless( {
                                                                                                         'description' => '/\\}/',
                                                                                                         'rdelim' => '/',
                                                                                                         'pattern' => '}',
                                                                                                         'hashname' => '__PATTERN3__',
                                                                                                         'lookahead' => 0,
                                                                                                         'ldelim' => '/',
                                                                                                         'mod' => '',
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Token' ),
                                                                                                bless( {
                                                                                                         'subrule' => 'COMMENT',
                                                                                                         'expected' => undef,
                                                                                                         'min' => 0,
                                                                                                         'argcode' => undef,
                                                                                                         'max' => 1,
                                                                                                         'matchrule' => 0,
                                                                                                         'repspec' => '?',
                                                                                                         'lookahead' => 0,
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Repetition' ),
                                                                                                bless( {
                                                                                                         'description' => '/\\\\n/',
                                                                                                         'rdelim' => '/',
                                                                                                         'pattern' => '\\n',
                                                                                                         'hashname' => '__PATTERN4__',
                                                                                                         'lookahead' => 0,
                                                                                                         'ldelim' => '/',
                                                                                                         'mod' => '',
                                                                                                         'line' => 149
                                                                                                       }, 'Parse::RecDescent::Token' ),
                                                                                                bless( {
                                                                                                         'hashname' => '__ACTION1__',
                                                                                                         'lookahead' => 0,
                                                                                                         'line' => 150,
                                                                                                         'code' => '{
    # Just return the list of problems/suites. These are objects.
    $return = $item{\'line(s?)\'};
  }'
                                                                                                       }, 'Parse::RecDescent::Action' )
                                                                                              ],
                                                                                   'line' => undef
                                                                                 }, 'Parse::RecDescent::Production' )
                                                                        ],
                                                             'name' => 'bracket_enclosed',
                                                             'vars' => '',
                                                             'changed' => 0,
                                                             'line' => 148
                                                           }, 'Parse::RecDescent::Rule' ),
                              'TIME_LIMIT' => bless( {
                                                       'impcount' => 0,
                                                       'calls' => [],
                                                       'opcount' => 0,
                                                       'prods' => [
                                                                    bless( {
                                                                             'number' => '0',
                                                                             'strcount' => 0,
                                                                             'dircount' => 0,
                                                                             'uncommit' => undef,
                                                                             'error' => undef,
                                                                             'patcount' => 1,
                                                                             'actcount' => 1,
                                                                             'items' => [
                                                                                          bless( {
                                                                                                   'description' => '/(?:(\\\\d+) :)? (\\\\d+)/x',
                                                                                                   'rdelim' => '/',
                                                                                                   'pattern' => '(?:(\\d+) :)? (\\d+)',
                                                                                                   'hashname' => '__PATTERN1__',
                                                                                                   'lookahead' => 0,
                                                                                                   'ldelim' => '/',
                                                                                                   'mod' => 'x',
                                                                                                   'line' => 232
                                                                                                 }, 'Parse::RecDescent::Token' ),
                                                                                          bless( {
                                                                                                   'hashname' => '__ACTION1__',
                                                                                                   'lookahead' => 0,
                                                                                                   'line' => 233,
                                                                                                   'code' => '{
   #print "TIME_LIMIT: ", Dumper(@item);
   $return = {$item[0] => $item[1]};
 }'
                                                                                                 }, 'Parse::RecDescent::Action' )
                                                                                        ],
                                                                             'line' => undef
                                                                           }, 'Parse::RecDescent::Production' )
                                                                  ],
                                                       'name' => 'TIME_LIMIT',
                                                       'vars' => '',
                                                       'changed' => 0,
                                                       'line' => 231
                                                     }, 'Parse::RecDescent::Rule' ),
                              'print_' => bless( {
                                                   'impcount' => 0,
                                                   'calls' => [
                                                                'STRING'
                                                              ],
                                                   'opcount' => 0,
                                                   'prods' => [
                                                                bless( {
                                                                         'number' => '0',
                                                                         'strcount' => 0,
                                                                         'dircount' => 0,
                                                                         'uncommit' => undef,
                                                                         'error' => undef,
                                                                         'patcount' => 2,
                                                                         'actcount' => 1,
                                                                         'items' => [
                                                                                      bless( {
                                                                                               'description' => '/print|error|warn(ing)?/i',
                                                                                               'rdelim' => '/',
                                                                                               'pattern' => 'print|error|warn(ing)?',
                                                                                               'hashname' => '__PATTERN1__',
                                                                                               'lookahead' => 0,
                                                                                               'ldelim' => '/',
                                                                                               'mod' => 'i',
                                                                                               'line' => 239
                                                                                             }, 'Parse::RecDescent::Token' ),
                                                                                      bless( {
                                                                                               'subrule' => 'STRING',
                                                                                               'matchrule' => 0,
                                                                                               'implicit' => undef,
                                                                                               'argcode' => undef,
                                                                                               'lookahead' => 0,
                                                                                               'line' => 239
                                                                                             }, 'Parse::RecDescent::Subrule' ),
                                                                                      bless( {
                                                                                               'description' => '/\\\\n/',
                                                                                               'rdelim' => '/',
                                                                                               'pattern' => '\\n',
                                                                                               'hashname' => '__PATTERN2__',
                                                                                               'lookahead' => 0,
                                                                                               'ldelim' => '/',
                                                                                               'mod' => '',
                                                                                               'line' => 239
                                                                                             }, 'Parse::RecDescent::Token' ),
                                                                                      bless( {
                                                                                               'hashname' => '__ACTION1__',
                                                                                               'lookahead' => 0,
                                                                                               'line' => 240,
                                                                                               'code' => '{
    $return = $suite->add_test(\'PRINT\', $item{STRING});
  }'
                                                                                             }, 'Parse::RecDescent::Action' )
                                                                                    ],
                                                                         'line' => undef
                                                                       }, 'Parse::RecDescent::Production' )
                                                              ],
                                                   'name' => 'print_',
                                                   'vars' => '',
                                                   'changed' => 0,
                                                   'line' => 238
                                                 }, 'Parse::RecDescent::Rule' ),
                              'paren_enclosed' => bless( {
                                                           'impcount' => 0,
                                                           'calls' => [],
                                                           'opcount' => 0,
                                                           'prods' => [
                                                                        bless( {
                                                                                 'number' => '0',
                                                                                 'strcount' => 0,
                                                                                 'dircount' => 1,
                                                                                 'uncommit' => undef,
                                                                                 'error' => undef,
                                                                                 'patcount' => 0,
                                                                                 'actcount' => 2,
                                                                                 'items' => [
                                                                                              bless( {
                                                                                                       'hashname' => '__DIRECTIVE1__',
                                                                                                       'name' => '<skip: "">',
                                                                                                       'lookahead' => 0,
                                                                                                       'line' => 119,
                                                                                                       'code' => 'my $oldskip = $skip; $skip= ""; $oldskip'
                                                                                                     }, 'Parse::RecDescent::Directive' ),
                                                                                              bless( {
                                                                                                       'hashname' => '__ACTION1__',
                                                                                                       'lookahead' => 0,
                                                                                                       'line' => 119,
                                                                                                       'code' => '{extract_bracketed($text,\'()\')}'
                                                                                                     }, 'Parse::RecDescent::Action' ),
                                                                                              bless( {
                                                                                                       'hashname' => '__ACTION2__',
                                                                                                       'lookahead' => 0,
                                                                                                       'line' => 120,
                                                                                                       'code' => '{
    # create a new Judge based on the condition/testproblem set contained in the
    # parens
    my $judge = TestHarness::Judge->new();
    my $string = $item[1];
    # $string should look lik this
    # (!noh, regression.suite, simple)

    # Remove surrounding parens.
    $string =~ s/\\s*[()]\\s*//g;
 
    # Split comma delimited list.
    my @tests = split /\\s*,\\s*/;

    # process each test/testsuite
    foreach my $test(@tests) {
      my $not = "true";
      my $name = "";

      # parse the test to look for condition/testname
      ($not, $name) = ($test =~ /^\\s* (\\-|!|not|fail|pass)? \\s* ([\\w.]+) \\s*$/ix);
      $judge->add_criterion($name, $not);
    }

    $return = $judge;
  }'
                                                                                                     }, 'Parse::RecDescent::Action' )
                                                                                            ],
                                                                                 'line' => undef
                                                                               }, 'Parse::RecDescent::Production' )
                                                                      ],
                                                           'name' => 'paren_enclosed',
                                                           'vars' => '',
                                                           'changed' => 0,
                                                           'line' => 119
                                                         }, 'Parse::RecDescent::Rule' ),
                              'if_else' => bless( {
                                                    'impcount' => 0,
                                                    'calls' => [
                                                                 'if_block',
                                                                 'elsif_block',
                                                                 'else_block'
                                                               ],
                                                    'opcount' => 0,
                                                    'prods' => [
                                                                 bless( {
                                                                          'number' => '0',
                                                                          'strcount' => 0,
                                                                          'dircount' => 2,
                                                                          'uncommit' => undef,
                                                                          'error' => undef,
                                                                          'patcount' => 0,
                                                                          'actcount' => 1,
                                                                          'items' => [
                                                                                       bless( {
                                                                                                'hashname' => '__DIRECTIVE1__',
                                                                                                'name' => '<skip: "">',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 37,
                                                                                                'code' => 'my $oldskip = $skip; $skip= ""; $oldskip'
                                                                                              }, 'Parse::RecDescent::Directive' ),
                                                                                       bless( {
                                                                                                'subrule' => 'if_block',
                                                                                                'matchrule' => 0,
                                                                                                'implicit' => undef,
                                                                                                'argcode' => undef,
                                                                                                'lookahead' => 0,
                                                                                                'line' => 37
                                                                                              }, 'Parse::RecDescent::Subrule' ),
                                                                                       bless( {
                                                                                                'hashname' => '__DIRECTIVE2__',
                                                                                                'name' => '<commit>',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 37,
                                                                                                'code' => '$commit = 1'
                                                                                              }, 'Parse::RecDescent::Directive' ),
                                                                                       bless( {
                                                                                                'subrule' => 'elsif_block',
                                                                                                'expected' => undef,
                                                                                                'min' => 0,
                                                                                                'argcode' => undef,
                                                                                                'max' => 100000000,
                                                                                                'matchrule' => 0,
                                                                                                'repspec' => 's?',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 37
                                                                                              }, 'Parse::RecDescent::Repetition' ),
                                                                                       bless( {
                                                                                                'subrule' => 'else_block',
                                                                                                'expected' => undef,
                                                                                                'min' => 0,
                                                                                                'argcode' => undef,
                                                                                                'max' => 1,
                                                                                                'matchrule' => 0,
                                                                                                'repspec' => '?',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 37
                                                                                              }, 'Parse::RecDescent::Repetition' ),
                                                                                       bless( {
                                                                                                'hashname' => '__ACTION1__',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 38,
                                                                                                'code' => '{
    print "\\nIF_ELSE: ", Dumper (%item), "\\n\\n" if $debug;
    # Foreach testproblem/suite from the elsif_blocks, assign the anti-judge
    # created in the if_block. 
    # 
    # Foreach testproblem/suite from the else_block, assign the anti-judge from
    # the if_block and from each else_block.
    #
    # Get a list of judges and files.
    #
    # If_block returns a judge.
    my @judges = ($item{if_block});

    # Elsif_block returns an array (judge, (lists of tests/testsuites)). We
    # must process each elsif_block. Each block has to worry about the judges
    # created in all previous elsif_blocks as well as the one created in the
    # if_block.
    foreach (@{$item{\'elsif_block(s?)\'}}) {
      my @tests = @{$_[1]};
      
      # apply the current set of anti-judges
      foreach my $test(@tests) {
	foreach my $judge(@judges) {
	  $test->add_anti_judge($judge);
	}
      }
      
      # Add the locally produced judge
      push @judges, $_[0];
    }

    # Else_block returns a list of tests/testsuites. There can only be one such
    # block but the parser will return it in an array if it exists.
    foreach (@{$item{\'else_block(?)\'}}) {
      # apply the current set of anti-judges
      foreach my $test(@$_) {
	foreach my $judge(@judges) {
	  $test->add_anti_judge($judge);
	}
      }
    }
    #
    $return = 1;
  }'
                                                                                              }, 'Parse::RecDescent::Action' )
                                                                                     ],
                                                                          'line' => undef
                                                                        }, 'Parse::RecDescent::Production' )
                                                               ],
                                                    'name' => 'if_else',
                                                    'vars' => '',
                                                    'changed' => 0,
                                                    'line' => 37
                                                  }, 'Parse::RecDescent::Rule' ),
                              '_alternation_1_of_production_2_of_rule_suitefile' => bless( {
                                                                                             'impcount' => 0,
                                                                                             'calls' => [
                                                                                                          'if_else',
                                                                                                          'line'
                                                                                                        ],
                                                                                             'opcount' => 0,
                                                                                             'prods' => [
                                                                                                          bless( {
                                                                                                                   'number' => '0',
                                                                                                                   'strcount' => 0,
                                                                                                                   'dircount' => 0,
                                                                                                                   'uncommit' => undef,
                                                                                                                   'error' => undef,
                                                                                                                   'patcount' => 0,
                                                                                                                   'actcount' => 0,
                                                                                                                   'items' => [
                                                                                                                                bless( {
                                                                                                                                         'subrule' => 'if_else',
                                                                                                                                         'matchrule' => 0,
                                                                                                                                         'implicit' => undef,
                                                                                                                                         'argcode' => undef,
                                                                                                                                         'lookahead' => 0,
                                                                                                                                         'line' => 263
                                                                                                                                       }, 'Parse::RecDescent::Subrule' )
                                                                                                                              ],
                                                                                                                   'line' => undef
                                                                                                                 }, 'Parse::RecDescent::Production' ),
                                                                                                          bless( {
                                                                                                                   'number' => '1',
                                                                                                                   'strcount' => 0,
                                                                                                                   'dircount' => 0,
                                                                                                                   'uncommit' => undef,
                                                                                                                   'error' => undef,
                                                                                                                   'patcount' => 0,
                                                                                                                   'actcount' => 0,
                                                                                                                   'items' => [
                                                                                                                                bless( {
                                                                                                                                         'subrule' => 'line',
                                                                                                                                         'matchrule' => 0,
                                                                                                                                         'implicit' => undef,
                                                                                                                                         'argcode' => undef,
                                                                                                                                         'lookahead' => 0,
                                                                                                                                         'line' => 263
                                                                                                                                       }, 'Parse::RecDescent::Subrule' )
                                                                                                                              ],
                                                                                                                   'line' => 263
                                                                                                                 }, 'Parse::RecDescent::Production' )
                                                                                                        ],
                                                                                             'name' => '_alternation_1_of_production_2_of_rule_suitefile',
                                                                                             'vars' => '',
                                                                                             'changed' => 0,
                                                                                             'line' => 263
                                                                                           }, 'Parse::RecDescent::Rule' ),
                              'BLANK_LINE' => bless( {
                                                       'impcount' => 0,
                                                       'calls' => [],
                                                       'opcount' => 0,
                                                       'prods' => [
                                                                    bless( {
                                                                             'number' => '0',
                                                                             'strcount' => 0,
                                                                             'dircount' => 0,
                                                                             'uncommit' => undef,
                                                                             'error' => undef,
                                                                             'patcount' => 1,
                                                                             'actcount' => 0,
                                                                             'items' => [
                                                                                          bless( {
                                                                                                   'description' => '/[ \\\\t]*\\\\n/',
                                                                                                   'rdelim' => '/',
                                                                                                   'pattern' => '[ \\t]*\\n',
                                                                                                   'hashname' => '__PATTERN1__',
                                                                                                   'lookahead' => 0,
                                                                                                   'ldelim' => '/',
                                                                                                   'mod' => '',
                                                                                                   'line' => 255
                                                                                                 }, 'Parse::RecDescent::Token' )
                                                                                        ],
                                                                             'line' => undef
                                                                           }, 'Parse::RecDescent::Production' )
                                                                  ],
                                                       'name' => 'BLANK_LINE',
                                                       'vars' => '',
                                                       'changed' => 0,
                                                       'line' => 255
                                                     }, 'Parse::RecDescent::Rule' ),
                              'elsif_block' => bless( {
                                                        'impcount' => 0,
                                                        'calls' => [
                                                                     'paren_enclosed',
                                                                     'bracket_enclosed'
                                                                   ],
                                                        'opcount' => 0,
                                                        'prods' => [
                                                                     bless( {
                                                                              'number' => '0',
                                                                              'strcount' => 0,
                                                                              'dircount' => 0,
                                                                              'uncommit' => undef,
                                                                              'error' => undef,
                                                                              'patcount' => 1,
                                                                              'actcount' => 1,
                                                                              'items' => [
                                                                                           bless( {
                                                                                                    'description' => '/else?\\\\s*if/i',
                                                                                                    'rdelim' => '/',
                                                                                                    'pattern' => 'else?\\s*if',
                                                                                                    'hashname' => '__PATTERN1__',
                                                                                                    'lookahead' => 0,
                                                                                                    'ldelim' => '/',
                                                                                                    'mod' => 'i',
                                                                                                    'line' => 97
                                                                                                  }, 'Parse::RecDescent::Token' ),
                                                                                           bless( {
                                                                                                    'subrule' => 'paren_enclosed',
                                                                                                    'matchrule' => 0,
                                                                                                    'implicit' => undef,
                                                                                                    'argcode' => undef,
                                                                                                    'lookahead' => 0,
                                                                                                    'line' => 97
                                                                                                  }, 'Parse::RecDescent::Subrule' ),
                                                                                           bless( {
                                                                                                    'subrule' => 'bracket_enclosed',
                                                                                                    'matchrule' => 0,
                                                                                                    'implicit' => undef,
                                                                                                    'argcode' => undef,
                                                                                                    'lookahead' => 0,
                                                                                                    'line' => 97
                                                                                                  }, 'Parse::RecDescent::Subrule' ),
                                                                                           bless( {
                                                                                                    'hashname' => '__ACTION1__',
                                                                                                    'lookahead' => 0,
                                                                                                    'line' => 98,
                                                                                                    'code' => '{
    print "\\nELSIF_BLOCK: ", Dumper (%item), "\\\\nn" if $debug;
    # Foreach bracket_enclosed assign the Judge created in the
    # paren_enclosed.
    my $judge = $item{paren_enclosed};
    foreach my $test (@{$item{bracket_enclosed}}) {
      $test->add_judge($judge);
    }

    # I need to return the judge and the included tests.
    my $return = [$judge, $item{bracket_enclosed}];
 }'
                                                                                                  }, 'Parse::RecDescent::Action' )
                                                                                         ],
                                                                              'line' => undef
                                                                            }, 'Parse::RecDescent::Production' )
                                                                   ],
                                                        'name' => 'elsif_block',
                                                        'vars' => '',
                                                        'changed' => 0,
                                                        'line' => 97
                                                      }, 'Parse::RecDescent::Rule' ),
                              'TESTNAME' => bless( {
                                                     'impcount' => 0,
                                                     'calls' => [],
                                                     'opcount' => 0,
                                                     'prods' => [
                                                                  bless( {
                                                                           'number' => '0',
                                                                           'strcount' => 0,
                                                                           'dircount' => 0,
                                                                           'uncommit' => undef,
                                                                           'error' => undef,
                                                                           'patcount' => 1,
                                                                           'actcount' => 1,
                                                                           'items' => [
                                                                                        bless( {
                                                                                                 'description' => '/[-\\\\/\\\\w0-9_]+(\\\\.test)?/i',
                                                                                                 'rdelim' => '/',
                                                                                                 'pattern' => '[-\\/\\w0-9_]+(\\.test)?',
                                                                                                 'hashname' => '__PATTERN1__',
                                                                                                 'lookahead' => 0,
                                                                                                 'ldelim' => '/',
                                                                                                 'mod' => 'i',
                                                                                                 'line' => 168
                                                                                               }, 'Parse::RecDescent::Token' ),
                                                                                        bless( {
                                                                                                 'hashname' => '__ACTION1__',
                                                                                                 'lookahead' => 0,
                                                                                                 'line' => 169,
                                                                                                 'code' => '{
    $return = $item[1];
  }'
                                                                                               }, 'Parse::RecDescent::Action' )
                                                                                      ],
                                                                           'line' => undef
                                                                         }, 'Parse::RecDescent::Production' )
                                                                ],
                                                     'name' => 'TESTNAME',
                                                     'vars' => '',
                                                     'changed' => 0,
                                                     'line' => 166
                                                   }, 'Parse::RecDescent::Rule' ),
                              'comment_line' => bless( {
                                                         'impcount' => 0,
                                                         'calls' => [
                                                                      'COMMENT'
                                                                    ],
                                                         'opcount' => 0,
                                                         'prods' => [
                                                                      bless( {
                                                                               'number' => '0',
                                                                               'strcount' => 0,
                                                                               'dircount' => 1,
                                                                               'uncommit' => undef,
                                                                               'error' => undef,
                                                                               'patcount' => 1,
                                                                               'actcount' => 0,
                                                                               'items' => [
                                                                                            bless( {
                                                                                                     'hashname' => '__DIRECTIVE1__',
                                                                                                     'name' => '<skip: "">',
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 250,
                                                                                                     'code' => 'my $oldskip = $skip; $skip= ""; $oldskip'
                                                                                                   }, 'Parse::RecDescent::Directive' ),
                                                                                            bless( {
                                                                                                     'subrule' => 'COMMENT',
                                                                                                     'matchrule' => 0,
                                                                                                     'implicit' => undef,
                                                                                                     'argcode' => undef,
                                                                                                     'lookahead' => 0,
                                                                                                     'line' => 250
                                                                                                   }, 'Parse::RecDescent::Subrule' ),
                                                                                            bless( {
                                                                                                     'description' => '/\\\\n/',
                                                                                                     'rdelim' => '/',
                                                                                                     'pattern' => '\\n',
                                                                                                     'hashname' => '__PATTERN1__',
                                                                                                     'lookahead' => 0,
                                                                                                     'ldelim' => '/',
                                                                                                     'mod' => '',
                                                                                                     'line' => 250
                                                                                                   }, 'Parse::RecDescent::Token' )
                                                                                          ],
                                                                               'line' => undef
                                                                             }, 'Parse::RecDescent::Production' )
                                                                    ],
                                                         'name' => 'comment_line',
                                                         'vars' => '',
                                                         'changed' => 0,
                                                         'line' => 249
                                                       }, 'Parse::RecDescent::Rule' ),
                              'exclude' => bless( {
                                                    'impcount' => 0,
                                                    'calls' => [
                                                                 'SUITENAME',
                                                                 'TESTNAME'
                                                               ],
                                                    'opcount' => 0,
                                                    'prods' => [
                                                                 bless( {
                                                                          'number' => '0',
                                                                          'strcount' => 0,
                                                                          'dircount' => 0,
                                                                          'uncommit' => undef,
                                                                          'error' => undef,
                                                                          'patcount' => 2,
                                                                          'actcount' => 1,
                                                                          'items' => [
                                                                                       bless( {
                                                                                                'description' => '/ \\\\- |  ! |  not |  exclude/ix',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => ' \\- |  ! |  not |  exclude',
                                                                                                'hashname' => '__PATTERN1__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => 'ix',
                                                                                                'line' => 174
                                                                                              }, 'Parse::RecDescent::Token' ),
                                                                                       bless( {
                                                                                                'subrule' => 'SUITENAME',
                                                                                                'matchrule' => 0,
                                                                                                'implicit' => undef,
                                                                                                'argcode' => undef,
                                                                                                'lookahead' => 0,
                                                                                                'line' => 174
                                                                                              }, 'Parse::RecDescent::Subrule' ),
                                                                                       bless( {
                                                                                                'description' => '/\\\\n/',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => '\\n',
                                                                                                'hashname' => '__PATTERN2__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => '',
                                                                                                'line' => 174
                                                                                              }, 'Parse::RecDescent::Token' ),
                                                                                       bless( {
                                                                                                'hashname' => '__ACTION1__',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 175,
                                                                                                'code' => '{
    $return = $suite->exclude_suite($item{SUITENAME});
  }'
                                                                                              }, 'Parse::RecDescent::Action' )
                                                                                     ],
                                                                          'line' => undef
                                                                        }, 'Parse::RecDescent::Production' ),
                                                                 bless( {
                                                                          'number' => '1',
                                                                          'strcount' => 0,
                                                                          'dircount' => 0,
                                                                          'uncommit' => undef,
                                                                          'error' => undef,
                                                                          'patcount' => 2,
                                                                          'actcount' => 1,
                                                                          'items' => [
                                                                                       bless( {
                                                                                                'description' => '/\\\\- | ! | not | exclude/ix',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => '\\- | ! | not | exclude',
                                                                                                'hashname' => '__PATTERN1__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => 'ix',
                                                                                                'line' => 178
                                                                                              }, 'Parse::RecDescent::Token' ),
                                                                                       bless( {
                                                                                                'subrule' => 'TESTNAME',
                                                                                                'matchrule' => 0,
                                                                                                'implicit' => undef,
                                                                                                'argcode' => undef,
                                                                                                'lookahead' => 0,
                                                                                                'line' => 178
                                                                                              }, 'Parse::RecDescent::Subrule' ),
                                                                                       bless( {
                                                                                                'description' => '/\\\\n/',
                                                                                                'rdelim' => '/',
                                                                                                'pattern' => '\\n',
                                                                                                'hashname' => '__PATTERN2__',
                                                                                                'lookahead' => 0,
                                                                                                'ldelim' => '/',
                                                                                                'mod' => '',
                                                                                                'line' => 178
                                                                                              }, 'Parse::RecDescent::Token' ),
                                                                                       bless( {
                                                                                                'hashname' => '__ACTION1__',
                                                                                                'lookahead' => 0,
                                                                                                'line' => 179,
                                                                                                'code' => '{
    $return = $suite->exclude_test($item{TESTNAME});
  }'
                                                                                              }, 'Parse::RecDescent::Action' )
                                                                                     ],
                                                                          'line' => 178
                                                                        }, 'Parse::RecDescent::Production' )
                                                               ],
                                                    'name' => 'exclude',
                                                    'vars' => '',
                                                    'changed' => 0,
                                                    'line' => 173
                                                  }, 'Parse::RecDescent::Rule' )
                            }
               }, 'Parse::RecDescent' );
}