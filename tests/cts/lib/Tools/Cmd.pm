package Tools::Cmd;

use Tools::Check    qw[check];
use Tools::Module   qw[can_load];

use ExtUtils::MakeMaker();
use File::Spec ();
use Config;

use strict;

BEGIN {
    use Exporter    ();
    use vars        qw[ @ISA $VERSION @EXPORT_OK $VERBOSE];

    $VERSION    =   0.02;
    $VERBOSE    =   0;

    @ISA        =   qw[ Exporter ];
    @EXPORT_OK  =   qw[can_run run];
}

### check if we can run some command ###
sub can_run {
    my $command = shift;

    for my $dir (split /$Config{path_sep}/, $ENV{PATH}) {
        my $abs = File::Spec->catfile($dir, $command);
        return $abs if $abs = MM->maybe_command($abs);
    }
}


### Execute a command: $cmd may be a scalar or an arrayref of cmd and args
### $bufout is an scalar ref to store outputs, $verbose can override conf
sub run {
    my %hash = @_;

    my $tmpl = {
        verbose => { default    => $VERBOSE },
        command => { required   => 1,
                     allow      => sub {!(ref $_[1]) or ref $_[1] eq 'ARRAY' }
                   },
    };

    my $args = check( $tmpl, \%hash, $VERBOSE )
                or ( warn(qq[Could not validate input!]), return );





    ### Kludge! This enables autoflushing for each perl process we launched.
    local $ENV{PERL5OPT} .= ' -MTools::Cmd::System=autoflush=1';

    my $verbose     = $args->{verbose};
    my $is_win98    = ($^O eq 'MSWin32' and !Win32::IsWinNT());

    my $err;                # error flag
    my $have_buffer;        # to indicate we executed via IPC::Run or IPC::Open3
                            # only then it makes sence to return the buffers

    my (@buffer,@buferr,@bufout);

    ### STDOUT message handler
    my $_out_handler = sub {
        my $buf = shift;
        return unless defined $buf;

        print STDOUT $buf if $verbose;
        push @buffer, $buf;
        push @bufout, $buf;
    };

    ### STDERR message handler
    my $_err_handler = sub {
        my $buf = shift;
        return unless defined $buf;

        print STDERR $buf if $verbose;
        push @buffer, $buf;
        push @buferr, $buf;
    };

    my $cmd = $args->{command};
    my @cmd = ref ($cmd) ? grep(length, @{$cmd}) : $cmd;

    print qq|Running [@cmd]...| if $verbose;

    ### First, we prefer Barrie Slaymaker's wonderful IPC::Run module.
    if (!$is_win98 and can_load(
                            modules => { 'IPC::Run' => '0.55' },
                            verbose => $verbose && ($^O eq 'MSWin32'),
    ) ) {
        STDOUT->autoflush(1); STDERR->autoflush(1);

        $have_buffer++;

         @cmd = ref($cmd) ? ( [ @cmd ] )
                         : map { /[<>|&]/
                                    ? $_
                                    : [ split / +/ ]
                               } split( /\s*([<>|&])\s*/, $cmd );
         IPC::Run::run(@cmd, \*STDIN, $_out_handler, $_err_handler) or $err++;


    ### Next, IPC::Open3 is know to fail on Win32, but works on Un*x.
    } elsif (   $^O !~ /^(?:MSWin32|cygwin)$/
                and can_load(
                    modules => { map{$_ => '0.0'} qw|IPC::Open3 IO::Select Symbol| },
                    verbose => $verbose
    ) ) {
        my $rv;
        ($rv,$err) = _open3_run(\@cmd, $_out_handler, $_err_handler);
        $have_buffer++;


    ### Abandon all hope; falls back to simple system() on verbose calls.
    } elsif ($verbose) {
        system(@cmd);
        $err = $? ? $? : 0;

    ### Non-verbose system() needs to have STDOUT and STDERR muted.
    } else {
        local *SAVEOUT; local *SAVEERR;

        open(SAVEOUT, ">&STDOUT")
            or warn "couldn't dup STDOUT: $!",      return;
        open(STDOUT, ">".File::Spec->devnull)
            or warn "couldn't reopen STDOUT: $!",   return;

        open(SAVEERR, ">&STDERR")
            or warn "couldn't dup STDERR: $!",      return;
        open(STDERR, ">".File::Spec->devnull)
            or warn "couldn't reopen STDERR: $!",   return;

        system(@cmd);

        open(STDOUT, ">&SAVEOUT")
            or warn "couldn't restore STDOUT: $!", return;
        open(STDERR, ">&SAVEERR")
            or warn "couldn't restore STDERR: $!", return;
    }

    ### unless $err has been set from _open3_run, set it to $? ###
    $err ||= $?;

    return wantarray
                ? $have_buffer
                    ? (!$err, $?, \@buffer, \@bufout, \@buferr)
                    : (!$err, $? )
                : !$err
}


### IPC::Run::run emulator, using IPC::Open3.
sub _open3_run {
    my ($cmdref, $_out_handler, $_err_handler, $verbose) = @_;
    my @cmd = @$cmdref;

    ### Following code are adapted from Friar 'abstracts' in the
    ### Perl Monastery (http://www.perlmonks.org/index.pl?node_id=151886).

    my ($infh, $outfh, $errfh); # open3 handles

    my $pid = eval {
        IPC::Open3::open3(
            $infh   = Symbol::gensym(),
            $outfh  = Symbol::gensym(),
            $errfh  = Symbol::gensym(),
            @cmd,
        )
    };


    return (undef, $@) if $@;

    my $sel = IO::Select->new; # create a select object
    $sel->add($outfh, $errfh); # and add the fhs

    STDOUT->autoflush(1); STDERR->autoflush(1);
    $outfh->autoflush(1) if UNIVERSAL::can($outfh, 'autoflush');
    $errfh->autoflush(1) if UNIVERSAL::can($errfh, 'autoflush');

    while (my @ready = $sel->can_read) {
        foreach my $fh (@ready) { # loop through buffered handles
            # read up to 4096 bytes from this fh.
            my $len = sysread $fh, my($buf), 4096;

            if (not defined $len){
                # There was an error reading
                warn "Error from child: $!";
                return(undef, $!);
            }
            elsif ($len == 0){
                $sel->remove($fh); # finished reading
                next;
            }
            elsif ($fh == $outfh) {
                $_out_handler->($buf);
            } elsif ($fh == $errfh) {
                $_err_handler->($buf);
            } else {
                warn "IO::Select error";
                return(undef, $!);
            }
        }
    }

    waitpid $pid, 0; # wait for it to die
    return 1;
}

1;

__END__

=pod

=head1 NAME

Tools::Cmd - finding and running system commands made easy

=head1 SYNOPSIS

    use Tools::Cmd qw[can_run run];

    my $full_path = can_run('wget') or warn 'wget is not installed!';


    ### commands can be arrayrefs or strings ###
    my $cmd = "$full_path -b theregister.co.uk";
    my $cmd = [$full_path, '-b', 'theregister.co.uk'];

    ### in scalar context ###
    if( run(command => $cmd, verbose => 0) ) {
        print "fetched webpage succesfully\n";
    }


    ### in list context ###
    my( $succes, $error_code, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );

    if( $success ) {
        print "this is what the command printed:\n";
        print join "", @$full_buf;
    }


    ### don't have Tools::Cmd be verbose, ie don't print to stdout or
    ### stderr when running commands -- default is '0'
    $Tools::Cmd::VERBOSE = 0;

=head1 DESCRIPTION

Tools::Cmd allows you to run commands, interactively if desisered,
platform independant but have them still work.

The C<can_run> function can tell you if a certain binary is installed
and if so where, whereas the C<run> function can actually execute any
of the commands you give it and give you a clear return value, as well
as adhere to your verbosity settings.

=head1 FUNCTIONS

=head2 can_run

C<can_run> takes but a single argument: the name of a binary you wish
to locate. C<can_run> works much like the unix binary C<which>, which
scans through your path, looking for the binary you asked for.

Unlike C<which> however, this function is platform independant and
will also work on, for example, Win32.

It will return the full path to the binary you asked for if it was
found, or C<undef> if it was not.

=head2 run

C<run> takes 2 arguments:

=over 4

=item command

This is the command to execute. It may be either a string or an array
reference.
This is a required argument.

=item verbose

This controls whether all output of a command should also be printed
to STDOUT/STDERR or should only be trapped in buffers (NOTE: buffers
require C<IPC::Run> to be installed or your system able to work with
C<IPC::Open3>).

It will default to the global setting of C<$Tools::Cmd::VERBOSE>,
which by default is 0.

=back

C<run> will return a simple C<true> or C<false> when called in scalar
context.
In list context, you will be returned a list of the following items:

=over 4

=item success

A simple boolean indicating if the command executed without errors or
not.

=item errorcode

If the first element of the return value (success) was 0, then some
error occurred. This second element is the error code the command
you requested exited with, if available.

=item full_buffer

This is an arrayreference containing all the output the command
generated.
Note that buffers are only available if you have C<IPC::Run> installed,
or if your system is able to work with C<IPC::Open3> -- See below).
This element will be C<undef> if this is not the case.

=item out_buffer

This is an arrayreference containing all the output sent to STDOUT the
command generated.
Note that buffers are only available if you have C<IPC::Run> installed,
or if your system is able to work with C<IPC::Open3> -- See below).
This element will be C<undef> if this is not the case.

=item error_buffer

This is an arrayreference containing all the output sent to STDERR the
command generated.
Note that buffers are only available if you have C<IPC::Run> installed,
or if your system is able to work with C<IPC::Open3> -- See below).
This element will be C<undef> if this is not the case.

=back

C<run> will try to execute your command using the following logic:

=over 4

=item *

If you are not on windows 98 and have C<IPC::Run> installed, use that
to execute the command. You will have the full output available in
buffers, interactive commands are sure to work  and you are guaranteed
to have your verbosity settings honored cleanly.

=item *

Otherwise, if you are not on MSWin32 or Cygwin, try to execute the
command by using C<IPC::Open3>. Buffers will be available, interactive
commands will still execute cleanly, and also your  verbosity settings
will be adhered to nicely;

=item *

Otherwise, if you have the verbose argument set to true, we fall back
to a simple system() call. We can not capture any buffers, but
interactive commands will still work.

=item *

Otherwise we will try and temporarily redirect STDERR and STDOUT, do a
system() call with your command and then re-open STDERR and STDOUT.
This is the method of last resort and will still allow you to execute
your commands cleanly. However, no buffers will be available.

=head1 Global Variables

The behaviour of Tools::Cmd can be altered by changing the following
global variables:

=head2 $Tools::Cmd::VERBOSE

This controls whether Tools::Cmd will print any output from the
commands to the screen or not. The default is 0;

=head1 See Also

C<IPC::Run>, C<IPC::Open3>

=head1 AUTHOR

This module by
Jos Boumans E<lt>kane@cpan.orgE<gt>.

=head1 COPYRIGHT

This module is
copyright (c) 2002 Jos Boumans E<lt>kane@cpan.orgE<gt>.
All rights reserved.

This library is free software;
you may redistribute and/or modify it under the same
terms as Perl itself.
