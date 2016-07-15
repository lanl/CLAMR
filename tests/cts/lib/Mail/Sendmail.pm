package Mail::Sendmail;
# Mail::Sendmail by Milivoj Ivkovic <mi\x40alma.ch>
# see embedded POD documentation after __END__
# or http://alma.ch/perl/mail.html

=head1 NAME

Mail::Sendmail v. 0.79 - Simple platform independent mailer

=cut

$VERSION = '0.79';

# *************** Configuration you may want to change *******************
# You probably want to set your SMTP server here (unless you specify it in
# every script), and leave the rest as is. See pod documentation for details

%mailcfg = (
    # List of SMTP servers:
    'smtp'    => [ qw( localhost ) ],
    #'smtp'    => [ qw( mail.mydomain.com ) ], # example

    'from'    => '', # default sender e-mail, used when no From header in mail

    'mime'    => 1, # use MIME encoding by default

    'retries' => 1, # number of retries on smtp connect failure
    'delay'   => 1, # delay in seconds between retries

    'tz'      => '', # only to override automatic detection
    'port'    => 25, # change it if you always use a non-standard port
    'debug'   => 0 # prints stuff to STDERR
);

# *******************************************************************

require Exporter;
use strict;
use vars qw(
            $VERSION
            @ISA
            @EXPORT
            @EXPORT_OK
            %mailcfg
            $address_rx
            $debug
            $log
            $error
            $retry_delay
            $connect_retries
           );

use Socket;
use Time::Local; # for automatic time zone detection
use Sys::Hostname; # for use of hostname in HELO

# use MIME::QuotedPrint if available and configured in %mailcfg
eval("use MIME::QuotedPrint");
$mailcfg{'mime'} &&= (!$@);

@ISA        = qw(Exporter);
@EXPORT     = qw(&sendmail);
@EXPORT_OK  = qw(
                 %mailcfg
                 time_to_date
                 $address_rx
                 $debug
                 $log
                 $error
                );

# regex for e-mail addresses where full=$1, user=$2, domain=$3
# see pod documentation about this regex

my $word_rx = '[\x21\x23-\x27\x2A-\x2B\x2D\x2F\w\x3D\x3F]+';
my $user_rx = $word_rx         # valid chars
             .'(?:\.' . $word_rx . ')*' # possibly more words preceded by a dot
             ;
my $dom_rx = '\w[-\w]*(?:\.\w[-\w]*)*'; # less valid chars in domain names
my $ip_rx = '\[\d{1,3}(?:\.\d{1,3}){3}\]';

$address_rx = '((' . $user_rx . ')\@(' . $dom_rx . '|' . $ip_rx . '))';
; # v. 0.61

sub time_to_date {
    # convert a time() value to a date-time string according to RFC 822

    my $time = $_[0] || time(); # default to now if no argument

    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @wdays  = qw(Sun Mon Tue Wed Thu Fri Sat);

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)
        = localtime($time);

    my $TZ = $mailcfg{'tz'};
    if ( $TZ eq "" ) {
        # offset in hours
        my $offset  = sprintf "%.1f", (timegm(localtime) - time) / 3600;
        my $minutes = sprintf "%02d", abs( $offset - int($offset) ) * 60;
        $TZ  = sprintf("%+03d", int($offset)) . $minutes;
    }
    return join(" ",
                    ($wdays[$wday] . ','),
                     $mday,
                     $months[$mon],
                     $year+1900,
                     sprintf("%02d:%02d:%02d", $hour, $min, $sec),
                     $TZ
               );
} # end sub time_to_date

sub sendmail {

    $error = '';
    $log = "Mail::Sendmail v. $VERSION - "    . scalar(localtime()) . "\n";

    my $CRLF = "\015\012";
    local $/ = $CRLF;
    local $\ = ''; # to protect us from outside settings
    local $_;

    my (%mail, $k,
        $smtp, $server, $port, $connected, $localhost,
        $fromaddr, $recip, @recipients, $to, $header,
       );

    # -------- a few internal subs ----------
    sub fail {
        # things to do before returning a sendmail failure
        print STDERR @_ if $^W;
        $error .= join(" ", @_) . "\n";
        close S;
        return 0;
    }

    sub socket_write {
        my $i;
        for $i (0..$#_) {
            # accept references, so we don't copy potentially big data
            my $data = ref($_[$i]) ? $_[$i] : \$_[$i];
            if ($mailcfg{'debug'} > 5) {
                if (length($$data) < 500) {
                    print ">", $$data;
                }
                else {
                    print "> [...", length($$data), " bytes sent ...]\n";
                }
            }
            print(S $$data) || return 0;
        }
        1;
    }

    sub socket_read {
        my $response; # for multi-line server responses
        do {
            chomp($_ = <S>);
            print "<$_\n" if $mailcfg{'debug'} > 5;
            if (/^[45]/ or !$_) {
                return; # return false
            }
            $response .= $_;
         } while (/^[\d]+-/);
         return $response;
    }
    # -------- end of internal subs ----------

    # all config keys to lowercase, to prevent typo errors
    foreach $k (keys %mailcfg) {
        if ($k =~ /[A-Z]/) {
            $mailcfg{lc($k)} = $mailcfg{$k};
        }
    }

    # redo mail hash, arranging keys case etc...
    while (@_) {
        $k = shift @_;
        if (!$k and $^W) {
            warn "Received false mail hash key: \'$k\'. Did you forget to put it in quotes?\n";
        }

        # arrange keys case
        $k = ucfirst lc($k);

        $k =~ s/\s*:\s*$//o; # kill colon (and possible spaces) at end, we add it later.
        # uppercase also after "-", so people don't complain that headers case is different
        # than in Outlook.
        $k =~ s/-(.)/"-" . uc($1)/ge;
        $mail{$k} = shift @_;
    }

    $smtp = $mail{'Smtp'} || $mail{'Server'};
    unshift @{$mailcfg{'smtp'}}, $smtp if ($smtp and $mailcfg{'smtp'}->[0] ne $smtp);

    # delete non-header keys, so we don't send them later as mail headers
    # I like this syntax, but it doesn't seem to work with AS port 5.003_07:
    # delete @mail{'Smtp', 'Server'};
    # so instead:
    delete $mail{'Smtp'}; delete $mail{'Server'};

    $mailcfg{'port'} = $mail{'Port'} || $mailcfg{'port'} || 25;
    delete $mail{'Port'};

    {    # don't warn for undefined values below
        local $^W = 0;
        $mail{'Message'} = join("", $mail{'Message'}, $mail{'Body'}, $mail{'Text'});
    }

    # delete @mail{'Body', 'Text'};
    delete $mail{'Body'}; delete $mail{'Text'};

    # Extract 'From:' e-mail address to use as envelope sender

    $fromaddr = $mail{'Sender'} || $mail{'From'} || $mailcfg{'from'};
    delete $mail{'Sender'};
    unless ($fromaddr =~ /$address_rx/) {
        return fail("Bad or missing From address: \'$fromaddr\'");
    }
    $fromaddr = $1;

    # add Date header if needed
    $mail{Date} ||= time_to_date() ;
    $log .= "Date: $mail{Date}\n";

    # cleanup message, and encode if needed
    $mail{'Message'} =~ s/\r\n/\n/go;     # normalize line endings, step 1 of 2 (next step after MIME encoding)

    $mail{'Mime-Version'} ||= '1.0';
    $mail{'Content-Type'} ||= 'text/plain; charset="iso-8859-1"';

    unless ( $mail{'Content-Transfer-Encoding'}
          || $mail{'Content-Type'} =~ /multipart/io )
    {
        if ($mailcfg{'mime'}) {
            $mail{'Content-Transfer-Encoding'} = 'quoted-printable';
            $mail{'Message'} = encode_qp($mail{'Message'});
        }
        else {
            $mail{'Content-Transfer-Encoding'} = '8bit';
            if ($mail{'Message'} =~ /[\x80-\xFF]/o) {
                $error .= "MIME::QuotedPrint not present!\nSending 8bit characters, hoping it will come across OK.\n";
                warn "MIME::QuotedPrint not present!\n",
                     "Sending 8bit characters without encoding, hoping it will come across OK.\n"
                     if $^W;
            }
        }
    }

    $mail{'Message'} =~ s/^\./\.\./gom;     # handle . as first character
    $mail{'Message'} =~ s/\n/$CRLF/go; # normalize line endings, step 2.

    # Get recipients
    {    # don't warn for undefined values below
        local $^W = 0;
        $recip = join(", ", $mail{To}, $mail{Cc}, $mail{Bcc});
    }

    delete $mail{'Bcc'};

    @recipients = ();
    while ($recip =~ /$address_rx/go) {
        push @recipients, $1;
    }
    unless (@recipients) {
        return fail("No recipient!")
    }

    # get local hostname for polite HELO
    $localhost = hostname() || 'localhost';

    foreach $server ( @{$mailcfg{'smtp'}} ) {
        # open socket needs to be inside this foreach loop on Linux,
        # otherwise all servers fail if 1st one fails !??! why?
        unless ( socket S, AF_INET, SOCK_STREAM, scalar(getprotobyname 'tcp') ) {
            return fail("socket failed ($!)")
        }

        print "- trying $server\n" if $mailcfg{'debug'} > 1;

        $server =~ s/\s+//go; # remove spaces just in case of a typo
        # extract port if server name like "mail.domain.com:2525"
        $port = ($server =~ s/:(\d+)$//o) ? $1 : $mailcfg{'port'};
        $smtp = $server; # save $server for use outside foreach loop

        my $smtpaddr = inet_aton $server;
        unless ($smtpaddr) {
            $error .= "$server not found\n";
            next; # next server
        }

        my $retried = 0; # reset retries for each server
        while ( ( not $connected = connect S, pack_sockaddr_in($port, $smtpaddr) )
            and ( $retried < $mailcfg{'retries'} )
              ) {
            $retried++;
            $error .= "connect to $server failed ($!)\n";
            print "- connect to $server failed ($!)\n" if $mailcfg{'debug'} > 1;
            print "retrying in $mailcfg{'delay'} seconds...\n" if $mailcfg{'debug'} > 1;
            sleep $mailcfg{'delay'};
        }

        if ( $connected ) {
            print "- connected to $server\n" if $mailcfg{'debug'} > 3;
            last;
        }
        else {
            $error .= "connect to $server failed\n";
            print "- connect to $server failed, next server...\n" if $mailcfg{'debug'} > 1;
            next; # next server
        }
    }

    unless ( $connected ) {
        return fail("connect to $smtp failed ($!) no (more) retries!")
    };

    {
        local $^W = 0; # don't warn on undefined variables
        # Add info to log variable
        $log .= "Server: $smtp Port: $port\n"
              . "From: $fromaddr\n"
              . "Subject: $mail{Subject}\n"
              . "To: ";
    }

    my($oldfh) = select(S); $| = 1; select($oldfh);

    socket_read()
        || return fail("Connection error from $smtp on port $port ($_)");
    socket_write("HELO $localhost$CRLF")
        || return fail("send HELO error");
    socket_read()
        || return fail("HELO error ($_)");
    socket_write("MAIL FROM: <$fromaddr>$CRLF")
        || return fail("send MAIL FROM: error");
    socket_read()
        || return fail("MAIL FROM: error ($_)");

    foreach $to (@recipients) {
        socket_write("RCPT TO: <$to>$CRLF")
            || return fail("send RCPT TO: error");
        socket_read()
            || return fail("RCPT TO: error ($_)");
        $log .= "$to\n    ";
    }

    # start data part

    socket_write("DATA$CRLF")
        || return fail("send DATA error");
    socket_read()
        || return fail("DATA error ($_)");

    # print headers
    foreach $header (keys %mail) {
        next if $header eq "Message";
        $mail{$header} =~ s/\s+$//o; # kill possible trailing garbage
        socket_write("$header: $mail{$header}$CRLF")
            || return fail("send $header: error");
    };

    #- test diconnecting from network here, to see what happens
    #- print STDERR "DISCONNECT NOW!\n";
    #- sleep 4;
    #- print STDERR "trying to continue, expecting an error... \n";

    # send message body (passed as a reference, in case it's big)
    socket_write($CRLF, \$mail{'Message'}, "$CRLF.$CRLF")
           || return fail("send message error");
    socket_read()
        || return fail("message transmission error ($_)");
    $log .= "\nResult: $_";

    # finish
    socket_write("QUIT$CRLF")
           || return fail("send QUIT error");
    socket_read();
    close S;

    return 1;
} # end sub sendmail

1;
__END__

=head1 SYNOPSIS

  use Mail::Sendmail;

  %mail = ( To      => 'you@there.com',
            From    => 'me@here.com',
            Message => "This is a very short message"
           );

  sendmail(%mail) or die $Mail::Sendmail::error;

  print "OK. Log says:\n", $Mail::Sendmail::log;

=head1 DESCRIPTION

Simple platform independent e-mail from your perl script. Only requires
Perl 5 and a network connection.

Mail::Sendmail contains mainly &sendmail, which takes a hash with the
message to send and sends it. It is intended to be very easy to setup and
use. See also L<"FEATURES"> below.

=head1 INSTALLATION

=over 4

=item Best

C<perl -MCPAN -e "install Mail::Sendmail">

=item Traditional

    perl Makefile.PL
    make
    make test
    make install

=item Manual

Copy Sendmail.pm to Mail/ in your Perl lib directory.

    (eg. c:\Perl\site\lib\Mail\
     or  /usr/lib/perl5/site_perl/Mail/
     or whatever it is on your system.
     They are listed when you type C< perl -V >)

=item ActivePerl's PPM

ppm install --location=http://alma.ch/perl/ppm Mail-Sendmail

But this way you don't get a chance to have a look at other files (Changes,
Todo, test.pl, ...).

=back

At the top of Sendmail.pm, set your default SMTP server(s), unless you specify
it with each message, or want to use the default (localhost).

Install MIME::QuotedPrint. This is not required but strongly recommended.

=head1 FEATURES

Automatic time zone detection, Date: header, MIME quoted-printable encoding
(if MIME::QuotedPrint installed), all of which can be overridden.

Bcc: and Cc: support.

Allows real names in From:, To: and Cc: fields

Doesn't send an X-Mailer: header (unless you do), and allows you to send any
header(s) you want.

Configurable retries and use of alternate servers if your mail server is
down

Good plain text error reporting

=head1 LIMITATIONS

Headers are not encoded, even if they have accented characters.

No suport for the SMTP AUTH extension.

Since the whole message is in memory, it's not suitable for
sending very big attached files.

The SMTP server has to be set manually in Sendmail.pm or in your script,
unless you have a mail server on localhost.

Doesn't work on OpenVMS, I was told. Cannot test this myself.

=head1 CONFIGURATION

=over 4

=item Default SMTP server(s)

This is probably all you want to configure. It is usually done through
I<$mailcfg{smtp}>, which you can edit at the top of the Sendmail.pm file.
This is a reference to a list of SMTP servers. You can also set it from
your script:

C<unshift @{$Mail::Sendmail::mailcfg{'smtp'}} , 'my.mail.server';>

Alternatively, you can specify the server in the I<%mail> hash you send
from your script, which will do the same thing:

C<$mail{smtp} = 'my.mail.server';>

A future version will (hopefully) try to set useful defaults for you
during the Makefile.PL.

=item Other configuration settings

See I<%mailcfg> under L<"DETAILS"> below for other configuration options.

=back

=head1 DETAILS

=head2 sendmail()

sendmail is the only thing exported to your namespace by default

C<sendmail(%mail) || print "Error sending mail: $Mail::Sendmail::error\n";>

It takes a hash containing the full message, with keys for all headers,
body, and optionally for another non-default SMTP server and/or port.

It returns 1 on success or 0 on error, and rewrites
C<$Mail::Sendmail::error> and C<$Mail::Sendmail::log>.

Keys are NOT case-sensitive.

The colon after headers is not necessary.

The Body part key can be called 'Body', 'Message' or 'Text'.

The SMTP server key can be called 'Smtp' or 'Server'. If the connection to
this one fails, the other ones in C<$mailcfg{smtp}> will still be tried.

The following headers are added unless you specify them yourself:

    Mime-Version: 1.0
    Content-Type: 'text/plain; charset="iso-8859-1"'

    Content-Transfer-Encoding: quoted-printable
    or (if MIME::QuotedPrint not installed)
    Content-Transfer-Encoding: 8bit

    Date: [string returned by time_to_date()]

If you wish to use an envelope sender address different than the
From: address, set C<$mail{Sender}> in your %mail hash.

The following are not exported by default, but you can still access them
with their full name, or request their export on the use line like in:
C<use Mail::Sendmail qw(sendmail $address_rx time_to_date);>

=head2 Mail::Sendmail::time_to_date()

convert time ( as from C<time()> ) to an RFC 822 compliant string for the
Date header. See also L<"%Mail::Sendmail::mailcfg">.

=head2 $Mail::Sendmail::error

When you don't run with the B<-w> flag, the module sends no errors to
STDERR, but puts anything it has to complain about in here. You should
probably always check if it says something.

=head2 $Mail::Sendmail::log

A summary that you could write to a log file after each send

=head2 $Mail::Sendmail::address_rx

A handy regex to recognize e-mail addresses.

A correct regex for valid e-mail addresses was written by one of the judges
in the obfuscated Perl contest... :-) It is quite big. This one is an
attempt to a reasonable compromise, and should accept all real-world
internet style addresses. The domain part is required and comments or
characters that would need to be quoted are not supported.

  Example:
    $rx = $Mail::Sendmail::address_rx;
    if (/$rx/) {
      $address=$1;
      $user=$2;
      $domain=$3;
    }

=head2 %Mail::Sendmail::mailcfg

This hash contains all configuration options. You normally edit it once (if
ever) in Sendmail.pm and forget about it, but you could also access it from
your scripts. For readability, I'll assume you have imported it
(with something like C<use Mail::Sendmail qw(sendmail %mailcfg)>).

The keys are not case-sensitive: they are all converted to lowercase before
use. Writing C<$mailcfg{Port} = 2525;> is OK: the default $mailcfg{port}
(25) will be deleted and replaced with your new value of 2525.

=over 4

=item $mailcfg{smtp}

C<$mailcfg{smtp} = [qw(localhost my.other.mail.server)];>

This is a reference to a list of smtp servers, so if your main server is
down, the module tries the next one. If one of your servers uses a special
port, add it to the server name with a colon in front, to override the
default port (like in my.special.server:2525).

Default: localhost.

=item $mailcfg{from}

C<$mailcfg{from} = 'Mailing script me@mydomain.com';>

From address used if you don't supply one in your script. Should not be of
type 'user@localhost' since that may not be valid on the recipient's
host.

Default: undefined.

=item $mailcfg{mime}

C<$mailcfg{mime} = 1;>

Set this to 0 if you don't want any automatic MIME encoding. You normally
don't need this, the module should 'Do the right thing' anyway.

Default: 1;

=item $mailcfg{retries}

C<$mailcfg{retries} = 1;>

How many times should the connection to the same SMTP server be retried in
case of a failure.

Default: 1;

=item $mailcfg{delay}

C<$mailcfg{delay} = 1;>

Number of seconds to wait between retries. This delay also happens before
trying the next server in the list, if the retries for the current server
have been exhausted. For CGI scripts, you want few retries and short delays
to return with a results page before the http connection times out. For
unattended scripts, you may want to use many retries and long delays to
have a good chance of your mail being sent even with temporary failures on
your network.

Default: 1 (second);

=item $mailcfg{tz}

C<$mailcfg{tz} = '+0800';>

Normally, your time zone is set automatically, from the difference between
C<time()> and C<gmtime()>. This allows you to override automatic detection
in cases where your system is confused (such as some Win32 systems in zones
which do not use daylight savings time: see Microsoft KB article Q148681)

Default: undefined (automatic detection at run-time).

=item $mailcfg{port}

C<$mailcfg{port} = 25;>

Port used when none is specified in the server name.

Default: 25.

=item $mailcfg{debug}

C<$mailcfg{debug} = 0;>

Prints stuff to STDERR. Current maximum is 6, which prints the whole SMTP
session, except data exceeding 500 bytes.

Default: 0;

=back

=head2 $Mail::Sendmail::VERSION

The package version number (you can not import this one)

=head2 Configuration variables from previous versions

The following global variables were used in version 0.74 for configuration.
As from version 0.78_1, they are not supported anymore.
Use the I<%mailcfg> hash if you need to access the configuration
from your scripts.

=over 4

=item $Mail::Sendmail::default_smtp_server

=item $Mail::Sendmail::default_smtp_port

=item $Mail::Sendmail::default_sender

=item $Mail::Sendmail::TZ

=item $Mail::Sendmail::connect_retries

=item $Mail::Sendmail::retry_delay

=item $Mail::Sendmail::use_MIME

=back

=head1 ANOTHER EXAMPLE

  use Mail::Sendmail;

  print "Testing Mail::Sendmail version $Mail::Sendmail::VERSION\n";
  print "Default server: $Mail::Sendmail::mailcfg{smtp}->[0]\n";
  print "Default sender: $Mail::Sendmail::mailcfg{from}\n";

  %mail = (
      #To      => 'No to field this time, only Bcc and Cc',
      #From    => 'not needed, use default',
      Bcc     => 'Someone <him@there.com>, Someone else her@there.com',
      # only addresses are extracted from Bcc, real names disregarded
      Cc      => 'Yet someone else <xz@whatever.com>',
      # Cc will appear in the header. (Bcc will not)
      Subject => 'Test message',
      'X-Mailer' => "Mail::Sendmail version $Mail::Sendmail::VERSION",
  );


  $mail{Smtp} = 'special_server.for-this-message-only.domain.com';
  $mail{'X-custom'} = 'My custom additionnal header';
  $mail{'mESSaGE : '} = "The message key looks terrible, but works.";
  # cheat on the date:
  $mail{Date} = Mail::Sendmail::time_to_date( time() - 86400 );

  if (sendmail %mail) { print "Mail sent OK.\n" }
  else { print "Error sending mail: $Mail::Sendmail::error \n" }

  print "\n\$Mail::Sendmail::log says:\n", $Mail::Sendmail::log;
 
Also see http://alma.ch/perl/Mail-Sendmail-FAQ.html for examples
of HTML mail and sending attachments. 

=head1 CHANGES

Main changes since version 0.78:

Added "/" (\x2F) as a valid character in mailbox part.

Removed old configuration variables which are not used anymore
since version 0.74.

Added support for different envelope sender (through C<$mail{Sender}>)

Changed case of headers: first character after "-" also uppercased

Support for multi-line server responses

Localized $\ and $_

Some internal rewrites and documentation updates

Fixed old bug of dot as 76th character on line disappearing.

Fixed very old bug where port number was not extracted from
stuff like 'my.server:2525'.

Fixed time_to_date bug with negative half-hour zones (only Newfoundland?)

Added seconds to date string

Now uses Sys::Hostname to get the hostname for HELO. (This may break the
module on some very old Win32 Perls where Sys::Hostname was broken)

Enable full session output for debugging

See the F<Changes> file for the full history. If you don't have it
because you installed through PPM, you can also find the latest
one on F<http://alma.ch/perl/scripts/Sendmail/Changes>.

=head1 AUTHOR

Milivoj Ivkovic <mi\x40alma.ch> ("\x40" is "@" of course)

=head1 NOTES

MIME::QuotedPrint is used by default on every message if available. It
allows reliable sending of accented characters, and also takes care of
too long lines (which can happen in HTML mails). It is available in the
MIME-Base64 package at http://www.perl.com/CPAN/modules/by-module/MIME/ or
through PPM.

Look at http://alma.ch/perl/Mail-Sendmail-FAQ.html for additional
info (CGI, examples of sending attachments, HTML mail etc...)

You can use this module freely. (Someone complained this is too vague.
So, more precisely: do whatever you want with it, but be warned that
terrible things will happen to you if you use it badly, like for sending
spam, or ...?)

Thanks to the many users who sent me feedback, bug reports, suggestions, etc.
And please excuse me if I forgot to answer your mail. I am not always reliabe
in answering mail. I intend to set up a mailing list soon.

Last revision: 06.02.2003. Latest version should be available on
CPAN: F<http://www.cpan.org/modules/by-authors/id/M/MI/MIVKOVIC/>.

=cut
