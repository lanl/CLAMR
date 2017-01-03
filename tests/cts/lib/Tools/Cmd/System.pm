package Tools::Cmd::System;

use strict;

sub import {
    my $class = shift;
    foreach my $flag (@_) {
        my ($key, $val) = split(/=/, $flag, 2);

        if ($key eq 'autoflush') {
            $| = $val;
        }
    }
}

1;