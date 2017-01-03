package read_output_files;

use     POSIX qw( strtod );
#use     diagnostics;
use     warnings;
use     Carp;
use     vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );
use     Exporter;
use     Cwd;
use     Time::Local;

use my_utils qw (
                 extrema
                 print_error
                 print_perl_obj
                 status_bar
                 which_exec
                );

$VERSION   = 1.00;

@ISA       = qw(
                Exporter
               );

@EXPORT    = qw(
               );

@EXPORT_OK = qw(
                parse_output_file
                parse_output_file_finish
               );


my(
   $include,
   $inc_dir,
   $extras,
   $GNUMBER_REGEXP,
   $HAS_EXTRAS,
   $require_file,
  );
$GNUMBER_REGEXP = '[+-]?\.?[0-9]+\.?[0-9]*([eE][+-]?\d+)?';

# include extras
#require read_output_files_extras;
#......................
#...global variables...
#......................
my(
   %G_RUN_COMMAND_SEEN,
  );

sub parse_output_file
  {
    my(
       $lines_ref,
       $first_cycle_ref,
       $nummat_ref,
       $data_ref,
       $found_fields_ref,
       $dim_ref,
       $verbosity_in
      ) = @_;
    my(
       $cycle,
       @dim_name,
       $dist,
       %field,
       @fields,
       $field_name,
       $hours,
       $ln,
       $line,
       $line_orig,
       $ln_max,
       $mday,
       $min,
       $mon,
       $sec,
       $time,
       $time_old,
       $year,
       @vals,
       $verbosity
      );

    # verbosity
    $verbosity = 0;
    if( defined( $verbosity_in ) ){
        $verbosity = $verbosity_in;
    }

    # strip leading/trailing whitespace and trash
    if( $verbosity > 0 ){
        print "  Cleaning up lines.\n";
    }
    grep( s/^\s*(.*?)\s*$/$1/, @{$lines_ref} );
    if( $#lines_ref > 0 && $$lines_ref[0] =~ /^0:/ ){
        grep( s/^0:\s+//, @$lines_ref );
    }

    if( $verbosity > 0 ){
        print "  Processing lines:\n    ";
    }
    $ln = 0;
    $ln_max = $#{$lines_ref};
    while ( $ln <= $ln_max ) {
      if( $verbosity > 0 ){
          &status_bar($ln+1, $ln_max+1);
      }
      $line_orig = $$lines_ref[$ln]; $ln++;
      $line = $line_orig;

      # nothing for blank lines
      if( $line !~ /\S/ ){
      }

      # input file parse (just simple for now)
      elsif ( $line =~ /^(tmax|ncmax)\s+=\s+(\S+)/ ) {
        $$data_ref{$1}[0] = $2;
      }

      # creating dump file (only in log files)
      elsif ( $line =~ /^(Create|Creating) PIO file:.*-dmp.*\s+cycle\s*=\s*(\d+)$/ ) {
        $cycle = $2;
        # get to dump time spent line
        $done = "false";
        $field_name = "dmp_write_time";
        while ( $done eq "false" ) {
          $line = $$lines_ref[$ln]; $ln++;
          if ( $line =~ /^Closed (Parallel IO|PIO) file:.*\s+sec\s+=\s+(\S+)/ ) {
            $$data_ref{$field_name}[$cycle] = $2;
            $$found_fields_ref{$field_name} = "true";
            $$found_fields_ref{ordering_diff}{"^dmp_write_time"} = "";
            $done = "true";
          }
          if ( $ln > $ln_max ) {
            $done = "true";
          }
        }
      }

      # creating hdf files
      elsif ( $line =~ /^HDFPLT called.*\s+(\S+)\s+(\S+)\s+(\S+)$/ ) {
        $cycle = $1;
        # get to dump time spent line
        $done = "false";
        $field_name = "hdf_write_time";
        while ( $done eq "false" ) {
          $line = $$lines_ref[$ln]; $ln++;
          if ( $line =~ /^HDFPLT finished.\s*cpu\s*=\s*(\S+)/ ) {
            $$data_ref{$field_name}[$cycle] = $1;
            $$found_fields_ref{$field_name} = "true";
            $$found_fields_ref{ordering_diff}{"^${field_name}"} = "";
            $done = "true";
          }
          if ( $ln > $ln_max ) {
            $done = "true";
          }
        }
      }

      #...### line
      elsif ( $line =~ /^\#\#\#\s+
                        (\S+)\s+    # 1
                        (\S+)\s+    # 2
                        (\S+)\s+    # 3
                        (\S+)\s+    # 4
                        (\S+)\s+    # 5
                        (\S+)\s+    # 6
                        (\S+)\s+    # 7
                        (\S+)\s+    # 8
                        (\S+)\s+    # 9
                        (\S+)\s+    # 10
                        ((\S+)\s+)? # 11
                        (\S+)\s+    # 13
                        (::)\s+     # 14
                        (\S+)\s+    # 15
                        (\S+.*)$    # 16
                        /x ) {
        undef(%field);
        $$found_fields_ref{ppp} = "true";
        $field{"cyc"}       = $1;
        $field{time}        = $2;
        $field{"dt"}        = $3;
        $field{"ncell"}     = $4;
        $field{"#pe"}       = $5;
        $field{"lvl"}       = $6;
        $field{"#ritr"}     = $7;
        $field{"#citr"}     = $8;
        $field{"sumcpu"}    = $9;
        $field{"cc/s/p"}    = $10;
        $field{"sumRSS_GB"} = $12;
        $field{"err"}       = $13;
        $field{"date"}      = $15;
        $field{"tstep"}     = $16;
        $field{"tstep"} =~ s/\s+/_/g;
        if ( $field{"cyc"} =~ /^($GNUMBER_REGEXP)$/ ) {
          $cycle = $field{cyc};
          $time = $field{time};
          foreach $field_name ((
                                "cyc",
                                "dt",
                                "ncell",
                                "#ritr",
                                "sumcpu",
                                "cc/s/p",
                                "sumRSS_GB",
                                "err",
                                "tstep"
                               )) {
            $$found_fields_ref{$field_name} = "true";
            $$found_fields_ref{ordering_diff}{"^cc\/s\/p"} = "";
            $$found_fields_ref{ordering_diff}{"^sumcpu"} = "";
            $$found_fields_ref{ordering_diff}{"^sumwall"} = "";
            $$found_fields_ref{ordering_diff}{"^sumRSS_"} = "";
            if( defined($field{$field_name})){
                $$data_ref{$field_name}[$cycle] = $field{$field_name};
            }
          }
          $field{"date"} =~ s/(\S{2})(\S{2}):(\S{2})(\S{2})(\S{2})/$1 $2 $3 $4 $5/;
          $year = "YEAR";
          ($mon, $mday, $hours, $min, $sec) = split( /\s/, $field{date});
          $mon = $mon - 1;
          $$data_ref{"secs"}[$cycle] = "$sec $min $hours $mday $mon $year";
        }
      }

      #...long cycle info...
      elsif ( $line =~ /^
                        cycle\s+          #  1
                        t\s+              #  2
                        dtnext\s+         #  3
                        timestep\s+       #  4
                        cstb\s+           #  5
                        tpct\s+           #  6
                        epct\s+           #  7
                        ritr\s+           #  8
                        hitr\s+           #  9
                        sumritr\s+        # 10
                        wallhr\s+         # 11
                        sumwallhr\s+      # 12
                        sumcpuhr\s+       # 13
                      date:time\s*      # 14
                        /x ) {
        $line_orig = $$lines_ref[$ln]; $ln++;
        $line = $line_orig;
        #...replace spaces in timestep with underscores
        @fields = split(/\s+/, $line );
        while ( $#fields > 13 ) {
          $fields[3] = "$fields[3]_$fields[4]";
          $fields[4] = "";
          $line = join( " ", @fields );
          @fields = split( /\s+/, $line );
        }
        if ( $line =~ /^
                       (\S+)\s+ # 1
                       (\S+)\s+ # 2
                       (\S+)\s+ # 3
                       (\S+)\s+ # 4
                       (\S+)\s+ # 5
                       (\S+)\s+ # 6
                       (\S+)\s+ # 7
                       (\S+)\s+ # 8
                       (\S+)\s+ # 9
                       (\S+)\s+ # 10
                       (\S+)\s+ # 11
                       (\S+)\s+ # 12
                       (\S+)\s+ # 13
                       (\S+)    # 14
                       $/x ) {
            undef(%field);
            $$found_fields_ref{long} = "true";
            $$found_fields_ref{ordering_diff}{"^sumcpu"} = "";
            $$found_fields_ref{ordering_diff}{"^sumwall"} = "";
            $field{"cyc"}       = $1;
            $field{time}        = $2;
            $field{"dt"}        = $3;
            $field{"tstep"}     = $4;
            $field{"date"}      = $14;
            $field{"#ritr"}     = $8;
            $field{"sumwallhr"} = $12;
            $field{"sumcpu"}    = $13;
            $field{"tstep"}     =~ s/\s+/_/g;
            if( $field{"cyc"} =~ /^($GNUMBER_REGEXP)$/ )
              {
                $time = $field{time};
                $cycle = $field{cyc};
                foreach $field_name ( keys %field )
                  {
                    if( $field_name eq "date" ||
                        $field_name eq "time" ){
                      next;
                    }
                    $$data_ref{$field_name}[$cycle] = $field{$field_name};
                  }
                $field{"date"} =~ s/(\S{4})(\S{2})(\S{2}):(\S{2})(\S{2})(\S{2})/$1 $2 $3 $4 $5 $6/;
                ($year, $mon, $mday, $hours, $min, $sec) = split( /\s/, $field{date});
                $year = $year - 1900;
                $mon = $mon - 1;
                $$data_ref{"secs"}[$cycle] = "$sec $min $hours $mday $mon $year";
              }
          }
      }

      #...cell info...
      elsif( $line =~ /^.*\ssum_cell\s.*\savg_cell\s.*\savg_top\s/ ){
          @fields = split( /\s+/, $line );
          $line_orig = $$lines_ref[$ln]; $ln++;
          $line = $line_orig;
          @vals = split( /\s+/, $line );
          for( $i = 0; $i <= $#fields; $i++ ){
              if( $fields[$i] eq "nummat" ){
                  $$nummat_ref = $vals[$i];
              }
              elsif( $fields[$i] eq "sum_cell"){
                  $$data_ref{ncell}[$cycle] = $vals[$i];                  
              }
              elsif( $fields[$i] eq "sum_top"){
                  $$data_ref{pct_top}[$cycle] = 100.0*$vals[$i]/$$data_ref{ncell}[$cycle];
              }
              elsif( $fields[$i] =~ /^(max_cell|avg_cell|min_cell)$/ ){
                  $$data_ref{$1}[$cycle] = $vals[$i];                  
              }
              elsif( $fields[$i] eq "free_mem"){
                  $$data_ref{procmon_machine_free_min}[$cycle] = $vals[$i];                  
              }
              elsif( $fields[$i] eq "pct"){
                  $$data_ref{procmon_free_mem_pct_min}[$cycle] = $vals[$i];                  
              }
          }
        $$found_fields_ref{cell} = "true";
        $$found_fields_ref{ordering_diff}{"^procmon"} = "";
      }

      #...what info...
      elsif( $line =~ /^
                       what\s+    #  1
                       max\s+     #  2
                       cell\s+    #  3
                       (xc)\s*    #  4
                       (yc)?\s*   #  5
                       (zc)?\s*   #  6
                       min\s+     #  7
                       cell\s+    #  8
                       xc\s*      #  9
                       (yc)?\s*   # 10
                       (zc)?\s*   # 11
                       /x ){
          $$found_fields_ref{what} = "";
          $$dim_ref = 1;
          $dim_name[1] = "1_$1";
          if ( defined($2) ) {
            $$dim_ref = 2;
            $dim_name[2] = "2_$2";
          }
          if ( defined($3) ) {
            $$dim_ref = 3;
            $dim_name[3] = "3_$3";
          }
          $done = "false";
          while ( $done eq "false" ) {
            $line_orig = $$lines_ref[$ln]; $ln++;
            if ( $line_orig !~ /^[A-Z]/ ||
                 $line_orig =~ /^FYI, / ) {
              $done = "true";
              next;
            }
            $line = $line_orig;
            #...boo! var name might have space in it
            $line =~ s/([A-Z])\s+([A-Z])/${1}_${2}/g;
            @fields = split( /\s+/, $line );
            $index = 0;
            $field_name = "what_$fields[$index]";
            $$found_fields_ref{what} .= " $field_name";
            $index++;
            $$data_ref{"${field_name}_max"}[$cycle] = $fields[$index];
            $index += 2;
            $dist = 0;
            for ( $i = 1; $i <= $$dim_ref; $i++ ) {
                $dist += $fields[$index]**2;
                $$data_ref{"${field_name}_max_$dim_name[$i]"}[$cycle] = $fields[$index];
              $index++;
            }
            $dist = $dist**.5;
            $$data_ref{"${field_name}_max_dist"}[$cycle] = $dist;
            
            $$data_ref{"${field_name}_min"}[$cycle] = $fields[$index];
            $index += 2;
            $dist = 0;
            for ( $i = 1; $i <= $$dim_ref; $i++ ) {
                $dist += $fields[$index]**2;
                $$data_ref{"${field_name}_min_$dim_name[$i]"}[$cycle] = $fields[$index];
              $index++;
            }
            $dist = $dist**.5;
            $$data_ref{"${field_name}_min_dist"}[$cycle] = $dist;
          }
          $$found_fields_ref{what} =~ s/^\s*//;
        }

      #...integrated state data...
      elsif( $line =~ /Integrated state data/ ||
          $line =~ /tmxd\s+=.*\s+tmass\s+=/ ){
          $name = "Integrated_state_data";
          if( $line =~ /Integrated state data for cycle number:\s*(\d+)\s*at time\s*(\S+)/ ){
              $ln += 1;
              $cycle = $1;
              $time = $2;
          }
          else{
              $ln -= 1;
          }
          $done = "false";
          while( $done eq "false" ){
              $line = $$lines_ref[$ln]; $ln++;
              if( $ln > $ln_max ){
                  $done = "true";
                  last;
              }
              # if a line has any <foo> = <bar>, then not done, otherwise done
              $done = "true";
              while( $line =~ /^(\S+)\s*=\s*(\S+)\s*(.*?)$/ ){
                  $field = $1;
                  $val = $2;
                  $line = $3;
                  # only do some fields requested
                  if( $field =~ /^(ein|eot|te|tie|tke|tpe|tre|tste|txe)$/ ){
                      $$found_fields_ref{"$name"}{$field} = "";
                      $$data_ref{"$name:$field"}[$cycle] = $val;
                  }
                  elsif( $field eq "error" ){
                      $$data_ref{"err"}[$cycle] = $val;
                  }
                  $done = "false";
              }
          }
          $ln--;
      }

      #...process info...
      elsif ( defined($cycle) &&
              $line =~ /^
                      procmon:
                        \s+(\S+)    # 1: field name
                        \s+(\S+)    # 2: min
                        \s+(\S+)    # 3: max
                        \s+(\S+)    # 4: avg
                        (\s+(\S+))? # 5: sum
                        /x ) {
        $$data_ref{"procmon_${1}_min"}[$cycle] = $2;
        $$data_ref{"procmon_${1}_max"}[$cycle] = $3;
        $$data_ref{"procmon_${1}_avg"}[$cycle] = $4;
        if ( defined($5) ) {
          $$data_ref{"procmon_${1}_sum"}[$cycle] = $5;
        }
        $$found_fields_ref{procmon} = "true";
      }

      #...cyc_ info...
      elsif ( $line =~ /^
                        (cyc_cc\/s)\s+     #  1
                        (cyc_sec)\s+       #  2
                        (cyc_cc\/s\/pe)\s+ #  3
                        /x ) {
        $field_1 = $1;
        $field_2 = $2;
        $field_3 = "cc/s/p";
        $line_orig = $$lines_ref[$ln]; $ln++;
        $line = $line_orig;
        if ( $line =~ /^
                       (\S+)\s+
                       (\S+)\s+
                       (\S+)\s+
                       /x ) {
            $$data_ref{$field_3}[$cycle] = $3;
            $$found_fields_ref{$field_3} = "true";
            $$found_fields_ref{ordering_diff}{"cc\/s\/p"} = "";
          }
      }

      #...matinfo...
      elsif ( defined($cycle) &&
              $line =~ /^
                        (mat)\s+    #  1
                        (md01)\s+   #  2
                        (md61)\s+   #  3
                        (mass)\s+   #  4
                        (eng)\s+    #  5
                        (rho)\s+    #  6
                        (sie)\s+    #  7
                        (ske)\s+    #  8
                        (spe)\s+    #  9
                        (sxe)\s+    #  10
                        (re)\s+     #  11
                        (matident)  #  12
                        /x ) {
        $done = "false";
        while ( $done eq "false" ) {
          $line_orig = $$lines_ref[$ln]; $ln++;
          if ( ! defined($line_orig) ) {
            $done = "true";
            next;
          }
          if ( $line_orig =~ /^
                              (\d+)\s+ # 1
                              (\S+)\s+ # 2
                              (\S+)\s+ # 3
                              (\S+)\s+ # 4
                              (\S+)\s+ # 5
                              (\S+)\s+ # 6
                              (\S+)\s+ # 7
                              (\S+)\s+ # 8
                              (\S+)\s+ # 9
                              (\S+)\s+ # 10
                              (\S+)\s+ # 11
                              (\S+.*?) # 12
                              $/x ) {
            undef( %fields );
            $fields{matnum} = $1;
            $fields{mat} = $12;
            $fields{mass} = $4;
            $fields{rho} = $6;
            if ( $fields{rho} != 0 ) {
              $fields{vol} = $fields{mass}/$fields{rho};
            } else {
              $fields{vol} = 0;
            }
            $fields{eng} = $5 * $fields{mass};
            $fields{ie} = $7 * $fields{mass};
            $fields{ke} = $8 * $fields{mass};
            $fields{pe} = $9 * $fields{mass};
            $fields{xe} = $10 * $fields{mass};
            $fields{re} = $11 * $fields{vol};
            $fields{mat} =~ s/\s+/_/g;
            # mat name needs to include matnum so that sorting prints in order
            $mat = sprintf( "%03d.%s", $fields{matnum},$fields{mat} );
            $$found_fields_ref{matinfo_num}{$fields{matnum}} = "$mat";
            # only found if at least one non-0 mass
            if ( $fields{mass} != 0 ) {
              $$found_fields_ref{matinfo}{$mat} = "$fields{matnum}";
            }
            $$data_ref{"matinfo_${mat}_mass"}[$cycle] = $fields{mass};
            $$data_ref{"matinfo_${mat}_rho"}[$cycle] = $fields{rho};
            $$data_ref{"matinfo_${mat}_vol"}[$cycle] = $fields{vol};
            $$data_ref{"matinfo_${mat}_eng"}[$cycle] = $fields{eng};
            $$data_ref{"matinfo_${mat}_ie"}[$cycle] = $fields{ie};
            $$data_ref{"matinfo_${mat}_ke"}[$cycle] = $fields{ke};
            $$data_ref{"matinfo_${mat}_pe"}[$cycle] = $fields{pe};
            $$data_ref{"matinfo_${mat}_xe"}[$cycle] = $fields{xe};
            $$data_ref{"matinfo_${mat}_re"}[$cycle] = $fields{re};
          } else {
            $done = "true";
            next;
          }
        }
      }

      #...mixed material cells
      elsif ( $line =~ /----- cells containing material ----/ ) {
        $done = "false";
        $line_orig = $$lines_ref[$ln]; $ln++;
        while ( $done eq "false" ) {
          $line_orig = $$lines_ref[$ln]; $ln++;
          if ( $ln > $ln_max ||
               !defined($line_orig) ||
               $line_orig !~ /\S/ ) {
            $done = "true";
            last;
          }
          $line = $line_orig;
          if ( $line =~ /^
                         (\d+)\s*:\s* # 1: material number
                         /x ) {
                $$nummat_ref = $1;
              }
        }
      }

      #...mixed material cells
      elsif ( $line =~ /----- cells with num materials -----/ ) {
        $done = "false";
        $line_orig = $$lines_ref[$ln]; $ln++;
        $ncell_all = 0;
        $ncell_top = 0;
        $ncell_all_mixed = 0;
        $ncell_top_mixed = 0;
        $ncell_all_mixed_mats = 0;
        $ncell_top_mixed_mats = 0;
        undef( %field );
        while ( $done eq "false" ) {
          $line_orig = $$lines_ref[$ln]; $ln++;
          if ( $ln > $ln_max ||
               !defined($line_orig) ||
               $line_orig !~ /\S/ ) {
            $done = "true";
            last;
          }
          $line = $line_orig;
          if ( $line =~ /^
                         (\d+)\s*:\s* # 1: number of materials
                         (\d+)\s*     # 2: number of cells total
                         \(\s*        # (
                         (\S+)        # 3: percent of tot
                         \s*[,%]\s*\)\s* # ,)
                         (\d+)\s*     # 4 : number of cells active
                         \(\s*        # (
                         (\S+)        # 3: percent of tot
                         \s*[,%]\s*\)\s* # ,)
                         /x ) {
                $num_mats = $1;
                $nall_num = $2;
                $ntop_num = $4;
                $field{tot}{$1} = $3;
                $field{top}{$1} = $5;
                $ncell_all += $nall_num;
                $ncell_top += $ntop_num;
                if( $num_mats > 1 )
                  {
                    $ncell_all_mixed += $nall_num;
                    $ncell_top_mixed += $ntop_num;
                  }
                $ncell_top_mixed_mats += $ntop_num*$num_mats;
                $ncell_all_mixed_mats += $nall_num*$num_mats;
              }
          }
        if( $ncell_all > 0 )
          {
            $$found_fields_ref{mixed_cells} = "true";
            $$data_ref{NcellAMix}[$cycle] = ($ncell_all_mixed/$ncell_all)*100.0;
            $$data_ref{NcellTMix}[$cycle] = ($ncell_top_mixed/$ncell_top)*100.0;
            $$data_ref{NcellTMixAvg}[$cycle] = $ncell_top_mixed_mats/$ncell_top;
            $$data_ref{NcellAMixAvg}[$cycle] = $ncell_all_mixed_mats/$ncell_all;
            $$data_ref{MatPctTot}[$cycle] = "";
            $$data_ref{MatPctTop}[$cycle] = "";
            for( $i = 1; $i <= $$nummat_ref; $i++ )
              {
                if( defined$field{tot}{$i} )
                  {
                    $$data_ref{MatPctTot}[$cycle] .= "$field{tot}{$i}:";
                    $$data_ref{MatPctTop}[$cycle] .= "$field{top}{$i}:";
                  }
                else
                  {
                    $$data_ref{MatPctTot}[$cycle] .= "-:";
                    $$data_ref{MatPctTop}[$cycle] .= "-:";
                  }
              }
            $$data_ref{MatPctTot}[$cycle] =~ s/:$//;
            $$data_ref{MatPctTop}[$cycle] =~ s/:$//;
          }
      }

      #...Resource
      elsif ( defined($cycle) &&
              $line =~ /Resources: start/ ) {
        while ( $ln <= $ln_max ) {
          $line_orig = $$lines_ref[$ln]; $ln++;
          if ( $line_orig =~ /^Type\s+Min\%\s+Mean\%\s+Max\%$/ ) {
            $$found_fields_ref{"Resources:memory"} = "true";
            while ( $ln <= $ln_max ) {
              $line_orig = $$lines_ref[$ln]; $ln++;
              if ( $line_orig =~ /(Estimate|RSS|Virt|RSS_MAX)\s+
                                  (\S+)\s+(\S+)\s+(\S+)/x ) {
                $field_name = "Resources:memory:$1:";
                $$data_ref{"${field_name}Min%"}[$cycle]  = $2;
                $$data_ref{"${field_name}Mean%"}[$cycle] = $3;
                $$data_ref{"${field_name}Max%"}[$cycle]  = $4;
                $$found_fields_ref{ordering_diff}{"^Resources:"} = "";
              } elsif ( $line_orig !~ /\S/ ) {
                last;
              }
            }
          }
          if ( $line_orig =~ /Resources: end/ ) {
            last;
          }
        }
      }

      #...restart_block...
      elsif( $line =~ /Echo restart block information/ ){
          $done = "false";
          while( $done eq "false" ){
              $line_orig = $$lines_ref[$ln]; $ln++;
              if ( ! defined($line_orig) ) {
                  $done = "true";
                  next;
              }
              if( $line_orig =~ /End of echo restart block information/ ){
                  $done = "true";
                  next;
              }
              if( $line_orig =~ /Echo restart block info, restart block name = (\S+.*)$/ ){
                  my($rs_block) = $1;
                  $rs_block =~ s/\s+/_/g;
                  $line_orig = $$lines_ref[$ln]; $ln++;
                  if( $line_orig =~ /Active flag = (\S+)/ ) {
                      my( $flag ) = $1;
                      if( ! defined($$found_fields_ref{restart_block}) ){
                          $$found_fields_ref{restart_block} = 1;
                      }
                      if( ! defined($$found_fields_ref{"restart_block:$rs_block"} ) ){
                          $$found_fields_ref{"restart_block:$rs_block"} =
                              $$found_fields_ref{restart_block};
                          $$found_fields_ref{restart_block} += 1;
                      }
                      if( $flag eq "true" ){
                          $flag = $$found_fields_ref{"restart_block:$rs_block"};
                      }
                      else{
                          $flag = 0;
                      }
                      $$data_ref{"restart_block:$rs_block"}[$cycle] = $flag;
                  }
              }
          }
      }

      #...other cycle lines: ncycle from io - use only if not defined since
      #...time is not set (restarts_block will use it)
      elsif ( ! defined($cycle) &&
              $line =~ /^ncycle\s*=\s*(\d+)$/ ) {
        undef(%field);
        $field{"cyc"}    = $1;
        $done = "false";
        while( $done eq "false" ){
            $line_orig = $$lines_ref[$ln]; $ln++;
            if ( ! defined($line_orig) ) {
                $done = "true";
                next;
            }
            if( $line_orig =~ /^time\s*=\s*(${GNUMBER_REGEXP})/ ){
                $field{time} = $1;
                $done = "true";
                next;
            }
        }
        if ( $field{"cyc"} =~ /^($GNUMBER_REGEXP)$/ ) {
          $cycle = $field{cyc};
          $time = $field{time};
        }
      }

      # this line also says what the cycle is, but do not know time
      # so cannot make first field
      elsif ( $line =~ /RESIZE: cycle, oldsize,.*=\s*(\d+)/ ){
          $cycle = $1;
      }

      # new cycle and time info
      elsif ( $line =~ /for cycle number:\s*(\d+)\s*at time\s*(${GNUMBER_REGEXP})/ ){
          $cycle = $1;
          $time = $2;
      }

      #...other cycle lines
      elsif ( $line =~ /^cycle\s*=\s*(\d+)\s*,?\s*time\s*=\s*(${GNUMBER_REGEXP})/ ) {
        $cycle = $1;
        $time = $2;
      }

      #...other cycle lines
      elsif ( $line =~ /^RMEDT:\s*(\d+)\s*(${GNUMBER_REGEXP})/ ) {
        $cycle = $1;
        $time = $2;
      }

      #...metric line...
      elsif ( $line =~ /^
                      metric:\s+
                        (contiguity)\s+(\S+)\s+(\S+) # contiguity
                        /x ) {
        $field_1 = $1;
        $field_2 = $2;
        $field_3 = $3;
        $field_name = "metric:$field_1";
        $$found_fields_ref{$field_name} = "true";
        if ( $field_1 eq "contiguity" ) {
          $$data_ref{"$field_name"}[$cycle]      =  $3;
        }
      }

      # crude start at reading other code output file
      elsif ( $line =~ /^Cycle:\s*(\S+)\s*Time:\s*(\S+)\s*dt:\s*(\S+)\s*/ ) {
        undef(%field);
        undef(%field);
        $field{"cyc"}    = $1;
        $field{time}   = $2;
        $field{"dt"}     = $3;
        if ( $field{"cyc"} =~ /^($GNUMBER_REGEXP)$/ ) {
          $cycle = $field{cyc};
          $time = $field{time};
          $$data_ref{dt}[$cycle] = $field{$field_name};
        }
      }

      # ISO: Isotope Summary
      elsif( $line =~ /ISO:\s+Isotope\s+Summary/i ){
          $done = "false";
          while( $done eq "false" ){
              $line_orig = $$lines_ref[$ln]; $ln++;
              if( $line_orig =~ /End\s+ISO:\s+Isotope\s+Summary/i ){
                  $done = "true";
              }
              elsif( $line_orig =~ /Time:\s*(\S+)\s*Cycle:\s*(\d+)/ ){
                  $time = $1;
                  $cycle = $2;
              }
              # Could read in column headers but just assume will not change.
              # If it does, simply add in reader and parse fields
              elsif( $line_orig =~ /^\s*
                                    (isosum)_(\S+)\s+  # 1: TAG
                                    (\S+)\s+           # 3: TIME (already have)
                                    (\S+)\s+           # 4: CYCLE (already have)
                                    (\S+)\s+           # 5: ISOTOPE
                                    (\S+)\s+           # 6: NUMBER
                                    (\S+)\s+           # 7: MOLES
                                    (\S+)\s*           # 8: MASS
                                   /x ){
                  $field_1 = sprintf( "%s_%03d_%s", ${1}, ${2}, ${5} );
                  $$data_ref{"${field_1}_NUMBER"}[$cycle] = $6;
                  $$data_ref{"${field_1}_MOLES"}[$cycle]  = $7;
                  $$data_ref{"${field_1}_MASS"}[$cycle]   = $8;
              }
              if( $ln > $ln_max ){
                  $done = "true";
              }
          }
      }

      # put parse_output_file_extras here if doing internal
      elsif( defined($HAS_EXTRAS) ){
          # parse_output_file_extras
          &parse_output_file_extras( $lines_ref, \$ln, \$cycle, \$time, $first_cycle_ref,
                                     $nummat_ref, $data_ref, $found_fields_ref,
                                     $dim_ref );
      }
      
      # set time if more precision
      if( defined($time) ){
          $time_old = $$data_ref{time}[$cycle];
          if( ! defined($time_old) ||
              length($time) > length($time_old) ){
              $$data_ref{time}[$cycle] = $time;
          }
          undef( $time );
      }

      # set cycle if defined
      if( defined($cycle) ){
          $$data_ref{cyc}[$cycle] = $cycle;
      }

      # if cycle set for first time for this file
      if( defined($cycle) && $$first_cycle_ref eq "true" ){
          $$first_cycle_ref = "false";
          $lost_cycles = $#{$$data_ref{time}} - $cycle;
          if( $lost_cycles < 0 ){
              $lost_cycles = 0;
          }
          $$data_ref{lost_cycles}[$cycle] = $lost_cycles;
          if( $lost_cycles > 0 ) {
              foreach $key ( keys %{$data_ref} ) {
                  $#{$$data_ref{$key}} = $cycle;
              }
          }
      }
    }

    undef( $ln );
    undef( $cycle );
    undef( $time );
    # parse_output_file_extras
    
  }

#...get "delta" stats...
sub parse_output_file_finish
  {
    my(
       $data_ref,
       $found_fields_ref,
       $cmd_ref,
      ) = @_;
    my(
       $cycle,
       $cycle_min,
       $cycle_previous,
       $cycle_start,
       $cycle_stop,
       $num,
       $secs,
       $secs_ideal,
       $secs_start,
      );

    # cycle_stop and cycle_min
    $$data_ref{cycle_min}[0] = 0;
    $cycle_stop = $#{$$data_ref{time}};
    for ( $cycle = 0; $cycle <= $cycle_stop; $cycle++ ) {
        if ( defined($$data_ref{secs}[$cycle]) ) {
            $$data_ref{cycle_min}[0] = $cycle;
            $cycle_min = $cycle;
            last;
        }
    }

    if ( defined($$found_fields_ref{ppp}) || defined($$found_fields_ref{long})) {
      $cycle_previous = 0;
      $secs_start = 0;
      $cycle_min = $$data_ref{cycle_min}[0];
      #...replace YEAR in secs field with actual year
      #...get to first real YEAR value
      for ( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ) {
        $secs = $$data_ref{secs}[$cycle];
        if ( defined($secs) && $secs =~ /(\d{3})$/) {
          $year = $1;
          @time_fields = split( / /, $secs );
          $secs_base = timegm(@time_fields);
          $cycle_base = $cycle;
          last;
        }
      }

      # just use current year if not found
      if( ! defined($year) ){
        $year = (gmtime)[5];
        $secs_base = time();
        $cycle_base = -1;
      }

      #...put total secs into secs field
      $secs_previous = 0;
      for ( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ) {
        $secs = $$data_ref{secs}[$cycle];
        if ( defined($secs) ) {
          @time_fields = split( / /, $secs );
          if( $time_fields[-1] eq "YEAR" ){
            $time_fields[-1] = $year;
            $secs = timegm(@time_fields);
            # must have been previous year, so decrease it
            while( $cycle < $cycle_base && $secs > $secs_base ){
              $time_fields[-1]--;
              $secs = timegm(@time_fields);
            }
            while( $cycle > $cycle_base && $secs < $secs_base ){
              $time_fields[-1]++;
              $secs = timegm(@time_fields);
            }
          }
          else{
            $year = $time_fields[-1];
            $secs_base = timegm(@time_fields);
            $cycle_base = $cycle;
          }
          $secs = timegm(@time_fields);
          $$data_ref{secs}[$cycle] = $secs;
          $secs_previous = $secs;
        }
      }

      #...get to first cycle with data and init
      for ( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ) {
        if ( defined($$data_ref{secs}[$cycle]) ) {
          $secs_start = $$data_ref{secs}[$cycle];
          $$data_ref{secs}[$cycle] = 0;
          $$data_ref{secs_ideal}[$cycle] = 0;
          $cycle_previous = $cycle;
          last;
        }
      }

      #...get to next cycle with data
      $cycle_start = $cycle_previous+1;
      for ( $cycle = $cycle_start; $cycle <= $cycle_stop; $cycle++ ) {
        if ( defined($$data_ref{secs}[$cycle]) ) {
          $cycle_start = $cycle;
          last;
        }
      }

      #...process remaining cycles
      $secs_ideal = 0;
      for ( $cycle = $cycle_start; $cycle <= $cycle_stop; $cycle++ ) {
        if ( ! defined($$data_ref{secs}[$cycle]) ) {
          next;
        }
        #...normalize given start at time 0
        $$data_ref{secs}[$cycle] -= $secs_start;
        #...do not count first cycle after restart
        if ( ! defined($$data_ref{lost_cycles}[$cycle]) ) {
          $secs = $$data_ref{secs}[$cycle] - $$data_ref{secs}[$cycle_previous];
          $secs_ideal += $secs;
          $$data_ref{"secs/cycle"}[$cycle] = $secs /($cycle-$cycle_previous);
          $$data_ref{"secs_ideal"}[$cycle] = $secs_ideal;
        }
        $cycle_previous = $cycle;
      }

      # secs_dmp
      if( defined($$data_ref{dmp_write_time}) ){
        $secs = 0;
        $num = 0;
        for ( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ) {
          $val = $$data_ref{dmp_write_time}[$cycle];
          if ( ! defined($val) ) {
            next;
          }
          $secs += $val;
          $num++;
          $$data_ref{secs_dmp}[$cycle] = $secs;
          if( ! defined($$data_ref{secs_io_write}[$cycle]) ){
              $$data_ref{secs_io_write}[$cycle] = 0;
          }
          $$data_ref{secs_io_write}[$cycle] += $val;
        }
        $$data_ref{dmp_write_num}[0] = $num;
        if( $num > 0 ){
            $$data_ref{"secs/dmp_write"}[0] = $secs/$num;
        }
        else{
            $$data_ref{"secs/dmp_write"}[0] = 0;
        }
      }

      # secs_hdf
      if( defined($$data_ref{hdf_write_time}) ){
        $secs = 0;
        for ( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ) {
          $val = $$data_ref{hdf_write_time}[$cycle];
          if ( ! defined($val) ) {
            next;
          }
          $secs += $val;
          $$data_ref{secs_hdf}[$cycle] = $secs;
          if( ! defined($$data_ref{secs_io_write}[$cycle]) ){
              $$data_ref{secs_io_write}[$cycle] = 0;
          }
          $$data_ref{secs_io_write}[$cycle] += $val;
        }
      }
      $secs = 0;
      for ( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ) {
          $val = $$data_ref{secs_io_write}[$cycle];
          if ( ! defined($val) ) {
              next;
          }
          $secs += $val;
          $$data_ref{secs_io_write}[$cycle] = $secs;
      }
      # for printing last
      $$found_fields_ref{ordering_diff}{"^secs"} = "";
      $$found_fields_ref{ordering_diff}{"^lost_cycles"} = "";
    }

    #...fill in last value
    my( $i, $last, $restart_block, @restart_blocks );
    @restart_blocks = sort grep( /^restart_block:/, keys %{$data_ref} );
    foreach $restart_block ( @restart_blocks ){
        $last = $$data_ref{$restart_block}[-1];
        $$data_ref{$restart_block}[$cycle_stop] = $last;
    }

    #...other finishing
    if( defined( $HAS_EXTRAS ) ){
        &parse_output_file_finish_extras( $data_ref, $found_fields_ref );
    }

    #...normalize time if defined
    if( ! defined($$cmd_ref{no_tshift}) && defined($$data_ref{time_start}) ){
      for( $cycle = $cycle_min; $cycle <= $cycle_stop; $cycle++ ){
        if( defined($$data_ref{time}[$cycle]) ){
          $$data_ref{time}[$cycle] -= $$data_ref{time_start};
        }
      }
    }
  }
1;

