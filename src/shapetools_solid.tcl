package provide shapetools_liquid 0.2

namespace eval shapetools_liquid {
    namespace export \
    zero_ampl \
    rectangular \
    wasp \
    wurst \
    cawurst \
    swurst \
    supergaussian \
    tanhpulse \
    hspulse \
    hahn_echo \
    shape_add \
    shape_sub \
    shape_mult \
    shape_integr \
    shape_loadbruker \
    shape_freqprofile


    ###########################################################################
    # Proc for empty, zero shape
    # dur and stepsize are in us
    # Version 1.0 Dec 2020 by MRH
    ###########################################################################
    proc zero_ampl {dur {stepsize 0.05}} {
        global par
        set i 1
        set nsteps [expr round($dur/$stepsize)]
        set checkduration [expr $nsteps*$stepsize]
        if {$checkduration != $dur } {
            puts "Error: zero_ampl, duration ($dur us) cannot be digitized with the selected stepsize ($stepsize us)"
            exit
        } else {
            while {$i <= $nsteps} {
                lappend wavelist [format "0.0 0.0"]
                incr i
            }
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for RECTANGULAR shape calculation
    # 
    # Changed 09.12.2020, MRH:
    #   - Added option for offset
    # 09.07.2019 by Max Busskamp
    ###########################################################################
    proc rectangular {dur rfmax {phase 0.0} {offset 0.0} {stepsize 0.05}} {
        global par

        set nsteps [expr round($dur/$stepsize)]
        
        for {set i 1} {$i <= $nsteps} {incr i} {
            set t [expr $stepsize*$i]
            set ampl $rfmax
            set ph [expr fmod(360.0e-6*$offset*$t+$phase,360)]
            set ph [expr fmod($ph+360.0,360)]
            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for WASP shape calculation
    #
    # 20.01.2020 by Max Busskamp:
    ###########################################################################
    proc wasp {dur rfmax N {phase 0.0} {offset 0.0} {stepsize 0.05}} {
        global par

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$stepsize)]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set t [expr $stepsize*$i]
            set ampl [expr $rfmax*(1.0-abs(pow(cos(($pi*$t)/($dur)),$N)))]
            set ph [expr fmod(360.0e-6*$offset*$t+$phase,360)]
            set ph [expr fmod($ph+360.0,360)]
            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for WURST shape calculation
	#
    # Changed 20.01.2020 by Max Busskamp:
    #   - Added Option for Phasecycle
    #   - Added rfmax to input variables
    #   - Added default values for stepsize, direction and offset
    ###########################################################################
    proc wurst {dur rfmax Delta N {direction 0} {phase 0.0} {offset 0.0} {stepsize 0.05}} {
        global par

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$stepsize)]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set ampl [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($dur)),$N)))]
            set ph [expr ((180.0/$pi)*2.0*$pi*(($offset*1e3+($Delta*1e3/2.0))*$stepsize*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($stepsize*1e-6*$i,2)))+$phase]
            if {$direction} {
                set ph [expr fmod($ph,360)]
            } else {
                set ph [expr fmod(-$ph,360)+360.0]
            }
            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for caWURST shape calculation 
	# Version 1.0 by JK
	#
    # Changed 16.08.2020 by Max Busskamp:
	#   - Fixed changin sweepwidth
	#   - Checked Direction
	#   - Added constant phase offset to enable phasecycling
    #   - Added rfmax to input variables
    #   - Added default values for stepsize, direction and offset
    ###########################################################################
    proc cawurst {dur rfmax Delta N {direction -1} {phase 0.0} {offset 0.0} {stepsize 0.05}} {
        global par

        #Variables
        set pi2 		[expr (4*atan(1))/2]

        # sweep / MHz
        set sweep		[expr $Delta/1000.0]

        set rate 		[expr 1.0e6*($sweep/$dur)]

        #On-resonance q-factor
        set q0			[expr (pow($rfmax,2))/(1e6*$rate)]

        set nsp			[expr round($dur/$stepsize)]

        #ampllitude
        set b 			[expr (2*$pi2)/$dur]

        #Frequency
        set w0 			0
        set k 			[expr 1.0e-6*$rate]

        #Calculate frequency correction factor
        set freq_prep_b 0.0
        set phase_prep1	0.0
        set freq_max 0.0
        for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
            set time				[expr $stepsize*$i]
            set ampl_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
            lappend ampllist	$ampl_prep
            
            set freq_prep 	[expr (1.0e-6*$stepsize*($ampl_prep**2/$q0)+$freq_prep_b)]
            lappend freqlist $freq_prep
            if {$freq_prep > $freq_max} {
                set freq_max $freq_prep
            }

            set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$dur)))]
            set phase_prep 	[expr fmod($phase_prep1*360+$phase,360)]
            lappend phaselist $phase_prep
            
            set freq_prep_b $freq_prep
        }

        set sweepfactor [expr $sweep*1e6/2.0/$freq_max]

        unset time
        unset ampl_prep
        unset ampllist
        unset freq_prep
        unset freqlist
        unset phase_prep1
        unset phase_prep
        unset phaselist
        unset freq_prep_b
        
        #Calculate ampllitude and frequency for positive offset values 
        set freq_prep_b 0.0
        set phase_prep1	0.0
        for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
            set time				[expr $stepsize*$i]
            set ampl_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
            lappend ampllist	$ampl_prep
            
            set freq_prep 	[expr (1.0e-6*$stepsize*$sweepfactor*($ampl_prep**2/$q0)+$freq_prep_b)]
            lappend freqlist $freq_prep

            set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$dur)))]
            set phase_prep 	[expr fmod($phase_prep1*360+$phase,360)]
            lappend phaselist $phase_prep
            
            set freq_prep_b $freq_prep
        }

        set freqlist_length [llength $freqlist]

        if {$direction == 1} {
            set direc 1.0
        } elseif {$direction == -1} {
            set direc -1.0
        } else {
            puts "direction error"
        }

        for {set j 1} {$j < $nsp} {incr j} {
            set index [expr $freqlist_length-$j]
            if {$index < 0} {
                set index [expr abs($index)]
                set help $direc
            } else {
                set help [expr -1.0*$direc]
            }
            set time  		[expr $stepsize*$j]
            set ampl			[lindex $ampllist $index]
            set freq_help 	[lindex $freqlist $index]
            set freq 		[expr $help*$freq_help]
            set phase		[lindex $phaselist $index]

            lappend wavelist [format "%6.2f %6.2f" $ampl $phase]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for saltire WURST shape calculation
	#
    # Changed 19.11.2020 by Max Busskamp:
    #   - Initial Version
    ###########################################################################
    proc swurst {dur rfmax Delta N {direction 0} {phase 0.0} {offset 0.0} {stepsize 0.05}} {
        global par

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$stepsize)]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set t [expr $stepsize*$i]
            set ampl [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($dur)),$N)))]
            set ampl [expr $ampl*cos(2.0e-6*$pi*(-1.0*$Delta*1e3/2.0)*$t-(2.0*$pi*(-1.0*$Delta*1e3/2.0)-2.0*$pi*(1.0*$Delta*1e3/2.0))*($t*$t/(2.0e6*$dur)))]
            
            if {$direction == 0} {
                if {$ampl < 0.0} {
                    set ampl [expr abs($ampl)]
                    set ph [expr 0.0+$phase]
                    set ph [expr fmod($ph,360)]
                } else {
                    set ph [expr 180.0+$phase]
                    set ph [expr fmod($ph,360)]
                }
            } else {
                if {$ampl < 0.0} {
                    set ampl [expr abs($ampl)]
                    set ph [expr 180.0+$phase]
                    set ph [expr fmod($ph,360)]
                } else {
                    set ph [expr 0.0+$phase]
                    set ph [expr fmod($ph,360)]
                }
            }
            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for supergaussian shape calculation
    # Version 1 16.09.2020 by Max Busskamp:
    ###########################################################################
    proc supergaussian {dur rfmax Delta N G {direction 1} {phase 0.0} {offset 0.0} {stepsize 0.05}} {
        global par

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$stepsize)]

        # t = $stepsize*$i
        for {set i 1} {$i <= $nsteps} {incr i} {

            set ampl [expr $rfmax*exp(-1.0*pow(2,($N+2))*pow(((($stepsize*$i)-$G)/($dur)),$N))]

            set ph [expr ((180.0/$pi)*2.0*$pi*(($offset*1e3+($Delta*1e3/2.0))*$stepsize*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($stepsize*1e-6*$i,2)))+$phase]
            if {$direction} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                set ph [expr fmod($ph,360)]
            }
            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for tanh/tan shape calculation
    # Version 1.0 Max Busskampl 21.09.2019
    ###########################################################################
    proc tanhpulse {dur rfmax Delta zeta tan_kappa {direction 0} {phase 0.0} {stepsize 0.05}} {
        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr int(round($dur/$stepsize))]
        set kappa [expr atan($tan_kappa)]
        set R [expr $Delta*1000*$dur/1000000.0]
        set A [expr ($R*$pi)/($dur/1000000.0)]

        set phi_max [expr -(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

        for {set i 1} {$i <= $nsteps} {incr i} {
            if {$i <= $nsteps/2} {
                set ampl [expr $rfmax*tanh((2*$i*$stepsize*$zeta)/$dur)]
            } else {
                set ampl [expr $rfmax*tanh((2*$zeta*(1.0-(($i*$stepsize)/$dur))))]
            }
            set ph [expr ($phi_max-(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$stepsize/1000000.0)/($dur/1000000.0))/cos($kappa)))*180/$pi+$phase]
            if {$direction} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                set ph [expr fmod($ph,360)]
            }

            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for HS shape calculation
    # Version 1.1 MRH Sept 2016
    ###########################################################################
    proc hspulse {dur rfmax Delta beta {direction 0} {phase 0.0} {offset 0.0} {stepsize 0.05}} {

        set nsteps [expr round($dur/$stepsize)]
        set phi0 [expr 180.0*$Delta*1000*$dur/10.6e6]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set x [expr cosh($beta*(2*$i*$stepsize/$dur-1))]
            set ampl [expr $rfmax/$x]
            set ph [expr $offset+$phi0*log($x)+$phase]

            if {$direction} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                set ph [expr fmod($ph,360)]
            }

            lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for a simple Hahn echo; dur is total duration (t-p180-t) with p180 centered; phase and offset applies to p180.
    # dur and stepsize are in us; offset in Hz; rfmax in kHz
    # Version 1.0 Dec 2020 by MRH
    # Version 1.1 Nov 2021 by JB fixed rf power (Hz instead of kHz)
    ###########################################################################
    proc hahn_echo {dur rfmax {phase 0.0} {stepsize 0.05} {offset 0.0} {verb 0}} {
        global par
        set rfmax [expr $rfmax*1e3]
        set i 1
        set j 1
        set nsteps [expr round($dur/$stepsize)]
        set p180 [expr 1e6/(2.0*$rfmax)]
        set p180steps [expr round($p180/$stepsize)]
        set delsteps [expr ($nsteps-$p180steps)/2]
        set checkduration [expr (2*$delsteps+$p180steps)*$stepsize]

        # check if total duration is possible with the selected rfmax
        if {$checkduration != $dur} {
            puts "Error: hahn_echo, p180 ($p180 us) cannot centered "
            puts "be digitized with the chosen stepsize ($stepsize us)"
            puts "Please adjust dur ($dur us), stepsize ($stepsize us) or rfmax ($rfmax Hz)!"
            puts "check: $checkduration us"
            exit
        }

        # first half echo
        while {$j <= $delsteps} {
            set ampl 0.0
            set ph 0.0
            lappend wavelist [format "%.6f %.6f" $ampl $ph]
            incr j
        }
        set j 1

        # make 180 pulse
        while {$i <= $p180steps} {
            set t [expr ($delsteps+$i)*$stepsize]
            set ampl $rfmax
            set ph [expr fmod(360.0e-6*$offset*$t+$phase,360)]
            set ph [expr fmod($ph+360.0,360)]
            lappend wavelist [format "%.6f %.6f" $ampl $ph]
            incr i
        }

        # second half echo
        while {$j <= $delsteps} {
            set ampl 0.0
            set ph 0.0
            lappend wavelist [format "%.6f %.6f" $ampl $ph]
            incr j
        }
        return $wavelist
    }


    ###########################################################################
    # Proc to bidurise add duro shape lists of the same length
	#
    # Changed 01.12.2020 by Max Busskamp:
    #   - Version 1.0
    #   - Version 1.1 Added safety check for list length
    #   - Version 1.2 Fixed Phase Angle calculations
    ###########################################################################
    proc shape_add {shape_list1 shape_list2 {return_complex 0}} {
        global par

        set rf_shape1_list []
        set rf_shape2_list []

        foreach pulse $shape_list1 {
            set rf_shape1_list [lappend rf_shape1_list {*}$pulse]
        }

        foreach pulse $shape_list2 {
            set rf_shape2_list [lappend rf_shape2_list {*}$pulse]
        }

        set pi [expr 4.0*atan(1.0)]

        set list_len1 [llength $rf_shape1_list]
        set list_len2 [llength $rf_shape2_list]

        if { $list_len1 == $list_len2 } {
            for {set i 0} {$i < $list_len1} {incr i} {
                set phase1 [lindex $rf_shape1_list $i 1]
                set phase2 [lindex $rf_shape2_list $i 1]
                set ampll1  [lindex $rf_shape1_list $i 0]
                set ampll2  [lindex $rf_shape2_list $i 0]

                set compl1_real [expr $ampll1*cos($phase1*2.0*$pi/360.0)]
                set compl1_imag [expr $ampll1*sin($phase1*2.0*$pi/360.0)]
                set compl2_real [expr $ampll2*cos($phase2*2.0*$pi/360.0)]
                set compl2_imag [expr $ampll2*sin($phase2*2.0*$pi/360.0)]

                set compl_real_sum [expr $compl1_real + $compl2_real]
                set compl_imag_sum [expr $compl1_imag + $compl2_imag]

                if {$return_complex == 0} {
                    set ampl [expr hypot($compl_real_sum, $compl_imag_sum)]
                    # We added 360degrees to get a correct WASP*MODLA Pulse. Please Check here if something seems wrong!
                    set ph  [expr atan2($compl_imag_sum, $compl_real_sum)/$pi*180.0+360.0]
                    set ph [expr fmod($ph,360)]

                    lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
                } else {
                    lappend wavelist [format "%6.2f %6.2f" $compl_real_sum $compl_imag_sum]
                }
            }
            return $wavelist
        } else {
            puts "ERROR: Shapes have not the same size!"
        }
    }


    ###########################################################################
    # Proc to bidurise subtract duro shape lists of the same length
	#
    # Changed 01.12.2020 by Max Busskamp:
    #   - Version 1.0
    #   - Version 1.1 Added safety check for list length
    ###########################################################################
    proc shape_sub {shape_list1 shape_list2 {return_complex 0}} {
        global par

        set rf_shape1_list []
        set rf_shape2_list []

        foreach pulse $shape_list1 {
            set rf_shape1_list [lappend rf_shape1_list {*}$pulse]
        }

        foreach pulse $shape_list2 {
            set rf_shape2_list [lappend rf_shape2_list {*}$pulse]
        }

        set pi [expr 4.0*atan(1.0)]

        set list_len1 [llength $rf_shape1_list]
        set list_len2 [llength $rf_shape2_list]

        if { $list_len1 == $list_len2 } {
            for {set i 0} {$i < $list_len1} {incr i} {
                set phase1 [lindex $rf_shape1_list $i 1]
                set phase2 [lindex $rf_shape2_list $i 1]
                set ampll1  [lindex $rf_shape1_list $i 0]
                set ampll2  [lindex $rf_shape2_list $i 0]

                set compl1_real [expr $ampll1*cos($phase1*2.0*$pi/360.0)]
                set compl1_imag [expr $ampll1*sin($phase1*2.0*$pi/360.0)]
                set compl2_real [expr $ampll2*cos($phase2*2.0*$pi/360.0)]
                set compl2_imag [expr $ampll2*sin($phase2*2.0*$pi/360.0)]

                set compl_real_sub [expr $compl1_real - $compl2_real]
                set compl_imag_sub [expr $compl1_imag - $compl2_imag]

                if {$return_complex == 0} {
                    set ampl [expr hypot($compl_real_sub, $compl_imag_sub)]
                    # We added 360degrees to get a correct WASP*MODLA Pulse. Please Check here if something seems wrong!
                    set ph  [expr atan2($compl_imag_sub, $compl_real_sub)/$pi*180.0+360.0]
                    set ph [expr fmod($ph,360)]

                    lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
                } else {
                    lappend wavelist [format "%6.2f %6.2f" $compl_real_sum $compl_imag_sum]
                }
            }
            return $wavelist
        } else {
            puts "ERROR: Shapes have not the same size!"
        }
    }


    ###########################################################################
    # Proc to bidurise multiplicate duro shape lists of the same length
	#
    # Changed 04.12.2020 by Max Busskamp:
    #   - Version 1.0
    ###########################################################################
    proc shape_mult {shape_list1 shape_list2 {rfmax 1.0} {return_complex 0}} {
        global par

        set rf_shape1_list []
        set rf_shape2_list []

        foreach pulse $shape_list1 {
            set rf_shape1_list [lappend rf_shape1_list {*}$pulse]
        }

        foreach pulse $shape_list2 {
            set rf_shape2_list [lappend rf_shape2_list {*}$pulse]
        }

        set pi [expr 4.0*atan(1.0)]

        set list_len1 [llength $rf_shape1_list]
        set list_len2 [llength $rf_shape2_list]

        set max1 0.0
        set max2 0.0

        for {set i 0} {$i < $list_len1} {incr i} {
            set max1_temp [tcl::mathfunc::max {*}[lindex $rf_shape1_list $i 0]]
            set max2_temp [tcl::mathfunc::max {*}[lindex $rf_shape2_list $i 0]]

            if { $max1_temp > $max1 } {
                set max1 $max1_temp
            }
            if { $max2_temp > $max2 } {
                set max2 $max2_temp
            }
        }

        if { $list_len1 == $list_len2 } {
            for {set i 0} {$i < $list_len1} {incr i} {
                set phase1 [lindex $rf_shape1_list $i 1]
                set phase2 [lindex $rf_shape2_list $i 1]
                set ampll1  [lindex $rf_shape1_list $i 0]
                set ampll2  [lindex $rf_shape2_list $i 0]

                set compl1_real [expr $ampll1*cos($phase1*2.0*$pi/360.0)]
                set compl1_imag [expr $ampll1*sin($phase1*2.0*$pi/360.0)]
                set compl2_real [expr $ampll2*cos($phase2*2.0*$pi/360.0)]
                set compl2_imag [expr $ampll2*sin($phase2*2.0*$pi/360.0)]

                set compl_real_mult [expr $compl1_real*$compl2_real - $compl1_imag*$compl2_imag]
                set compl_imag_mult [expr $compl1_real*$compl2_imag + $compl1_imag*$compl2_real]

                if {$return_complex == 0} {
                    set ampl [expr hypot($compl_real_mult, $compl_imag_mult)/($max1*$max2)*$rfmax]
                    # We added 360degrees to get a correct WASP*MODLA Pulse. Please Check here if something seems wrong!
                    set ph  [expr atan2($compl_imag_mult, $compl_real_mult)/$pi*180.0+360.0]
                    set ph [expr fmod($ph,360)]

                    lappend wavelist [format "%6.2f %6.2f" $ampl $ph]
                } else {
                    lappend wavelist [format "%6.2f %6.2f" $compl_real_sum $compl_imag_sum]
                }
            }
            return $wavelist
        } else {
            puts "ERROR: Shapes have not the same size!"
        }
    }


    ###########################################################################
    # Returns integral of given shape
	#
    # Changed 01.12.2020 by Max Busskamp:
    #   - Version 1.0
    ###########################################################################
    proc shape_integr {rf_shape1 {stepsize 0.05}} {
        global par

        set rf_shape1_list [shape2list $rf_shape1]
        set list_len1 [llength $rf_shape1_list]
        set rf_shape1_integr 0.0
        set integr_norm 0.0
        for {set i 0} {$i < $list_len1} {incr i} {
            set ampll1  [lindex $rf_shape1_list $i 0]
            if {$ampll1 > $integr_norm} {
                set integr_norm $ampll1
            }
        }

        for {set i 0} {$i < $list_len1} {incr i} {
            set ampll1  [lindex $rf_shape1_list $i 0]
            set rf_shape1_integr [expr $rf_shape1_integr + $ampll1*$stepsize/$integr_norm]
        }
        return $rf_shape1_integr
    }


    ###########################################################################
    # Proc to load Bruker shapefiles into Simpson
    # Version 1.0, MRH Dec 2020
    ###########################################################################
    proc shape_loadbruker {file {rfmax 100.0}} {
        global par
        set np 1e9

        if {[file exists $file]} {
            set fp [open $file r]
        } else {
            puts "Error: $file does not exist!"
            exit
        }
        set i 0
        while {[eof $fp] != 1} {
            gets $fp buf
            # read number of shape points 
            if {[lindex $buf 0]== "\#\#NPOINTS="} {
                set np [lindex $buf 1]
            }
            # Bruker shape files uses "," as a delimiter bedureen ampllitude and phase values
            if {[string index [lindex $buf 0] end] == ","} {
            set ampl [expr ($rfmax/100.00)*[string range [lindex $buf 0] 0 end-1]]
            set phi [lindex $buf 1]
            if {![string is double -strict $ampl] || ![string is double -strict $phi]} {
                    puts "Error: ampllitude or phase is not a number (check line $i of shape data)!"
                    exit
                }
                incr i
                lappend wavelist [list $ampl $phi]
            }
            # exit while loop if the number of read shape points ($i) is equal to the number of points ($np)
            if {$i == $np} {
                break
            }
        }
        return $wavelist
    }


    ###########################################################################
    # Proc to calculate the frequency profile of shapefiles
    # Version 1.0, MRH Dec 2020
    ###########################################################################
    proc shape_freqprofile {rf_shape1 {name "freq_rf_shape1"} {stepsize 0.05} {rfmax 100.0}} {
        global par

        set pi [expr 4.0*atan(1.0)]

        set rf_shape1_list [shape2list $rf_shape1]
        set list_len1 [llength $rf_shape1_list]
        set dur [expr $list_len1*$stepsize]

        set np [expr int(pow(2,ceil(log($dur/$stepsize)/log(2.0))))]
        set f [fcreate -type fid -np $np -sw [expr 1.0e6/$stepsize] -ref 0]
        fzero $f


        set i 1
        for {set i 0} {$i < $list_len1} {incr i} {
            set ampl [lindex $rf_shape1_list $i 0]
            set ph [lindex $rf_shape1_list $i 1]
            fsetindex $f [expr $i+1] [expr $ampl*cos($ph*$pi/180.0)] [expr $ampl*sin($ph*$pi/180.0)]
            incr i
        }
        # fsave $f ${name}_profile.shape -xreim
        fft $f
        # convert to magnitude data
        fexpr $f [list sqrt(pow(\$re,2)+pow(\$im,2))] 0.0

        # and save data
        fsave $f ${name}_profile.spe -binary
        fsave $f ${name}_profile.dat -xreim
        funload $f
        # puts "Pulse shape saved to ${name}_profile.shape"
        puts "Inversion profile saved to: ${name}_profile.spe and .dat"
    }
}
