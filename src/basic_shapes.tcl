package provide basic_shapes 1.0

namespace eval basic_shapes {

    namespace export \
    zero_ampl \
    rectangular \
    gaussian \
    sinc \
    wurst \
    cawurst \
    supergaussian \
    tanhtan \
    hs \
    hahn_echo \
    duration_check \
    shape_add


    ###########################################################################
    # Proc for empty, zero shape
    # dur and stepsize are in us
    # Version 1.0 Dec 2020 by MRH
    ###########################################################################
    proc zero_ampl {dur args} {
        array set options { -stepsize 0.05 -verbose 0 }
        array set options $args

        set i 1
        duration_check $dur $options(-stepsize) $options(-verbose)
        set nsteps [expr round($dur/$options(-stepsize))]
        set checkduration [expr $nsteps*$options(-stepsize)]
        if {$checkduration != $dur } {
            puts "Error: zero_ampl, duration ($dur us) cannot be digitized with the selected stepsize ($options(-stepsize) us)"
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
    proc rectangular {dur rfmax args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        set nsteps [expr round($dur/$options(-stepsize))]
        
        for {set i 1} {$i <= $nsteps} {incr i} {
            set t [expr $options(-stepsize)*$i]
            set ampl $rfmax
            set ph [expr fmod(360.0e-6*$options(-offset)*$t+$options(-phase),360)]
            set ph [expr fmod($ph+360.0,360)]
            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
        }
        return $wavelist
    }

    ###########################################################################
    # Proc for GAUSSIAN amplitude shape calculation
    # dur is in us, rfmax in Hz
    # fhwm (full width half maximum) in us
    # Version 1.0 Mar 2022 by MRH
    ###########################################################################
    proc gaussian {dur rfmax args} {
        set fwhm [expr $dur/2.0]
        array set options { -fwhm $fwhm -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)
        # Nedeed: option to specify the truncation level and recalculate fwhm
        set nsteps [expr round($dur/$options(-stepsize))]
        set mid [expr ($nsteps/2)*$options(-stepsize)]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set t [expr $i*$options(-stepsize)]
            set ampl [expr $rfmax*exp(-(4.0*log(2.0)*pow(($t-$mid),2)/(pow(-$options(-fwhm),2))))]
            set ph [expr fmod(360.0e-6*$options(-offset)*$t+$options(-phase),360)]
            set ph [expr fmod($ph+360.0,360)]
            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
        }
        return $wavelist
    }

    ###########################################################################
    # Proc for SINC amplitude shape calculation
    # dur is in us, rfmax in Hz
    # N is the number of zero crossings (sinc3 or sinc5)
    # Version 1.0 Mar 2022 by MRH
    ###########################################################################
    proc sinc {dur rfmax N args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)
        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$options(-stepsize))]
        set mid [expr ($nsteps/2)*$options(-stepsize)]
        
        for {set i 1} {$i <= $nsteps} {incr i} {
            set t [expr $i*$options(-stepsize)]
            if {$t == $mid} {
                set ampl $rfmax
            } else {
                set x [expr 2.0*$pi*($t-$mid)/$N]
                set ampl [expr $rfmax*sin($x)/$x]
            }
            if {$ampl < 0} {
                #set ampl [expr abs($ampl)]
                set phadd 180.0
            } else {
                set phadd 0.0
            }
            set ph [expr fmod(360.0e-6*$options(-offset)*$t+$options(-phase)+$phadd,360)]
            set ph [expr fmod($ph+360.0,360)]
            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
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
    proc wurst {dur rfmax Delta N args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$options(-stepsize))]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set ampl [expr $rfmax*(1.0-abs(pow(cos(($pi*$options(-stepsize)*$i)/($dur)),$N)))]
            set ph [expr ((180.0/$pi)*2.0*$pi*(($options(-offset)*1e3+($Delta*1e3/2.0))*$options(-stepsize)*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($options(-stepsize)*1e-6*$i,2)))+$options(-phase)]
            if {$options(-direction) == 1} {
                set ph [expr fmod($ph,360)]
            } elseif {$options(-direction) == 0} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                puts stderr "Direction has to be 1 or 0."
                exit 2
            }
            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
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
    proc cawurst {dur rfmax Delta N args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        #Variables
        set pi2 		[expr (4*atan(1))/2]

        # Sweep / MHz
        set sweep		[expr $Delta/1000.0]

        set rate 		[expr 1.0e6*($sweep/$dur)]

        #On-resonance q-factor
        set q0			[expr (pow($rfmax,2))/(1e6*$rate)]

        set nsp			[expr round($dur/$options(-stepsize))]

        #Amplitude
        set b 			[expr (2*$pi2)/$dur]

        #Frequency
        set w0 			0
        set k 			[expr 1.0e-6*$rate]

        #Calculate frequency correction factor
        set freq_prep_b 0.0
        set phase_prep1	0.0
        set freq_max 0.0
        for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
            set time				[expr $options(-stepsize)*$i]
            set ampl_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
            lappend ampllist	$ampl_prep
            
            set freq_prep 	[expr (1.0e-6*$options(-stepsize)*($ampl_prep**2/$q0)+$freq_prep_b)]
            lappend freqlist $freq_prep
            if {$freq_prep > $freq_max} {
                set freq_max $freq_prep
            }

            set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$dur)))]
            set phase_prep 	[expr fmod($phase_prep1*360+$options(-phase),360)]
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
        
        set freq_prep_b 0.0
        set phase_prep1	0.0
        for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
            set time				[expr $options(-stepsize)*$i]
            set ampl_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
            lappend ampllist	$ampl_prep
            
            set freq_prep 	[expr (1.0e-6*$options(-stepsize)*$sweepfactor*($ampl_prep**2/$q0)+$freq_prep_b)]
            lappend freqlist $freq_prep

            set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$dur)))]
            set phase_prep 	[expr fmod($phase_prep1*360+$options(-phase),360)]
            lappend phaselist $phase_prep
            
            set freq_prep_b $freq_prep
        }

        set freqlist_length [llength $freqlist]

        if {$options(-direction) == 1} {
            set direc 1.0
        } elseif {$options(-direction) == 0} {
            set direc -1.0
        } else {
            puts stderr "Direction has to be 1 or 0."
            exit 2
        }

        for {set j 1} {$j < $nsp} {incr j} {
            set index [expr $freqlist_length-$j]
            if {$index < 0} {
                set index [expr abs($index)]
                set help $direc
            } else {
                set help [expr -1.0*$direc]
            }
            set time  		[expr $options(-stepsize)*$j]
            set ampl		[lindex $ampllist $index]
            set freq_help 	[lindex $freqlist $index]
            set freq 		[expr $help*$freq_help]
            set phase		[lindex $phaselist $index]

            lappend wavelist [format "%6.4f %6.4f" $ampl $phase]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for supergaussian shape calculation
    # Version 1 16.09.2020 by Max Busskamp:
    ###########################################################################
    proc supergaussian {dur rfmax Delta N G args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)
        
        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$options(-stepsize))]

        # t = $options(-stepsize)*$i
        for {set i 1} {$i <= $nsteps} {incr i} {

            set ampl [expr $rfmax*exp(-1.0*pow(2,($N+2))*pow(((($options(-stepsize)*$i)-$G)/($dur)),$N))]

            set ph [expr ((180.0/$pi)*2.0*$pi*(($options(-offset)*1e3+($Delta*1e3/2.0))*$options(-stepsize)*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($options(-stepsize)*1e-6*$i,2)))+$options(-phase)]
            if {$options(-direction) == 1} {
                set ph [expr fmod($ph,360)]
            } elseif {$options(-direction) == 0} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                puts stderr "Direction has to be 1 or 0."
                exit 2
            }
            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for tanh/tan shape calculation
    # Version 1.0 Max Busskampl 21.09.2019
    ###########################################################################
    proc tanhtan {dur rfmax Delta zeta tan_kappa args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr int(round($dur/$options(-stepsize)))]
        set kappa [expr atan($tan_kappa)]
        set R [expr $Delta*1000*$dur/1000000.0]
        set A [expr ($R*$pi)/($dur/1000000.0)]

        set phi_max [expr -(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

        for {set i 1} {$i <= $nsteps} {incr i} {
            if {$i <= $nsteps/2} {
                set ampl [expr $rfmax*tanh((2*$i*$options(-stepsize)*$zeta)/$dur)]
            } else {
                set ampl [expr $rfmax*tanh((2*$zeta*(1.0-(($i*$options(-stepsize))/$dur))))]
            }
            set ph [expr ($phi_max-(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$options(-stepsize)/1000000.0)/($dur/1000000.0))/cos($kappa)))*180/$pi+$options(-phase)]
            if {$options(-direction) == 1} {
                set ph [expr fmod($ph,360)]
            } elseif {$options(-direction) == 0} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                puts stderr "Direction has to be 1 or 0."
                exit 2
            }

            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for HS shape calculation
    # Version 1.1 MRH Sept 2016
    ###########################################################################
    proc hs {dur rfmax Delta beta args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        set nsteps [expr round($dur/$options(-stepsize))]
        set phi0 [expr 180.0*$Delta*1000*$dur/10.6e6]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set x [expr cosh($beta*(2*$i*$options(-stepsize)/$dur-1))]
            set ampl [expr $rfmax/$x]
            set ph [expr $options(-offset)+$phi0*log($x)+$options(-phase)]

            if {$options(-direction) == 1} {
                set ph [expr fmod($ph,360)]
            } elseif {$options(-direction) == 0} {
                set ph [expr fmod(-$ph,360)+360.0]
            } else {
                puts stderr "Direction has to be 1 or 0."
                exit 2
            }

            lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for a simple Hahn echo; dur is total duration (t-p180-t) with p180 centered; phase and offset applies to p180.
    # dur and stepsize are in us; offset in Hz; rfmax in kHz
    # Version 1.0 Dec 2020 by MRH
    # Version 1.1 Nov 2021 by JB fixed rf power (Hz instead of kHz)
    ###########################################################################
    proc hahn_echo {dur rfmax args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)
        
        set i 1
        set j 1
        set nsteps [expr round($dur/$options(-stepsize))]
        set p180 [expr 1e6/(2.0*$rfmax)]
        set p180steps [expr round($p180/$options(-stepsize))]
        set delsteps [expr ($nsteps-$p180steps)/2]
        set checkduration [expr (2*$delsteps+$p180steps)*$options(-stepsize)]

        # check if total duration is possible with the selected rfmax
        if {$checkduration != $dur} {
            puts "Error: hahn_echo, p180 ($p180 us) cannot centered "
            puts "be digitized with the chosen stepsize ($options(-stepsize) us)"
            puts "Please adjust dur ($dur us), stepsize ($options(-stepsize) us) or rfmax ($rfmax Hz)!"
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
            set t [expr ($delsteps+$i)*$options(-stepsize)]
            set ampl $rfmax
            set ph [expr fmod(360.0e-6*$options(-offset)*$t+$options(-phase),360)]
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
    # Proc to calculate the frequency profile of shapefiles
    # Version 1.0, MRH Dec 2020
    ###########################################################################
    proc duration_check {dur stepsize {verbose 0}} {
        set check [expr $dur/$stepsize]
        if {$check == floor($check)} {
            if {$verbose == 1} {
                # should the number of points not be an integer?
                puts [format "Shape can be resolved. Resulting number of points: %s" $check]
            }
        } else {
            if {$verbose == 1} {
                puts stderr "Shape cannot be resolved using this stepsize. Resulting number of points would be: $check"
                exit 2
            }
        }
    }


    ###########################################################################
    # Proc to bitwise add two shape lists of the same length
	#
    # Changed 01.12.2020 by Max Busskamp:
    #   - Version 1.0
    #   - Version 1.1 Added safety check for list length
    #   - Version 1.2 Fixed Phase Angle calculations
    ###########################################################################
    proc shape_add {shape_list1 shape_list2 args} {
        array set options { -return_complex 0 -verbose 0 }
        array set options $args
        
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

                if {$options(-return_complex) == 0} {
                    set ampl [expr hypot($compl_real_sum, $compl_imag_sum)]
                    # We added 360degrees to get a correct WASP*MODLA Pulse. Please Check here if something seems wrong!
                    set ph  [expr atan2($compl_imag_sum, $compl_real_sum)/$pi*180.0+360.0]
                    set ph [expr fmod($ph,360)]

                    lappend wavelist [format "%6.4f %6.4f" $ampl $ph]
                } else {
                    lappend wavelist [format "%6.4f %6.4f" $compl_real_sum $compl_imag_sum]
                }
            }
            return $wavelist
        } else {
            puts "ERROR: Shapes have not the same size!"
        }
    }
}