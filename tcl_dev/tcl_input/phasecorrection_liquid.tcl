spinsys {
        channels 1H
        nuclei 1H
        shift 1 0 0 0.0 0 0 0
}

par {
    method           direct
    proton_frequency 500e6
    spin_rate        0
    crystal_file     alpha0beta0
    gamma_angles     1
    np               2
    start_operator   Inz
    detect_operator  Inp
    sw               2000e3
    variable tsw     1e6/sw
    verbose          0000
    conjugate_fid    false
}

proc pulseq {} {
    global par rfsh1 rfsh2 rfsh3
    reset
    matrix set 1 totalcoherence { -1 }
    matrix set 2 totalcoherence { 1 }
    matrix set 3 totalcoherence { -1 1 }

    # Experiment selection done by scanning for the type string
    if {[string equal $par(type) "double_echo"] || [string equal $par(type) "CHORUS"] || [string equal $par(type) "loadshape_double_echo"]} {
        pulse_shaped $par(tw1) $rfsh1
        filter 1
        delay $par(tau1)
        pulse_shaped $par(tw2) $rfsh2
        filter 2
        delay $par(tau2)
        pulse_shaped $par(tw3) $rfsh3
        filter 1
        delay $par(tau3)
        acq
    } elseif {[string equal $par(type) "double_echo_cycled"] || [string equal $par(type) "CHORUS_cycled"]} {
        pulse_shaped $par(tw1) $rfsh1
        delay $par(tau1)
        pulse_shaped $par(tw2) $rfsh2
        delay $par(tau2)
        pulse_shaped $par(tw3) $rfsh3
        delay $par(tau3)
        store 1
        acq 2 1 $par(ph31)
    } elseif {[string equal $par(type) "double_echo_zerophase"]} {
        pulse_shaped $par(tw1) $rfsh1
        delay $par(tau1)
        pulse_shaped $par(tw2) $rfsh2
        delay $par(tau2)
        pulse_shaped $par(tw3) $rfsh3
        delay $par(tau3)
        acq
    } elseif {[string equal $par(type) "double_chirp"]} {
        pulse_shaped $par(tw1) $rfsh1
        pulse_shaped $par(tw2) $rfsh2
        delay $par(tau1)
        store 1
        acq 2 1 $par(ph31)
    } elseif {[string equal $par(type) "SHAPEsingle"]} {
        pulse_shaped $par(tw1) $rfsh1
        store 1
        acq 2 1 $par(ph31)
    } elseif {[string equal $par(type) "create_shapes_SHAPEsingle"] || [string equal $par(type) "create_shapes_CHORUS"] || [string equal $par(type) "create_shapes_double_echo"] || [string equal $par(type) "create_shapes_ABSTRUSE"] || [string equal $par(type) "create_shapes_CHORUS_cycled"] || [string equal $par(type) "create_shapes_double_echo_cycled"] || [string equal $par(type) "create_shapes_double_echo_zerophase"] || [string equal $par(type) "create_shapes_loadshape_double_echo"]} {
    } else {
        puts "Please select excitation mode in main!"
        exit
    }
}

proc main {} {
    global par rfsh1 rfsh2 rfsh3 argc argv

    # Read Arguments from commandline
    if { $argc != 30 } {
        puts "Wrong number of Inputs"
        puts "Please try again."
    } else {
        set par(filename)                   [lindex $argv 1]
        set par(tw1)                        [lindex $argv 2]
        set par(tw2)                        [lindex $argv 3]
        set par(tw3)                        [lindex $argv 4]
        set par(type)                       [lindex $argv 5]
        set par(Delta1)                     [lindex $argv 6]
        set par(Delta2)                     [lindex $argv 7]
        set par(Delta3)                     [lindex $argv 8]
        set par(var11)                      [lindex $argv 9]
        set par(var12)                      [lindex $argv 10]
        set par(var13)                      [lindex $argv 11]
        set par(rf_factor1)                 [lindex $argv 12]
        set par(rf_factor2)                 [lindex $argv 13]
        set par(rf_factor3)                 [lindex $argv 14]
        set par(tau1)                       [lindex $argv 15]
        set par(tau2)                       [lindex $argv 16]
        set par(tau3)                       [lindex $argv 17]
        set par(filename_phasecorrect)      [lindex $argv 18]
        set par(ss_offset)                  [lindex $argv 19]
        set par(var21)                      [lindex $argv 20]
        set par(var22)                      [lindex $argv 21]
        set par(var23)                      [lindex $argv 22]
        set par(shape_type)                 [lindex $argv 23]
        set par(var31)                      [lindex $argv 24]
        set par(var32)                      [lindex $argv 25]
        set par(var33)                      [lindex $argv 26]
        set par(phaseoff1)                  [lindex $argv 27]
        set par(phaseoff2)                  [lindex $argv 28]
        set par(phaseoff3)                  [lindex $argv 29]}

    set par(stepsize)   0.05

    set par(np_tau1)    [expr round($par(tau1)/$par(stepsize))]
    set par(np_tau2)    [expr round($par(tau2)/$par(stepsize))]
    set par(np_tau3)    [expr round($par(tau3)/$par(stepsize))]

    #Define start and end values of constant list in Hz
    set par(start_offset)   [expr -1*($par(Delta1)*1000/2.0)]
    set par(end_offset)     [expr    ($par(Delta1)*1000/2.0)]
    set offset_value_list   [list_offset]

    set results {}

    set par(phasecycles) 16
    set par(ph1_list)  { 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 }
    set par(ph2_list)  { 0 90 180 270 0 90 180 270 0 90 180 270 0 90 180 270 }
    set par(ph3_list)  { 0 0 0 0 90 90 90 90 180 180 180 180 270 270 270 270 }
    set par(ph31_list) { 0 180 0 180 180 0 180 0 0 180 0 180 180 0 180 0 }

    if {[string first "create_shapes_" $par(type)] == 0} {
        set offset_value_list { 0 }
        set par(phasecycles) 1
    } elseif {[string equal $par(type) "double_chirp"]} {
        set par(phasecycles) 4
    } elseif {[string equal $par(type) "CHORUS_cycled"]} {
        set par(phasecycles) 16
    } elseif {[string equal $par(type) "SHAPEsingle"]} {
        set par(phasecycles) 1
    } elseif {[string equal $par(type) "CHORUS"]} {
        set par(phasecycles) 1
    } elseif {[string equal $par(type) "double_echo"]} {
        set par(phasecycles) 1
    } elseif {[string equal $par(type) "loadshape_double_echo"]} {
        set par(phasecycles) 1
    } elseif {[string equal $par(type) "double_echo_cycled"]} {
        set par(phasecycles) 16
    } elseif {[string equal $par(type) "double_echo_zerophase"]} {
        set par(phasecycles) 1
    }

for {set index 0} {$index<$par(phasecycles)} {incr index} {
    set par(ph1) [expr $par(phaseoff1) + [lindex $par(ph1_list) $index]]
    set par(ph2) [expr $par(phaseoff2) + [lindex $par(ph2_list) $index]]
    set par(ph3) [expr $par(phaseoff3) + [lindex $par(ph3_list) $index]]
    set par(ph31) [lindex $par(ph31_list) $index]

    # Set shapes for WURST type experiments
    if {[string equal $par(type) "double_echo"] || [string equal $par(type) "CHORUS"] || [string equal $par(type) "double_echo_cycled"] || [string equal $par(type) "CHORUS_cycled"] || [string equal $par(type) "double_echo_zerophase"]} {
        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
        set rfsh2 [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]

        # Set third WURST pulse (refocussing)
        set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
        set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
        set rfsh3 [pulsegen $par(shape_type) $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(ph3) $par(stepsize)]
    } elseif {[string equal $par(type) "loadshape_double_echo"]} {

        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [load_shape $par(var13)]
        shape_manipulate $rfsh1 -ampl {$ampl*$par(rf1)}
        shape_manipulate $rfsh1 -phase {$phase+$par(ph1)}

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
        set rfsh2 [load_shape $par(var23)]
        shape_manipulate $rfsh2 -ampl {$ampl*$par(rf2)}
        shape_manipulate $rfsh2 -phase {$phase+$par(ph2)}

        # Set third WURST pulse (refocussing)
        set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
        set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
        set rfsh3 [load_shape $par(var33)]
        shape_manipulate $rfsh3 -ampl {$ampl*$par(rf3)}
        shape_manipulate $rfsh3 -phase {$phase+$par(ph3)}
    } elseif {[string equal $par(type) "double_chirp"]} {
        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
        set rfsh2 [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]

    } elseif {[string equal $par(type) "SHAPEsingle"]} {
        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]
    } elseif {[string equal $par(type) "create_shapes_SHAPEsingle"]} {
        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [pulsegen_corr $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) $par(filename_phasecorrect)]
        set rfsh_shape1 [pulsegen_corr_list $par(shape_type) $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) $par(filename_phasecorrect)]

        printwave $rfsh_shape1 1
        save_shape $rfsh1 $par(filename).simpson1
    } elseif {[string equal $par(type) "create_shapes_ABSTRUSE"]} {
        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [pulsegen_corr $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) $par(filename_phasecorrect)]
        set rfsh_shape1 [pulsegen_corr_list $par(shape_type) $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) $par(filename_phasecorrect)]

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
        set rfsh2 [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2) $par(stepsize)]
        set rfsh2 [shape2list [pulsegen $par(shape_type) $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2) $par(stepsize)]]

        printwave $rfsh_shape1 1
        printwave $rfsh_shape2 2

        save_shape $rfsh1 $par(filename).simpson1
        save_shape $rfsh2 $par(filename).simpson2
    } elseif {[string equal $par(type) "create_shapes_loadshape_double_echo"]} {

        set phasecorr_list [listFromFile $par(filename_phasecorrect)]

        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]

        set rfsh1 [load_shape $par(var13)]
        shape_manipulate $rfsh1 -ampl {$ampl*$par(rf1)}
        shape_manipulate $rfsh1 -phase {$phase+$par(phaseoff1)}

        set rfsh_shape1 [load_shape $par(var13)]
        shape_manipulate $rfsh_shape1 -ampl {$ampl*100}
        shape_manipulate $rfsh1 -phase {$phase+$par(phaseoff1)}
        shape_manipulate $rfsh1 -phase {$phase-[lindex $phasecorr_list $i-1]}
        shape_manipulate $rfsh_shape1 -phase {$phase-[lindex $phasecorr_list $i-1]}
        set rfsh_shape1 [shape2list $rfsh_shape1]

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]

        set rfsh2 [load_shape $par(var23)]
        shape_manipulate $rfsh2 -ampl {$ampl*$par(rf2)}
        shape_manipulate $rfsh2 -phase {$phase+$par(phaseoff2)}

        set rfsh_shape2 [load_shape $par(var23)]
        shape_manipulate $rfsh_shape2 -ampl {$ampl*100}
        shape_manipulate $rfsh2 -phase {$phase+$par(phaseoff2)}
        set rfsh_shape2 [shape2list $rfsh_shape2]

        # Set third WURST pulse (refocussing)
        set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
        set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
        set rfsh3 [load_shape $par(var33)]
        shape_manipulate $rfsh3 -ampl {$ampl*$par(rf3)}
        shape_manipulate $rfsh2 -phase {$phase+$par(phaseoff2)}

        set rfsh_shape3 [load_shape $par(var33)]
        shape_manipulate $rfsh_shape3 -ampl {$ampl*100}
        shape_manipulate $rfsh2 -phase {$phase+$par(phaseoff2)}
        set rfsh_shape3 [shape2list $rfsh_shape3]

        printwave $rfsh_shape1 1
        printwave $rfsh_shape2 2
        printwave $rfsh_shape3 3

        save_shape $rfsh1 $par(filename).simpson1
        save_shape $rfsh2 $par(filename).simpson2
        save_shape $rfsh3 $par(filename).simpson3
    } elseif {[string equal $par(type) "create_shapes_CHORUS"] || [string equal $par(type) "create_shapes_double_echo"] || [string equal $par(type) "create_shapes_double_echo_cycled"] || [string equal $par(type) "create_shapes_CHORUS__cycled"] || [string equal $par(type) "create_shapes_double_echo_zerophase"]} {
        # Set first WURST pulse (excitation)
        set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
        set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
        set rfsh1 [pulsegen_corr $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) $par(filename_phasecorrect)]
        set rfsh_shape1 [pulsegen_corr_list $par(shape_type) $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) $par(filename_phasecorrect)]

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
        set rfsh2 [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2) $par(stepsize)]
        set rfsh_shape2 [shape2list [pulsegen $par(shape_type) $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2) $par(stepsize)]]

        # Set third WURST pulse (refocussing)
        set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
        set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
        set rfsh3 [pulsegen $par(shape_type) $par(tw3) $par(Delta3) $par(rf2) $par(var31) $par(var32) $par(phaseoff3) $par(stepsize)]
        set rfsh_shape3 [shape2list [pulsegen $par(shape_type) $par(tw3) $par(Delta3) 100 $par(var31) $par(var32) $par(phaseoff3) $par(stepsize)]]

        if {[string equal $par(shape_type) "WURST"]} {
            set rfsh3 [list2shape [wurst $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(phaseoff3)]]
            set rfsh_shape3 [wurst $par(tw3) $par(Delta3) 100 $par(var31) $par(phaseoff3)]
        } elseif {[string equal $par(shape_type) "tanhpulse"]} {
            set rfsh3 [list2shape [tanhpulse $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(phaseoff3)]]
            set rfsh_shape3 [tanhpulse $par(tw3) $par(Delta3) 100 $par(var31) $par(var32) $par(phaseoff3)]
        } elseif {[string equal $par(shape_type) "hspulse"]} {
            set rfsh3 [list2shape [hspulse $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(phaseoff3)]]
            set rfsh_shape3 [hspulse $par(tw3) $par(Delta3) 100 $par(var31) $par(phaseoff3)]
        } elseif {[string equal $par(shape_type) "caWURST"]} {
            set rfsh3 [list2shape [cawurst $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(phaseoff3)]]
            set rfsh_shape3 [cawurst $par(tw3) $par(Delta3) 100 $par(var31) $par(phaseoff3)]
        } elseif {[string equal $par(shape_type) "supergaussian"]} {
            set rfsh3 [list2shape [supergaussian $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(phaseoff3)]]
            set rfsh_shape3 [supergaussian $par(tw3) $par(Delta3) 100 $par(var31) $par(var32) $par(phaseoff3)]
        }


        printwave $rfsh_shape1 1
        printwave $rfsh_shape2 2
        printwave $rfsh_shape3 3

        save_shape $rfsh1 $par(filename).simpson1
        save_shape $rfsh2 $par(filename).simpson2
        save_shape $rfsh3 $par(filename).simpson3
    }

        foreach offset_value        $offset_value_list {
            set par(offset)         $offset_value

            # Make simulation
            set f [fsimpson [list [list shift_1_iso $par(offset)]]]
            # Save final spectra and move them to the corresponding output directories
            set re [findex $f 1 -re]
            set im [findex $f 1 -im]
            # Calculate magnitude of xy projection (of I1p)
            set xyproj [expr sqrt(($re)**2+($im)**2)]

            lappend results [format "%s %s %s %s" $par(offset) $re $im $xyproj]

            funload $f
        }
free_all_shapes
}

# Extraction and summation
set results_sum {}
set n [llength $offset_value_list]

for {set i 0} {$i < $n} {incr i} {
    set results_sum_temp {}
    set results_sum_temp [lmap [lreplace [lrepeat $n b] $i $i a] $results {set a}]

    set re 0.0
    set im 0.0

    foreach elem $results_sum_temp {
        set offs     [lindex $elem 0]
        set re       [expr $re+[lindex $elem 1]]
        set im       [expr $im+[lindex $elem 2]]
        set xyproj   [expr sqrt(($re)**2+($im)**2)]
    }
    if {[string equal $par(type) "double_chirp"]} {
        lappend results_sum [format "%s %s %s %s" $offs [expr $re/2.0] [expr $im/2.0] [expr $xyproj/2.0]]
    } elseif {[string equal $par(type) "double_echo_cycled"] || [string equal $par(type) "CHORUS_cycled"]} {
        lappend results_sum [format "%s %s %s %s" $offs [expr $re/8.0] [expr $im/8.0] [expr $xyproj/8.0]]
    } else {
        lappend results_sum [format "%s %s %s %s" $offs [expr $re*2.0] [expr $im*2.0] [expr $xyproj*2.0]]
    }
}

# Write output
    if {[string first "create_shapes_" $par(type)] != 0} {
        set fileID [open $par(filename).out "w"]
        foreach l $results_sum {
            puts $fileID $l
        }
        close $fileID
    }
}


###########################################################################
# Wrapper for shape generation
# Changed 02.03.2022 by Max Bußkamp:
#   - Initial version
###########################################################################
proc pulsegen {shape_type dur sweepwidth rfmax var1 var2 phase stepsize} {
    # var1: N/Zeta/Beta
    # var2: none/tan(kappa)/mean

    if {[string equal $shape_type "WURST"]} {
        set shape [list2shape [wurst $dur $sweepwidth $rfmax $var1 $phase $stepsize]]
    } elseif {[string equal $shape_type "tanhpulse"]} {
        set shape [list2shape [tanhpulse $dur $sweepwidth $rfmax $var1 $var2 $phase $stepsize]]
    } elseif {[string equal $shape_type "hspulse"]} {
        set shape [list2shape [hspulse $dur $sweepwidth $rfmax $var1 $phase $stepsize]]
    } elseif {[string equal $shape_type "caWURST"]} {
        set shape [list2shape [cawurst $dur $sweepwidth $rfmax $var1 $phase $stepsize]]
    } elseif {[string equal $shape_type "supergaussian"]} {
        set shape [list2shape [supergaussian $dur $sweepwidth $rfmax $var1 $var2 $phase $stepsize]]
    }

    return $shape
}


###########################################################################
# Wrapper for shape generation with phasecorrection
# Changed 02.03.2022 by Max Bußkamp:
#   - Initial version
###########################################################################
proc pulsegen_corr {shape_type dur sweepwidth rfmax var1 var2 phase stepsize filename_phasecorrect} {
    # var1: N/Zeta/Beta
    # var2: none/tan(kappa)/mean

    if {[string equal $par(shape_type) "WURST"]} {
        set shape [list2shape [wurst_phasecorr $dur $sweepwidth $rfmax $var1 $filename_phasecorrect $phase $stepsize]]
    } elseif {[string equal $par(shape_type) "tanhpulse"]} {
        set shape [list2shape [tanhpulse_phasecorr $dur $sweepwidth $rfmax $var1 $var2 $filename_phasecorrect $phase $stepsize]]
    } elseif {[string equal $par(shape_type) "hspulse"]} {
        set shape [list2shape [hspulse_phasecorr $dur $sweepwidth $rfmax $var1 $filename_phasecorrect $phase $stepsize]]
    } elseif {[string equal $par(shape_type) "caWURST"]} {
        set shape [list2shape [cawurst_phasecorr $dur $sweepwidth $rfmax $var1 $filename_phasecorrect $phase $stepsize]]
    } elseif {[string equal $par(shape_type) "supergaussian"]} {
        set shape [list2shape [supergaussian_phasecorr $dur $sweepwidth $rfmax $var1 $var2 $filename_phasecorrect $phase $stepsize]]
    }

    return $shape
}


proc pulsegen_corr_list {shape_type dur sweepwidth rfmax var1 var2 phase stepsize filename_phasecorrect} {
    # var1: N/Zeta/Beta
    # var2: none/tan(kappa)/mean

    if {[string equal $par(shape_type) "WURST"]} {
        set shape_list       [wurst_phasecorr $dur $sweepwidth $rfmax $var1 $filename_phasecorrect $phase $stepsize]
    } elseif {[string equal $par(shape_type) "tanhpulse"]} {
        set shape_list       [tanhpulse_phasecorr $dur $sweepwidth $rfmax $var1 $var2 $filename_phasecorrect $phase $stepsize]
    } elseif {[string equal $par(shape_type) "hspulse"]} {
        set shape_list       [hspulse_phasecorr $dur $sweepwidth $rfmax $var1 $filename_phasecorrect $phase $stepsize]
    } elseif {[string equal $par(shape_type) "caWURST"]} {
        set shape_list       [cawurst_phasecorr $dur $sweepwidth $rfmax $var1 $filename_phasecorrect $phase $stepsize]
    } elseif {[string equal $par(shape_type) "supergaussian"]} {
        set shape_list       [supergaussian_phasecorr $dur $sweepwidth $rfmax $var1 $var2 $filename_phasecorrect $phase $stepsize]
    }

    return $shape_list
}


###########################################################################
# Proc to read file into list
# Changed 08.06.2020 by Max Bußkamp:
###########################################################################
proc listFromFile {filename} {
    set f [open $filename r]
    set data [split [string trim [read $f]]]
    close $f
    return $data
}


###########################################################################
# Proc for caWURST shape calculation 
# Only use sweep direction 1
# Changed 16.08.2020 by Max Bußkamp:
#   - Added rfmax to input variables
#   - Added default values for stepsize, direction and offset
###########################################################################
proc cawurst {tw Delta rfmax N {phaseoffset 0.0} {stepsize 0.05} {direction -1} {offset 0.0}} {
    global par

    #Variables
    set pi2 		[expr (4*atan(1))/2]

    # sweep / MHz
    set sweep		[expr $Delta/1000.0]

    set rate 		[expr 1.0e6*($sweep/$tw)]

    #On-resonance q-factor
    set q0			[expr (pow($rfmax,2))/(1e6*$rate)]

    set nsp			[expr round($tw/$stepsize)]

    #Amplitude
    set b 			[expr (2*$pi2)/$tw]

    #Frequency
    set w0 			0
    set k 			[expr 1.0e-6*$rate]

    #Calculate frequency correction factor
    set freq_prep_b 0.0
    set phase_prep1	0.0
    set freq_max 0.0
    for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
        set time				[expr $stepsize*$i]
        set amp_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
        lappend amplist	$amp_prep
        
        set freq_prep 	[expr (1.0e-6*$stepsize*($amp_prep**2/$q0)+$freq_prep_b)]
        lappend freqlist $freq_prep
        if {$freq_prep > $freq_max} {
            set freq_max $freq_prep
        }

        set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$tw)))]
        set phase_prep 	[expr fmod($phase_prep1*360+$phaseoffset,360)]
        lappend phaselist $phase_prep
        
        set freq_prep_b $freq_prep
    }

    set sweepfactor [expr $sweep*1e6/2.0/$freq_max]

    unset time
    unset amp_prep
    unset amplist
    unset freq_prep
    unset freqlist
    unset phase_prep1
    unset phase_prep
    unset phaselist
    unset freq_prep_b
    
    #Calculate amplitude and frequency for positive offset values 
    set freq_prep_b 0.0
    set phase_prep1	0.0
    for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
        set time				[expr $stepsize*$i]
        set amp_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
        lappend amplist	$amp_prep
        
        set freq_prep 	[expr (1.0e-6*$stepsize*$sweepfactor*($amp_prep**2/$q0)+$freq_prep_b)]
        lappend freqlist $freq_prep

        set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$tw)))]
        set phase_prep 	[expr fmod($phase_prep1*360+$phaseoffset,360)]
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
        set amp			[lindex $amplist $index]
        set freq_help 	[lindex $freqlist $index]
        set freq 		[expr $help*$freq_help]
        set phase		[lindex $phaselist $index]

        lappend wavelist [format "%6.2f %6.2f" $amp $phase]
    }
    return $wavelist
}


###########################################################################
# Proc for WURST shape calculation
# Version 1 16.09.2020 by Max Bußkamp:
###########################################################################
proc supergaussian {tw Delta rfmax N G {phaseoffset 0.0} {stepsize 0.05} {direction 1} {offset 0.0}} {
    global par

    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]

    # t = $stepsize*$i
    for {set i 1} {$i <= $nsteps} {incr i} {

        set amp [expr $rfmax*exp(-1.0*pow(2,($N+2))*pow(((($stepsize*$i)-$G)/($tw)),$N))]

        set ph [expr ((180.0/$pi)*2.0*$pi*(($offset*1e3+($Delta*1e3/2.0))*$stepsize*1e-6*$i-($Delta*1e3/(2.0*$tw*1e-6))*pow($stepsize*1e-6*$i,2)))+$phaseoffset]
        if {$direction} {
            set ph [expr fmod(-$ph,360)+360.0]
        } else {
            set ph [expr fmod($ph,360)]
        }
        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for WURST shape calculation
# Changed 20.01.2020 by Max Bußkamp:
#   - Added Option for Phasecycle
#   - Added rfmax to input variables
#   - Added default values for stepsize, direction and offset
###########################################################################
proc wurst {tw Delta rfmax N {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {
    global par

    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set amp [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($tw)),$N)))]
        set ph [expr ((180.0/$pi)*2.0*$pi*(($offset*1e3+($Delta*1e3/2.0))*$stepsize*1e-6*$i-($Delta*1e3/(2.0*$tw*1e-6))*pow($stepsize*1e-6*$i,2)))+$phaseoffset]
        if {$direction} {
            set ph [expr fmod($ph,360)]
        } else {
            set ph [expr fmod(-$ph,360)+360.0]
        }
        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for tanh/tan shape calculation
# Version 1.0 Max Busskamp 21.09.2019
###########################################################################
proc tanhpulse {tw Delta rfmax zeta tan_kappa {phaseoffset 0.0} {stepsize 0.05} {direction 0}} {
    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr int(round($tw/$stepsize))]
    set kappa [expr atan($tan_kappa)]
    set R [expr $Delta*1000*$tw/1000000.0]
    set A [expr ($R*$pi)/($tw/1000000.0)]

    set phi_max [expr -(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

    for {set i 1} {$i <= $nsteps} {incr i} {
        if {$i <= $nsteps/2} {
            set amp [expr $rfmax*tanh((2*$i*$stepsize*$zeta)/$tw)]
        } else {
            set amp [expr $rfmax*tanh((2*$zeta*(1.0-(($i*$stepsize)/$tw))))]
        }
        set ph [expr ($phi_max-(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$stepsize/1000000.0)/($tw/1000000.0))/cos($kappa)))*180/$pi+$phaseoffset]
        if {$direction} {
            set ph [expr fmod(-$ph,360)+360.0]
        } else {
            set ph [expr fmod($ph,360)]
        }

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for HS shape calculation
# Version 1.1 MRH Sept 2016
###########################################################################
proc hspulse {tw Delta rfmax beta {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {

    set nsteps [expr round($tw/$stepsize)]
    set phi0 [expr 180.0*$Delta*1000*$tw/10.6e6]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set x [expr cosh($beta*(2*$i*$stepsize/$tw-1))]
        set amp [expr $rfmax/$x]
        set ph [expr $offset+$phi0*log($x)+$phaseoffset]

        if {$direction} {
            set ph [expr fmod(-$ph,360)+360.0]
        } else {
            set ph [expr fmod($ph,360)]
        }

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for caWURST shape calculation 
# Only use sweep direction 1
# Changed 16.08.2020 by Max Bußkamp:
#   - Added rfmax to input variables
#   - Added default values for stepsize, direction and offset
###########################################################################
proc cawurst_phasecorr {tw Delta rfmax N filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 1} {offset 0.0}} {
    global par

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]

    #Variables
    set pi2 		[expr (4*atan(1))/2]

    # sweep / MHz
    set sweep		[expr $Delta/1000.0]

    set rate 		[expr 1.0e6*($sweep/$tw)]

    #On-resonance q-factor
    set q0			[expr (pow($rfmax,2))/(1e6*$rate)]

    set nsp			[expr round($tw/$stepsize)]

    #Amplitude
    set b 			[expr (2*$pi2)/$tw]

    #Frequency
    set w0 			0
    set k 			[expr 1.0e-6*$rate]

    #Calculate frequency correction factor
    set freq_prep_b 0.0
    set phase_prep1	0.0
    set freq_max 0.0
    for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
        set time				[expr $stepsize*$i]
        set amp_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
        lappend amplist	$amp_prep
        
        set freq_prep 	[expr (1.0e-6*$stepsize*($amp_prep**2/$q0)+$freq_prep_b)]
        lappend freqlist $freq_prep
        if {$freq_prep > $freq_max} {
            set freq_max $freq_prep
        }

        set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$tw)))]
        set phase_prep 	[expr fmod($phase_prep1*360+$phaseoffset,360)]
        lappend phaselist $phase_prep
        
        set freq_prep_b $freq_prep
    }

    set sweepfactor [expr $sweep*1e6/2.0/$freq_max]

    unset time
    unset amp_prep
    unset amplist
    unset freq_prep
    unset freqlist
    unset phase_prep1
    unset phase_prep
    unset phaselist
    unset freq_prep_b

    #Calculate amplitude and frequency for positive offset values 
    set freq_prep_b 0.0
    set phase_prep1	0.0
    for {set i [expr round($nsp/2)+1]} {$i <= $nsp} {incr i} {
        set time				[expr $stepsize*$i]
        set amp_prep			[expr $rfmax*(1-abs((cos($b*$time))**($N)))]
        lappend amplist	$amp_prep
        
        set freq_prep 	[expr (1.0e-6*$stepsize*$sweepfactor*($amp_prep**2/$q0)+$freq_prep_b)]
        lappend freqlist $freq_prep
        
        set phase_prep1 [expr $phase_prep1+($freq_prep/(round($nsp)/(1.0e-6*$tw)))]
        set phase_prep 	[expr $phase_prep1*360+$phaseoffset]
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

    set index2 0
    for {set j 1} {$j < $nsp} {incr j} {
        set index [expr $freqlist_length-$j]
        if {$index < 0} {
            set index [expr abs($index)]
            set help $direc
        } else {
            set help [expr -1.0*$direc]
        }
        set time  		[expr $stepsize*$j]
        set amp			[lindex $amplist $index]
        set freq_help 	[lindex $freqlist $index]
        set freq 		[expr $help*$freq_help]
        set phase		[expr fmod([lindex $phaselist $index]-[lindex $phasecorr_list $index2 0],360)]

        lappend wavelist [format "%6.2f %6.2f" $amp $phase]
        incr index2
    }
    return $wavelist
}


###########################################################################
# Proc for WURST shape calculation
# Version 1 16.09.2020 by Max Bußkamp:
###########################################################################
proc supergaussian_phasecorr {tw Delta rfmax N G filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 1} {offset 0.0}} {
    global par

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]
    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]

    # t = $stepsize*$i
    for {set i 1} {$i <= $nsteps} {incr i} {
        set j [expr $i-1]

        set amp [expr $rfmax*exp(-1.0*pow(2,($N+2))*pow(((($stepsize*$i)-$G)/($tw)),$N))]

        set ph [expr ((180.0/$pi)*2.0*$pi*(($offset*1e3+($Delta*1e3/2.0))*$stepsize*1e-6*$i-($Delta*1e3/(2.0*$tw*1e-6))*pow($stepsize*1e-6*$i,2)))+$phaseoffset+[lindex $phasecorr_list $j 0]]
        if {$direction} {
            set ph [expr fmod(-$ph,360)+360.0]
        } else {
            set ph [expr fmod($ph,360)]
        }
        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for phasecorrected WURST shape calculation
# Changed 09.07.2019 by Max Bußkamp:
###########################################################################
proc wurst_phasecorr {tw Delta rfmax N filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {
    global par

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]
    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set j [expr $i-1]
        set amp [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($tw)),$N)))]
        set ph [expr ((180.0/$pi)*2.0*$pi*(($offset*1e3+($Delta*1e3/2.0))*$stepsize*1e-6*$i-($Delta*1e3/(2.0*$tw*1e-6))*pow($stepsize*1e-6*$i,2)))+$phaseoffset+[lindex $phasecorr_list $j 0]]
        if {$direction} {
            set ph [expr fmod($ph,360)]
        } else {
            set ph [expr fmod(-$ph,360)+360.0]
        }
        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for tanh/tan shape calculation
# Version 1.0 Max Busskamp 21.09.2019
###########################################################################
proc tanhpulse_phasecorr {tw Delta rfmax zeta tan_kappa filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 0}} {

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]
    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr int(round($tw/$stepsize))]
    set kappa [expr atan($tan_kappa)]
    set R [expr $Delta*1000*$tw/1000000.0]
    set A [expr ($R*$pi)/($tw/1000000.0)]

    set phi_max [expr -(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set j [expr $i-1]
        if {$i <= $nsteps/2} {
            set amp [expr $rfmax*tanh((2*$i*$stepsize*$zeta)/$tw)]
        } else {
            set amp [expr $rfmax*tanh((2*$zeta*(1.0-(($i*$stepsize)/$tw))))]
        }
        set ph [expr ($phi_max-(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$stepsize/1000000.0)/($tw/1000000.0))/cos($kappa)))*180/$pi+$phaseoffset+[lindex $phasecorr_list $j 0]]
        if {$direction} {
            set ph [expr fmod($ph,360)]
        } else {
            set ph [expr fmod(-$ph,360)+360.0]
        }

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for HS shape calculation
# Version 1.1 MRH Sept 2016
###########################################################################
proc hspulse_phasecorr {tw Delta rfmax beta filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]
    set nsteps [expr round($tw/$stepsize)]
    set phi0 [expr 180.0*$Delta*1000*$tw/10.6e6]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set j [expr $i-1]
        set x [expr cosh($beta*(2*$i*$stepsize/$tw-1))]
        set amp [expr $rfmax/$x]
        set ph [expr $offset+$phi0*log($x)+$phaseoffset+[lindex $phasecorr_list $j 0]]

        if {$direction} {
            set ph [expr fmod(-$ph,360)+360.0]
        } else {
            set ph [expr fmod($ph,360)]
        }

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for sigmoid shape calculation
# Changed 29.04.2020 by Max Bußkamp:
#
# Version 1.0 Max Busskamp 29.04.2020
###########################################################################
proc sigmoid {tw Delta rfmax k d N {phaseoffset 0.0} {stepsize 0.05}} {
    global par

    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set amp [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($tw)),$N)))]
        set ph [expr (($d*exp($k*pow((2.0*$i/$nsteps)-1.0,$d))*$Delta*1e3*$k*pow((2.0*$i/$nsteps)-1.0,$d-1))/pow((1+exp($k*pow((2.0*$i/$nsteps)-1.0,$d))),2))*180/$pi+$phaseoffset]
        
        set ph [expr fmod($ph,360)]

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for WURST shape with TANH phase calculation
# Changed 14.11.2019 by Max Bußkamp:
#   - Version 1.0
###########################################################################
proc wurst_amp_tanh_phase {tw Delta rfmax N tan_kappa {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {
    global par

    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]
    set kappa [expr atan($tan_kappa)]
    set R [expr $Delta*1000*$tw/1000000.0]
    set A [expr ($R*$pi)/($tw/1000000.0)]

    set phi_max [expr -(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set amp [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($tw)),$N)))]

        set ph [expr ($phi_max-(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$stepsize/1000000.0)/($tw/1000000.0))/cos($kappa)))*180/$pi+$phaseoffset]
        if {$direction} {
            set ph [expr fmod($ph,360)]
        } else {
            set ph [expr fmod(-$ph,360)+360.0]
        }

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for sigmoid shape calculation
# Changed 29.04.2020 by Max Bußkamp:
#
# Version 1.0 Max Busskamp 29.04.2020
###########################################################################
proc sigmoid_phasecorr {tw Delta rfmax k d N filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05}} {
    global par

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]
    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set j [expr $i-1]
        set amp [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($tw)),$N)))]
        set ph [expr (($d*exp($k*pow((2.0*$i/$nsteps)-1.0,$d))*$Delta*1e3*$k*pow((2.0*$i/$nsteps)-1.0,$d-1))/pow((1+exp($k*pow((2.0*$i/$nsteps)-1.0,$d))),2))*180/$pi+$phaseoffset+[lindex $phasecorr_list $j 0]]
        
        set ph [expr fmod($ph,360)]

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


###########################################################################
# Proc for WURST shape with TANH phase calculation
# Changed 14.11.2019 by Max Bußkamp:
#   - Version 1.0
###########################################################################
proc wurst_amp_tanh_phase_phasecorr {tw Delta rfmax N tan_kappa filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {
    global par

    set phasecorr_file  [open $filename_phasecorrect]
    set phasecorr       [read $phasecorr_file]
    set phasecorr_list  [split $phasecorr ""]
    set pi [expr 4.0*atan(1.0)]
    set nsteps [expr round($tw/$stepsize)]
    set kappa [expr atan($tan_kappa)]
    set R [expr $Delta*1000*$tw/1000000.0]
    set A [expr ($R*$pi)/($tw/1000000.0)]

    set phi_max [expr -(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

    for {set i 1} {$i <= $nsteps} {incr i} {
        set j [expr $i-1]
        set amp [expr $rfmax*(1.0-abs(pow(cos(($pi*$stepsize*$i)/($tw)),$N)))]

        set ph [expr ($phi_max-(($A*$tw/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$stepsize/1000000.0)/($tw/1000000.0))/cos($kappa)))*180/$pi+$phaseoffset+[lindex $phasecorr_list $j 0]]
        if {$direction} {
            set ph [expr fmod($ph,360)]
        } else {
            set ph [expr fmod(-$ph,360)+360.0]
        }

        lappend wavelist [format "%6.2f %6.2f" $amp $ph]
    }
    return $wavelist
}


proc list_offset {} {
        global par

        for {set i $par(start_offset)} {$i<=$par(end_offset)} {set i [expr $i+$par(ss_offset)]} {
                set offs [expr -1*$i]
                lappend list_offs $offs
        }
        return $list_offs
}


proc list_gen { list_start list_end list_step} {

        for {set i $list_start} {$i<=$list_end} {set i [expr $i+$list_step]} {
                set list_entry [expr $i]
                set list_entry [format "%.2f" $list_entry]
                lappend list $list_entry
        }
        return $list
}


###########################################################################
# Proc for generating an output file from shapes
# 09.07.2019 by Max Bußkamp
###########################################################################
proc printwave {wave counter} {
    global par
    set filename $par(filename).shape$counter
    if {[file exists $filename]} {
    puts "Warning: $filename exists and will be overwritten!"
    }
    set fp [open $filename w]
    set np [llength $wave]
    
    puts $fp "##TITLE= WURST-${par(var11)} shape (duration $par(tw1) us, sweep width $par(Delta1) kHz, step size: 0.05 us, sweep rate: $par(sweep_rate1) MHz/ms)"
    puts $fp "##USAGE= WURST pulse for inversion and excitation"
    puts $fp "##JCAMP-DX= 5.00 \$\$ Bruker JCAMP library"
    puts $fp "##DATA TYPE= Shape Data"
    puts $fp "##ORIGIN= Generated from wurst program"
    puts $fp "##DATE= "
    puts $fp "##TIME= "
    puts $fp "##\$SHAPE_PARAMETERS= Type: Wurst ; Total Sweep-Width \[Hz\] [expr $par(Delta1)*1000.0] ; Length of Pulse \[usec\] [expr 1.0*$par(tw1)] ; Amplitude Power Index ${par(var11)}.0 ; 1=High to low field, -1=Low to high field 1"
    puts $fp "##MINX= 0.000000e+00"
    puts $fp "##MAXX= 1.000000e+02"
    puts $fp "##MINY= 0.000000e+00"
    puts $fp "##MAXY= 3.600000e+02"
    puts $fp "##\$SHAPE_EXMODE= Adiabatic"
    puts $fp "##\$SHAPE_TOTROT= 1.800000e+02"
    puts $fp "##\$SHAPE_TYPE= Inversion"
    puts $fp "##\$SHAPE_BWFAC= 0.000000e+00"
    puts $fp "##\$SHAPE_INTEGFAC= 0.000000e+00"
    puts $fp "##\$SHAPE_MODE= 1"
    puts $fp "##NPOINTS= $np"
    puts $fp "##XYPOINTS= (XY..XY)"

    foreach l $wave {
    puts $fp [format "%.6e, %.6e" [lindex $l 0] [lindex $l 1]]
    }
    puts $fp "##END= "
    close $fp
}