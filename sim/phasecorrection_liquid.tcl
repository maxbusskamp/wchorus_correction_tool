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
        global par rfsh1 rfsh2 rfsh3 rfsh_combined
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
        } elseif {[string equal $par(type) "compressedCHORUS_cycled"]} {
            pulse_shaped $par(seq_duration) $rfsh_combined
            store 1
            acq 2 1 $par(ph31)
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
        } elseif {[string equal $par(type) "create_shapes_SHAPEsingle"] || [string equal $par(type) "create_shapes_CHORUS"] || [string equal $par(type) "create_shapes_double_echo"] || [string equal $par(type) "create_shapes_ABSTRUSE"] || [string equal $par(type) "create_shapes_CHORUS_cycled"] || [string equal $par(type) "create_shapes_compressedCHORUS_cycled"] || [string equal $par(type) "create_shapes_double_echo_cycled"] || [string equal $par(type) "create_shapes_double_echo_zerophase"] || [string equal $par(type) "create_shapes_loadshape_double_echo"]} {
        } else {
            puts "Please select excitation mode in main!"
            puts "You selected: $par(type)"
            exit
        }
    }

    proc main {} {
        global par rfsh1 rfsh2 rfsh3 rfsh_combined argc argv

        # Read Arguments from commandline
        if { $argc != 33 } {
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
            set par(phaseoff3)                  [lindex $argv 29]
            set par(compression)                [lindex $argv 30]
            set par(phasecyclesteps)            [lindex $argv 31]
            set par(deadtime)                   [lindex $argv 32]}
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
            if {[string equal $par(type) "create_shapes_compressedCHORUS_cycled"]} {
                set par(phasecycles) $par(phasecyclesteps)
            }
        } elseif {[string equal $par(type) "double_chirp"]} {
            set par(phasecycles) 4
        } elseif {[string equal $par(type) "CHORUS_cycled"]} {
            set par(phasecycles) 16
        } elseif {[string equal $par(type) "compressedCHORUS_cycled"]} {
            set par(phasecycles) $par(phasecyclesteps)
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
        } elseif {[string equal $par(type) "compressedCHORUS_cycled"]} {
            # Set first WURST pulse (excitation)
            set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
            set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
            set rfsh1 [shape2list [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]

            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
            set rfsh2 [shape2list [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]]

            # Set third WURST pulse (refocussing)
            set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
            set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
            set rfsh3 [shape2list [pulsegen $par(shape_type) $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(ph3) $par(stepsize)]]

            if {$par(compression) < $par(tau2)} {
                set delay1_length       [expr $par(tw1)-$par(compression)]
                set delay2_length       [expr $par(tw2)]
                set delay2_short_length [expr $par(tw2)-$par(compression)]
                set delay3_length       [expr $par(tau2)-$par(compression)]
                set delay4_length       [expr $par(tw3)]
                set delay1              [zero_ampl $delay1_length]
                set delay2              [zero_ampl $delay2_length]
                set delay2_short        [zero_ampl $delay2_short_length]
                set delay3              [zero_ampl $delay3_length]
                set delay4              [zero_ampl $delay4_length]

                # puts "Delay 1:$delay1_length"
                # puts "Delay 2:$delay2_length"
                # puts "Delay 2 short:$delay2_short_length"
                # puts "Delay 3:$delay3_length"
                # puts "Delay 4:$delay4_length"

                set par(seq_duration) [expr $delay1_length+$delay2_length+$delay3_length+$par(tw3)]
                # puts "Sequence duration: $par(seq_duration)"

                set rfsh_combined [shape_add [list $rfsh1 $delay2_short $delay3 $delay4]                                            [list $delay1 $rfsh2 $delay3 $delay4]]
                # puts "First Combined Shape: [expr [llength $rfsh_combined]*0.05]"
                set rfsh_combined [list2shape [shape_add [list $rfsh_combined]                                                        [list $delay1 $delay2 $delay3 $rfsh3]]]
                # puts "Second Combined Shape: [expr [llength [shape2list $rfsh_combined]]*0.05]"
            } else {
                set delay1_length       [expr $par(tw1)-$par(compression)]
                set delay2_length       [expr $par(tw2)-($par(compression)-$par(tau2))]
                set delay2_short_length [expr $par(tw2)-$par(compression)-($par(compression)-$par(tau2))]
                set delay4_length       [expr $par(tw3)]
                set delay4_short_length [expr $par(tw3)-($par(compression)-$par(tau2))]
                set delay1              [zero_ampl $delay1_length]
                set delay2              [zero_ampl $delay2_length]
                set delay2_short        [zero_ampl $delay2_short_length]
                set delay4              [zero_ampl $delay4_length]
                set delay4_short        [zero_ampl $delay4_short_length]

                # puts "Delay 1: $delay1_length"
                # puts "Delay 2: $delay2_length"
                # puts "Delay 2 short: $delay2_short_length"
                # puts "Delay 3: 0"
                # puts "Delay 4: $delay4_length"
                # puts "Delay 4 short: $delay4_short_length"

                set par(seq_duration) [expr $delay1_length+$delay2_length+$delay4_length]
                # puts "Sequence duration: $par(seq_duration)"
                
                set rfsh_combined [shape_add [list $rfsh1 $delay2_short $delay4]                                            [list $delay1 $rfsh2 $delay4_short]]
                # puts "First Combined Shape: [expr [llength $rfsh_combined]*0.05]"
                set rfsh_combined [list2shape [shape_add [list $rfsh_combined]                                                        [list $delay1 $delay2 $rfsh3]]]
                # puts "Second Combined Shape: [expr [llength [shape2list $rfsh_combined]]*0.05]"
            }
            printwave [shape2list $rfsh_combined] _combined_before_$index
            save_shape $rfsh_combined $par(filename).simpson_combined_before_$index
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
            set rfsh1 [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]
            set rfsh_shape1 [shape2list [pulsegen $par(shape_type) $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]]

            printwave $rfsh_shape1 1
            save_shape $rfsh1 $par(filename).simpson1
        } elseif {[string equal $par(type) "create_shapes_ABSTRUSE"]} {
            # Set first WURST pulse (excitation)
            set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
            set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
            set rfsh1 [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]
            set rfsh_shape1 [shape2list [pulsegen $par(shape_type) $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]]

            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
            set rfsh2 [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2) $par(stepsize)]
            set rfsh_shape2 [shape2list [pulsegen $par(shape_type) $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2) $par(stepsize)]]

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
            set rfsh1 [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]
            set rfsh_shape1 [shape2list [pulsegen $par(shape_type) $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(phaseoff1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]]

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

            printwave $rfsh_shape1 1
            printwave $rfsh_shape2 2
            printwave $rfsh_shape3 3

            save_shape $rfsh1 $par(filename).simpson1
            save_shape $rfsh2 $par(filename).simpson2
            save_shape $rfsh3 $par(filename).simpson3
        } elseif {[string equal $par(type) "create_shapes_compressedCHORUS_cycled"]} {
            # Set first WURST pulse (excitation)
            set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
            set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
            set rfsh1 [shape2list [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize) -filename_phasecorrect $par(filename_phasecorrect)]]

            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
            set rfsh2 [shape2list [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]]

            # Set third WURST pulse (refocussing)
            set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
            set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
            set rfsh3 [shape2list [pulsegen $par(shape_type) $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(ph3) $par(stepsize)]]

            # Add deadtime before third pulse
            set par(tau2)       [expr [lindex $argv 16]+$par(deadtime)]

            if {$par(compression) < $par(tau2)} {
                set delay1_length       [expr $par(tw1)-$par(compression)]
                set delay2_length       [expr $par(tw2)]
                set delay2_short_length [expr $par(tw2)-$par(compression)]
                set delay3_length       [expr $par(tau2)-$par(compression)]
                set delay4_length       [expr $par(tw3)]
                set delay1              [zero_ampl $delay1_length]
                set delay2              [zero_ampl $delay2_length]
                set delay2_short        [zero_ampl $delay2_short_length]
                set delay3              [zero_ampl $delay3_length]
                set delay4              [zero_ampl $delay4_length]

                # puts "Delay 1:$delay1_length"
                # puts "Delay 2:$delay2_length"
                # puts "Delay 2 short:$delay2_short_length"
                # puts "Delay 3:$delay3_length"
                # puts "Delay 4:$delay4_length"

                set par(seq_duration) [expr $delay1_length+$delay2_length+$delay3_length+$par(tw3)]
                # puts "Sequence duration: $par(seq_duration)"

                set rfsh_combined [shape_add [list $rfsh1 $delay2_short $delay3 $delay4]                                            [list $delay1 $rfsh2 $delay3 $delay4]]
                # puts "First Combined Shape: [expr [llength $rfsh_combined]*0.05]"
                set rfsh_combined [list2shape [shape_add [list $rfsh_combined]                                                        [list $delay1 $delay2 $delay3 $rfsh3]]]
                # puts "Second Combined Shape: [expr [llength [shape2list $rfsh_combined]]*0.05]"
            } else {
                set delay1_length       [expr $par(tw1)-$par(compression)]
                set delay2_length       [expr $par(tw2)-($par(compression)-$par(tau2))]
                set delay2_short_length [expr $par(tw2)-$par(compression)-($par(compression)-$par(tau2))]
                set delay4_length       [expr $par(tw3)]
                set delay4_short_length [expr $par(tw3)-($par(compression)-$par(tau2))]
                set delay1              [zero_ampl $delay1_length]
                set delay2              [zero_ampl $delay2_length]
                set delay2_short        [zero_ampl $delay2_short_length]
                set delay4              [zero_ampl $delay4_length]
                set delay4_short        [zero_ampl $delay4_short_length]

                # puts "Delay 1: $delay1_length"
                # puts "Delay 2: $delay2_length"
                # puts "Delay 2 short: $delay2_short_length"
                # puts "Delay 3: 0"
                # puts "Delay 4: $delay4_length"
                # puts "Delay 4 short: $delay4_short_length"

                set par(seq_duration) [expr $delay1_length+$delay2_length+$delay4_length]
                # puts "Sequence duration: $par(seq_duration)"
                
                set rfsh_combined [shape_add [list $rfsh1 $delay2_short $delay4]                                            [list $delay1 $rfsh2 $delay4_short]]
                # puts "First Combined Shape: [expr [llength $rfsh_combined]*0.05]"
                set rfsh_combined [list2shape [shape_add [list $rfsh_combined]                                                        [list $delay1 $delay2 $rfsh3]]]
                # puts "Second Combined Shape: [expr [llength [shape2list $rfsh_combined]]*0.05]"
            }

            printwave [shape2list $rfsh_combined] _combined_$index
            save_shape $rfsh_combined $par(filename).simpson_combined_$index
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
        if {[string equal $par(type) "double_chirp"] || [string equal $par(type) "compressedCHORUS_cycled"]} {
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
    proc pulsegen {shape_type dur sweepwidth rfmax var1 var2 phase stepsize args} {
        array set options { -filename_phasecorrect 'none' }
        array set options $args
        # var1: N/Zeta/Beta
        # var2: none/tan(kappa)/mean

        if {[string equal $shape_type "WURST"]} {
            set shape [list2shape [wurst $dur $rfmax $sweepwidth $var1 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
        } elseif {[string equal $shape_type "tanhtan"]} {
            set shape [list2shape [tanhtan $dur $rfmax $sweepwidth $var1 $var2 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
        } elseif {[string equal $shape_type "hspulse"]} {
            set shape [list2shape [hs $dur $rfmax $sweepwidth $var1 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
        } elseif {[string equal $shape_type "caWURST"]} {
            set shape [list2shape [cawurst $dur $rfmax $sweepwidth $var1 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
        } elseif {[string equal $shape_type "supergaussian"]} {
            set shape [list2shape [supergaussian $dur $rfmax $sweepwidth $var1 $var2 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
        }

        return $shape
    }


    ###########################################################################
    # Proc for empty, zero shape
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
    # Proc for WURST shape calculation
    ###########################################################################
    proc wurst {dur rfmax Delta N args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 -filename_phasecorrect none }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$options(-stepsize))]

        if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
        } else {
            set phasecorr_file  [open $options(-filename_phasecorrect)]
            set phasecorr       [read $phasecorr_file]
            set phasecorr_list  [split $phasecorr "
"]
            # puts $phasecorr_list
        }

        for {set i 1} {$i <= $nsteps} {incr i} {
            set j [expr $i-1]

            set ampl [expr $rfmax*(1.0-abs(pow(cos(($pi*$options(-stepsize)*$i)/($dur)),$N)))]
            if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
                set ph [expr ((180.0/$pi)*2.0*$pi*(($options(-offset)*1e3+($Delta*1e3/2.0))*$options(-stepsize)*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($options(-stepsize)*1e-6*$i,2)))+$options(-phase)]
            } else {
                set ph [expr ((180.0/$pi)*2.0*$pi*(($options(-offset)*1e3+($Delta*1e3/2.0))*$options(-stepsize)*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($options(-stepsize)*1e-6*$i,2)))+$options(-phase)+[lindex $phasecorr_list $j 0]]
            }
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
    ###########################################################################
    proc cawurst {dur rfmax Delta N args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 -filename_phasecorrect none }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
        } else {
            set phasecorr_file  [open $options(-filename_phasecorrect)]
            set phasecorr       [read $phasecorr_file]
            set phasecorr_list  [split $phasecorr "
"]
        }

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

        set index2 0
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
            if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
                set phase	[lindex $phaselist $index]
            } else {
                set phase	[expr fmod([lindex $phaselist $index]-[lindex $phasecorr_list $index2 0],360)]
            }

            lappend wavelist [format "%6.4f %6.4f" $ampl $phase]
            incr index2
        }
        return $wavelist
    }


    ###########################################################################
    # Proc for supergaussian shape calculation
    ###########################################################################
    proc supergaussian {dur rfmax Delta N G args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 -filename_phasecorrect none }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)
        
        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr round($dur/$options(-stepsize))]

        if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
        } else {
            set phasecorr_file  [open $options(-filename_phasecorrect)]
            set phasecorr       [read $phasecorr_file]
            set phasecorr_list  [split $phasecorr "
"]
        }

        # t = $options(-stepsize)*$i
        for {set i 1} {$i <= $nsteps} {incr i} {
            set j [expr $i-1]

            set ampl [expr $rfmax*exp(-1.0*pow(2,($N+2))*pow(((($options(-stepsize)*$i)-$G)/($dur)),$N))]

            if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
                set ph [expr ((180.0/$pi)*2.0*$pi*(($options(-offset)*1e3+($Delta*1e3/2.0))*$options(-stepsize)*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($options(-stepsize)*1e-6*$i,2)))+$options(-phase)]
            } else {
                set ph [expr ((180.0/$pi)*2.0*$pi*(($options(-offset)*1e3+($Delta*1e3/2.0))*$options(-stepsize)*1e-6*$i-($Delta*1e3/(2.0*$dur*1e-6))*pow($options(-stepsize)*1e-6*$i,2)))+$options(-phase)+[lindex $phasecorr_list $j 0]]
            }
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
    ###########################################################################
    proc tanhtan {dur rfmax Delta zeta tan_kappa args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -direction 0 -filename_phasecorrect none }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
        } else {
            set phasecorr_file  [open $options(-filename_phasecorrect)]
            set phasecorr       [read $phasecorr_file]
            set phasecorr_list  [split $phasecorr "
"]
        }

        set pi [expr 4.0*atan(1.0)]
        set nsteps [expr int(round($dur/$options(-stepsize)))]
        set kappa [expr atan($tan_kappa)]
        set R [expr $Delta*1000*$dur/1000000.0]
        set A [expr ($R*$pi)/($dur/1000000.0)]

        set phi_max [expr -(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos($kappa))]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set j [expr $i-1]
            if {$i <= $nsteps/2} {
                set ampl [expr $rfmax*tanh((2*$i*$options(-stepsize)*$zeta)/$dur)]
            } else {
                set ampl [expr $rfmax*tanh((2*$zeta*(1.0-(($i*$options(-stepsize))/$dur))))]
            }
            if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
                set ph [expr ($phi_max-(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$options(-stepsize)/1000000.0)/($dur/1000000.0))/cos($kappa)))*180/$pi+$options(-phase)]
            } else {
                set ph [expr ($phi_max-(($A*$dur/1000000.0)/(2*$kappa*$tan_kappa))*log(cos((2*$kappa*($i-($nsteps/2))*$options(-stepsize)/1000000.0)/($dur/1000000.0))/cos($kappa)))*180/$pi+$options(-phase)+[lindex $phasecorr_list $j 0]]
            }
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
    ###########################################################################
    proc hs {dur rfmax Delta beta args} {
        array set options { -stepsize 0.05 -verbose 0 -phase 0.0 -offset 0.0 -direction 0 -filename_phasecorrect none }
        array set options $args
        duration_check $dur $options(-stepsize) $options(-verbose)

        if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
        } else {
            set phasecorr_file  [open $options(-filename_phasecorrect)]
            set phasecorr       [read $phasecorr_file]
            set phasecorr_list  [split $phasecorr "
"]
        }

        set nsteps [expr round($dur/$options(-stepsize))]
        set phi0 [expr 180.0*$Delta*1000*$dur/10.6e6]

        for {set i 1} {$i <= $nsteps} {incr i} {
            set j [expr $i-1]

            set x [expr cosh($beta*(2*$i*$options(-stepsize)/$dur-1))]
            set ampl [expr $rfmax/$x]
            if {[string equal $options(-filename_phasecorrect) none] || [string equal $options(-filename_phasecorrect) 'none']} {
                set ph [expr $options(-offset)+$phi0*log($x)+$options(-phase)]
            } else {
                set ph [expr $options(-offset)+$phi0*log($x)+$options(-phase)+[lindex $phasecorr_list $j 0]]
            }
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
    # Proc to check if the given duration is an even divider of the stepsize
    ###########################################################################
    proc duration_check {dur stepsize {verbose 0}} {
        set check [expr $dur/$stepsize]
        if {$check == floor($check)} {
            if {$verbose == 1} {
                puts "Shape can be resolved. Resulting number of points: $check"
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


    ###########################################################################
    # Proc to read file into list
    ###########################################################################
    proc listFromFile {filename} {
        set f [open $filename r]
        set data [split [string trim [read $f]]]
        close $f
        return $data
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