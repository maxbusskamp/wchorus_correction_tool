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

    lappend ::auto_path ../src/
    if {![namespace exists shapetools_liquid]} {
        package require shapetools_liquid
        namespace import shapetools_liquid::*
    }

    # Read Arguments from commandline
    if { $argc != 31 } {
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
        set par(compression)                [lindex $argv 30]}
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
    } elseif {[string equal $par(type) "compressedCHORUS_cycled"]} {
        set par(phasecycles) 1
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

        if {[string equal $par(shape_type) "WURST"]} {
            set rfsh1 [shape2list [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]
            
            set testshape [pulsegen $par(shape_type) $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]
            save_shape $testshape original.shape1
        } else {
            puts "Only useable with WURST shape"
            exit
        }

        # Set second WURST pulse (refocussing)
        set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
        set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]

        if {[string equal $par(shape_type) "WURST"]} {
            set rfsh2 [shape2list [pulsegen $par(shape_type) $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]]
        } else {
            puts "Only useable with WURST shape"
            exit
        }


        # Set third WURST pulse (refocussing)
        set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
        set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]

        if {[string equal $par(shape_type) "WURST"]} {
            set rfsh3 [shape2list [pulsegen $par(shape_type) $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(ph3) $par(stepsize)]]
        } else {
            puts "Only useable with WURST shape"
            exit
        }

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

            set rfsh_combined [shape_add [list $rfsh1 $delay2_short $delay3 $delay4]\
                                         [list $delay1 $rfsh2 $delay3 $delay4]]
            # puts "First Combined Shape: [expr [llength $rfsh_combined]*0.05]"
            set rfsh_combined [list2shape [shape_add [list $rfsh_combined]\
                                                     [list $delay1 $delay2 $delay3 $rfsh3]]]
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
            
            set rfsh_combined [shape_add [list $rfsh1 $delay2_short $delay4]\
                                         [list $delay1 $rfsh2 $delay4_short]]
            # puts "First Combined Shape: [expr [llength $rfsh_combined]*0.05]"
            set rfsh_combined [list2shape [shape_add [list $rfsh_combined]\
                                                     [list $delay1 $delay2 $rfsh3]]]
            # puts "Second Combined Shape: [expr [llength [shape2list $rfsh_combined]]*0.05]"
        }
        save_shape $rfsh_combined combined.shape
        
        if {$par(compression) < $par(tau2)} {
            set delay1_length [expr $par(tw1)-$par(compression)]
            set delay2_length [expr $par(tw2)-$par(compression)]
            set delay3_length [expr $par(tau2)-$par(compression)]
            set delay4_length [expr $par(tw3)-$par(compression)]
            set delay1 [zero_ampl $delay1_length]
            set delay2 [zero_ampl $delay2_length]
            set delay3 [zero_ampl $delay3_length]
            set delay4 [zero_ampl $delay4_length]

            set par(seq_duration) [expr $delay1_length+$delay2_length+$delay3_length+$par(tw3)]

            set rfsh_combined [shape_add [list [shape2list $rfsh1] $delay2 $delay3 $delay4] \
                                        [list $delay1 [shape2list $rfsh2] $delay3 $delay4]]
            set rfsh_combined [list2shape [shape_add [list $rfsh_combined] \
                                                    [list $delay1 $delay2 $delay3 [shape2list $rfsh3]]]]
        } else {
            set delay1_length_extra [expr $par(tw1)-$par(compression)]
            set delay1_length [expr $par(tw1)-$par(compression)-($par(compression)-$par(tau2))]
            set delay2_length [expr $par(tw2)-$par(compression)-($par(compression)-$par(tau2))]
            set delay4_length [expr $par(tw3)-$par(compression)]
            set delay1_extra [zero_ampl [expr $par(tw1)-$par(compression)]]
            set delay1 [zero_ampl $delay1_length]
            set delay2 [zero_ampl $delay2_length]
            set delay4 [zero_ampl $delay4_length]

            set par(seq_duration) [expr $delay1_length+$delay2_length+$par(tw3)+($par(compression)-$par(tau2))]
            set seq_duration1 [expr $delay2_length+$delay4_length+$par(tw1)]
            set seq_duration2 [expr $delay1_length+$delay4_length+$par(tw2)]
            set seq_duration3 [expr $delay1_length_extra+$delay2_length+$par(tw3)]
            
            set rfsh_combined [shape_add [list [shape2list $rfsh1] $delay2 $delay4] \
                                        [list $delay1 [shape2list $rfsh2] $delay4]]
            set rfsh_combined [list2shape [shape_add [list $rfsh_combined] \
                                                    [list $delay1_extra $delay2 [shape2list $rfsh3]]]]
        }








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

        if {$par(compression) < $par(tau2)} {
            set delay1_length [expr $par(tw1)-$par(compression)]
            set delay2_length [expr $par(tw2)-$par(compression)]
            set delay3_length [expr $par(tau2)-$par(compression)]
            set delay4_length [expr $par(tw3)-$par(compression)]
            set delay1 [zero_ampl $delay1_length]
            set delay2 [zero_ampl $delay2_length]
            set delay3 [zero_ampl $delay3_length]
            set delay4 [zero_ampl $delay4_length]

            set par(seq_duration) [expr $delay1_length+$delay2_length+$delay3_length+$par(tw3)]

            set rfsh_combined [shape_add [list [shape2list $rfsh1] $delay2 $delay3 $delay4]                                             [list $delay1 [shape2list $rfsh2] $delay3 $delay4]]
            set rfsh_combined [list2shape [shape_add [list $rfsh_combined]                                                         [list $delay1 $delay2 $delay3 [shape2list $rfsh3]]]]
        } else {
            set delay1_length_extra [expr $par(tw1)-$par(compression)]
            set delay1_length [expr $par(tw1)-$par(compression)-($par(compression)-$par(tau2))]
            set delay2_length [expr $par(tw2)-$par(compression)-($par(compression)-$par(tau2))]
            set delay4_length [expr $par(tw3)-$par(compression)]
            set delay1_extra [zero_ampl [expr $par(tw1)-$par(compression)]]
            set delay1 [zero_ampl $delay1_length]
            set delay2 [zero_ampl $delay2_length]
            set delay4 [zero_ampl $delay4_length]

            set par(seq_duration) [expr $delay1_length+$delay2_length+$par(tw3)+($par(compression)-$par(tau2))]
            set seq_duration1 [expr $delay2_length+$delay4_length+$par(tw1)]
            set seq_duration2 [expr $delay1_length+$delay4_length+$par(tw2)]
            set seq_duration3 [expr $delay1_length_extra+$delay2_length+$par(tw3)]
            
            set rfsh_combined [shape_add [list [shape2list $rfsh1] $delay2 $delay4]                                             [list $delay1 [shape2list $rfsh2] $delay4]]
            set rfsh_combined [list2shape [shape_add [list $rfsh_combined]                                                         [list $delay1_extra $delay2 [shape2list $rfsh3]]]]
        }

        printwave $rfsh_shape1 1
        printwave $rfsh_shape2 2
        printwave $rfsh_shape3 3
        printwave [shape2list $rfsh_combined] _combined

        save_shape $rfsh1 $par(filename).simpson1
        save_shape $rfsh2 $par(filename).simpson2
        save_shape $rfsh3 $par(filename).simpson3
        save_shape $rfsh_combined $par(filename).simpson_combined
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
        lappend results_sum [format "%s %s %s %s" $offs [expr $re*2.0] [expr $im*2.0] [expr $xyproj*W2.0]]
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
    } elseif {[string equal $shape_type "tanhpulse"]} {
        set shape [list2shape [tanhpulse $dur $rfmax $sweepwidth $var1 $var2 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
    } elseif {[string equal $shape_type "hspulse"]} {
        set shape [list2shape [hspulse $dur $rfmax $sweepwidth $var1 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
    } elseif {[string equal $shape_type "caWURST"]} {
        set shape [list2shape [cawurst $dur $rfmax $sweepwidth $var1 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
    } elseif {[string equal $shape_type "supergaussian"]} {
        set shape [list2shape [supergaussian $dur $rfmax $sweepwidth $var1 $var2 -phase $phase -stepsize $stepsize -filename_phasecorrect $options(-filename_phasecorrect)]]
    }

    return $shape
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