#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:42:31 2019
@author: Max Busskamp
"""
import os
from subprocess import call
from subprocess import run
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy import interpolate
import PySimpleGUI as sg
import glob
from datetime import datetime
import base64

plt.rc('lines', linewidth=1)
plt.rc('axes', titlesize=18, labelsize=12)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)

plt.rcParams['figure.dpi'] = 140

sg.theme('DarkBlue')


def create_simpson():  # Write simpson input files
    simpson_inputfile = open('phasecorrection_liquid.tcl', 'w')
    simpson_input = """spinsys {
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

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh1 [list2shape [wurst $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh1 [list2shape [tanhpulse $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh1 [list2shape [hspulse $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh1 [list2shape [cawurst $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh1 [list2shape [supergaussian $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1)  $par(stepsize)]]
            }


            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh2 [list2shape [wurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh2 [list2shape [tanhpulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh2 [list2shape [hspulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh2 [list2shape [cawurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh2 [list2shape [supergaussian $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2)  $par(stepsize)]]
            }


            # Set third WURST pulse (refocussing)
            set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
            set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh3 [list2shape [wurst $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(ph3) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh3 [list2shape [tanhpulse $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(ph3) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh3 [list2shape [hspulse $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(ph3) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh3 [list2shape [cawurst $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(ph3) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh3 [list2shape [supergaussian $par(tw3) $par(Delta3) $par(rf3) $par(var31) $par(var32) $par(ph3) $par(stepsize)]]
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

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh1 [list2shape [wurst $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh1 [list2shape [tanhpulse $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh1 [list2shape [hspulse $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh1 [list2shape [cawurst $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh1 [list2shape [supergaussian $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]
            }


            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh2 [list2shape [wurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh2 [list2shape [tanhpulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh2 [list2shape [hspulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh2 [list2shape [cawurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(ph2) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh2 [list2shape [supergaussian $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(ph2) $par(stepsize)]]
            }
        } elseif {[string equal $par(type) "SHAPEsingle"]} {
            # Set first WURST pulse (excitation)
            set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
            set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh1 [list2shape [wurst $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh1 [list2shape [cawurst $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh1 [list2shape [tanhpulse $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh1 [list2shape [supergaussian $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(ph1) $par(stepsize)]]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh1 [list2shape [hspulse $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(ph1) $par(stepsize)]]
            }
        } elseif {[string equal $par(type) "create_shapes_SHAPEsingle"]} {
            # Set first WURST pulse (excitation)
            set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
            set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh1 [list2shape [wurst_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [wurst_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh1 [list2shape [tanhpulse_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [tanhpulse_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh1 [list2shape [hspulse_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [hspulse_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh1 [list2shape [cawurst_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [cawurst_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh1 [list2shape [supergaussian_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [supergaussian_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]
            }


            printwave $rfsh_shape1 1
            save_shape $rfsh1 $par(filename).simpson1
        } elseif {[string equal $par(type) "create_shapes_ABSTRUSE"]} {
            # Set first WURST pulse (excitation)
            set par(sweep_rate1) [expr ($par(Delta1)*1e3)/($par(tw1)*1e-6)]
            set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh1 [list2shape [wurst_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [wurst_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh1 [list2shape [tanhpulse_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [tanhpulse_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh1 [list2shape [hspulse_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [hspulse_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh1 [list2shape [cawurst_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [cawurst_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh1 [list2shape [supergaussian_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [supergaussian_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]
            }


            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh2 [list2shape [wurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(phaseoff2)]]
                set rfsh_shape2 [wurst $par(tw2) $par(Delta2) 100 $par(var21) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh2 [list2shape [tanhpulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2)]]
                set rfsh_shape2 [tanhpulse $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh2 [list2shape [hspulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(phaseoff2)]]
                set rfsh_shape2 [hspulse $par(tw2) $par(Delta2) 100 $par(var21) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh2 [list2shape [cawurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(phaseoff2)]]
                set rfsh_shape2 [cawurst $par(tw2) $par(Delta2) 100 $par(var21) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh2 [list2shape [supergaussian $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2)]]
                set rfsh_shape2 [supergaussian $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2)]
            }

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

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh1 [list2shape [wurst_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [wurst_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh1 [list2shape [tanhpulse_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [tanhpulse_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh1 [list2shape [hspulse_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [hspulse_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh1 [list2shape [cawurst_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [cawurst_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(filename_phasecorrect) $par(phaseoff1)]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh1 [list2shape [supergaussian_phasecorr $par(tw1) $par(Delta1) $par(rf1) $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]]
                set rfsh_shape1 [supergaussian_phasecorr $par(tw1) $par(Delta1) 100 $par(var11) $par(var12) $par(filename_phasecorrect) $par(phaseoff1)]
            }


            # Set second WURST pulse (refocussing)
            set par(sweep_rate2) [expr ($par(Delta2)*1e3)/($par(tw2)*1e-6)]
            set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]

            if {[string equal $par(shape_type) "WURST"]} {
                set rfsh2 [list2shape [wurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(phaseoff2)]]
                set rfsh_shape2 [wurst $par(tw2) $par(Delta2) 100 $par(var21) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "tanhpulse"]} {
                set rfsh2 [list2shape [tanhpulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2)]]
                set rfsh_shape2 [tanhpulse $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "hspulse"]} {
                set rfsh2 [list2shape [hspulse $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(phaseoff2)]]
                set rfsh_shape2 [hspulse $par(tw2) $par(Delta2) 100 $par(var21) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "caWURST"]} {
                set rfsh2 [list2shape [cawurst $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(phaseoff2)]]
                set rfsh_shape2 [cawurst $par(tw2) $par(Delta2) 100 $par(var21) $par(phaseoff2)]
            } elseif {[string equal $par(shape_type) "supergaussian"]} {
                set rfsh2 [list2shape [supergaussian $par(tw2) $par(Delta2) $par(rf2) $par(var21) $par(var22) $par(phaseoff2)]]
                set rfsh_shape2 [supergaussian $par(tw2) $par(Delta2) 100 $par(var21) $par(var22) $par(phaseoff2)]
            }


            # Set third WURST pulse (refocussing)
            set par(sweep_rate3) [expr ($par(Delta3)*1e3)/($par(tw3)*1e-6)]
            set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]

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
        set phasecorr_list  [split $phasecorr "\n"]

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
        set phasecorr_list  [split $phasecorr "\n"]
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
        set phasecorr_list  [split $phasecorr "\n"]
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
        set phasecorr_list  [split $phasecorr "\n"]
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
        set phasecorr_list  [split $phasecorr "\n"]
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
        set phasecorr_list  [split $phasecorr "\n"]
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
        set phasecorr_list  [split $phasecorr "\n"]
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
    }"""
    simpson_inputfile = open('phasecorrection_liquid.tcl', 'w')
    simpson_inputfile.write(simpson_input)
    simpson_inputfile.close()


def simulate_sequence(exp_type, shape_type):  # Start simulation with set parameter
    global delta
    global tw1
    global tw2
    global tw3
    global rffactor1
    global rffactor2
    global rffactor3
    global tau1
    global tau2
    global tau3
    global var11
    global var12
    global var13
    global ss_offset
    global var21
    global var22
    global var23
    global filename
    global var31
    global var32
    global var33
    global phaseoff1
    global phaseoff2
    global phaseoff3

    run(['simpson',
            'phasecorrection_liquid.tcl',
            filename,
            tw1,
            tw2,
            tw3,
            exp_type,
            delta1,
            delta2,
            delta3,
            var11,
            var12,
            var13,
            rffactor1,
            rffactor2,
            rffactor3,
            tau1,
            tau2,
            tau3,
            'none',
            ss_offset,
            var21,
            var22,
            var23,
            shape_type,
            var31,
            var32,
            var33,
            phaseoff1,
            phaseoff2,
            phaseoff3])

    # filename_phasecorr = glob.glob('*.out')[0]
    filename_phasecorr = filename + '.out'

    fig = plt.figure(figsize=(6.5, 6))
    ax1 = plt.subplot2grid((2, 1), (0, 0))
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    fig.tight_layout()

    dat_phasecorr = np.genfromtxt(filename_phasecorr, delimiter=' ')

    phase = np.zeros((len(dat_phasecorr), 2))
    for i in range(len(dat_phasecorr)):
        phase[i, 0] = dat_phasecorr[i, 0]
        phase[i, 1] = np.angle(complex(dat_phasecorr[i, 1], dat_phasecorr[i, 2]))
    phase[:,1] = np.rad2deg(np.unwrap(phase[:,1]))

    ###############################################################################
    # Interpolation
    ###############################################################################

    np_tw = round(float(tw1)/0.05)
    pulselength = np.linspace(phase[0, 0], phase[-1, 0], np_tw)


    interpol = interpolate.interp1d(phase[:, 0], phase[:, 1])
    interpol = interpol(pulselength)


    ###############################################################################
    # Plotting
    ###############################################################################
    np.savetxt(filename_phasecorr.replace('.out', '.phasecorr'), interpol, delimiter=' ')

    ax1.plot(dat_phasecorr[:, 0]/1000, dat_phasecorr[:, 1], label=r'Real')
    ax1.plot(dat_phasecorr[:, 0]/1000, dat_phasecorr[:, 2], label=r'Imag')
    ax1.plot(dat_phasecorr[:, 0]/1000, np.sqrt(dat_phasecorr[:, 1] ** 2 + dat_phasecorr[:, 2] ** 2), label=r'Magnitude')
    ax1.legend(fontsize=6)
    ax1.invert_xaxis()
    ax1.set_ylabel('Magnetization / a.u.')
    ax2.plot(phase[:, 0]/1000, phase[:, 1], ls='-', c='k', label='Phase')
    ax2.plot(pulselength/1000, interpol, ".", ms=0.5, c='r', label='Interpolation')
    ax2.legend(fontsize=6)
    ax2.invert_xaxis()
    ax2.set_xlabel('Frequency Offset / kHz')
    ax2.set_ylabel('Phase / degree')
    plt.tight_layout()
    plt.savefig(filename_phasecorr.replace('.out', '.png'))

    ##############################################################################
    # Create corrected Shape Files
    ##############################################################################
    run(['simpson',
        'phasecorrection_liquid.tcl',
        filename,
        tw1,
        tw2,
        tw3,
        'create_shapes_'+exp_type,
        delta1,
        delta2,
        delta3,
        var11,
        var12,
        var13,
        rffactor1,
        rffactor2,
        rffactor3,
        tau1,
        tau2,
        tau3,
        filename_phasecorr.replace('.out', '.phasecorr'),
        ss_offset,
        var21,
        var22,
        var23,
        shape_type,
        var31,
        var32,
        var33,
        phaseoff1,
        phaseoff2,
        phaseoff3])

    run(['rm', '-f', 'phasecorrection_liquid.tcl'])
    # run(['rm', '-f', filename_phasecorr])
    # run(['rm', '-f', glob.glob('*create_shapes_*.out')[0]])

    sg.popup('Simulation finished! Close to show Plot!')

    plt.show()


def sequence_diagram(exp_type):  # Generate base64 image of pulsesequence
    if(exp_type == 'CHORUS') or (exp_type == 'CHORUS_cycled'):
        diagram = b'iVBORw0KGgoAAAANSUhEUgAABQEAAADwCAYAAACwjKJ3AAAACXBIWXMAAC4jAAAuIwF4pT92AAAF7GlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUgWE1QIENvcmUgNS42LWMxNDIgNzkuMTYwOTI0LCAyMDE3LzA3LzEzLTAxOjA2OjM5ICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIgeG1sbnM6cGhvdG9zaG9wPSJodHRwOi8vbnMuYWRvYmUuY29tL3Bob3Rvc2hvcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RFdnQ9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZUV2ZW50IyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ0MgMjAxOCAoV2luZG93cykiIHhtcDpDcmVhdGVEYXRlPSIyMDIwLTA1LTE4VDE1OjI4OjUxKzAyOjAwIiB4bXA6TW9kaWZ5RGF0ZT0iMjAyMC0wNS0xOFQxNjoxOToxNyswMjowMCIgeG1wOk1ldGFkYXRhRGF0ZT0iMjAyMC0wNS0xOFQxNjoxOToxNyswMjowMCIgZGM6Zm9ybWF0PSJpbWFnZS9wbmciIHBob3Rvc2hvcDpDb2xvck1vZGU9IjMiIHBob3Rvc2hvcDpJQ0NQcm9maWxlPSJzUkdCIElFQzYxOTY2LTIuMSIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDo4MGMxMTExMy01NjFjLWE0NGItYTAxMC05OWZmY2UyNTlkYTMiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6ZWJkOTZkM2EtYWYyMi1lYzQ3LTg1ZWUtN2JiNmUzNDk0MTcxIiB4bXBNTTpPcmlnaW5hbERvY3VtZW50SUQ9InhtcC5kaWQ6ZWJkOTZkM2EtYWYyMi1lYzQ3LTg1ZWUtN2JiNmUzNDk0MTcxIj4gPHhtcE1NOkhpc3Rvcnk+IDxyZGY6U2VxPiA8cmRmOmxpIHN0RXZ0OmFjdGlvbj0iY3JlYXRlZCIgc3RFdnQ6aW5zdGFuY2VJRD0ieG1wLmlpZDplYmQ5NmQzYS1hZjIyLWVjNDctODVlZS03YmI2ZTM0OTQxNzEiIHN0RXZ0OndoZW49IjIwMjAtMDUtMThUMTU6Mjg6NTErMDI6MDAiIHN0RXZ0OnNvZnR3YXJlQWdlbnQ9IkFkb2JlIFBob3Rvc2hvcCBDQyAyMDE4IChXaW5kb3dzKSIvPiA8cmRmOmxpIHN0RXZ0OmFjdGlvbj0ic2F2ZWQiIHN0RXZ0Omluc3RhbmNlSUQ9InhtcC5paWQ6ODBjMTExMTMtNTYxYy1hNDRiLWEwMTAtOTlmZmNlMjU5ZGEzIiBzdEV2dDp3aGVuPSIyMDIwLTA1LTE4VDE2OjE5OjE3KzAyOjAwIiBzdEV2dDpzb2Z0d2FyZUFnZW50PSJBZG9iZSBQaG90b3Nob3AgQ0MgMjAxOCAoV2luZG93cykiIHN0RXZ0OmNoYW5nZWQ9Ii8iLz4gPC9yZGY6U2VxPiA8L3htcE1NOkhpc3Rvcnk+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+T4/djAAAYT5JREFUeJzt3Xe8HFX9//HXvWmkkEIKIaGEUEMJHWkqXSkKiBSliPqlV6WLioBYaIIIIkWK0kQUBQEBAamKCAiC9F5CC0kgvdzfH58zv3N27szubJ0t7+fjsY89Mzt399x7d2dnPvM559PV09ODZLIFcCKwNjAK6Mr4c8cDZ9SpTyKdZBzwQ2AbYAzQP+PPPQxsUq9OiXSgocApwBeBpYCBGX/uPWDJenVKRNrGjsAxwOrYMXdW+wOX1qVHIiIibaJv3h1oEV8GrgP65N0RkQ41HvinuxeR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3Z/MsyMiUl9dPT09efdBRERERERERERE6kjDgUVERERERERERNqcgoAiIiIiIiIiIiJtTkFAERERERERERGRNqcgoIiIiIiIiIiISJurpjrwccA41z4RmF19d6SNTQCOcu3HgStz64mISO19C1jOtU8GpufYFxFJtgWwk2v/Abgvx76ISLIlgO+79kvA+Tn2pZNsDOzh2rcCd+TYFxGpo2qCgPsAa7j2KSgIKMWNA4507RtQEFBE2suewIaufSYKAoo0o/XxxyKvoCCgSDMahv+cPoCCgI2yJv7vPhUFAUXaloYDi4iIiIiIiIiItDkFAUVERERERERERNqcgoAiIiIiIiIiIiJtTkFAERERERERERGRNqcgoIiIiIiIiIiISJtTEFBERERERERERKTNKQgoIiIiIiIiIiLS5hQEFBERERERERERaXMKAoqIiIiIiIiIiLQ5BQFFRERERERERETaXFdPT0+lP9sXH0ScV5vuSBvrAvq59kJ3ExFpF/pOFGl+3dhnFWABsCjHvohIuv7ufhH2WZX60/5RpENUEwQUERERERERERGRFqDhwCIiIiIiIiIiIm1OQUAREREREREREZE2pyCgiIiIiIiIiACsCkwF5gPL5twXEamxvqU3SbUyMNC1/4sKPUhxg4CVXHsa8Fp+XZEW0QcYix18jAJGuvvRrj0UWBx7bw0GhmMFaMD2TYsFzzULmOvaC4EZwMdu/UzgI+AD4MPg/l3gTXevyVOllBWx9yHA02gic5FmNBoY59pvYft7aV79sWDECsBywPLAkvjjgMWx7/3hbvuZWGGmOdj3+IfAe9gx52vAy8AzwPRG/QJSkf7AJNf+BHgpx750kiWAZVx7JWCEa08AXs+jQyJSH9UEAW8E1nDtJbCTaJE0awMPuvYNwO75dUWayBgseLKSu18RC/otiwUAq9lHhUaU3iTVPOBtLCD4KvCCu73obtr3CcDVwIauvTQWYBCR5rIfcIZrHwWcl1tPJK4PsCawCbAxdty4CtCvjOfI+l3/KvAU8E/gIeARLIAozWE88IRrPwB8Or+udJQvA79y7d8H6ycA9zW8NyJSN7U6wRYRKWYQsBZ2gL8mdgFhMnYBodn1xw6AJgCbJTw+BTuZeNLdR7f5jemeiIhIS1oG2A74HLAVMKxBrzvB3b7glhdgwcC/ArcDj6MRANLZBgft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'
    elif(exp_type == 'double_echo') or (exp_type == 'double_echo_cycled') or (exp_type == 'double_echo_zerophase') or (exp_type == 'loadshape_double_echo'):
        diagram = b'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'
    elif(exp_type == 'double_chirp'):
        diagram = b'iVBORw0KGgoAAAANSUhEUgAAAwYAAADWCAYAAABrEJVmAAAACXBIWXMAAC4jAAAuIwF4pT92AAAFHGlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUgWE1QIENvcmUgNS42LWMxNDIgNzkuMTYwOTI0LCAyMDE3LzA3LzEzLTAxOjA2OjM5ICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIgeG1sbnM6cGhvdG9zaG9wPSJodHRwOi8vbnMuYWRvYmUuY29tL3Bob3Rvc2hvcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RFdnQ9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZUV2ZW50IyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ0MgMjAxOCAoV2luZG93cykiIHhtcDpDcmVhdGVEYXRlPSIyMDIwLTA1LTIwVDExOjQxOjM2KzAyOjAwIiB4bXA6TW9kaWZ5RGF0ZT0iMjAyMC0wNS0yMFQxMTo1MDo0NSswMjowMCIgeG1wOk1ldGFkYXRhRGF0ZT0iMjAyMC0wNS0yMFQxMTo1MDo0NSswMjowMCIgZGM6Zm9ybWF0PSJpbWFnZS9wbmciIHBob3Rvc2hvcDpDb2xvck1vZGU9IjMiIHBob3Rvc2hvcDpJQ0NQcm9maWxlPSJzUkdCIElFQzYxOTY2LTIuMSIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDoxZjU5YWI0Yi0xOTBjLWFlNDMtYjhjMi1mZDJmYjliNzQyZmUiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6MWY1OWFiNGItMTkwYy1hZTQzLWI4YzItZmQyZmI5Yjc0MmZlIiB4bXBNTTpPcmlnaW5hbERvY3VtZW50SUQ9InhtcC5kaWQ6MWY1OWFiNGItMTkwYy1hZTQzLWI4YzItZmQyZmI5Yjc0MmZlIj4gPHhtcE1NOkhpc3Rvcnk+IDxyZGY6U2VxPiA8cmRmOmxpIHN0RXZ0OmFjdGlvbj0iY3JlYXRlZCIgc3RFdnQ6aW5zdGFuY2VJRD0ieG1wLmlpZDoxZjU5YWI0Yi0xOTBjLWFlNDMtYjhjMi1mZDJmYjliNzQyZmUiIHN0RXZ0OndoZW49IjIwMjAtMDUtMjBUMTE6NDE6MzYrMDI6MDAiIHN0RXZ0OnNvZnR3YXJlQWdlbnQ9IkFkb2JlIFBob3Rvc2hvcCBDQyAyMDE4IChXaW5kb3dzKSIvPiA8L3JkZjpTZXE+IDwveG1wTU06SGlzdG9yeT4gPC9yZGY6RGVzY3JpcHRpb24+IDwvcmRmOlJERj4gPC94OnhtcG1ldGE+IDw/eHBhY2tldCBlbmQ9InIiPz5I64SMAABRp0lEQVR4nO3deXgsR3kv/u9bM9J0j5Yjb8cYMNjYZjMYg8ExNgZicFjDvhMuEAhLwMANv5sFQiCXkAA3CYGwLzfXEMJmEiABAyFsxgZsDDYxqw02YPDuc7RNVWtm6v39oapRaTQzks6R1Brp+3me80xXT3VP6Yxm1G+/tQiI1kFVcwCvXaVaC8AMgJsBXA7gchFprfN17gHgVQCOFJGHHEhbiQ6Wqj4WwP1XqTaLxd/3XwK4SERuWudrZAB+H8CfAHiGiFx4IG2lA6eqJwN42hqrT2PxPf8VgJ8C+ImI+IN47SMBvBLA80TksAM9DxER0ZZTVVHVQ1T1JFX9jK70G1X9nqrenOzbr6r/W1Un13D+k1X1k6raDsd+ZSt+LqJeVLWuqrdT1Wer6mzX77pV1StU9aeq2kz2f1VVH7DGc79CVX+dHHvGVvxctJyqjqrqXlU9W1Wv6nqfvapeqapfU9UvqOqlqjqTPD+tqv+83vcu/F79g6rOh/MUm/XzERERbTpVParrD+gDk+eMqj5CVa9Jnv++qt6mz7kyVf1/qnqeql6bHMPAgLYFVf1A8nt5vqpK8txRqvqu5Pm2qr5wwLmepKqfCheanoHB9qGqz03ej1lV3dOjzoiqnqWqn+v6DvyIru0GyOtU9d908SZKxMCAiIiGm6o2kj9sp/V4/naqemtS51JVHelzrtHweF8GBrTdqOqrk9/Lz/Sp8/dJHa+qj+lTbzTZ/hQDg+1DVc9M3o/pNdR/ki7d9ddwsT+1yjHxu25El7KrDAyIqHSm7AbQ0JsZ9KSI/BrAG5Nd9wHwlD51F8LmTzamaUQbanYNdf4cS58JAfA3vSolv+vA4tgE2j4Gfqd1E5HzADwcQHxPTwbwEU0ySj2OWQiPTQDXHVgziYg2HgMDOli6hjr/2lV+2Cr15wAc8GA+ok2y6u+6iDQAnJ/suruq3mGVw9Z1IUqbbi3facuIyAUA/izZ9XAAz17j4Xz/iWjbYGBAW+EaAO2k3HOcQSQiCmB+MxtEtIl+3lU+cpX6a8lE0Pb3j1ie/fkLVa2s4Ti+/0S0bTAwoK0gWP67tm8Nx6xrelOibaT7e3W13/f2Ks/TEAjdgt6Z7DoWwIPWcCi/64ho22BgQFvhZCwGB9G3SmoH0Va4d7J9E1ZmEGjn+q+u8iNKaQUR0QFiYEBb4RXJ9iyAD5fUDqJNpar3BJAuyPeeg1n8iobOZQCaSfnkcppBRHRgGBjQplLV/w/AM2IRwItE5MYSm0S0KVT1GADnAYj9yi9Bn1mJaGcKK7zfmuxabXwJEdG2Ui27AbSj/JWqXgzgegBjAB4DIK5t8GsAfygiPed/Jxoyp6rqewBcjcUs2IkAngVgHIszan0YwEvDLEW0u0xjKSDIy2wIEdF6MTCgjfRpLF4UHQ/gEAA/BPAlAN8G8MWuuduJhtmVWOxPfncAx2AxIPh/AH4E4PMiwnEFu1e6UvItpbWCiOgAMDCgjXSJiHBgMe0G+0Tk42U3grYXVTVYvCkS3VBWW4iIDgTHGBAREW2MkwCMJuULy2oIEdGBYGBARES0Mc7qKn+xlFYQER0gBgZEREQHKaxy/KJk11dF5LKSmkNEdEAYGBARER285wE4IWx7AK8psS1ERAeEgQEdLOmzvVHn3chzEh2Mzfpd7/caVI51vweqei8Ab0l2/bWIfGOdr8f3nohKx8CADtZYsj2xEScMM3tMbuQ5iTZA+rs+2bfW+tWSbf6+l298PZVV9cEAvgKgHnZ9CMDr1nGKOFh5RFWz9bw2EdFGY2BAB0xV92L5xdLxG3TqQ7D0uzm1QeckOljHptshgN0IU322qRzp95ioas87+ap6V1X9ABbXszgEwD4ArwDwbBFpr+P10ulNp9bXVCIiohKpqqjqIap6iqp+Spf7hao+QlX3qmpt9bOtOPeUqp6oqu9KzulV9RxVPUFVuYoobSlVravq0ar6P1R1ruv3/S3h93LP6mdacd4RVb29qj5cVa9NzvkNVT1NVY/ajJ+HelPV0fC99TBVvarrfZ5R1W+p6sdV9QOq+u+q+tPk++kKVX21qq45i6Squaoep6ovUNVm8lrvVtV7qurUJv64RER9sU8jrUu4OH/tGqp+UkQuWee5X4TFVWT7+aqIfH495yQ6GKr6WAD3X6XaPhF50zrPexcAz12l2qtExK/nvHRgVPVkAE9bY/VpADcDuA7At0XkpgN4vUcAeNCAKteIyLvXe14iIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIioqHlnLuLtfYia+1Fzrk3lt0eItp4RVHcI37OrbV/VXZ7aPMVRXFS8p7/Zdntoc03Nzd3cvKer2ViCiLaxqplvGi73R43xtwfALz3657RgYi2v3a7PSEi9wcAEfl12e2hzee9n8TSLE6/LLMttDUqlcoeLL3nV5fZFiI6eFzgjIiIiIiIGBgQEREREREDAyIiIiIiAgMDIiIiIiICAwMiIiIiIgIDAyIiIiIiAgMDIiIiIiICAwMiIiIiIgIDAyIiIiIiQkkrH1er1YaqXgYAIvLzMtpARJurUqnMx885uCLqrtBut+eq1eplAKCqfM93Ae/9bKVSuQwAVPWacltDRERERERERERERERERERERERERERERERERERERERERBtCynjRmZmZw0dHR58QitdkWfbFMtpBRJtndnZ278jIyONC8edZln2pzPbQ5pubmzuyWq0+FgBU9ao8z79cdptoc83Pzx9VqVR+FwBU9co8z79SdpuIaMjMz8+fYq1Va602Go1Pl90eItp4jUbj/vFz7pz7RNntoc1nrX1AfM+ttR8tuz20+ay1D0re8w+X3R4iOjhc4IyIiIiIiBgYEBERERERAwMiIiIiIgIDAyIiIiIiAgMDIiIiIiICAwMiIiIiIgIDAyIiIiIiAgMDIiIiIiICAwMiIiIiIgIDAyIiIiIiAiBlvKiq5s1m884A0G63p/M8v6aMdhDR5lHVerPZPAEA2u32/jzPf1F2m2hzqepYs9k8HgDa7fa+PM9/WXabaHPxPSciIiIiIiIiIiIiIiIiIiIiIiIiIiIiom6qWrXWfr3RaFzWaDTuX3Z7iGj9qmW8qKoKgNFQbItIq4x2ENHm2U2fc1WdKIriwap6XwDHArgNFn/2pqrOA7jGGPNDAF/JsuzKMtu6mXbTe06L0vfcWrtXRM4UEQC4M4Bvltk2Ilq/UgKDRqNxH2PMdwBAVT8D4LFltIOINo+19jQRuQgAROQ8AE8uuUkbam5u7shKpfI8AA93zp0GYKRXvXCRBFUFAFhrrwHwJQAfzvP8q1vR1q3inDsDwAWh+DEATyuxObQFnHMPBPBVABCR85OnjimjPUR0cEoJDIiIhpVz7njv/Tki8gcA8gM4xTEAng/g+c65y1T1LVmW/QvvrtMOUEu2jymrEUR04BgYEBGtwczMzOEjIyN/r6q/JzENsGRaRL4C4CsAfiIi16nqHIAR7/1hInKsqp4mIg8FcNd4kKqeDOBc59xrnHMvzrLsS1v2AxFtvDQwOKq0VhDRAWNgQES0ikaj8XQReSuAI5LdLQCfBPCOLMu+ucod/4sAfDic62gAzxCRcwDcLjx/vKp+0Vp77sLCwiv37Nlz6yb8GESbSlXzJGY+kGwaEZXMlN0AIqLtSlX3NBqNT4vIv2ApKFhQ1bcCOD7P86fleX7BeroB1ev1X9Xr9TdlWXYnAM8G8NPwlAB4zujo6A+dc2dv6A9CtAVEJB1nUy+tIUR0wBgYEBH14Jy7k3PuIhF5TNwnIpd57+9fr9dfkef5Lw7m/CKykOf5B7MsuxeAvwSwEJ46UlXPd86dczDnJypBJdlmxoBoCDEwICLqYq09U1W/DeDuYVehqq+s1Wr3HRsb++5GvpaIuDzPX1epVE4VkcvC7oqqvq3RaLxVVSuDjifaRtLfVWYMiIYQAwMiokRRFI/F4nSih4ddN6rqQ+r1+t+LSHuzXnd0dPTyWq12hoh8Mu4TkZcVRXGeqvacCpVom2HGgGjIMTAgIgqstWd57z+KpUW6rgBwWr1ev3ArXl9EGrVa7clY7FqkAKCqjyuK4iPMHNB2JyLpNQUzBkRDiIEBERGARqNxOoDPAMjCrguyLDs9z/Ort7IdIqJ5nr8OwAuwFBw80Tn3rrDKLNF2xYwB0ZBjYEBEu97CwsK9AHwWwFjY9Z0syx4tIrNltSnP8/eLyMuSXX9QFMXfltUeotV0ZbVGVJVTohMNmVI+tCJyg4i8JRSvKKMNRLS5ROQ38XOuqpeV3Jy+9u/ff0i73f43EZkKu65oNpsPz/N8psx2AUCWZW9vNBojIvL3AKCqf2St/UGe5/+37Lb18evkPf9e2Y2hzScivwIQ3/OHArh98nQdQOmfIyIiIqJVqaqx1p5vrdXw7+fz8/PbbsVW59ybkjY25ubmTi67TUTdrLVfS35PdW5u7siy20RE68OuRES0aznnXgfg4bHovX/y2NjYdSU2qadarfZnAM4PxbxSqXxqZmbmsDLbRNRDlhYqlQoHIBMNGQYGRLQrOeceCeDVya4/HBsbu7Ss9gwiIr7ZbD4LQFxU7Y4jIyPncjAybTPLBhwbYzgAmWjIMDAgol1n3759U6r6Xix9B743z/N/KrNNq5mcnLyl3W4/DoANux5VFMULSmwSUbdaWmi1WgwMiIYMAwMi2nWyLPsHALcDABG5LMuylw0+YnsYHx+/TEReGcuq+uZGo3F0mW0iSiwLBERktF9FItqeSpmVyFp7LIC/BgBV/U69Xv+7MtpBRJvHOXeCqv5vABCRb2VZ9tay2wQA1tqHAPgfodhqt9vPF5GizDatR61We7dz7vEAzgYwKSLvAfDIkpsFAHDO3UVVXwcAInJRlmX/WHKTaJMVRXE37/1fAICIHKKqnecYGBANn1IyBt77QwE8Lfx7YBltIKLN5b0/HEuf8weU3BwAgKqOAXgvgNg3/83bdVxBPyKiWFz8bC7seoS19lklNqlDVY9AeM9V9Yyy20Obz3u/F0vv+bLBxwwMiIYPuxIR0a5hrX0DgDuF4k+yLHt9me05UHmeX6Oqr0l2vWVmZubw0hpEtKj7moKBAdGQYWBARLuCc+7OIvKHoegBPF9EXJltOhh5nv8jgG+G4mHVavU1g+oTbYFl1xTtdnukrIYQ0YFhYEBEu8WbAMQLlQ/mef6NMhtzsESkXalUXozFIAci8mLn3AklN4uog12JiIYPAwMi2vGstQ9U1cfFoqr+RZnt2Sijo6OXA/hwKI4gTOpAtB0wMCAaPgwMiGhHC4uAvTHZ9bf1ev1XZbVno6nqnwFohO0nWWu3xUBvInCMAdHQYWBARDuatfapAO4fijdkWfZ/ymzPRqvX678G8LZk15vKagtRUADMGBANIwYGRLRjqaqIyJ/Hsoj8hYjMltmmzZBl2V8DuCEUT3fOnV1me2jXmw+PDAyIhgwDAyLasay1jwdwYij+vFar/d8y27NZRGRWVf82llX1VWW2h3a9eQDw3nNWIqIhw8CAiHYsEfmzZPvNItIqsz2bKc/zdwG4ORQfzLEGVKJ5gF2JiIZRKYHByMhIC8CtAG4VkZky2kBEm0tVmwifc1Xd8u47zrmHA7hvKF5fq9XO3eo2bCURmQfwjmTXn251G9L3HMCO67JFK8X3XESm477wuwiwKxERERFtB9bar1tr1VqrjUbjj8puz1aYnp4+1Fo7E3/u+fn5U8puE+0ORVGcFH/vrLVfC49/VXa7iGh92JWIiHacRqNxBoAzQ/GWPM/fU2Z7tsqePXtuFZF3x7Ix5n+V2R7aPVqtVpodYFcioiHFwICIdhwReVlSfFvStWHHa7fbbwHgQvGJjUbjdmW2h3YHERlJtuPg49pqxxVFcU/n3MOstXfczPYR0dowMCCiHWV+fv4oAI8PxaLVar17UP2dZmxs7DoAHwvFqoi8sMz20O6QBAZtVXVh36oZA+/9F1X18wCevpntI6K1YWBARDuKMeZFAOJFykcnJiZuLLM9ZfDepwuevVBVV71zS3SQ4mduIfxL9/WkqqMA9obiXTepXUS0DtUyXlRVJ4uiOB0A2u32DWNjY98rox1EtHn27ds3lef5aQDQarWuHx8fv2yzX1NVR5xzz49l7/07BtXfqcbGxr5rrb0YwKkA9lprnwDgI5v9uul7boy5bnR09PLNfk0q1/79+w/Jsuy3VDXOANbEUmAwMGPgnLstlm5Q3m2TmkhE61BKYNBoNE4wxpwPACLyGQCPLaMdRLR5arXa3VT1fACoVqvnAXjyZr+mc+7JAG4bit8eGxu7ZLNfcxt7BxYDA4jIS7AFgUGWZfeI73m73f4YgKdt9mtSuWq12knxPQ+aqtoUEWD16Upvn2wfueGNI6J1Y1ciItpJ/jDZfmdprdgGsiz7OICbQvGM+fn5e5fZHto1msaYBWD1MQaqepukOL6prSKiNWFgQEQ7QlEUdwdwRijeFC6Mdy0RcQDen5SfV2JzaPdYCIuewXs/sFeCiKTBAAMDom2AgQER7Qiq+py4LSLnhgvj3e59ABQAROTpqpqV3B7a+ZoAWgAgIqsFBmNJsaaqAwcrE9HmY2BARENPVauq+nuxLCIfKrM920We51cD+HooHmqtfUyZ7aFdoRMYYJVZibz3aWCA6elpZg2ISsbAgIiGXlEUDwNwVCh+p1arfb/M9mwz58YNEXl2mQ2hXaEZ/gGrTHAiIvW0XKvVJjarUUS0NgwMiGgn6Fzwisi5gyruNmGsxVwoPqzRaNx+UH2igyEiTVVdU8agOzAwxjBjQFQyBgZENNSmp6cPVdXYRWZhYWHho6U2aJsRkXkAnwzFCoBnltgc2uFUdQFrzBh0dyVqt9sMDIhKxsCAiIZarVZ7GoAaAIjIv09OTt5ccpO2o3+KG+xORJusZYxpAoCIrCtjICLsSkRUMgYGRDTUVPXpSZHdiHrIsuzrAK4OxbvNzc2dXGJzaGfrDD5ew3SlyzIG7EpEVD4GBkQ0tEJ/+dMBQET21Wq1L5TcpG1JRBRAp4tVpVJ5aonNoZ2tHdcxWC1joKrLMgbdXYuIaOsNjOY3S71ev6ooiscAgPf++jLaQESba2Fh4cdZlj0GANrt9m826WWegnCDQ1U/KSILm/Q6Q6/dbn+8Uqn8WSg+VVVfFQKGDbOwsPDDWq0W3/Nfb+S5aXtaWFj471qt9hhVfSaAp2L5dKWrXWN0BwJcZ4OoZKUEBiIyDeDfy3htItoaU1NT+7DJn3MReWqy/bHNfK1hNz4+fpm19kcA7gbg2EajcQqA72zka+zZs+dW8Lt9V4nvubX2/gCgqi0RWdOsROgKDESktglNJKJ1YFciIhpK1to7ALhfKN5Uq9W+WmJzhsV5caNSqTylzIbQzhJXOTbGdLoSYfWbj3la8N4zY0BUMgYGRDSUVPVpACQUz0vuUlIfxph/iduq+jRVlUH1idYhBgFrXvkYXV2HmDEgKh8DAyIaSsaYdAAtuxGtQa1W+zGAH4Ti0dba3yqzPbRzxBmIVLUzXSlWzxh0BwLMGBCVjIEBEQ0da+0xqnqfULwuy7JvlNqg4fLxuGGMeWKZDaGdI3YlAtDG0gJnq2UMRrvOwcCAqGSlBAZFUZxkrb0h/OO840Q70Pz8/CnJ5/x9G3luVX1cUvyUiLQ38vw7mTEmroLc/f940BqNxqnJe/7ujTw3bU+NRuP+1tobROQ5YdeaZyVS1ZgxcOFxtF9dItoapQQGrVZrBMBeAHtVdaqMNhDR5hKRUYTPuYhMbfC5H5tsf3ojz73T1Wq1HwC4MhSPL4ri7ht17vQ9BzC1Ueel7Su+56oaBxK3VHVNYwySMQUzAKCqzBgQlYxdiYhoqExPTx8K4AGhOMvZiNZPRDpTinrvHzuoLtE6tdYxK1HMEMyERw4+JioZAwMiGiqjo6OPxtIFx+dEpCizPcNIVdMsCwMD2jAhWxAzBqKqlT71qgAqACAiMTBgxoCoZAwMiGiopN2Iui5waY2yLLsQwE2heGqj0bh9me2hncMYk2YMgP7diTrjCVR1GuB0pUTbAQMDIhoaqlpT1bNDsbmwsPD5Uhs0pMJg7c/GojHmUWW2h3aUdrVaTdcU6RkYTE9Pd4IAVY1jDPJedVNzc3NHzs/P38dae+xBt5SIVmBgQERDoyiKhwCYCMWvT01N7SuzPcPMe9/JtqgquxPRhghdiToZg+np6Z7jDKrVaidjkHQlWjVjUKlUPm6MuRTAaw6yqUTUAwMDIhoa3vtHxm3ORnRw6vX6FwHYUPxtVa2X2R7aMZrJrESoVCo9A4NKpbIiY4C1DT4+Pjze74BbSER9MTAgoqEhIg9PiueX1pAdQEQaAL4eillRFA8usTm0c7TTMQbGmJ5didLxBCIyHTYHDj4OwetRoXg3VeU1DNEG44eKiIaCc+4uAI4LxSuzLLuqzPbsBCLSCa68948osy20Y6TrGAwKDNLBx2vKGCwsLBwLQEKxMj09vecg20pEXRgYENFQSC9cVZXZgo3xubghIhyATAdNVVve+07GQER6diVqtVqdIMAYMxs2B6583G63p9JyrVZjYEC0wVZbfGRTVCqVfQA+CQCqenEZbSCizWWMuQXhc+69/+bBnk9EOoGBMYaBwQbIsuxKa+1VWOy3faxz7oQsy65c7bh+jDE3q2p8z7+1Ue2k7csYc5OqflJVHwzgMAAt732rUuksX9DzOiPtSqSqc2FzYGBgjJlU1U65O1AgooNXSmCQZdnPATypjNcmoq2RZdlPsUGfc1WtO+ceGIq2Vqt9feABtGaq+nkReWkoPgLAAQcGtVrtx+B3+65Sq9V+COBJ1tqLARxmjGlNTEw0nXMAABFZrStRU1WdiACrBAbe+4lQDwBQqVSmDv4nIKIUuxIR0bZXFMVZWBqY+JUwcJY2QJp9UVWOM6ADknQZagFoJ/t7rnwsIvHzvICl6U37LYYWj5lMy8YYdiUi2mAMDIho2/Ped2YjEhEuaraBarXaVwG4UHwwpy2lA6GqaWDQGXzcbDb79UyI2YHCGFOEcwzMGKjqRFr23k8dUGOJqC8GBkS07YnI2UmRgcEGCtmXr4ViVhTFmWW2h4ZWFegscNYJDPoNPvbexyCgkzHo1+0oOdeervLUgTeXiHphYEBE25q19g4A7hyKvziYwbHUm6p+KSk+pLSG0DDrZAxEpA1Au/b3q99U1YWwPTBjICLjaVlVJ/vVJaIDU8rg40ajcTtjzMsAQFV/kOf5B8toBxFtHmvtHUTkJQCgqt/P8/zDB3iqhybb/3nwLaNu3vsvJbPIPHRQ3UGstceIyIvDOS+r1+sf2Yj20fZlrT1WRF6kqnvDrpgtaAOoDhhjkA4+boZBxSOqKiKivY7pEQiw2xvRBislY6Cqt1HVPw7/nlhGG4hoc6nq7eLnXEQedxCn6tzBVtX/OviWUbexsbHLAdwIAKp68uzs7BEHeKrbJ+/5YzeuhbSN3UFV/xjAHgBQ1TjwOAYI/boHxf1pxgAYnDVYNsaA42GINh67EhHRtqWqAuCsWPTef6XM9uxU4Q7tl2OxUqmcNag+0QCt9HFAxqATGFSr1TQw6DvOQFXHus7BwIBogzEwIKJta2Fh4Z4AbgMAInL5+Pj4DSU3aSfrZGNEhOMM6IBUKpU49WjMHAyclUhE0ulKMTMz0zdjICJ5167uMhEdJAYGRLRttdvttL/7l/pWpI2Qjt/4ndJaQUMt6UoUL/b7zUoUZzFa1pWoUqkM6kqUpQVmDIg2XimDj2n4zM7OHjE6Onqk936v937SGDMpIpOqGh8PATApIrWuuaYn4nR13vsREclEZDY+qaqzWEo9N0VkTlXnVHVGRGZUdcYYM+O93xceZ0Tk11mW3SAixdb9D1AZuu5cc3zBJsrz/BfW2qsAHA/gjs6547Is+1nZ7aKhs6wrUQwAeuh0JfLeL8QVjQcFBiKSqyoAzAKY4BgDoo3HwIAwNzd3m0qlcjyA4wAcA+A2qnpbEdkL4PYA9gIY9d4DAIxZTDSFL+jOY/d29774xd+rTr+6IgJVXfYIAM45WGtvAXA9gOsAXCciNwC4FsDPRORno6OjP2fwMLxUteqce2AoLtRqtQtKbdDu8CUsBgbA4qBvBga0XumsRACw6qxE3vtmnBVr0FoGqhq7Dt2KxYHI7EpEtMEYGOwSqjphrb2nMeaeWAwAjvfeHycixwEY664fL8C3ucPCvxOBlQGKc85ba6/F4sXNz1T1KlX9iTHm8jzPry6lxbRm1tpTknnLLxaR+VIbtAuo6ldF5EVh+0wA7y25STRkwgJnwFKAsFrGYKHdbi8kgcGgrkR5qHOrqt4RnK6UaMMxMNhhVFWKojjWe38vETlJRE5S1Xs55+4kIpJePK/j4v9WLN2ZvxHAjIjsV9X9sbuPiEyr6oz3fg4AKpXKDMIdI1VtFUUxm56wWq2OVqvVTkBijJloNptVEakaY/Z47/eIyFToqrQn6bJ0mIgcCeB2CJmMAe02AO4Q/v22iHR+ZmvtDIDvq+r3jTGXe+8vz/P8Cl58bisPTLa/XlordhHv/deS9QweXGJTaHgtm66038rHWAoMWpOTkwvOOYT6q44xUNVbQpmBAdEGY2Aw5FQ1c87dD8ADAJxRFMXpqnrIWrrtBBbAVSLyMyzeWb9aRK5tt9s3DkNf/jj2QVWPUtWjAByNxa4Qx4XHo/ocOgngASLygNhFyTnXds79t/f+GwAuAvCNer3+qy35QaiXM+OGiLAb0RYYHx+/3lp7JYATANzeWntMnufXlNwsGiI91jHoeZ0hIiPhu7cJoDP4uNls9u1KJCJZ+Jt2a9i1amBgrX0NgHMAfCfP80eu/hMQ7W4MDIaMqo4XRXEWgDNV9XTn3H2R3DUfEAj8BsD3ReRyVf0JFrvW/Kxer/9681u9eSYmJm4CcBOAK3o9r6r1hYWF49rt9vEAjheRe4Qsyt2xMttQUdWTReRkAC8FAGvtrwBcICIXich/1Wq1H2/eT0ORqhrn3Bmh2K7VaheV2qDd5etYDAyAxazNNeU1hYZQu+ux36xEI2Hs2LLpSgdlDLrGGABrG2NwFoAjADxQVUfD9KhE1AcDgyFQFMWJ7Xb7ESLyCOfcAzC4+0xLRK5Q1ctV9fsicnmz2bx8cnLy5q1q73YiIg0A/x3+dajqyMLCwt2897G71ckA7oPFMQupowE8Q1Wfoaqw1l4N4PMicn6tVvsyux5tjmazeU8Ah4bi90Rkpsz27DIXAHhe2D4TwAdLbAsNmTjGQFVboftmz8HHSGYlEhG11jYBjPQLDFRVnHO1UFxPV6KTwuPYwsLCcQB+tIZjiHatUgKDer3+PSwtbd4aVHc3UtVaURS/o6qPAvBw7/0dB4wHmAPwLQAXisg3arXat0RkbssaO6RC+vr74R+AxT88CwsLd/Xen4HFrlmnY+nOaXQsgBer6oudc4W19uuqer6I/Bu7XCyX5/m3cYCfc+99Z3yBiHxtI9tFq0rHczywb60esiy7CEvveXNQXdoZsiy7AMCRzrkbgKWuRMaYVshg9+1KFDbj78kCFgODfl2JMgASXuPW8DdxYMZAVSedc/EGA7z3R4CBAdFApQQGIuKxeEFLQeg68WAAzyiK4glhXYBe5lX1y8aYL7Xb7W/U6/XviwiDqw0gIorFPxo/AvB+AJibmzvSGHMGgAeJyMMB3Dk5pAbgbBE5G8DfWWu/KSIfaTabH5+YmLhxq9u/3Rzk55zjC0qS5/nV1tpfYnHQ/p3n5+dvOzY29pu1HMvv9t1HRLyqulgeGRmJGYM2MHAdg5gZWAj1m+Fiv2fGYGZmJh8dHY2vuT/srqpqRUTavY5Jg4LwGoev4Uci2tXYlahkRVGcqKrPcs49C8BtgZ7jBH4O4Esi8h+1Wu0/RcR1V6DNMT4+fgOAfw3/YK09NgQCD1XVh2Pp7qgAOF1VT69Wq28NQcIHa7XaR9kNZv1U9QFxsyiKb5TamN3pAgDPBAAReQCAj5fbHNrOpqenK7Va7OXTGVswcOVjhK5EsetR7Pvvve8ZGIyMjHSyA2FGvFjMAPTs0um9PySuuwMAxhgGBkSrMKtXoY2mqhPOuRdba6/w3l+hqn+CEBQEHsDXROQlAI7N8/y4PM9fmGXZvzMoKFee51dnWfbeLMuekmXZ3pBFeBeAdAyHAXCGqr7HOXe9tfb9c3NzJ5fS4CHknLsLlmaTumJycvKWQfVp43Vlac7sW5EIQKVS6Vz8J+sYrLbA2QgAGGPiYOCBsxgZYzqBQZIxwPT0dN8xd8YYZgyI1okZgy1UFMVd2+32Hzrnno3F6TKXEZHvee//BcBH6/X6tVvfQlqPEKR9AcAXVPXlRVGcrapPB/A4AHFhrhzA8yqVyvOstReq6tvzPP9kGONAvT0obqgq1y8ogYh8PclcrmucAe0+XQOM17uOQTrGoO+sRCJSi7+T3vv9MRMwMjKS9WuXqk6l4/NEhIEB0SpKCQxUNVtYWLgTALTb7dmdPFd8WHDskar6cu/9Q2XlKOJrAHzQGPMRToU5vMKF/ucAfE5V69baxwB4hog8Ekt3zM4QkTOcc9dZa9/TbDbfsZNni1LVfGFh4VgAaLfbM2sNdsOKuxHHF5RgdHT0x865GwHsFZF7zszMHLaWzE3Xez497NMh0+pUte6cu1tSXtM6BqoauxLFwGBg16NWq1WLwUClUpmJQYIxpm9gICLLxuqpavesc0TUpZSuRI1G40Tv/Q+89z8A8PYy2rAVnHMPdc59W1X/A8DZCDMqBBeq6lOyLDshz/PXMijYOUSkUa/XP1qv1x/jvb8DgL/E4loL0VEAXjcyMnJNo9F469zc3JHltHRzWWtPjp9zY8xb1nFo5w61qjIwKIGIaNKdSEZHR88YeEDgnDslvuci8neb2ETaJsICm1+O5bV2JcLSIOMYEMTjes5KJCKdQQzNZnMm2T9o3YNlgYCIjPerS0SLOMZgEzjnHmqt/baq/ieA+yVPzQJ4rzHmHnmeP6Ber3+CMwrtbGNjY7/J8/x1WZYdrapPweKKyp2nReRllUrlqp0cIKyHtfYYLM6GAwBXrnU2HNoUHGdA6+a9bwPLMgH9ph+NmYFmqB+7EvULDGIA0BofH2/E/atkDPakZVWd6FeXiBYxMNhAISC4OAQEpyZPXScir8iy7LZ5nr+wVqv9oKw2UjlEpKjX65/I8/wMVT0di92OovEQIFxprX39Lv/jlfZn5/oFJWq3253xHarKcQa0Jt77FgAYY9rAivEHqVFgKSAwxqw1Y1CEfwAWuxj1a4uqdmcImDEgWgUDgw3QaDRub639YI8MwU2q+qdZlh2fZdlbufAYAUC9Xv9mnuePqlQqJ4vIJwDEUZ4TAP7cOfcT59wLVHU3fj45vmCbqNfrl6vq/lC8zy4PWGmNYsYAq4wxwFIAENc9iNOV9sswjAKL05qGMV0x8BiUMYgrI6ffsatyznUvbEm0a+zGC48NEwaZvk5ErgTwrOSpGBAcU6/X3yQijX7noN1rdHT08izLntJut+/TFSAcFaY6/Vaj0TitzDaWIL0zzRmJSiQiXkQuDMVqURT3L7VBNBT27NkTuxK1gP4LnMXZiowxy8YY9OtK5L2vhfPGbEER6vfNGAAYC4/7wuNag9s3qWq/TAfRjsbA4AA1Go0nO+d+DOC1WFxgBQAsgL9kQEDrMT4+flmWZU/B4urKlyVP3U9ELrTWvmdmZmbHz6YxNzd3GyytLP2rPM+vKbE5hOXTxbI7Ea3RmgYfx1mJsDT4eOCsRF1diToZBiyuQN9T0pXohvC4alciVa2p6iOdcw9erS7RTsTAYJ1mZ2ePcM6dJyIfB3B03B9mHjoxz/PXMSCgA5Hn+QW1Wu0UAM8GcGPYbQC8YGRk5IeNRuPx5bVu8xljTk+KzBZsD+n7cHrfWkRLWunjausYJLMYrTZYOQ4+joOUHbCUSehFRGLG4PrwuGpg4Jw7FYvBxlNXq0u0EzEwWAfn3KOq1erlqvrEZPdPROQR9Xr9d/M8v7q0xtGOICI+z/MPFkVxV1V9G5b+yO4VkX91zn18enr60EHnGFZpYCAiFw2qS1sjz/PvAoirrZ+qqlwUk1bjw+NqYwwqab3VZjFKZiWKmQIX9vcdY4ClrkQxY1BT1b7Tm4Z2xO+hJyRZDaJdg4HBGuzfv/8Qa+0/h6zAUWF3oap/kmXZPbIs+3yZ7aOdZ2pqal+9Xn95pVK5b9q9SFWfPDo6erlz7mElNm9TqGpnPEW73f5WmW2hRSKyAOB7oTg2Pz9/jzLbQ0PhgAYfx7EG/cYYIHQZEpFlYwwwoCsRQmAgIjFjgJmZmdWyBnHNjsOccxxXQ7sOA4NVWGsfUKvV/hvAM+M+EbnUGHNKvV5/M9choM00Ojp6ea1W+y0Af4WlP7S3V9XzG43G21a7+zUswp25+4Rio16vf7/M9tASEflm3K5Wq7ttMDytj4qIAsu6CPULDKqhXjN97Dcrkfd+2fSmqrrmjIGqxowBRkdHBw5AFpGTk+IdB9Ul2okYGAzgnHsBFld0vF3Y1RKRN9VqtdO5FgFtFRFZyPP8NcaYU1T18qXdco5z7sKwKNhQazQaJwPIQ/ESBtzbh/e+k71JszpEPWiyHTMHAwMDLN3wWK0r0bLBxzFzMGiMAZa6EsUxWzDG1PvUhaoKgL1J+bYDzk20I5XSX9R7/8tKpfJSABCRn5fRhkFUdbIoin9S1Scku3/svX/G2NjY9/oeSLSJarXa91X1tKIo/lpVXwFAANwXwMXOuWdkWfalclu4nPf+59Vq9aWheNWgupVK5f6qi9cUIsJuRNtL+n4M7Frhvf9Z/G4HcOXmNYm2C+/9lZVK5R9V9RwsXdzDGNNS1UELnI3EeqG8WleimB1d93SlInJz/H5ptVp5v8rT09NTtVotPR8DA9p1SgkMJiYmbgLwjjJeezVFUZzknDsPQLrAycezLHu+iMyW1S4ioDMTxx/Nz89/3Rjz/wDsAXCEqn7eWvu6LMveEFP5ZRsfH78Ba/ycp3ei064rVL56vf4ra+21AG4P4ISZmZnDJycnb+5Vd2xs7Dps0+922hxjY2O/cc59FsA5WBoYDISMQb91DLB0/dE9XenAMQbJNKXxsWd3SlU1zrk8tOGW8L0ogzIGWZbtjQEEABhjGBjQrsOuRAnn3Nne+wuwFBS0wkJlT2NQQNvJ2NjYp0TkvknXogqA1xdF8QlV7XtHbBvrBAbNZvPbZTaEVkqyODI6OnpqqY2h7agKAGkXwGSMQb+MwbKuRGuYlagGAMaY7sHH/cZZ5VjMqkJVG1iaxajv96Oq7u0qH9WvLtFOxcAgsNb+vqp+FsBk2PVrVX1QWKhsW9yBJUplWXZVnuenAXhf3KeqT3TOfXl2dvaIEpu2LrOzs3sBHBuKV4+Pj18/qD5tPY4zoFVUAEBV28m+NrCmdQyWrXyMPoFB9+DjMGNW365Es7OznUHJ1WrVAmiE+oPGGBzZtWtHTg1NNMiuDwxUVay1fwPgA1j6Qvp2u90+pV6vcy512tZExOV5/gIReRmWBvudVq1WL3DO3anMtq1V18Jm7Ea0Pa15nAHtPkl3oXTSgHVlDFabrrTfyscxYFhx8mq1kxnw3jsAFgDa7XbfjIExpnuF+al+dYl2ql0dGKjqqHPuwwD+NO4TkU9lWXZW6B9NNBSyLPtHY8wTEe6KAbiLqn6z0Whs+24flUolHV/AgcfbUJ7nl2KpT/dvqWq/iz3aneLvw4qMAXqMZVRVg3D9EbscrTZdKbpWPsbS4OOegYExppMxaLfbFiEwGNSVyHs/2bVrql/dlHPuuLXUIxoGpQQGRVHc1Vp7ibX2Eufc35bRBlUdtdZ+AsDTk31vq9VqTxKRxoBDibalWq32ae/9g7G0yudeEfmytfasMtpTFMU94+c8ZOV64sJm218Y9H5ZKE4sLCzcrVe9hYWFeyXv+Ru2rIFUmvn5+XsbY94IACKSLh42KGPQufhPxiKsaVYiVY1jC2KA0LMrURoAjI2NuTDOYGBXIhGJgcFMeMzXslaMqr6+KIqenwmiYVNKYNBqtcawOM3ifb33J6xWf6Opat059+8i8piwy4vIK+r1+stFpD3wYKJtbGxs7BIsdvX4cdwF4DNlBAftdnsc4XMuIsf3qqOqVQCnhKKr1+uX96pH5VPVTtDmve/Znajdbk8gvOcAeBd1FzDGTKpq7LbYCQKS8Qa9AoNOFmFkZKR7jEHPMQkxM2CM6Z6VqGcg0Wq10oXPrIhYAPDeDxp8HAODX8R9c3NzU/3qh2OqqvoI7/05g+oRDYtd15UoBAWfAfA7YVcbwO9nWfbWEptFtGHyPL+6KIrTAXwn7BoD8Dnn3O+W2KyeGo3GSQDiXcbvxAGFtC2l2RwOQKYVNJ3rc0BXov3796cZgzVNV9pj5eOB05UaY2IA4MP3yqpdiZKMwS/jvpGRkal+9QHAOXemiEwBePbMzEz3GAWiobOrAgNVnXDOfQHAQ8Kulqr+Xp7n55bZLqKNNjU1ta/ZbD5cRC4Nu2qq+gnn3KNLbVgXji8YHl3rSzAwoF46gUGycNmKjEGlUukEC8kYg4EZgKSLURysHAch9+tKFDMGLjzGLsJr6Ur0q7jPez/Vr3445pHxvCMjI48dVJdoGOyawCAMNP4kgAeEXU1VfWq9Xv9ome0i2iyTk5O31Gq1h2DpTm9NVT/pnHtYme1KpeML0ikxafvJ8/waAL8Jxbvt37//kBKbQ9uQiPik2Ar7VmQMjDHV7npxViL0X+Bs2eDjZNrSfisfx8AgZgpsOG4tXYluRAgojDF7+tUPx9wvKW5512iijbYrAgNVrVprPw7g7LCrEJEn1uv1fy2zXUSbTUSmsyx7GIALw65RVT1vG81WlPZVZ2CwzYnIxXEzy7Lt8jtE28eKrkRhHNEyaWCQDD6OmYN+gUFc92BNKx8nYwlcOG7VwccI6xip6gzCAGTv/fiA+gBwdLJ9bN9aRENixwcGqirOufeISEzxtVX1WVmW/XupDSPaIiIyk2XZI7E05mAcwBcWFhbuVWKzEPrjxgGqv6zX678usz20Ou99pzsRFzqjbukYg0ErH6czD3nvl01XOmBBtBgANEP9gSsfJ12JbNdj1qN6NBmOnQYwF7b7BgaqKgBul+wairVjiAbZ8YFBURT/B8Dvh6ICeFG9Xv9EiU0i2nIiMtNsNh8B4EehPNVutz9nrS3tDtfo6OhpACS059tltYPWrmscCAMD6tZr8HGvwKBz8e+9X9PKxzGY6DEr0WqBQcwYuK79vYyHurOqumpgMD8/vxfLp0s9ZsC5iYbCjg4MGo3G/6eqr4xlVf2TPM/fX2abiMoyOTl5s6o+DEszbtwWwOenp6cPLaM9ad9c7z0DgyGQZdl3sHQBd79wx5QoSgODvtOPphmDdru9bB2DXvWBxXGCab0YIPQbfNzdlcgY48L+voGBiIyFuvMiMhfqTwyof3TXrsO4+B8Nux0bGDjnHiUib0x2vbFer/+f0hpEtA3U6/VfichDAFwfdt15dHT039ayiM8m6PRRT/qu0zYmIg1V/WEoHlYUBftUUyrtSrSmdQz27NkTuwYNzBiga4xBMvh4TV2J1pIxUNV6eJzHGroSVSqVw7t2mbm5uVVvtMzMzBzOoJq2qx0ZGMzPz99bVT+GpS+kj2VZ9qoy20S0XWRZdpX3/tFYmr7vgc65d5fQlLiwWTvLsu+V8Pp0AETkkrjtvb/foLq066wpY9BqtVbMSoRVMgbompUIwMAxBlgaS+C6HnsGBmGQdC1sdzIGq4wxiEFAPDdGR0e7g4UVRkZGXmGtfdpq9YjKUEpgUK1W5wFcAuASY8xPN/Lc8/PzRxljPo3FRZ0A4MIsy54tIjroOKLdZGxs7FJVfTaAOL3gcxuNxv/ayNeoVCqzCJ9zVb0yfS6Mbdgbij+Kf4Rp+0sDA2PMssAgfc8BXLXFTaMSeO9nAMSJA+aTpwaNMUizAk0AqFQqA1c+RggA4iDlZGG0gRmDmClIBiv3m940XjOgUqnMxzEGqjrWpz5EJAYBVyF8l3rvj+hXP5zPAPgfIvLmmKEg2k76fQA3Va1W+zGSbgQbRVXzoij+Q1Vjv7+ft1qtx4tIMfBAol2oXq+fZ619LYDXA4CIvNE595Msyz6zEeev1WpXoM/nXFXvJ9LJpF/Sqw5tT+12+xJjFu8pdc3hjtHR0cuxCd/ttH2NjY19z1r7fgCvBdC5AWCMaYVJinpdZ8R9KiJxWtN4wT+wK1GlUunOGPS70K+FdhTh0Yb29FzHoNFojMXfa+99Q1Vnw3dU3zEGqjoVNm8GsB/Aoao6MGNQFMXvIExxWhTF7wN4+6D6RFttR3Ulcs69Q1XvE4ozxpjHTkxM3FRqo4i2sSzL3gDgQ6FoVPWfi6K422a/bnqnOb0DTdtfvV7/byxN/XgfDrYkLGUF2sm+NgCIyKBZiVrJ7rgg2sAFzmIAkQw+7pkxiIOSY6YgZg7QpytRpVLpZAZarda8MSZmDAZ1JTostPkWLAYHMMYMDAy8949Kjr//oLpEZdgxgYFz7kUAnhuKbRF5crhjSUR9iIhmWfYHWFpcbMJ7/7FB6fONkN5pbrfbDAyGiIg0AXw/FMcXFhbuWmZ7aFuoAMsGHANLC5YNmpWoExgkXYP69WSIxyzLGPQbfIylTMKaAoNWq9X5zhsfH+90JRo0xkBEDg3nvhXALQDgvR+4IriInJBs331QXaIy7IjAoNFonKqq/5Ds+ossy75YVnuIhomIFN77JwC4Luy6p3PufZv1eqGP7b1DsQh3oGmIqCoHIFNHkhVYkTHA4FmJmsm+1WYlWtcYAywFAGsKDIwxsb+/ArBxpWSEtQ36OAQARORWLHYlgojsGVAfADqBgarelRk32m6GPjCYnp4+VEQ+hqXZBP4jy7I3rnIYESXGxsauA/BMLP1xfrpz7pzNeK1wh3kyFC/nGKDhk3b/6h5nQLtSBQCMMZ3AIJl+dNAYg14ZAxNuHiB5zsRjuqcrBVDtrg8AIlILjw4AKpXKaoFBzBhYEfHGmDiQetAA4T2hLfsBTIftvoFB6PZ0x2RXZq09asD5ibZcKYOPZ2dnjxgZGXlKKP48y7LzD+Q8qirW2nOxtNrgVUVRPKter/sBhxFRD3mef6XRaPx5XP9DVf+u0WhcXK/XD2jxsbm5uSOr1eqTQvGqLMu+ACy/w5zeeabhYYy5xPvFr1kR6byf8/PzR1UqlSeE4pXM3O588/PztwVwHwDoGnjbN2PQbrerYaBvGhi0kgkJRrA0uDiWY71OYJDUH0UyZSgAeO9rIrLmMQbtdjsOPp4P9WPGoG9goKqTIgIRmQ3BAURkql9959xRWPn/cTsA1/Y7Bli8AVqpVEbGx8dvGFSPaCOUkjEwxtxBVd+uqm/33r/oQM9jrf2fIvLoWKxUKk865JBD9m9MK4l2nzzP3ywinwrFERH50KDBd4MYY+4UP+cAnh/3p3eYOfB4OI2Ojv4EwEwo3isO9DTGHBffc1X9/fJaSFvFGHOCqv52KHbuhicZg0r3Yl69Bh+PjIyk3Yq6b1p2ugtVq9VmeFxInl/R/ShmDLCyK9FIr+47ScagOzAYNF3pZKg7IyLTYfegjMFtkmJsz2371Y9GR0dfUalUzuWiaLQVhrYr0fz8/H1E5G9iWUT+Z5gqj4gOkIiotfa5AK4Ju05wzv3jBr9GJzAwxjAwGEIi4gF8NxRHG43GPctsD20b6XpB6XiD7gvxXl2JWj2eBwDMzMykF/4xIOgEErOzs73GGSwLDLz3aUahV9YgZgYsABhjVs0YIHSJFJFpVY2BwVS/ysaYI2OTEb5jjTEDAwNV3aOq5wB4WFEULx1Ul2gjDGVgoKpjxph/QbiLICL/mmXZe0puFtGOELJuz8LSH/bnWGufuRHnDnOUnxSKc+HOMw2hNNtTqVQ4zoCAJDAYdKHfK2OQjDHA7OzssgxApVIZ7a6XjDFY9nwijjt0ADAyMtIJDGZnZ1cEBiKSh0cbjuuMMRhwpz5mDGZFZH/YN2jwccwY3ICwKJyq3m5AfTjnnhS7J6nqEwfVJdoIQxkYOOfeDuAuofhL59zzB9UnovXJ8/wbCAufBe8MqxUflEajcRKW7tZdGhc3ouHjvecAZFomZJIAACMjI+vKGKTbxphlgYQxphMoeO8X0sfw/IrAQFWz8FzMGHTGLFSr1RWBgfc+D8fZUI4ZA4Mei6iFVYuroe7MWjIGqhozBtdjaRa4vX2qx2NOToqncBYj2mxDFxgURfFYAM8JRQ/gOVNTU/vKaxHRzpRl2esBfCUUJwF86GD/KKV3ljm+YLh1vX8MDAhYY8YAPQID730nY5AGAqE82l2v3W53AoNei6LFMQZx8HHalcgYs+JCP2YMELoSVavVGBhgdnZ2xTiD+fn5OLMaKpVKOsZg0ErJe8Nr3SAiN4d9AxdEE5F7JcXxhYWFTV+Akna3oQoMZmdnj/DevzfZ9YY8z7/S9wAiOmDh7t9zEKbhA3CGtfYVB3PO9M5yeseZhk+e59cAuDEU766qfS+IaNfoOcZgZmZm2Q2FpCtRr+lN0+djOc0YNAGg1WqlgcGqYwxarVaR1O/Vlage2tEIrxO7EqFara4YZ1CtVjuBQavVmlHVOBh/ol/Xo2RBtJtV9eawb2BgAODEtNBut0/oVzFSVWOtPWa1ekS9DFVgMDIy8k6EtJuIfC/c0SSiTZLn+S8BvDyWReQNRVEczEBTZgx2lkvDo3HO3XtgTdrxVLVnxqC7axBCxkBEOnXa7XYnY9AjA9ApT0xMLADA1NRUp36z2Rw0xqAAgMnJyU7GoNVqrcgYxK5EcYxBq9XqZAySxc862u12JzAYGxubabfbs7E6+s9kdGh4vDVmDAAc1qcuQrAdj4lB1KpdOouieCGAr1tr77BaXaJuQxMYWGt/T1XjnOiFiDxbRJoDDyKig5bn+bki8slQrKnquWEQ8bqEPrkxDX5rlmXXbFQbqTQcZ0CpNDDoZAOMMcsyBt77uFhZGhj0zRi0Wq30+yb+3R84xgBhLJMxJgYEq2UMlnUlGh8f7wQG7XZ7xYV+Mr1pS0RstVqNgQEajUa/7FnMGNzqvY+BQd+MwcLCQufCXkTirIt36lcfAGZmZg5X1b8BcDSATVvBnnauoQgM5ufnjxKRt8Wyqr62Vqv9d5ltItpNWq3WSwDcBACqem/n3J+u9xzOuftgqW/xJSKig+rT9pdmfUTkvmW2hcqXfqYHdQ1CjzEGaQagO2PQVW6GfU0sBSK9blQs60oUJjqIx67IGABYFhiERwV6ZwxEZBwAVHUOAJrNZicwqFQqk931g0PDsWnG4NB+Y7dUNa4LseC9vygce0yfcwMAqtXqWViaGelB4YYM0ZoNRWBQqVTepqqHhOI38zz/21IbRLTLjI+P32CM+YNk158XRXH39ZxDRE5JiuxGtAO0Wq3vJEUGBtQzYyAi/cYY9JyVqNls9htj4LtmMlsIz6/alSiIQUKvwcf18NgIjx5hEbIkm9DhvR8Pz80BwPj4eBxjgHa73S9jcEg49lZjTAwMzMzMTL8pTo8Oj9cCuDq83sDpTQGcmWzXiqI4s29Noh62fWDgnHtk2oXIGPN8TnFItPVqtdqnReQToTjqvf+Aqq75O0RVOxeOxpjvDKpLw2F8fPx6hPnYARznve93p5R2h05g4L1fV8YAyYJlPcYYjHbXScvdgUG4Ax+7K6WBQbzQX9GVSFWXTVcaxCChb8YAwFwou9geY8yKwCB8V06F529ttVqd2RRrtdoh3fXDMXHdg2sRpjdV1aN61U3ataxLn/f+Xv3qJq+TOefOWc/3Oe1cpfwS1Ov1HxhjTjTGnAig70p+qjqhqu9Odr2hVqv9cPNbSES9tFqtl4lI/IN2WlEUf9Cvbp7nl8XPuff+fyK5o9xuty/tdxwNF1WNQZ5UKhXE91xVX1lqw2hLZFl2CYCvAYD3/r/i/jRjgDWsYxDu0Puw3S9j0B0YxHEG3RmDdEG0hWR/AQDtdnstXYk62+12e0XGIAYGIjKf7J4N+1YEBtPT03sQ/h+MMfvq9XonMGi32/0Cg3R609+E3UeoaneglR4T13iK2ZR79KsbOefeqapvc869f7W6tPOVEhiIiKvVaj+s1Wo/rNfrv+pXz1obB9AAwI+zLHvz1rSQiHoZHx+/XlX/OJZV9c2NRqNnaltEbPyc53m+H8Cdw1M31Ov1a7egubQFRKQT5Hnv75F8t/960HG0M4SuNzZs74/7JyYm1psxAJYu/PvNSrTQtT8udras/v79+zsX/tVqdU0ZA4TAIE5XGvTNGKjqWHicS3bPhH0rMme1Wq3TXch7Px1mP7IAUKlUegYGxpgjQ/0bRCQuiGastbfpVX92dnZvXCUZYQ0aERnY5bPRaNweS2tDPWtubu7IAdVpF9i2aaNGo3GqiLwoFL2qPl9EioEHEdGmy7LsAwDincFJY8xbVzsmDDyO3zfMFuwgaWDQNY6Edo9eF/p9MwYxUOhaBK1zfHcgkVz4ryljUK1WB2YM0GMlYwB1ADDGdDIGsVtRrzEGAJZ1JQr6ZgwqlUonMMiybH/Y3AcA3vtDu+uH1z8ynO+GWq32m2R/z8CgWq3Gmy8QkU+Fusf3qhsZYx4DIK67UK1Wq08YVJ92vm0ZGKhq1RjzHoQvE1V9Z71ev7DkZhERFmceEZEXI9x9U9Unzs/PP27QMaqaXjAyMNhBms1mOpCcA5B3p3jhnwYDnYv+VqvVM2NgjFkWGKjqahmDNY0xqFQqnQv/dIxBGAfQL2OQhfo2qR8XO+uVMegVGMyF+uPd9dvt9lTY9AgBBEJgYIzpmTEAEAODG8Mg5/nw8+3tUz/OYnSj9/67YXvP9PR0z8Aj/BwP6iqf0a9u5Jy7S6PRuLzRaDxttbo0fLZlYOCc+1NVPTkUf5Pn+Z+X2R4iWi7LsisB/FUsG2PeuW/fvql+9dM7ycYYBgY7yMTExE0AYpfQY2dmZvou2EQ7Vq/AoG/GIK5jgK6uRMmCZ91jDPoNPu6ZMUjXNWi3253AIAYJ3vu+YwxU1SX7+mYMugcfp9vJc2n9PeH8M8m0rvvCvn6BQWeMQSjfGOr3Cwziuge/bLfbV8ed1Wp10KJocWHCGKz81oC6CK9/noicJCLn9utKSsOrlMBAVY2qjod/yyJ359ydAbw6lo0xfygi01veSCIaKMuyN4nIZaF4VK1We0P6fPo5x/KBx5yRaIdJuhPJ6Ojo6eF973XxRTtMmMkmXohL8tSgMQaxN0DPrkTdYwbQJ2MQuwn1WPeg87vXarVWdCUalDHo6kq06qxEqjqf7JtNn+tq657w3HSyb1/YtyIwCItITgJAu92OU5uuFhgcHc73qxCwx0DlmF6Vw1iI40Lxn8PjcWHF5Z4ajcZpAOKA5lEReWG/ul2v1W9KVtpmSgkMGo3GvZ1zs865WWvtx+J+VRVVfRfCB1REPlGr1T5dRhuJaDARaXnvX4hwZ1BEXmStfUB83lr7W/FzLiLpwGMOSt1hVPXSZPsz4X0/t8w20dZwzp0J4DQAEJE4tXhc7CzOMrSWdQyApQv/Nc1KJCI9MwatVqsTGOzZs2fF4GP0HmOQd9XpBAne+xUZgzj4OK5jEPbNhcdBGYPpZN902DfVXX9+fv4whECrUqncEurdFNrVMzBQ1aNDe38Zdl0T6t++V/1Go3EXhOtA7/0HYrOstX0HLIvIY7rKj+5XN7LWPsRae4219kuq2mvNCdpGtlVXIufccwCcFYrTYYpDItqm6vX6xSGYBxa/T94d7nQto6rxTiKzBTuQiPB9JSBZxyDo2TUI/WclioOPe2YMwmrHqSYAeO+XXWx2jTnoZAxiV6LujEH4fooLorlkf9+MAYCx8Dif1I9BwqCMwf5k3/6wb6q7/sjISKdLXlEUN4d6AzMGxpjYrefarseegYExJk5tOj02NnZpnIp6lZmMHhwe4wxMJ8/Ozh7Rr3JYU+Jd4Wd8iLX2JQPOTdvAtgkM9u/ffwiAN8ayiPwx7ywSbX95nr8KS3+ATrTWnjOgOscX7EALCwsMDAiq2h0YxHEGq65jEPTMGMSuRV0zDAF9Vj5OuhK1k3ELncHHPcYYjGLpzvmKdQywlE1IxZWSO4GBMSZmDFZ0xemTMdgf9q3oZuO9PzxuTkxMxAv2m8K+nhficfGzZM2Da8P+noEBgBPC409DvbhO1J17VQ53++N4sb/DYiAoIyMjfcclFEXx8OR1ICIvSW4U9eScu5O19tvW2s+yC9LW2zaBQeifHKPg79RqNS60QTQERGRWVf8oKb9uwNoGDAx2oMnJyZsB/HLVirTT9cwYrGMdg1he66xEC13PI7zeaNfzUcwYLAsM9u/f37nwr1arnYxBnJVoUMYgHWMQMwZ9xhhMhudmkn0xSJjqUT8GBvtEpA0A3vsYGBzeo34VSwHDssAAfTIGWBpf8NPw+JPQxhN6VW40GicidNsyxnwcwM/Ca5/a5/xQ1cd3v6a19n49Ky8d808ATgXwSOfc2wbVjebm5m5jrR00yJrWaFsEBvPz86cAeEEoelV9SVgFkYiGQL1e/wSA80NxQkT+rlc97z0Dgx2K3YkIXYFBsvpxz4xB9+DjOF3pOlY+juXufuuxW1D32kfxon9ZV6KRkZFOOZ2uNG6raq+MQexK1Ejqx5l9VgQGCAOJEbrgAMu6FU11VzbGxK5ENyf14/aKjEGj0diL8P+cLIYWA4N+MwcdEx5/Fh6vAgDvfc+1D4wxJ4fNudHR0Z8CuBgAVPXeveoHvxMeP4alWZ4e1q+ytfZMAA9Mdj3TOXenAeeHtfasSqVyNYCfOedeMKgurW47BAZijHk7lmYpeFe9Xr+45DYR0TqJyMuw9If3qQgDEhMceLyDpQOQaddaFhgMmH605zoGSbk7Y9BzutI4+Li7K1Ecc9C9KGoSKCzLGFQqlc6Fv/c+na40XvT3zRgYY9JZifqOMUAIDPpkDHp1JYozFd2SnD8GBisyBgCOihu1Wu26UD9mDo7q030n3mG/JrTnqnDccb3qi8iJYfP7IuLjZ15E7tXj3Ah38ONMSR9U1S+Gpx7Sq37w3K5yRVWf3a9yCNo+hMVgT1T1rasFEjMzM4dbaz9rrf15o9F46qC6u9F2CAzuiKULiBsXFhZeU2ZjiOjAZFl2FRb7nQIARGTZ5AGqesmKg2jHYDcx6jfGwHu/pjEGMWOQrHMQDZyutMf0pvHCf1lXImNMz65ExphOxqDZbK7IGKD3GIO+XYkwIDBQ1U5gYIzZHzanuisbYw4N9W+N+7z3MTAY685iVCqVuBryfAw+jDExMMhmZ2eXLXIWphOOx1wT9sXMwfj8/PyKAc6qereweUU4//dD+egwTrTb6eGxVavVLjDG/Gcon9prdqIwUPl3Q/FjAL4ZtvtevBdF8TwAt012Zar6Z/3qq6oZGRk5D8AjARwrIv9srX1gv/rA4hhYa+1HrbW/sNa+PkzRO5Cqyvz8/FHhZxoqpQcGInLXpPjHU1NT+0prDBEdlCzL3gAgLqxzdPocLxx3Ng5Aph5dgHt2JYoX/v3WMUCf6UqTlZGj1QYfd2cM4mrty7oStVqtzkX25ORkOl1pz4xBuDDMws+ypsHHvcYYtNvt/WGz1n2hr6qHhvqdayJjTBxjAGvt4V31jwybcTE0FEURuxShVqulF88oiuJoLF0DXgMAzWazsyiaMaZXf/27hdf6UagfAwPUarUTuyur6v3Cz3C5iMyKyAXhqdxau6L7URh7cHg45v8CeF946i5hjasVVDV2HZoB8Iuw/Xt9AhU4554GIF3tuQrgnf0u4FW1GqbNfyoWF5D786Io/rZX3agoihOdc1cYY37jnPtR6C4/UKPRONVa++pGo/H0MF5kVZsVdJQeGGApRXhhlmUfLLUlRHRQRMSKyCv6PMfAYAebnJy8RVWvj+Xuu7K0K6xputLV1jHoN12pMaZ7MHG/MQYDBx+jqytR1/Sl6XSl/TIGdSytMbCmwccisiJjUKlUOjMUzc/P7+mqvyJjUKvVOuMNksHJ0YrAYHx8/CaE/4M4Y1FyfFwl2WdZdi0ATE1N7UsCkWO76udY7OEBY8yPAGBiYuJGADeF9q4IDETkfuHYiwFgdHT0ijhFqzFmxUxGInJ22Jyp1Wpfcc79G5Z+Rx7RXb8oinsCuGcovllEXhq2syzLnt5dP3SPSrMJMcNzYggYVrDWvhzAmV3neYW1tmd3qEajcbT3/ksA4pSvJxhj/qsoinv0qh9e4w0i8i0AfyUi/+Kcu2h+fv6ofvWLori7tfZzzrnCOXerc+7NcV2NXlS10mg0nmatfb+19r2NRuOJg7Ie2yEwAIB2u91+qSwtE05EQyrLss+o6n9072+32wwMdjgR+UlS7HnHjna0nl2Juhc4w+qzEvWbrnStsxLVwuOaxhgkgUGRZj3iOgboyhjMz8+PJXU6gYH3Pl5oVruzEghdibz3ncCg2Wzuj9sjIyNTXW2NGYNOYBCyDQsAUKlUemYMVPXGpL4CuD7sX5YxwFJG9/qu9SGuDscuCwwWFhaOR7hmVNX0c/6D8HMtW/sgXHjGsQeXhHN6EYkDlntNcfrg8NxXRKR5yCGH7AfwjbBvxYBl731cUK/ZarXeV6vVPoulqVef2V3fOffbWFq1+X3GmFMRgkVV/ePucRUzMzOHiUine3sMarAYFL6juzuUqhoR+RCWumhFe7z3/xqzRilr7esAvArLVw2/nzHma3Nzc93ngXPu0d77i7EYKFVU9RBV/V/OuUucc3fprl8UxUlFUXxHRD4C4HkA/kBEziuK4tJGo9Fzmlmx1t7Y64nNJCJVVU3/aFgsRW5ENPwqWLwwFGDxD5Sq3jz4EBp2qjoRL7JEpKWq7Bq6w4nIiC6t3NsEsD95+lAsfhfMIrkTj8U+9SNYXByskezfg8W7/Q0kC4cBmACQiYhLZv4BFvvy51i8uEtn+8lVdVxEmsnFXGc/FgOQfcn+UVXdo6qaDPDtuz/8TLHP/i0IKzyn+1X1lq4g4wgRgYjsTwIcQeg6o6r7uy7QDwFQFZG5dKYkAIcBMCIy0zXr0iQWAx6Hxf/v+DNM6eKik8v+r0WkrqpjPf6PJsP4g2XnCeeexGLwl/5fjAPIRWQhGUzd/X+xLxmIPobFQKsNIA16YhZEun7mejhGsfh/nQafh2AxiFwAMN1VH+H8sTtb+rNpeM6LyHjsxiUi05qslSEiY6oag0KLxfEbe8L/J1T1lfV6/e9jfefcOaoap1dtAvgigIdiKRD9WJ7nncyEtfZZAM7F8qAg9QPv/dljY2PXhfrPBfBerFwwEKE9+0XkpVmWfQTAqLX2ZSLyl+iahSvhsThw+0MA2qr6eBH5nlhrS79Lr6oQGbjeBRERERFR2dpYCnjvmef5L5xzd1bV7yIEMar61Hq9/onQ5eh8LGW13phl2auttU8K2YU4g9alAP639/5sEXkJloKF6wG8V0ROVNUnJm1QLGZiboPFsQ+pm7EYjKwY67IG7xVrbXcqj4hoQ4S7QGU3g4jKYbB4geOR3OkVERO6bSzbn9RXLN2F7+wPd+5X7O9RX8JzPfeHm5HtHvWB5A7z4svJmveH77vYbar7Z4trDLS7vhMH1u+xf+D/UZhGtFM/dG/pW797vyymd1f8bANet+f/XdIeVdXuQem9frZ+50nrp689qH78fUnf5059VfVdXdfjz5a2qde+xROJ/FTDDE0i8l0Ar1bVfwAQu/K8O8/zF8f6zrmXh+ejeSxlNQDgO1mWnSUis6H+S0Pmod8d81sAPDXP8/9S1apz7k8A/CVWrhcStVX1HZVK5QOtVqsiIs8XkRegd+bhkv8fa5lXz3cEozcAAAAASUVORK5CYII='
    elif(exp_type == 'SHAPEsingle'):
        diagram = b'iVBORw0KGgoAAAANSUhEUgAAAfQAAADPCAYAAAADHCNQAAAACXBIWXMAAC4jAAAuIwF4pT92AAAFHGlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUgWE1QIENvcmUgNS42LWMxNDIgNzkuMTYwOTI0LCAyMDE3LzA3LzEzLTAxOjA2OjM5ICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIgeG1sbnM6cGhvdG9zaG9wPSJodHRwOi8vbnMuYWRvYmUuY29tL3Bob3Rvc2hvcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RFdnQ9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZUV2ZW50IyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ0MgMjAxOCAoV2luZG93cykiIHhtcDpDcmVhdGVEYXRlPSIyMDIwLTA1LTIyVDA5OjM1OjAzKzAyOjAwIiB4bXA6TW9kaWZ5RGF0ZT0iMjAyMC0wNS0yMlQwOTozNjoyNSswMjowMCIgeG1wOk1ldGFkYXRhRGF0ZT0iMjAyMC0wNS0yMlQwOTozNjoyNSswMjowMCIgZGM6Zm9ybWF0PSJpbWFnZS9wbmciIHBob3Rvc2hvcDpDb2xvck1vZGU9IjMiIHBob3Rvc2hvcDpJQ0NQcm9maWxlPSJzUkdCIElFQzYxOTY2LTIuMSIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDo2YWRjMjU1Yi0xODhl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'

    diagram = base64.b64decode(diagram)
    filename = 'sequence.png'
    with open(filename, 'wb') as f:
        f.write(diagram)

    return()


def exp_layout(exp_type, shape_type):  #Generate a CHORUS Layout
    sequence_diagram(exp_type)
    if(exp_type == 'CHORUS') or (exp_type == 'CHORUS_cycled'):  # define CHORUS layout
        pulse1 = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                [sg.Input(key='tw1', size=(10,1))],
                [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                [sg.Input(default_text='0.265', key='rffactor1', size=(10,1))],
                [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse1')],
                [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var11_out')],
                [sg.Input(key='var11', size=(10,1))],
                [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var12_out')],
                [sg.Input(key='var12', size=(10,1))],
                [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff1_out')],
                [sg.Input(default_text='0.0', key='phaseoff1', size=(10,1))]]
        if(shape_type == 'supergaussian'):
            pulse2 = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                    [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                    [sg.Text('Parameter 1 of P2:'), sg.Text(size=(10,1), key='var21_out')],
                    [sg.Text('Parameter 2 of P2:'), sg.Text(size=(10,1), key='var22_out')],
                    [sg.Input(key='var22', size=(10,1))],
                    [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff2_out')],
                    [sg.Input(default_text='0.0', key='phaseoff2', size=(10,1))]]
        else:
            pulse2 = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                    [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                    [sg.Text('Parameter 1 of P2:'), sg.Text(size=(10,1), key='var21_out')],
                    [sg.Text('Parameter 2 of P2:'), sg.Text(size=(10,1), key='var22_out')],
                    [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff2_out')],
                    [sg.Input(default_text='0.0', key='phaseoff2', size=(10,1))]]
        delay  = [[sg.Text('Duration of D (us):'), sg.Text(size=(10,1), key='tau2_out')]]
        pulse3 = [[sg.Text('Duration of P3 (us):'), sg.Text(size=(10,1), key='tw3_out')],
                [sg.Input(key='tw3', size=(15,1))],
                [sg.Text('RF-Factor of P3:'), sg.Text(size=(10,1), key='rffactor3_out')],
                [sg.Input(default_text='1.0', key='rffactor3', size=(15,1))],
                [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse3')],
                [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var31_out')],
                [sg.Input(key='var31', size=(10,1))],
                [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var32_out')],
                [sg.Input(key='var32', size=(10,1))],
                [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff3_out')],
                [sg.Input(default_text='0.0', key='phaseoff3', size=(10,1))]]

        layout = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                [sg.Image(filename='./sequence.png', key='sequence')],
                # [sg.Button('', image_filename=sequence,
                #  button_color=(sg.theme_background_color(),sg.theme_background_color()),
                #  border_width=0, key='sequence')],
                [sg.Column(pulse1), sg.Column(pulse2), sg.Column(delay), sg.Column(pulse3)],
                [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta1_out')],
                [sg.Input(key='delta1', size=(15,1))],
                [sg.Text('Offset Stepsize (Hz):'), sg.Text(size=(10,1), key='ss_offset_out')],
                [sg.Input(key='ss_offset', size=(10,1))],
                [sg.Button('Start Simulation'), sg.Button('Exit')]]
    elif(exp_type == 'double_echo') or (exp_type == 'double_echo_cycled') or (exp_type == 'double_echo_zerophase'):
        pulse1  = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                   [sg.Input(key='tw1', size=(10,1))],
                   [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                   [sg.Input(key='rffactor1', size=(10,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse1')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var11_out')],
                   [sg.Input(key='var11', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var12_out')],
                   [sg.Input(key='var12', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta1_out')],
                   [sg.Input(key='delta1', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff1_out')],
                   [sg.Input(default_text='0.0', key='phaseoff1', size=(10,1))]]
        delay1  = [[sg.Text('Duration of D1 (us):'), sg.Text(size=(10,1), key='tau1_out')],
                   [sg.Input(default_text='0', key='tau1', size=(15,1))]]
        pulse2  = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                   [sg.Input(key='tw2', size=(15,1))],
                   [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                   [sg.Input(key='rffactor2', size=(15,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse2')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var21_out')],
                   [sg.Input(key='var21', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var22_out')],
                   [sg.Input(key='var22', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta2_out')],
                   [sg.Input(key='delta2', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff2_out')],
                   [sg.Input(default_text='0.0', key='phaseoff2', size=(10,1))]]
        delay2  = [[sg.Text('Duration of D2 (us):'), sg.Text(size=(10,1), key='tau2_out')],
                   [sg.Input(default_text='0', key='tau2', size=(15,1))]]
        pulse3  = [[sg.Text('Duration of P3 (us):'), sg.Text(size=(10,1), key='tw3_out')],
                   [sg.Input(key='tw3', size=(15,1))],
                   [sg.Text('RF-Factor of P3:'), sg.Text(size=(10,1), key='rffactor3_out')],
                   [sg.Input(key='rffactor3', size=(15,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse3')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var31_out')],
                   [sg.Input(key='var31', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var32_out')],
                   [sg.Input(key='var32', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta3_out')],
                   [sg.Input(key='delta3', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff2_out')],
                   [sg.Input(default_text='0.0', key='phaseoff2', size=(10,1))]]
        delay3  = [[sg.Text('Duration of D3 (us):'), sg.Text(size=(10,1), key='tau3_out')],
                   [sg.Input(default_text='0', key='tau3', size=(15,1))]]

        layout  = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                   [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                   [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                   [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                   [sg.Image(filename='./sequence.png', key='sequence')],
                #    [sg.Button('', image_data=sequence,
                #    button_color=(sg.theme_background_color(),sg.theme_background_color()),
                #    border_width=0, key='sequence')],
                   [sg.Column(pulse1), sg.Column(delay1), sg.Column(pulse2), sg.Column(delay2), sg.Column(pulse3), sg.Column(delay3)],
                   [sg.Text('Offset Stepsize (Hz):'), sg.Text(size=(10,1), key='ss_offset_out')],
                   [sg.Input(key='ss_offset', size=(10,1))],
                   [sg.Button('Start Simulation'), sg.Button('Exit')]
                  ]
    elif(exp_type == 'loadshape_double_echo'):
        pulse1  = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                   [sg.Input(key='tw1', size=(10,1))],
                   [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                   [sg.Input(key='rffactor1', size=(10,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse1')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var11_out')],
                   [sg.Input(key='var11', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var12_out')],
                   [sg.Input(key='var12', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta1_out')],
                   [sg.Input(key='delta1', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff1_out')],
                   [sg.Input(default_text='0.0', key='phaseoff1', size=(10,1))],
                   [sg.Text('loadshape:'), sg.Text(size=(10,1), key='var13_out')],
                   [sg.Input(default_text='./shape1', key='var13', size=(10,1))]]
        delay1  = [[sg.Text('Duration of D1 (us):'), sg.Text(size=(10,1), key='tau1_out')],
                   [sg.Input(default_text='0', key='tau1', size=(15,1))]]
        pulse2  = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                   [sg.Input(key='tw2', size=(15,1))],
                   [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                   [sg.Input(key='rffactor2', size=(15,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse2')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var21_out')],
                   [sg.Input(key='var21', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var22_out')],
                   [sg.Input(key='var22', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta2_out')],
                   [sg.Input(key='delta2', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff2_out')],
                   [sg.Input(default_text='0.0', key='phaseoff2', size=(10,1))],
                   [sg.Text('loadshape:'), sg.Text(size=(10,1), key='var23_out')],
                   [sg.Input(default_text='./shape2', key='var23', size=(10,1))]]
        delay2  = [[sg.Text('Duration of D2 (us):'), sg.Text(size=(10,1), key='tau2_out')],
                   [sg.Input(default_text='0', key='tau2', size=(15,1))]]
        pulse3  = [[sg.Text('Duration of P3 (us):'), sg.Text(size=(10,1), key='tw3_out')],
                   [sg.Input(key='tw3', size=(15,1))],
                   [sg.Text('RF-Factor of P3:'), sg.Text(size=(10,1), key='rffactor3_out')],
                   [sg.Input(key='rffactor3', size=(15,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='pulse3')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var31_out')],
                   [sg.Input(key='var31', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var32_out')],
                   [sg.Input(key='var32', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta3_out')],
                   [sg.Input(key='delta3', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff3_out')],
                   [sg.Input(default_text='0.0', key='phaseoff3', size=(10,1))],
                   [sg.Text('loadshape:'), sg.Text(size=(10,1), key='var33_out')],
                   [sg.Input(default_text='./shape3', key='var33', size=(10,1))]]
        delay3  = [[sg.Text('Duration of D3 (us):'), sg.Text(size=(10,1), key='tau3_out')],
                   [sg.Input(default_text='0', key='tau3', size=(15,1))]]

        layout  = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                   [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                   [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                   [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                   [sg.Image(filename='./sequence.png', key='sequence')],
                #    [sg.Button('', image_data=sequence,
                #    button_color=(sg.theme_background_color(),sg.theme_background_color()),
                #    border_width=0, key='sequence')],
                   [sg.Column(pulse1), sg.Column(delay1), sg.Column(pulse2), sg.Column(delay2), sg.Column(pulse3), sg.Column(delay3)],
                   [sg.Text('Offset Stepsize (Hz):'), sg.Text(size=(10,1), key='ss_offset_out')],
                   [sg.Input(key='ss_offset', size=(10,1))],
                   [sg.Button('Start Simulation'), sg.Button('Exit')]
                  ]
    elif(exp_type == 'double_chirp'):
        pulse1  = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                   [sg.Input(key='tw1', size=(10,1))],
                   [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                   [sg.Input(default_text='0.265', key='rffactor1', size=(10,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='shapes')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var11_out')],
                   [sg.Input(key='var11', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var12_out')],
                   [sg.Input(key='var12', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta1_out')],
                   [sg.Input(key='delta1', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff1_out')],
                   [sg.Input(default_text='0.0', key='phaseoff1', size=(10,1))]]
        pulse2  = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                   [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                   [sg.Input(default_text='1.0', key='rffactor2', size=(10,1))],
                   [sg.Text('Parameter 1 of P2:'), sg.Text(size=(10,1), key='var21_out')],
                   [sg.Text('Parameter 2 of P2:'), sg.Text(size=(10,1), key='var22_out')],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta2_out')],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff2_out')],
                   [sg.Input(default_text='0.0', key='phaseoff2', size=(10,1))]]
        delay1  = [[sg.Text('Duration of D1 (us):'), sg.Text(size=(10,1), key='tau1_out')]]

        layout  = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                   [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                   [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                   [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                   [sg.Image(filename='./sequence.png', key='sequence')],
                #    [sg.Button('', image_data=sequence,
                #    button_color=(sg.theme_background_color(),sg.theme_background_color()),
                #    border_width=0, key='sequence')],
                   [sg.Column(pulse1), sg.Column(pulse2), sg.Column(delay1)],
                   [sg.Text('Offset Stepsize (Hz):'), sg.Text(size=(10,1), key='ss_offset_out')],
                   [sg.Input(key='ss_offset', size=(10,1))],
                   [sg.Button('Start Simulation'), sg.Button('Exit')]
                  ]
    elif(exp_type == 'SHAPEsingle'):
        pulse1  = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                   [sg.Input(key='tw1', size=(10,1))],
                   [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                   [sg.Input(key='rffactor1', size=(10,1))],
                   [sg.Text('Parameter for WURST/TANH/HS:'), sg.Text(size=(10,1), key='shapes')],
                   [sg.Text('N/Zeta/Beta:'), sg.Text(size=(10,1), key='var11_out')],
                   [sg.Input(key='var11', size=(10,1))],
                   [sg.Text('none/tan_kappa/none:'), sg.Text(size=(10,1), key='var12_out')],
                   [sg.Input(key='var12', size=(10,1))],
                   [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta1_out')],
                   [sg.Input(key='delta1', size=(15,1))],
                   [sg.Text('Additional Phase Offset:'), sg.Text(size=(10,1), key='phaseoff1_out')],
                   [sg.Input(default_text='0.0', key='phaseoff1', size=(10,1))]]

        layout  = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                   [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                   [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                   [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                   [sg.Image(filename='./sequence.png', key='sequence')],
                #    [sg.Button('', image_data=sequence,
                #    button_color=(sg.theme_background_color(),sg.theme_background_color()),
                #    border_width=0, key='sequence')],
                   [sg.Column(pulse1)],
                   [sg.Text('Offset Stepsize (Hz):'), sg.Text(size=(10,1), key='ss_offset_out')],
                   [sg.Input(key='ss_offset', size=(10,1))],
                   [sg.Button('Start Simulation'), sg.Button('Exit')]
                  ]
    return layout


def selection_layout():  # Generates a selection Layout
    layout = [[sg.Text('Select Shape:'), sg.Combo(['WURST', 'caWURST', 'supergaussian', 'tanhpulse', 'hspulse', 'loadshape'], key='shape_type')],
              [sg.Text('Select Sequence:'), sg.Combo(['SHAPEsingle', 'double_echo', 'CHORUS', 'double_chirp', 'double_echo_zerophase', 'double_echo_cycled', 'CHORUS_cycled', 'loadshape_double_echo'], key='exp_type')],
              [sg.Button('Continue'), sg.Button('Exit')]]
    return(layout)


def read_parameter(exp_type, shape_type):  # Update the "output" text element to be the value of "input" element and set input to parameter
    global simulation_window
    global delta1
    global delta2
    global delta3
    global tw1
    global tw2
    global tw3
    global rffactor1
    global rffactor2
    global rffactor3
    global tau1
    global tau2
    global tau3
    global var11
    global var12
    global var13
    global ss_offset
    global var21
    global var22
    global var23
    global filename
    global var31
    global var32
    global var33
    global phaseoff1
    global phaseoff2
    global phaseoff3

    if(exp_type == 'CHORUS') or (exp_type == 'CHORUS_cycled'):  # define CHORUS parameter

        delta1       = simulation_values['delta1']
        delta2       = simulation_values['delta1']
        delta3       = simulation_values['delta1']
        tw1          = simulation_values['tw1']
        tw3          = simulation_values['tw3']
        rffactor1    = simulation_values['rffactor1']
        rffactor3    = simulation_values['rffactor3']
        var11        = simulation_values['var11']
        var12        = simulation_values['var12']
        var31        = simulation_values['var31']
        var32        = simulation_values['var32']
        ss_offset    = simulation_values['ss_offset']
        tau2         = '{:.1f}'.format(float(tw1)/2.0)
        tw2          = '{:.1f}'.format(float(tau2)+float(tw3))
        rffactor2    = rffactor3
        if(shape_type == 'supergaussian'):
            var21           = var31
            var22           = simulation_values['var22']
        else:
            var21           = var31
            var22           = var32
        phaseoff1       = simulation_values['phaseoff1']
        phaseoff2       = simulation_values['phaseoff2']
        phaseoff3       = simulation_values['phaseoff3']

        simulation_window['phaseoff1_out'].update(simulation_values['phaseoff1'])
        simulation_window['phaseoff2_out'].update(simulation_values['phaseoff2'])
        simulation_window['phaseoff3_out'].update(simulation_values['phaseoff3'])
        simulation_window['delta1_out'].update(simulation_values['delta1'])
        simulation_window['tw1_out'].update(simulation_values['tw1'])
        simulation_window['tw2_out'].update(tw2)
        simulation_window['tw3_out'].update(simulation_values['tw3'])
        simulation_window['rffactor1_out'].update(simulation_values['rffactor1'])
        simulation_window['rffactor2_out'].update(rffactor2)
        simulation_window['rffactor3_out'].update(simulation_values['rffactor3'])
        simulation_window['tau2_out'].update(tau2)
        simulation_window['var11_out'].update(simulation_values['var11'])
        simulation_window['var12_out'].update(simulation_values['var12'])
        if(shape_type == 'supergaussian'):
            simulation_window['var21_out'].update(var21)
            simulation_window['var22_out'].update(simulation_values['var22'])
        else:
            simulation_window['var21_out'].update(var21)
            simulation_window['var22_out'].update(var22)
        simulation_window['var31_out'].update(simulation_values['var31'])
        simulation_window['var32_out'].update(simulation_values['var32'])
        simulation_window['ss_offset_out'].update(simulation_values['ss_offset'])

        simpson_info = (f"Experiment = {exp_type} \n"
                        f"Shape = {shape_type} \n"
                        f"Offset Stepsize = {ss_offset} \n"
                        f"delta1 = {delta1} \n"
                        f"delta2 = {delta2} \n"
                        f"delta3 = {delta3} \n"
                        f"tw1 = {tw1} \n"
                        f"tw2 = {tw2} \n"
                        f"tw3 = {tw3} \n"
                        f"rffactor1 = {rffactor1} \n"
                        f"rffactor2 = {rffactor2} \n"
                        f"rffactor3 = {rffactor3} \n"
                        f"phaseoff1 = {phaseoff1} \n"
                        f"phaseoff2 = {phaseoff2} \n"
                        f"phaseoff3 = {phaseoff3} \n"
                        f"tau2 = {tau2} \n"
                        f"Parameter: \n"
                        f"If shape_type is WURST: \n"
                        f"varX1 is N and varX2 is none \n"
                        f"If shape_type is tanhpulse: \n"
                        f"varX1 is Zeta and varX2 is tan_kappa \n"
                        f"var11 = {var11} \n"
                        f"var12 = {var12} \n"
                        f"var21 = {var21} \n"
                        f"var22 = {var22} \n"
                        f"var31 = {var31} \n"
                        f"var32 = {var32} \n")
    elif(exp_type == 'double_echo') or (exp_type == 'double_echo_cycled') or (exp_type == 'double_echo_zerophase'):  # define double_echo parameter
        simulation_window['delta1_out'].update(simulation_values['delta1'])
        simulation_window['delta2_out'].update(simulation_values['delta2'])
        simulation_window['delta3_out'].update(simulation_values['delta3'])
        simulation_window['tw1_out'].update(simulation_values['tw1'])
        simulation_window['tw2_out'].update(simulation_values['tw2'])
        simulation_window['tw3_out'].update(simulation_values['tw3'])
        simulation_window['rffactor1_out'].update(simulation_values['rffactor1'])
        simulation_window['rffactor2_out'].update(simulation_values['rffactor2'])
        simulation_window['rffactor3_out'].update(simulation_values['rffactor3'])
        simulation_window['tau1_out'].update(simulation_values['tau1'])
        simulation_window['tau2_out'].update(simulation_values['tau2'])
        simulation_window['tau3_out'].update(simulation_values['tau3'])
        simulation_window['var11_out'].update(simulation_values['var11'])
        simulation_window['var12_out'].update(simulation_values['var12'])
        simulation_window['var21_out'].update(simulation_values['var21'])
        simulation_window['var22_out'].update(simulation_values['var22'])
        simulation_window['var31_out'].update(simulation_values['var31'])
        simulation_window['var32_out'].update(simulation_values['var32'])
        simulation_window['ss_offset_out'].update(simulation_values['ss_offset'])
        simulation_window['phaseoff1_out'].update(simulation_values['phaseoff1'])
        simulation_window['phaseoff2_out'].update(simulation_values['phaseoff2'])
        simulation_window['phaseoff3_out'].update(simulation_values['phaseoff3'])

        delta1      = simulation_values['delta1']
        delta2      = simulation_values['delta2']
        delta3      = simulation_values['delta3']
        tw1         = simulation_values['tw1']
        tw2         = simulation_values['tw2']
        tw3         = simulation_values['tw3']
        rffactor1   = simulation_values['rffactor1']
        rffactor2   = simulation_values['rffactor2']
        rffactor3   = simulation_values['rffactor3']
        tau1        = simulation_values['tau1']
        tau2        = simulation_values['tau2']
        tau3        = simulation_values['tau3']
        var11           = simulation_values['var11']
        var12           = simulation_values['var12']
        var21           = simulation_values['var21']
        var22           = simulation_values['var22']
        var31           = simulation_values['var31']
        var32           = simulation_values['var32']
        ss_offset    = simulation_values['ss_offset']
        phaseoff1       = simulation_values['phaseoff1']
        phaseoff2       = simulation_values['phaseoff2']
        phaseoff3       = simulation_values['phaseoff3']

        simpson_info = (f"Experiment = {exp_type} \n"
                        f"Shape = {shape_type} \n"
                        f"Offset Stepsize = {ss_offset} \n"
                        f"delta1 = {delta1} \n"
                        f"delta2 = {delta2} \n"
                        f"delta3 = {delta3} \n"
                        f"tw1 = {tw1} \n"
                        f"tw2 = {tw2} \n"
                        f"tw3 = {tw3} \n"
                        f"rffactor1 = {rffactor1} \n"
                        f"rffactor2 = {rffactor2} \n"
                        f"rffactor3 = {rffactor3} \n"
                        f"phaseoff1 = {phaseoff1} \n"
                        f"phaseoff2 = {phaseoff2} \n"
                        f"phaseoff3 = {phaseoff3} \n"
                        f"tau1 = {tau1} \n"
                        f"tau2 = {tau2} \n"
                        f"tau3 = {tau3} \n"
                        f"Parameter: \n"
                        f"If shape_type is WURST: \n"
                        f"varX1 is N and varX2 is none \n"
                        f"If shape_type is tanhpulse: \n"
                        f"varX1 is Zeta and varX2 is tan_kappa \n"
                        f"var11 = {var11} \n"
                        f"var12 = {var12} \n"
                        f"var21 = {var21} \n"
                        f"var22 = {var22} \n"
                        f"var31 = {var31} \n"
                        f"var32 = {var32} \n")
    elif(exp_type == 'loadshape_double_echo'):  # define double_echo parameter
        simulation_window['delta1_out'].update(simulation_values['delta1'])
        simulation_window['delta2_out'].update(simulation_values['delta2'])
        simulation_window['delta3_out'].update(simulation_values['delta3'])
        simulation_window['tw1_out'].update(simulation_values['tw1'])
        simulation_window['tw2_out'].update(simulation_values['tw2'])
        simulation_window['tw3_out'].update(simulation_values['tw3'])
        simulation_window['rffactor1_out'].update(simulation_values['rffactor1'])
        simulation_window['rffactor2_out'].update(simulation_values['rffactor2'])
        simulation_window['rffactor3_out'].update(simulation_values['rffactor3'])
        simulation_window['tau1_out'].update(simulation_values['tau1'])
        simulation_window['tau2_out'].update(simulation_values['tau2'])
        simulation_window['tau3_out'].update(simulation_values['tau3'])
        simulation_window['var11_out'].update(simulation_values['var11'])
        simulation_window['var12_out'].update(simulation_values['var12'])
        simulation_window['var13_out'].update(simulation_values['var13'])
        simulation_window['var21_out'].update(simulation_values['var21'])
        simulation_window['var22_out'].update(simulation_values['var22'])
        simulation_window['var23_out'].update(simulation_values['var23'])
        simulation_window['var31_out'].update(simulation_values['var31'])
        simulation_window['var32_out'].update(simulation_values['var32'])
        simulation_window['var33_out'].update(simulation_values['var33'])
        simulation_window['ss_offset_out'].update(simulation_values['ss_offset'])
        simulation_window['phaseoff1_out'].update(simulation_values['phaseoff1'])
        simulation_window['phaseoff2_out'].update(simulation_values['phaseoff2'])
        simulation_window['phaseoff3_out'].update(simulation_values['phaseoff3'])

        delta1      = simulation_values['delta1']
        delta2      = simulation_values['delta2']
        delta3      = simulation_values['delta3']
        tw1         = simulation_values['tw1']
        tw2         = simulation_values['tw2']
        tw3         = simulation_values['tw3']
        rffactor1   = simulation_values['rffactor1']
        rffactor2   = simulation_values['rffactor2']
        rffactor3   = simulation_values['rffactor3']
        tau1        = simulation_values['tau1']
        tau2        = simulation_values['tau2']
        tau3        = simulation_values['tau3']
        var11           = simulation_values['var11']
        var12           = simulation_values['var12']
        var13           = simulation_values['var13']
        var21           = simulation_values['var21']
        var22           = simulation_values['var22']
        var23           = simulation_values['var23']
        var31           = simulation_values['var31']
        var32           = simulation_values['var32']
        var33           = simulation_values['var33']
        ss_offset    = simulation_values['ss_offset']
        phaseoff1       = simulation_values['phaseoff1']
        phaseoff2       = simulation_values['phaseoff2']
        phaseoff3       = simulation_values['phaseoff3']

        simpson_info = (f"Experiment = {exp_type} \n"
                        f"Shape = {shape_type} \n"
                        f"Offset Stepsize = {ss_offset} \n"
                        f"delta1 = {delta1} \n"
                        f"delta2 = {delta2} \n"
                        f"delta3 = {delta3} \n"
                        f"tw1 = {tw1} \n"
                        f"tw2 = {tw2} \n"
                        f"tw3 = {tw3} \n"
                        f"rffactor1 = {rffactor1} \n"
                        f"rffactor2 = {rffactor2} \n"
                        f"rffactor3 = {rffactor3} \n"
                        f"phaseoff1 = {phaseoff1} \n"
                        f"phaseoff2 = {phaseoff2} \n"
                        f"phaseoff3 = {phaseoff3} \n"
                        f"tau1 = {tau1} \n"
                        f"tau2 = {tau2} \n"
                        f"tau3 = {tau3} \n"
                        f"Parameter: \n"
                        f"If shape_type is WURST: \n"
                        f"varX1 is N and varX2 is none \n"
                        f"If shape_type is tanhpulse: \n"
                        f"varX1 is Zeta and varX2 is tan_kappa \n"
                        f"var11 = {var11} \n"
                        f"var12 = {var12} \n"
                        f"var12 = {var13} \n"
                        f"var21 = {var21} \n"
                        f"var22 = {var22} \n"
                        f"var22 = {var23} \n"
                        f"var31 = {var31} \n"
                        f"var32 = {var32} \n"
                        f"var32 = {var33} \n")
    elif(exp_type == 'SHAPEsingle'):  # define SHAPEsingle parameter
        simulation_window['delta1_out'].update(simulation_values['delta1'])
        simulation_window['tw1_out'].update(simulation_values['tw1'])
        simulation_window['rffactor1_out'].update(simulation_values['rffactor1'])
        simulation_window['var11_out'].update(simulation_values['var11'])
        simulation_window['var12_out'].update(simulation_values['var12'])
        simulation_window['ss_offset_out'].update(simulation_values['ss_offset'])
        simulation_window['phaseoff1_out'].update(simulation_values['phaseoff1'])

        delta1      = simulation_values['delta1']
        tw1         = simulation_values['tw1']
        rffactor1   = simulation_values['rffactor1']
        var11           = simulation_values['var11']
        var12           = simulation_values['var12']
        ss_offset    = simulation_values['ss_offset']
        phaseoff1       = simulation_values['phaseoff1']

        simpson_info = (f"Experiment = {exp_type} \n"
                        f"Shape = {shape_type} \n"
                        f"Offset Stepsize = {ss_offset} \n"
                        f"delta1 = {delta1} \n"
                        f"tw1 = {tw1} \n"
                        f"rffactor1 = {rffactor1} \n"
                        f"phaseoff1 = {phaseoff1} \n"
                        f"Parameter: \n"
                        f"If shape_type is WURST: \n"
                        f"varX1 is N and varX2 is none \n"
                        f"If shape_type is tanhpulse: \n"
                        f"varX1 is Zeta and varX2 is tan_kappa \n"
                        f"var11 = {var11} \n"
                        f"var12 = {var12} \n")
    elif(exp_type == 'double_chirp'):  # define double_chirp parameter
        delta1      = simulation_values['delta1']
        delta2      = delta1
        tw1         = simulation_values['tw1']
        tw2         = '{:.1f}'.format(float(tw1)/2.0)
        rffactor1   = simulation_values['rffactor1']
        rffactor2   = simulation_values['rffactor2']
        tau1        = tw2
        var11           = simulation_values['var11']
        var12           = simulation_values['var12']
        var21           = var11
        var22           = var12
        ss_offset    = simulation_values['ss_offset']
        phaseoff1       = simulation_values['phaseoff1']
        phaseoff2       = simulation_values['phaseoff2']

        simulation_window['delta1_out'].update(simulation_values['delta1'])
        simulation_window['delta2_out'].update(delta2)
        simulation_window['tw1_out'].update(simulation_values['tw1'])
        simulation_window['tw2_out'].update(tw2)
        simulation_window['rffactor1_out'].update(simulation_values['rffactor1'])
        simulation_window['rffactor2_out'].update(simulation_values['rffactor2'])
        simulation_window['tau1_out'].update(tau1)
        simulation_window['var11_out'].update(simulation_values['var11'])
        simulation_window['var12_out'].update(simulation_values['var12'])
        simulation_window['var21_out'].update(var21)
        simulation_window['var22_out'].update(var22)
        simulation_window['ss_offset_out'].update(simulation_window['ss_offset'])
        simulation_window['phaseoff1_out'].update(simulation_values['phaseoff1'])
        simulation_window['phaseoff2_out'].update(simulation_values['phaseoff2'])

        simpson_info = (f"Experiment = {exp_type} \n"
                        f"Shape = {shape_type} \n"
                        f"Offset Stepsize = {ss_offset} \n"
                        f"delta1 = {delta1} \n"
                        f"delta2 = {delta2} \n"
                        f"tw1 = {tw1} \n"
                        f"tw2 = {tw2} \n"
                        f"rffactor1 = {rffactor1} \n"
                        f"rffactor2 = {rffactor2} \n"
                        f"phaseoff1 = {phaseoff1} \n"
                        f"phaseoff2 = {phaseoff2} \n"
                        f"tau1 = {tau1} \n"
                        f"Parameter: \n"
                        f"If shape_type is WURST: \n"
                        f"varX1 is N and varX2 is none \n"
                        f"If shape_type is tanhpulse: \n"
                        f"varX1 is Zeta and varX2 is tan_kappa \n"
                        f"var11 = {var11} \n"
                        f"var12 = {var12} \n"
                        f"var21 = {var21} \n"
                        f"var22 = {var22} \n")

    filename = str(datetime.today()).replace('-','').replace(':','').replace(' ','').split('.')[0] + '_' + exp_type + '_' + shape_type
    simpson_infofile = open(filename + '.info', 'w')
    simpson_infofile.write(simpson_info)
    simpson_infofile.close()

# Declare Parameter
filename     = '0'
delta1       = '0'
delta2       = '0'
delta3       = '0'
tw1          = '0'
tw2          = '0'
tw3          = '0'
rffactor1    = '0'
rffactor2    = '0'
rffactor3    = '0'
tau1         = '0'
tau2         = '0'
tau3         = '0'
ss_offset    = '0'
var11        = '0'
var12        = '0'
var13        = '0'
var21        = '0'
var22        = '0'
var23        = '0'
var31        = '0'
var32        = '0'
var33        = '0'
phaseoff1    = '0'
phaseoff2    = '0'
phaseoff3    = '0'
# Create Sequence Selection Window
sequence_selection = sg.Window('Phasecorrection Tool', selection_layout(), grab_anywhere=True)

while True:  # Event loop creating windows
    event, values = sequence_selection.read()   # Read the event that happened and the values dictionary
    print(event, values)
    if event == sg.WIN_CLOSED or event == 'Exit':  # If user closed window with X or if user clicked "Exit" button then exit
        break
    if event == 'Continue':  # Starting new window with selected experiment
        simulation_window = sg.Window('Simulation Parameter', exp_layout(values['exp_type'], values['shape_type']), grab_anywhere=True)
        while True:
            simulation_event, simulation_values = simulation_window.read()   # Read the event that happened and the values dictionary
            print(simulation_event, simulation_values)
            if simulation_event == sg.WIN_CLOSED or simulation_event == 'Exit':     # If user closed window with X or if user clicked "Exit" button then exit
                run(['rm', '-f', 'sequence.png'])
                break
            if simulation_event == 'Start Simulation':  # Starting simulations with input parameter
                read_parameter(values['exp_type'], values['shape_type'])
                create_simpson()
                simulate_sequence(values['exp_type'], values['shape_type'])
    simulation_window.close()
sequence_selection.close()