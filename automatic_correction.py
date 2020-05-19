#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:42:31 2019

@author: m_buss13
"""

#%%
###############################################################################
# Standard parameter
###############################################################################
import os
from subprocess import call
from subprocess import run
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy import interpolate
import PySimpleGUI as sg

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

            # Experiment selection done by scanning for the type string
            if {[string equal $par(type) "double_echo"]} {
                pulse_shaped $par(tw1) $rfsh1
                delay $par(tau1)
                pulse_shaped $par(tw2) $rfsh2
                delay $par(tau2)
                pulse_shaped $par(tw3) $rfsh3
                delay $par(tau3)
            } elseif {[string equal $par(type) "double_echo_phasecorr"]} {
                pulse_shaped $par(tw1) $rfsh1
                delay $par(tau1)
                pulse_shaped $par(tw2) $rfsh2
                delay $par(tau2)
                pulse_shaped $par(tw3) $rfsh3
                delay $par(tau3)
            } elseif {[string equal $par(type) "create_shapes"]} {
            } else {
                puts "Please select excitation mode in main!"
                exit
            }

            store 1
            acq 2 1 $par(ph31)
        }

        proc main {} {
            global par rfsh1 rfsh2 rfsh3 argc argv

            # Read Arguments from commandline
            if { $argc != 15 } {
                puts "Wrong number of Inputs"
                puts "Please try again."
            } else {
                set par(tw1)                        [lindex $argv 1]
                set par(tw2)                        [lindex $argv 2]
                set par(tw3)                        [lindex $argv 3]
                set par(type)                       [lindex $argv 4]
                set par(Delta)                      [lindex $argv 5]
                set par(N)                          [lindex $argv 6]
                set par(rf_factor1)                 [lindex $argv 7]
                set par(rf_factor2)                 [lindex $argv 8]
                set par(rf_factor3)                 [lindex $argv 9]
                set par(tau1)                       [lindex $argv 10]
                set par(tau2)                       [lindex $argv 11]
                set par(tau3)                       [lindex $argv 12]
                set par(filename_phasecorrect)      [lindex $argv 13]
                set par(ss_offset)                  [lindex $argv 14]
            }

            set par(stepsize)   0.05

            set par(np_tau1)    [expr round($par(tau1)/$par(stepsize))]
            set par(np_tau2)    [expr round($par(tau2)/$par(stepsize))]
            set par(np_tau3)    [expr round($par(tau3)/$par(stepsize))]

            #Define start and end values of constant list in Hz
            set par(start_offset)   [expr -1*($par(Delta)*1000/2.0)]
            set par(end_offset)     [expr    ($par(Delta)*1000/2.0)]
            set offset_value_list   [list_offset]

            set results {}

            set par(phasecycles) 16
            set par(ph1_list)  { 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 }
            set par(ph2_list)  { 0 90 180 270 0 90 180 270 0 90 180 270 0 90 180 270 }
            set par(ph3_list)  { 0 0 0 0 90 90 90 90 180 180 180 180 270 270 270 270 }
            set par(ph31_list) { 0 180 0 180 180 0 180 0 0 180 0 180 180 0 180 0 }

            if {[string equal $par(type) "create_shapes"]} {
                set offset_value_list { 0 }
                set par(phasecycles) 1
            }

        for {set index 0} {$index<$par(phasecycles)} {incr index} {
            set par(ph1) [lindex $par(ph1_list) $index]
            set par(ph2) [lindex $par(ph2_list) $index]
            set par(ph3) [lindex $par(ph3_list) $index]
            set par(ph31) [lindex $par(ph31_list) $index]

            # Set shapes for WURST type experiments
            if {[string equal $par(type) "double_echo"]} {
                # Set first WURST pulse (excitation)
                set par(sweep_rate1) [expr ($par(Delta)*1e3)/($par(tw1)*1e-6)]
                set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
                set rfsh1 [list2shape [wurst $par(N) $par(tw1) $par(Delta) $par(rf1) $par(ph1) $par(stepsize)]]

                # Set second WURST pulse (refocussing)
                set par(sweep_rate2) [expr ($par(Delta)*1e3)/($par(tw2)*1e-6)]
                set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
                set rfsh2 [list2shape [wurst $par(N) $par(tw2) $par(Delta) $par(rf2) $par(ph2) $par(stepsize)]]

                # Set third WURST pulse (refocussing)
                set par(sweep_rate3) [expr ($par(Delta)*1e3)/($par(tw3)*1e-6)]
                set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
                set rfsh3 [list2shape [wurst $par(N) $par(tw3) $par(Delta) $par(rf3) $par(ph3) $par(stepsize)]]

                # Set filenames for WURST type experiments
                set par(filename) $par(name)_rffactor_$par(rf_factor1)_$par(rf_factor2)_$par(rf_factor3)_tw_$par(tw1)_$par(tw2)_$par(tw3)_delays_$par(tau1)_$par(tau2)_$par(tau3)_delta_$par(Delta)_N_$par(N)
            } elseif {[string equal $par(type) "create_shapes"]} {
                # Set first WURST pulse (excitation)
                set par(sweep_rate1) [expr ($par(Delta)*1e3)/($par(tw1)*1e-6)]
                set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
                set rfsh1 [list2shape [wurst_phasecorr $par(N) $par(tw1) $par(Delta) $par(rf1) $par(filename_phasecorrect) 0.0]]

                # Set second WURST pulse (refocussing)
                set par(sweep_rate2) [expr ($par(Delta)*1e3)/($par(tw2)*1e-6)]
                set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
                set rfsh2 [list2shape [wurst $par(N) $par(tw2) $par(Delta) $par(rf2) 0.0]]

                # Set third WURST pulse (refocussing)
                set par(sweep_rate3) [expr ($par(Delta)*1e3)/($par(tw3)*1e-6)]
                set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
                set rfsh3 [list2shape [wurst $par(N) $par(tw3) $par(Delta) $par(rf3) 0.0]]

                # Set filenames for WURST type experiments
                set par(filename) $par(name)_$par(type)_rffactor_$par(rf_factor1)_$par(rf_factor2)_$par(rf_factor3)_tw_$par(tw1)_$par(tw2)_$par(tw3)_delays_$par(tau1)_$par(tau2)_$par(tau3)_delta_$par(Delta)_N_$par(N)

                set rfsh_shape1 [wurst_phasecorr $par(N) $par(tw1) $par(Delta) 100 $par(filename_phasecorrect)]
                set rfsh_shape2 [wurst $par(N) $par(tw2) $par(Delta) 100]
                set rfsh_shape3 [wurst $par(N) $par(tw3) $par(Delta) 100]

                printwave $rfsh_shape1 1
                printwave $rfsh_shape2 2
                printwave $rfsh_shape3 3

                set rfsh11 [list2shape [wurst_phasecorr $par(N) $par(tw1) $par(Delta) 1 $par(filename_phasecorrect)]]
                set rfsh22 [list2shape [wurst $par(N) $par(tw2) $par(Delta) 1]]
                set rfsh33 [list2shape [wurst $par(N) $par(tw3) $par(Delta) 1]]

                save_shape $rfsh1 $par(filename).simpson1
                save_shape $rfsh2 $par(filename).simpson2
                save_shape $rfsh3 $par(filename).simpson3
            } elseif {[string equal $par(type) "double_echo_phasecorr"]} {
                # Set first WURST pulse (excitation)
                set par(sweep_rate1) [expr ($par(Delta)*1e3)/($par(tw1)*1e-6)]
                set par(rf1) [format "%.2f" [expr $par(rf_factor1)*sqrt($par(sweep_rate1))]]
                set rfsh1 [list2shape [wurst_phasecorr $par(N) $par(tw1) $par(Delta) $par(rf1) $par(filename_phasecorrect) $par(ph1) 0.05]]

                # Set second WURST pulse (refocussing)
                set par(sweep_rate2) [expr ($par(Delta)*1e3)/($par(tw2)*1e-6)]
                set par(rf2) [format "%.2f" [expr $par(rf_factor2)*sqrt($par(sweep_rate2))]]
                set rfsh2 [list2shape [wurst $par(N) $par(tw2) $par(Delta) $par(rf2) $par(ph2) 0.05]]

                # Set third WURST pulse (refocussing)
                set par(sweep_rate3) [expr ($par(Delta)*1e3)/($par(tw3)*1e-6)]
                set par(rf3) [format "%.2f" [expr $par(rf_factor3)*sqrt($par(sweep_rate3))]]
                set rfsh3 [list2shape [wurst $par(N) $par(tw3) $par(Delta) $par(rf3) $par(ph3) 0.05]]

                # Set filenames for WURST type experiments
                set par(filename) $par(name)_$par(type)_rffactor_$par(rf_factor1)_$par(rf_factor2)_$par(rf_factor3)_tw_$par(tw1)_$par(tw2)_$par(tw3)_delays_$par(tau1)_$par(tau2)_$par(tau3)_delta_$par(Delta)_N_$par(N)
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
            set xyproj 0.0

            foreach elem $results_sum_temp {
                set offs     [lindex $elem 0]
                set re       [expr $re+[lindex $elem 1]]
                set im       [expr $im+[lindex $elem 2]]
                set xyproj   [expr $xyproj+[lindex $elem 3]]
            }
            lappend results_sum [format "%s %s %s %s" $offs [expr $re/8.0] [expr $im/8.0] [expr $xyproj/8.0]]
        }

        # Write output
                set fileID [open $par(filename).out "w"]
                foreach l $results_sum {
                    puts $fileID $l
                }
                close $fileID
        }



        ###########################################################################
        # Proc for WURST shape calculation
        # Changed 20.01.2020 by Max Bußkamp:
        #   - Added Option for Phasecycle
        #   - Added rfmax to input variables
        #   - Added default values for stepsize, direction and offset
        ###########################################################################
        proc wurst {N tw Delta rfmax {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {
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
        # Proc for WURST shape calculation
        # Changed 09.07.2019 by Max Bußkamp:
        #   - Added rfmax to input variables
        #   - Added default values for stepsize, direction and offset
        ###########################################################################
        proc wurst_phasecorr {N tw Delta rfmax filename_phasecorrect {phaseoffset 0.0} {stepsize 0.05} {direction 0} {offset 0.0}} {
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
            
            puts $fp "##TITLE= WURST-${par(N)} shape (duration $par(tw1) us, sweep width $par(Delta) kHz, step size: 0.05 us, sweep rate: $par(sweep_rate1) MHz/ms)"
            puts $fp "##USAGE= WURST pulse for inversion and excitation"
            puts $fp "##JCAMP-DX= 5.00 \$\$ Bruker JCAMP library"
            puts $fp "##DATA TYPE= Shape Data"
            puts $fp "##ORIGIN= Generated from wurst program"
            puts $fp "##DATE= "
            puts $fp "##TIME= "
            puts $fp "##\$SHAPE_PARAMETERS= Type: Wurst ; Total Sweep-Width \[Hz\] [expr $par(Delta)*1000.0] ; Length of Pulse \[usec\] [expr 1.0*$par(tw1)] ; Amplitude Power Index ${par(N)}.0 ; 1=High to low field, -1=Low to high field 1"
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


def simulate_sequence():  # Start simulation with set parameter
    #######################################################################
        # Calculate Phasecorrection for given pulselength from simulation
        ###############################################################################
        run(['simpson',
                'phasecorrection_liquid.tcl',
                tw1,
                tw2,
                tw3,
                'double_echo',
                delta,
                N,
                rffactor1,
                rffactor2,
                rffactor3,
                tau1,
                tau2,
                tau3,
                'none',
                ss_offset])

        filename_phasecorr = ('phasecorrection_liquid_' +
                            'rffactor_' + rffactor1 + '_' + rffactor2 + '_' + rffactor3 +
                            '_tw_' + tw1 + '_' + tw2 + '_' + tw3 +
                            '_delays_' + tau1 + '_' + tau2 + '_' + tau3 +
                            '_delta_' + delta +
                            '_N_' + N +
                            '.out')

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
        weight = []
        for i in range(len(phase[:, 0])):
            if(i < len(phase[:, 0])*start):
                weight.append(outer_weight)
            elif(i > len(phase[:, 0])*end):
                weight.append(outer_weight)
            else:
                weight.append(1.0)

        coefs = poly.polyfit(phase[:, 0], phase[:, 1], poly_order, w=weight)
        ffit = poly.Polynomial(coefs)

        np_tw = round(float(tw1)/0.05)
        pulselength = np.linspace(phase[0, 0], phase[-1, 0], np_tw)

        phase_interpol = ffit(pulselength)

        phase_interpol2 = []
        phase_interpol2.append(pulselength)
        phase_interpol2.append(phase_interpol)

        weight = []
        for i in range(len(phase_interpol2[0])):
            if(i < len(phase_interpol2[0])*start):
                weight.append(outer_weight)
            elif(i > len(phase_interpol2[0])*end):
                weight.append(outer_weight)
            else:
                weight.append(1.0)

        phase_interpol2[1]=phase_interpol2[1]*weight

        interpol_org = interpolate.interp1d(phase[:, 0], phase[:, 1])
        interpol_org = interpol_org(phase_interpol2[0])

        rms_start = 0.0
        for i in range(len(phase_interpol2[0])):
            if(i < len(phase_interpol2[0])*start):
                rms_start = rms_start + 0.0
            elif(i > len(phase_interpol2[0])*end):
                rms_start = rms_start + 0.0
            else:
                rms_start = rms_start + abs(phase_interpol2[1][i] - interpol_org[i])
        rms = rms_start

        while abs(rms_start-rms) < rms_limit:
            rms_start = rms
            start = start - rmssteps
            end = end + rmssteps
            weight = []
            for i in range(len(phase[:, 0])):
                if(i < len(phase[:, 0])*start):
                    weight.append(outer_weight)
                elif(i > len(phase[:, 0])*end):
                    weight.append(outer_weight)
                else:
                    weight.append(1.0)

            coefs = poly.polyfit(phase[:, 0], phase[:, 1], poly_order, w=weight)
            ffit = poly.Polynomial(coefs)

            pulselength = np.linspace(phase[0, 0], phase[-1, 0], np_tw)

            phase_interpol = ffit(pulselength)

            phase_interpol2 = []
            phase_interpol2.append(pulselength)
            phase_interpol2.append(phase_interpol)

            weight = []
            for i in range(len(phase_interpol2[0])):
                if(i < len(phase_interpol2[0])*start):
                    weight.append(outer_weight)
                elif(i > len(phase_interpol2[0])*end):
                    weight.append(outer_weight)
                else:
                    weight.append(1.0)

            phase_interpol2[1]=phase_interpol2[1]*weight

            rms = 0.0
            for i in range(len(phase_interpol2[0])):
                if(i < len(phase_interpol2[0])*start):
                    rms = rms + 0.0
                elif(i > len(phase_interpol2[0])*end):
                    rms = rms + 0.0
                else:
                    rms = rms + abs(phase_interpol2[1][i] - interpol_org[i])
            if (start < 0.0-rmssteps):
                break

        start = start + rmssteps
        end = end - rmssteps
        weight = []
        for i in range(len(phase[:, 0])):
            if(i < len(phase[:, 0])*start):
                weight.append(outer_weight)
            elif(i > len(phase[:, 0])*end):
                weight.append(outer_weight)
            else:
                weight.append(1.0)

        coefs = poly.polyfit(phase[:, 0], phase[:, 1], poly_order, w=weight)
        ffit = poly.Polynomial(coefs)

        pulselength = np.linspace(phase[0, 0], phase[-1, 0], np_tw)

        phase_interpol = ffit(pulselength)

        phase_interpol2 = []
        phase_interpol2.append(pulselength)
        phase_interpol2.append(phase_interpol)

        weight = []
        for i in range(len(phase_interpol2[0])):
            if(i < len(phase_interpol2[0])*start):
                weight.append(outer_weight)
            elif(i > len(phase_interpol2[0])*end):
                weight.append(outer_weight)
            else:
                weight.append(1.0)

        phase_interpol2[1]=phase_interpol2[1]*weight

        rms = 0.0
        for i in range(len(phase_interpol2[0])):
            if(i < len(phase_interpol2[0])*start):
                rms = rms + 0.0
            elif(i > len(phase_interpol2[0])*end):
                rms = rms + 0.0
            else:
                rms = rms + abs(phase_interpol2[1][i] - interpol_org[i])

        ###############################################################################
        # Plotting
        ###############################################################################
        np.savetxt(filename_phasecorr.replace('.out', '.phasecorr'), phase_interpol2[1], delimiter=' ')

        ax1.plot(dat_phasecorr[:, 0]/1000, dat_phasecorr[:, 1], label=r'Real')
        ax1.plot(dat_phasecorr[:, 0]/1000, dat_phasecorr[:, 2], label=r'Imag')
        ax1.plot(dat_phasecorr[:, 0]/1000, np.sqrt(dat_phasecorr[:, 1] ** 2 + dat_phasecorr[:, 2] ** 2), label=r'Magnitude')
        ax1.legend(fontsize=6)
        ax1.invert_xaxis()
        ax1.set_ylabel('Magnetization / a.u.')
        ax2.plot(phase[:, 0]/1000, phase[:, 1], ls='-', c='k', label='Phase')
        ax2.plot(phase_interpol2[0][:]/1000, phase_interpol2[1][:], ".", ms=0.5, c='r', label='Interpolation')
        ax2.legend(fontsize=6)
        ax2.invert_xaxis()
        ax2.set_xlabel('Frequency Offset / kHz')
        ax2.set_ylabel('Phase / degree')
        plt.tight_layout()
        plt.savefig(filename_phasecorr.replace('.out', '_phase.png'))

        ##############################################################################
        # Create corrected Shape Files
        ##############################################################################
        run(['simpson',
            'phasecorrection_liquid.tcl',
            tw1,
            tw2,
            tw3,
            'create_shapes',
            delta,
            N,
            rffactor1,
            rffactor2,
            rffactor3,
            tau1,
            tau2,
            tau3,
            filename_phasecorr.replace('.out', '.phasecorr'),
            ss_offset])

        run(['rm', '-f', 'phasecorrected.spinsys'])
        run(['rm', '-f', 'offset.spinsys'])
        run(['rm', '-f', 'phasecorrection_liquid.tcl'])
        run(['rm', '-f', filename_phasecorr])
        run(['rm', '-f', filename_phasecorr.replace('phasecorrection_liquid', 'phasecorrection_liquid_create_shapes')])

        sg.popup('Simulation finished! Close to show Plot!')

        plt.show()


def sequence_diagram():  # Generate base64 image of pulsesequence
    diagram = b'iVBORw0KGgoAAAANSUhEUgAABQEAAADwCAYAAACwjKJ3AAAACXBIWXMAAC4jAAAuIwF4pT92AAAF7GlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUgWE1QIENvcmUgNS42LWMxNDIgNzkuMTYwOTI0LCAyMDE3LzA3LzEzLTAxOjA2OjM5ICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIgeG1sbnM6cGhvdG9zaG9wPSJodHRwOi8vbnMuYWRvYmUuY29tL3Bob3Rvc2hvcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RFdnQ9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZUV2ZW50IyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ0MgMjAxOCAoV2luZG93cykiIHhtcDpDcmVhdGVEYXRlPSIyMDIwLTA1LTE4VDE1OjI4OjUxKzAyOjAwIiB4bXA6TW9kaWZ5RGF0ZT0iMjAyMC0wNS0xOFQxNjoxOToxNyswMjowMCIgeG1wOk1ldGFkYXRhRGF0ZT0iMjAyMC0wNS0xOFQxNjoxOToxNyswMjowMCIgZGM6Zm9ybWF0PSJpbWFnZS9wbmciIHBob3Rvc2hvcDpDb2xvck1vZGU9IjMiIHBob3Rvc2hvcDpJQ0NQcm9maWxlPSJzUkdCIElFQzYxOTY2LTIuMSIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDo4MGMxMTExMy01NjFjLWE0NGItYTAxMC05OWZmY2UyNTlkYTMiIHhtcE1NOkRvY3VtZW50SUQ9InhtcC5kaWQ6ZWJkOTZkM2EtYWYyMi1lYzQ3LTg1ZWUtN2JiNmUzNDk0MTcxIiB4bXBNTTpPcmlnaW5hbERvY3VtZW50SUQ9InhtcC5kaWQ6ZWJkOTZkM2EtYWYyMi1lYzQ3LTg1ZWUtN2JiNmUzNDk0MTcxIj4gPHhtcE1NOkhpc3Rvcnk+IDxyZGY6U2VxPiA8cmRmOmxpIHN0RXZ0OmFjdGlvbj0iY3JlYXRlZCIgc3RFdnQ6aW5zdGFuY2VJRD0ieG1wLmlpZDplYmQ5NmQzYS1hZjIyLWVjNDctODVlZS03YmI2ZTM0OTQxNzEiIHN0RXZ0OndoZW49IjIwMjAtMDUtMThUMTU6Mjg6NTErMDI6MDAiIHN0RXZ0OnNvZnR3YXJlQWdlbnQ9IkFkb2JlIFBob3Rvc2hvcCBDQyAyMDE4IChXaW5kb3dzKSIvPiA8cmRmOmxpIHN0RXZ0OmFjdGlvbj0ic2F2ZWQiIHN0RXZ0Omluc3RhbmNlSUQ9InhtcC5paWQ6ODBjMTExMTMtNTYxYy1hNDRiLWEwMTAtOTlmZmNlMjU5ZGEzIiBzdEV2dDp3aGVuPSIyMDIwLTA1LTE4VDE2OjE5OjE3KzAyOjAwIiBzdEV2dDpzb2Z0d2FyZUFnZW50PSJBZG9iZSBQaG90b3Nob3AgQ0MgMjAxOCAoV2luZG93cykiIHN0RXZ0OmNoYW5nZWQ9Ii8iLz4gPC9yZGY6U2VxPiA8L3htcE1NOkhpc3Rvcnk+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+T4/djAAAYT5JREFUeJzt3Xe8HFX9//HXvWmkkEIKIaGEUEMJHWkqXSkKiBSliPqlV6WLioBYaIIIIkWK0kQUBQEBAamKCAiC9F5CC0kgvdzfH58zv3N27szubJ0t7+fjsY89Mzt399x7d2dnPvM559PV09ODZLIFcCKwNjAK6Mr4c8cDZ9SpTyKdZBzwQ2AbYAzQP+PPPQxsUq9OiXSgocApwBeBpYCBGX/uPWDJenVKRNrGjsAxwOrYMXdW+wOX1qVHIiIibaJv3h1oEV8GrgP65N0RkQ41HvinuxeR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3Z/MsyMiUl9dPT09efdBRERERERERERE6kjDgUVERERERERERNqcgoAiIiIiIiIiIiJtTkFAERERERERERGRNqcgoIiIiIiIiIiISJurpjrwccA41z4RmF19d6SNTQCOcu3HgStz64mISO19C1jOtU8GpufYFxFJtgWwk2v/Abgvx76ISLIlgO+79kvA+Tn2pZNsDOzh2rcCd+TYFxGpo2qCgPsAa7j2KSgIKMWNA4507RtQEFBE2suewIaufSYKAoo0o/XxxyKvoCCgSDMahv+cPoCCgI2yJv7vPhUFAUXaloYDi4iIiIiIiIiItDkFAUVERERERERERNqcgoAiIiIiIiIiIiJtTkFAERERERERERGRNqcgoIiIiIiIiIiISJtTEFBERERERERERKTNKQgoIiIiIiIiIiLS5hQEFBERERERERERaXMKAoqIiIiIiIiIiLQ5BQFFRERERERERETaXFdPT0+lP9sXH0ScV5vuSBvrAvq59kJ3ExFpF/pOFGl+3dhnFWABsCjHvohIuv7ufhH2WZX60/5RpENUEwQUERERERERERGRFqDhwCIiIiIiIiIiIm1OQUAREREREREREZE2pyCgiIiIiIiIiACsCkwF5gPL5twXEamxvqU3SbUyMNC1/4sKPUhxg4CVXHsa8Fp+XZEW0QcYix18jAJGuvvRrj0UWBx7bw0GhmMFaMD2TYsFzzULmOvaC4EZwMdu/UzgI+AD4MPg/l3gTXevyVOllBWx9yHA02gic5FmNBoY59pvYft7aV79sWDECsBywPLAkvjjgMWx7/3hbvuZWGGmOdj3+IfAe9gx52vAy8AzwPRG/QJSkf7AJNf+BHgpx750kiWAZVx7JWCEa08AXs+jQyJSH9UEAW8E1nDtJbCTaJE0awMPuvYNwO75dUWayBgseLKSu18RC/otiwUAq9lHhUaU3iTVPOBtLCD4KvCCu73obtr3CcDVwIauvTQWYBCR5rIfcIZrHwWcl1tPJK4PsCawCbAxdty4CtCvjOfI+l3/KvAU8E/gIeARLIAozWE88IRrPwB8Or+udJQvA79y7d8H6ycA9zW8NyJSN7U6wRYRKWYQsBZ2gL8mdgFhMnYBodn1xw6AJgCbJTw+BTuZeNLdR7f5jemeiIhIS1oG2A74HLAVMKxBrzvB3b7glhdgwcC/ArcDj6MRANLZBgftCXl1QkTqQ0FAEam1bizQt6G7bQCsTvn7m/ewbKo3sUDb+xQO2Z2BDeedhh/u+zHFh2H2xw5sBmKByWFueRiFw41HYUOOlsOGjhULVo51t22CdXOxq9iPAP9y98+V/pVFRETa2ggs42hvLMOrq8i207Aph57BhvO+CryBH+o7DT/VR6QbGx480t3G4oN+K2MXIZeJ/Uxf4DPudrp7rd9jI1f+Wc4vJ9ImBgXt0bn1QkTqQkFAEalWH2zYzmeBzbGD+uEZfm4BNpz2efzQ2heBV7DA35ya99SG9s6j/CG8g7AhystjQ5ej4csruXXxIksDgE+5W2QKNpzi7+72DMo0EBGR9jcAy/jbB9jBLcctBB7FMvIewi6gVTJ/9CKsoMFUbOqOJCOAdfFDjzfD5heMLAcc7W7PAZcBV2FzBIt0gjATcGDqViLSkhQEFJFKjAc+jw3h2ZrS8/B8hB3QP4Fd1X8KC4LNq18Xa2oW8Ky73RZ7bDCW6TgZy4BcGzu5GBLbbiw2F2Y0H+YUbOjRX4E70QT1IiLSXsYDhwAHYBn2cR8Bf8S+V/9G4+bY/ci93t/ccj9gU+yYZid8UQqweQnPwDIEfw+ci2X3i7SzMFA/KHUrEWlJCgKKSBZdwHrAl4AdsWBXmoXAY/gr+Y+QfjW+HczEfsfwpKAPsBp+SPSmbjkc9jQW+Jq7LcIyIP4E/AELNoqIiLSiScB3gD3oXdhjLnAr8Bt3Hx/Om4f5wL3udiI2ZHg3bMjyRLdNP+Ar7vYQ8CPgLw3up0ijhJ9bZQKKtBkFAUUkTRcWvNoV2AUbHpMkGsJzLzbc9X5sbr5OthBfIOQyt240Nt9QNGw6DKR24wOGp2NZkn/Esg6eaESHRUREqjQJ+D6W8R6fJuMh4Epsnr1GZfxV6r/u9gNgC+Cb2DyG/d3jmwC3YMc+p7i2SDsJYwQKAoq0GQUBRSRuEnb1fm9ghZRt3sPmtbvF3aY2pmst7X3gRncDGIMFBL+AzZEUFh9Zzd1OAv4H/A7LmnipUZ0VERHJaAQWMDuEwnOLeViG+8+Ahxvfrar1AHe727eB/YAjsIJhAOsDN2MXQY8Enmx4D0XqI/wcaziwSJuJX6UTkc40AjuwfRzLQjuZ3gHAZ4AfYsOCl8Su9F+FAoCVeg/LiNgXHxA8F3g9tt0k7P/xAhZ43Q9dlRURkfx1Y99hz2LHEFHgYB524WoSdqzQigHAuHeBn2LDgw8E3g4e2xw7froKVVKV9qBMQJE2piCgSGfbGLgceAs4DytqEXoa+C52IL868D1svj+prYXYUOpvAROADYCfUFgZsQsbTnw5dvLxc+x/IiIi0mifBv6NDfEd49YtAi7FLiLuC7ycT9fqai5wMbAiVj14mlvfjVU//h9wMDY3sEirUiagSBtTEFCk8wzAsskex+bo2Y/Cq3xvA+dgFW7XwOaoU6GKxunB5hk6EVgeyxC8BH+iATAcOBybs+jv2JyN2p+LiEi9DcMy3v5O4YXDh7B5bfcH3mx8txpuNnastDL2Hb3IrR8JXIh9jxcroibSzJQJKNLGunp6eir92d2wLzqAX2Op/yJplsQCFWDzmt2ZY1861SjgIOBQrDJtaCE2t9+vgDvcsjSXAVh15gOwwGBX7PGXsGzOy4FPGts1wfZvS7r2VcCsHPsiIsnWwjLgAR7EijdJdhsBv6VwupB3sPkAL8UHwjrROliG/mbBurnYdB5n0tl/m3ItDuzl2u9g80pK/U3Cji/BCvws5drv0vu8QURaWDVBQBFpDeOA47DgUfxq3pvYFezLsCHB0hpWwf6f+1FYUASs6uJ57jatob0SEZF21A+rgnscfpjrAuBsbK5gXXgyXcBXgbMoDJrciX1fv53wMyLN6F38MP8ZWAawiLQJBQFF2teS2BxzR9A7+Pc4VoTiWmB+Y7slNTQAq+R8PFZNOPQxNiTpDFS8RUREKrM8lv23SbDuFWz+uwdz6VHzG4VlRu4UrJuGVU++No8OiZRpOjDUtRdgFwJEpE0oCCjSfkZjxTwOxIJEkUXYkIpzgAdy6JfUTzewAzZJ+Wdjj32MZWucjbI1REQku72waUIGB+suxS4w6vuktAOwY67w7/dL4Eh0AVaa21ygf7DcH71nRdqGgoAi7WMQViziRArT9hcBt2Lz0qiyb/vbFDgB2DG2/gNsiNK52MGdiIhIki5s6O+P8fPPTsOq3l6XU59aVVIm5QPYHL/v59IjkeK66T03+FDsorKItAEFAUVaXzfwdeA0/CS+YMG/67H5ep7JoV+Sr02wwO+2sfXPA8cANze8RyIi0uwGAldiBQAjdwN7Y0UapHz9gFOxqTuioOoL2MW65/PqlEiKQcDM2LpRwIc59EVE6qC7ip+9Dvi3uy1em+5IG5uMf7/8JOe+tJN1gPux4TlhAPAuYH1sgmoFADvTQ8DnsMzAcN6mlYE/YxOVr5pDv9rVFfh93Oh8uyIiKfbBf073zLkvzWgp4O8UBgAvBj6PAoDVmI+N0tgdXzl+JeAfwFZ5daqJjcN/Ti/NuS+dZBfsb/5wwmP9E9aJSIuqJgi4OrCuu/WtTXekjQ3Bv18m5tyXdjACK/rwLwqHmDwGbONuj+fQL2k+DwGfxgqIvBis3xp7j5xK78IxUr5J+H2cDpZFmtNY/Od0yZz70mwmY0GpDdzyQmxqiQPRXGC18nvswtwbbnkEcDtwWG49ak4D8J/TVXLuSycZjf3NJyc8puMakTZSTRBQRPKxI/AUNjdPH7duKnAUsCGWBSgS6gF+h1UQPgo/r8tiwPeA/wJb5tIzERHJ26bYqIJl3fLHwM7AT/PqUBt7AtgIy7gCS6Q4H42SkeamIKBIG1EQUKR1LIkFcm4Gxrt1i4CLsGEl59F7Il+R0HzsfbIa9l6KTMSGB/8Cy9oVEZHOsDnwV2zif4CXsSDVLXl1qAO8DXwWuClYdzxwSi69ESmtX94dEJHaURBQpDXsDjxN4Tw9/wE+hWUETs2jU9Ky3sSGB28LvOLWdQOHYlmBm+fTLRERaaCtgL8Ag93yk8DGaC7hRpgJ7Ar8Klj3feBH+XRHpChlAoq0EQUBRZrbQCxz63pgpFs3HxuisyHwaE79kvZwJ7AG9n6KskiXwypBnocO+kRE2tXnsJEFg9zyE1hQ8L28OtSBFmEXcn8erDsROCuf7oik0vGgSBtREFCkeW2AXZU/Ilj3T6wi8AnAvDw6JW1nFvZ+2hx43q3rwt539wMr5tMtERGpk+2woahRUajHsGJRH+TVoQ7Wg83Ve26w7mis+FtXDv0RSaIgoEgbURBQpDkdCTyID8AsBE4DNsOGBYvU2gNYVbhLgnUbYieHuyX+hIiItJptsADgYm75n1gG4Id5dUjoAb4FnB2sOxg4I5/uiAD2vpzr2goCirQRBQFFmssQ4DrsinA0Ce/rWOXW7wML8umWdIiZwAHAl/AnhItjRUR+hSaGFhFpZWsCN+BP6B/E5oadlleHpMAx2AXfcPmofLoiwnz8qCMFAUXaiIKAIs1jVeARrGBD5EZgLeC+XHokneqPwNpYdmDkAOAuYGweHRIRkaosh1UBHuaW/wV8HpiRW48kyfcpzAg8G/hyTn2RzqYgoEibUhBQpDlsAzwMTHLLC7B52nZDV+glH28CW2BFQ3rcus8A/wbWz6tTIiJStiWA24Gl3PLLwBeAT3LrkRRzLPAb1+4GrsaGbIs0koKAIm2qbxU/+xQw27U1RFFK+Ri76gzwYp4daUIHAL/AD7V8H/gK8LfceiRiomD0w8CVWAbJOOBeYB8sY1DM0/hJ3FW0R6Q5vYM/FpmSZ0caaDHgT9hoA7DiH9sD7+bWIymlB/gmlnm/DRaAuRG7EPdkjv1qlLn4z+kzeXakw7yH/d1HYHOSz0NBQJG21NXT01N6KxGph77Az4DDgnWPAjsBb+fSI5F0qwN/Bia65UXAiWjichGRZtWNzem6q1uejVUBfii3Hkk5hmLTwazllt8CNgbeyK1H0gl2B67H3m+zgJWArwNX5NgnEakhDQcWyccA4FoKA4A3Ap9FAUBpTk9j1YL/7pa7saHCl1FdVrmIiNTH6fgA4EJgTxQAbCUzgB2x6TkAxgN/wFd2FqmHaGRSOBxYheFE2oiCgCKNNxy4k8KJnn+OXXmblUeHRDL6EPgcNjQ48g2s2qROSkREmseXgOOD5cOwbG5pLW8C2+Hnh14fuCC33kgnSAoCZhkOvBdwOBasFpEmpiCgSGONB+4HPu2WFwD/BxyJDa8UaXZzgf2AU4N1OwO3YkOXREQkX6sAv8bPU3o+cFF+3ZEq/Re7ULzQLX8Dm09apB4qCQKOBH6LJTVsVqd+iUiNKAgo0jjLAw8Aa7jl2dgwncty65FI5U4GDsUHr7fAitkskVuPRERkCDZkdJhbfhg4Jr/uSI3cSeHFt18Am+bUF2lvUcBvHhYIDNelWSZoT6p5j0SkpqqZx2lPLOoPcCmWHSKSZix+XpoXgb/m2Jc8LI8FSCa45Y+wAiD359UhkRq4EJvD8lpsOPD62Pt8W6zKdSf5MrCka18BzMyvKyKSYm184OR+2q/SaheWAbiaW34X2A1VLG8Xp2Hv4V2wbK0bgPWwqtftZHFgX9d+G/hjjn3pJKthF3Q3d8vlZAKGQ4BXTd1KRJpCNUHAk/AZTdegIKAUNxG7agl20NJJQcBVscDIOLf8LrAN8FRuPRKpnZuAL7r7QdgJyt1YBcp38+pUDo7FCqeA/S0UBBRpPtvgK5ofRfsFAY/Hgn5gJ/C7YxU+pT30YEOB1wRWBJYCrsYuvC3IsV+1Ngp/zvAACgI2ymb4vztUHgQcVctOiUjtaTiwSH2tAtyLDwC+g11lUwBQ2smdwA74wNcawD3AmNx6JCLSWbYAfhgsHwPcl1NfpH6mYSNrokJy8f+7SK2UEwRcMmgPqU93RKRWFAQUqZ9lsYzH6ItxCpYd9b/ceiRSP/dilYNnuOVJWHBQcwSKiNTXcGwagj5u+Rpsgn5pT09iReUixwJb5tQXaV/z8BmmfYptCAxOaYtIE1IQUKQ+lsaCIsu55bewNPtn8uqQSAM8CHwe+NgtT8aqBi+eW49ERNrfr7ALjwDPAvvn2BdpjGvxFZ+7sSDwiNx6I+1oPr4wSL9iG1IY+NMxn0iTUxBQpPbGAHdhxUAA3sMyAF/KrUcijfMwNkfgbLf8KeBmYGBuPRIRaV/7YnP/gZ2wfw0/VFTa29FY0BesOuslOfZF2k+YCViqjkAYBNRwYJEmpyCgSG0NAv6EzQUINnfLdviDNJFOcC9W/XqOW/4s8DuqK0YlIiKFlgfOD5a/BzySU1+k8WYBe+HnbdsV2Ce/7kibKScTcFDQVhBQpMkpCChSO/2AG4GN3PLH2Bxpj+XWI5H83Al8BX8VeUcKq86JiEjluoHLgaFu+X7grPy6Izl5DDg5WL4QqxwsUq35VJYJOIDSQUMRyZGCgCK10YXNyfN5tzwf2A1dkZfOdhPwdaDHLR9I4cmKiIhU5jtYljXAdCwDbGF+3ZEcnQHc49pDKCwSI1KpSjMBQfMCijQ1BQFFauP7WLADLODxTawysEin+y2Fgb8fYHNYiYhIZdbFjjsiBwOv5dQXyd8iYD9sChqATYHj8uqMtI0wE7DcIKCGBIs0MQUBRaq3C4VBjhOA3+TUF5FmdBrw82D5V/hh8yIikl1f4FL8SfnVWKVY6WyvY8HgyPfx81OLVGIBPhOwVGZpPOinIKBIE1MQUKQ6k4GrsOHAYJXZzsivOyJN69vA7a69GDZ/5rj8uiMi0pKOAdZx7SnA4Tn2RZrLdcAfXHsx4GL88alIuRZSeSbg4MStRKQpdPX09JTeKtkK2BcMwP+wVHSRNAOBia49HXgzx77Uykhszr/o93oQ2AqYm1uPRJrbUOBhYDW3/BjwaazCYatbHn8Q/Bz+wFlEmsdIYKxrvwNMzbEvlVgJ+A92TAXwZeyCikhkLPAMMMItH4Rl37eSfsDKrj0TeDW/rnSUEdjF2bOxwobnYsdn3wH+ghV4SzMFWDJY/jTwQF16KSJVK1Xpp5iXatYL6QSzgafz7kQN9QNuwAcAXwe+hAKAIsXMwD4n/wCGY/Na/Qqb0L7VvZJ3B0SkpA/drRV1Ab/EBwBvQQFA6W0KNh/gJW75DCyA00oX3+fTXucMreIjd4suzIbDgcvNBBxQw36JSI1pOLBIZX4ObOHas4Fdgffy645Iy3gO2ANfxXJv4Oj8uiMi0hIOwEYbgI2oOCjHvkhzuwy4y7WHAhfl2BdpPVGS0EKyBwEHxpYXS9xKRJqCgoAi5TsEf/Ddg1UFfjS/7oi0nDuw4SWRM4AdcuqLiEizWwr4SbB8LPBWTn2R5teDBY1nuuUdgN3y6460mCgIOB9/wbbY6ME+CY8rE1CkiSkIKFKezbA5MiI/BK7PpysiLe0MrKol2HfRb7G5ZkVEpNCF2BQKAPdi1YFFinkFODVY/jl+nkCRYpIyAYsFAZMCfsoEFGli1QQBB2Cpv/H0X5Ek3fj3S6mU8mY1HAtURP3/E/CDvDoj0gb+DyuuA/b5ug7on1tvqqPvRJHm1xf/Oe2Tc1+y2gbY2bXnAgdjmV4ipZwD/Nu1x9I6x6xd+M9pqx4TtKI+FP7Ns84JmPQ/UiagSBOrJgj4KDZx6Cx0ZUlK2wj/frm6xLbN6tfAcq79PFbMQFWxRSo3B6tuGU3Uvz6WXduK7sPv48bn3BcRSfYt/Of0sJz7kkV/4BfB8g+BZ3Pqi7SeBdj0NdGx6iHA5Py6k9kE/Of0b/l2paN8E/ubf9YtL3A3yJ4JGBUVUSagSBPTcGCRbA4GdnHt+VgA8OP8uiPSNt4A9sVnthyD5gcUEQH4NrCya78MnJVjX6Q1PYpdxAYL5FyAZdqJlFJJEHBGwjoRaTIKAoqUtjpwdrB8PH4Io4hU71Zsziuwk5NfYxPhi4h0qrHAicHyUVj2tEi5TsBn3G8G7JFjX6R1hEHAYsOBk4KAygQUaWIKAooUtxhwDX6er9spLAwiIrVxNPCEa48BrkDfUSLSuc4Bhrr2HcDNOfZFWtuHwCnB8lnAkJz6Iq0jnBOwWCZgOCegMgFFWoBOsESKOxc/f8q7wNfRhNwi9TAX2As/n8y2WGBQRKTTbArs6drzgCNy7Iu0hwuB/7j2eOCkHPsiraGSTMDpCetEpMkoCCiS7kvAga69CJsHcEp+3RFpe89gQ94ip2NFhUREOkUfCudtOwd4Lr/uSJtYiBXDiS5kfxtYJb/uSAuopDpwlAk4MGnDmA2x4nArlt81EamGgoAiyZYCLgmWzwTuzKkvIp3kEuAG1+4HXAUMyq87IiINdQCwlmu/hV0MEamFB4DrXLs/KjQjxWUdDhxl/S3Aj+bIkgl4K3a899WKeiciFVMQUCTZhcASrv0I8L0c+yLSaQ4AXnXtlYDT8uuKiEjDLA6cHCwfC3ySU1+kPYXvqR2BLXLsizS3cqsDz3W3cF2a4cBI116/ks6JSOUUBBTpbU9gZ9eei80DOD91axGptWkUzr95FDZHlohIOzsOWNK1H8FnbYnUylvA2cHyOeh8UJKVOxx4Lr6CeanqwCsF7VXL75qIVEM7fZFCo4DzguWTsXnKRKSx7gUudu1u4FJKH1SKiLSqccC3guVjUCEyqY+zgHdce218ERqRULmZgPPIngm4XNAeUX7XRKQaCgKKFDofGOPaj2NXSEUkH8cBb7j2qsB3c+yLiEg9nQYMdu2bgPvz64q0uU+AHwTLp6NqrtJb1kzApOHA/VO2jQwP2sPK7pmIVKVYVL+UK7HiCeBTf0XSvIUPqD2RYz+K2RF/NXQe8DU0DFgkTzOwCt23uuXjgT8C/86tR+muxiZdB/g4z46ISKp/4Y9FHsuzIzGTgH1deyFwUo59kc5wGXAEsDowAascfHaxH2ig6fjP6ct5dqTDPIn93Q/BRl4swPZHYFXLu0jOTg6DgFmChgBDg3Y/7ALIzPK7LCKV6Orp0UgDEewq1NPAeLd8MnBqft0RkcCV+BPkJ7FJpBWgF5F2cSuwnWv/EjsJF6m3LwB/du1pwIrAh7n1RprFdCxItx32fnjErR+AJUnEHQhcBPwXuB7Lar4H2LLIa5xMYTbq0ljCiIg0gIYDi5hz8QHAJ4Gf5NcVEYk5En9wOBk4Ice+iIjU0ub4AOAn6AKkNM7NwN9cezhwYn5dkSYSZfGFcwKCZQMmqTYTEDQkWKShFAQUga2B/Vx7AfANkq90iUg+pgGHB8snoWpyItL6uoAzg+UzgSk59UU604n4IZ6HUViwQTpTFOyLBwHTphELqwNnnRNw8djy8KydE5HqKQgona4/VgwkchbNOd+YSKf7I/A71x4A/DzHvoiI1MLO2PQGYNVam2VONukc/8KGcIJ9t34vx75IcwgzARcG69MyAaOA3zx80LDcTMDhWTsnItVTEFA63VH4jKLXgR/m1xURKeEIbK4agG2AXXPsi4hINbqxebEip6OJ8SUf38cHb/YDVsmvK5KzbixDGQqrA0N6YC8MGkaZgKWqTcczAeNBQRGpo2qCgGcA17nb4Np0R9rYyvj3y5E59yUylsIKfEeiA3CRZvYuNuF05Gc0z/fPafh93Iic+yIiyXbEf04/n3Nf9gDWcu3XgUtz7It0theAq1y7D/lXpx6N/5yeXGJbqZ0tsb95JGsmYBQEnEflcwIOytJBEamNaoKA22EHMHtQety/yCj8+2XTnPsSOQf/JXQncFN+XRGRjM4DnnLtZWieicy3xe/jdDAr0pwm4T+neWY79cGyryKn4TNoRPJwKv49+FVgtRz7MgT/Od06x350mhWB3YLlrHMCRgG/+fg51cudE1DHTSINpOHA0qk2A/Z07XkUFh0Qkea1ADgUP5H5sWjokoi0ln3wU5G8is/CEsnLa8CvXbsP8IP8uiJNotzCIGEQsFQmYDzoN7C8rolINRQElE7UF7gAP+fF2cBz+XVHRMp0P3CDa/dHRUJEpHX0o7D4wvfxJ84ieToNmO3aXwbWzq8r0gTiw4GzZAJGw4FLZQLGg37KBBRpIAUBpRMdBkx27TeBH+XYFxGpzNHAJ669LfDFHPsiIpLVN4CJrv08cG2OfREJvQNc5NpdwCk59kXyFy8MkhYEjNbPxw8pLxUEXMzdR6M6FAQUaaC0D7Mk6wOsBKyJzUW1LLA0MA6bW24gNkl9tOObDczBTlRnAG8Ab7n7V4En3b00zpIUDnEIAwki0jrexKp5/8Qtn4vN7Tk77QdERHI2gMKiC9+jcLidSN5+AhyAnc98EdgQeCTXHkleGpEJOAMYhoYDizSUgoDFjQG2Aj4DrIMF/8q5UpGlQuQ04D/A48A9wL3YDlHq4xTsywbgPvyQQhFpPT8D9sPm1loeOAL4aZ4dEhEp4hvYRWSA/wK/z7EvIkneA84HTnDL38eqakvnic8JmFYdOAr4zSP7nIBRJuBUFAQUaTgNBy7UBWwCnAk8AUwBrgEOAj5F9gDgx8CsjNsOBz4LHAX8CfgQeBDLVls943NINqsC33TthRQWFxCR1jMPy+aNnACMzKkvIiLF9AOOC5ZPBhbl1BeRYs7CzmUAtgfWzbEvkp+FlF8dOMoE7EPxoGH02FR3r+HAIg2kTECzJvAVd5tQZLvZ2JXbx4GXsKG9r2PD0qZhgb+5sZ8ZhA3/GAOMx4YPLwusDKwFTKLwaklfLBC5CXaA+CQ2X8y1WOUuqdxP8e/5y7H/pYi0tluBu4CtsYsq3wW+lWeHREQS7IM/xnwGuCm3nogU9yHwSyxo3QV8BysUIp2l3CDgAgrPg/uTPEVLmPX3obvPGgQc7PoxPeP2IpKgk4OAfYBdsJPFTVK2eRf4G3aC+Q9sAueFKdummeVuH5FcgbY/sAawKXYSuzk2v2BksrudDtyCDX+7t8w+CHwaXzhgNnBqjn0Rkdo6Fvg3lt1+CFb9+8VceyQi4vUBjg+Wf4iyAKW5nYMV0hsEfAlLmHgq1x5Joy10tx4sGFxOJiCkBwEXC9rlZALeCmyHBagPybC9iKToxOHAA7Ghty9g88HFA4D/w+a/WAtYCtgLyxr7H+UHALOYBzyGzb+xEzaUbTPgbCzDMNKNBbHucdvvTXqatRTqwhcPABvm8EZOfRGR2nsCuM61+wOn5dcVEZFe9sRGgIBdoNB8xNLs3gUude0u/ByB0jmi894oG7BUEDCcExDSi4NUmgm4gbvfKsO2IlJENUHA6dgQ2Gm0xtXMPth8cC9g2XTLB4+9j80DuA6wGnYC+ST5zBe3AJsT8BhgOSwz8BJgZrDNOsBvsBPfHRrbvYotwL9fGl2Nd3d8sPd9LAgoIu3lJPwwlD2AjRv8+jNore9EkU40B/85ndOg14wHUE5HFYGlNZxJ4ffqykW2raWF+M/px0W3lFqaS+HfO9pPRcHAtOSTMBNwXsL6uDAIODVhXZJxwCjXXqnIc4tIBtUEATfDqt+OoPnH5e+ABfUuxeblizwDHIDN0XccFlRrJouAv+P7eCKF2YFrYEOE/w6s3/DelecR/PvlGw183X7YsJvID1D1ZZF29CrwC9eOZ/82wjb4fdw7DX5tEcnmfPzn9FcNes1dseM1sHmkr2nQ64pU603gStfuQ+OyAV/Hf063b9Briv2vPxssl5sJGA8CpmUCJg0HLhUEXCZod6EicCJVaffhwKOxghq3YBl+kSewwOAaWJZdo64GV2MqdlI7EdgP+4KMfAabs/AsVGI97hBgRdd+Hvt/i0h7Oh1/QPkZYMcc+yIiAoWBkx9ReJIs0uxOx8/ztjeFI6mk/YSBvngmYFoQMAr2xecEzJIJOM3dDyjRryViy6MStxKRTNo5CLgb8DQ2D0vkDeBALGvuVvIZ7lut+diVmpWw3+U9t74PcDRW8XbLfLrWdIZilUIjJ1L45SQi7eUjCjMAf0J7f8+JSHPbEVjPtd8CrsivKyIVCbNX+2Ejp6R9hUN+o+BfdO6UFgTsG2xXTibgXKx4JpQOAo6ILY8usb2IFNGOJ0eDgd8Cv8PvIOYBJ2OBs4upT4GPRpuH/S6rYtltUUBzIlbN+HRUOOQI/JWih4A/5tgXEWmM84HXXHt1Ci8EiYg0UlgROJxfTaSV/Bh/7rQfsGR+XZE6iwJ6Pfj5jUvNCRgF++aRLRMwKgIyGx80XCxl24gyAUVqqN2CgCtjw2L3Ctb9C7sKeyrtefD1ETZn4NbAy25dF/Ad4HY6dyc5DPhWsHwSrZn5KSLlmUPveUDTrl6LiNTLp7D5s8EqYF5aZFuRZvYc/kL6YsDhOfZF6isK9IXFi6J2WlCvb7BdliBglPU3Bz8lV6lMwOGx5U49vxWpiWqCgBsD27pbM5xg7YQF/KLJlxdggZ+NsSGy7e5uYDJwUbBua+DfNEfRkGH498vkBrzet/BXje4G7m3Aa4pIc7gcmwMULAP8Kw14zQ3x+7hSB7Miko/l8J/TZev8Wt8J2ucDM+v8eiL1dEbQPgRYvI6vNRD/Od2gjq8jhcZjFy/AZwGCDwKWygScT2HwMC0IGGYOZg0CxjMBVRhEpArVBAEvBv7qbvX8IsjiYOAP2BxwAB8A22ETMLfD0N+sZmJ/i73xB5vLYgGw7XLqU2R1/PvluyW2rdZw4Mhg+ZQ6v56INJeF2PClyPep/8Wq8/H7OF2hFmlOu+M/p7vU8XVWwRcmmgVcUMfXEmmEfwH3uPYI4P/q+Fpj8Z/Tc+r4OlJoB+Cnrt0VrC9VHTgM6vVk2D4K+M3FDwfuT/G4RDwIOLjItiJSQjsMB/4BcCH+d3kIWBubF69TXQ1sCrzklgcDfwK+mluPGutofNr4HcB9+XVFRHLyG2wIE1iF8H1y7IuIdJbj8Mell2EXp0VaXZgN+G3SCz9I60vKBCxVGCTaLhoSnPb+iNbPxWcCQvFswGGx5aGJW4lIJq0cBOzCrqyeHKy7AauM+1YuPWou/wE2Av7plvthJ8XtPo/HSKwgSOS0vDoiIrlaSOHcgCejExYRqb9x2IgMsP3QeTn2RaSWbgced+2lUeGtdhbOox6NqsuSCQilg4ZhdeBwvv5iQcAhsWVlAopUoZWDgGdhc1JELscy3dqx+EelPgC2Am5zy93YweghqT/R+o7BXx26FXggx76ISL6uBf7n2ssB++bYFxHpDEfhT4qvx4/KEGkHZwftE2jtc0lJV86cgFGwL8oAjIKBWeYEDM/bi1UIjioKR/3KOhXZuhm3E+korbrj/gmWhh45FfgGhZORipkJ7IydDINlUJ5Pe54MjwIODZY1F6BIZ1uIfT9EvoeyAUWkfoYCBwTLZ+XVEZE6uR54zbUnkf+c41If5QwHjoJ9URCwVDXhcE7AcDhwseOzKBNwamy5mD4UBq1FxGnFIOB3geOD5dMoHBIsvc3D5sO6wS13A78Gvpxbj+rjWPyVoZuBR3Lsi4g0h98BT7n2ssDXc+yLiLS3A/FzV92BHzop0i4WUFis49i8OiJ1VUkQMD4nYJbhw1kzAaPhv1PcfZYg4JrA5sCqGbYV6SitFgTcl8I53n6OVX2U0hYCewG3uOU+WAGRz+bWo9oagVVGBpvH4gf5dUVEmsgiCrOCj6P+lYJFpPP0BQ4Lls9I21CkxV0GfOjanwXWybEvUh/lzAkYHw4c3WfJBMw6J2AUBHzX3WcZDrypu2+3pBeRqrVSEHAz4OJg+WJs3hXJbj6wG3C3W+6PZclMyKtDNXQ4/gvhVuCxHPsiIs3lD8B/XXsisEeOfRGR9vQlLNsYbH9zd5FtRVrZTArPyY7MqyNSN2EmYLHMvm58PCEKFkZzApaqDjyPwuHAWeYEjDIBswQBN3H3OuYTiWmVIOAE4Eb8FYK/YnO/9aT9gKSag80RGJ0QjwH+Qu/S661kMIVX33+SV0dEpCn1AGcGyyfROt9/ItIawkDIOegYVdrbL/DBoT2BsTn2RWovayZguC7rnIBhEHAePuBYTiZgluHAUSbgGlg1axFxWuEkaBA2v9sYt/w0sDsqAlKNj4GdsOrBAKsBV2FFQ1rRgcBo1/47qggsIr1dA7zq2pOAHfPrioi0mfXwWSfv44uxibSrt7EEDbDgzYE59kVqL2t14H4J25UzHBh85mBaEHCx4LWzBgEXB5YLlhUEFAm0QhDwAiyCDzb/xM7AjNx60z5eBnbB74C/CBydX3cq1o/Cq+8/zqsjItLUFlBYqfOkvDoiIm3nW0H7IgqHuIm0q7Dy6iEUH84prSVrYZCkTMCshUGic9Bof5n2/hkctN8Ptk0KSkaWjC2PK7KtSMepZnL0I/FDSD+pQV+SfA3Yz7UXYvPZvVin1+pED2AHrhe65R8BDwIP1+G1nsXmywF4s4bPux9+Dp4nsGp8IiJJLsOCf0sBGwJbAPfU6LmPA5Zw7Q+LbSgiubkJfxz5ZI2ecyns+BTs5PdXNXpekWb3KPBP4FPYiK09gCtr8Lzv4s8ZPii2odTUHcAlwP4UHsdkHQ5cbiZglAE4N7Y+LgwChu+HgaTHIMbElhUEFAlUEwSs94THKwPnB8unULuTNfF+ic2ZsBe2s74eq/JV65PYqcAfa/ycfYBjg+XT0Rw8IpJuDlZVPsoYPpHafa/8vUbPIyL184K71dIh+MyW64G3avz8Is3sPGy6DbAEkVoEAWdR+3MGKe1V/Jzxs4L1tR4OHM8ELFVIREFAkRpr1uHAUdXaqPLP37AAj9THIfgr48tggcFWsBuwkms/h1UAFREp5kJgmmtvA2yUX1dEpMUNwLJmIuenbSjSpm7Aj/BZB/hMjn2R6kWBvoXBunILg5SbCVgqCDgoaH+Ysj4uPhxYhWtEAs0aBPwusJZrv4tlqS1K31yqNANL4Y+uyOxG85dT7wJOCJZ/it4jIlLaDGyu2chxeXVERFre3viTzQeBR3Lsi0geFlCYPHBk2obSEqKgXliAs1gmYNJw4FLVgeOFQaL7tCBgNFdgD/BRsH5gyvbQOxNwicStRDpUMwYB16EwuHMIvhKQ1M9jwKnB8gU091WTbfCB4jeAq3Psi4i0lvPwQ112wmcUi4iU4/CgfV5uvRDJ10UUfqdOzLEvUp1imYDlDgcuVRgknglYrDowWLBwZrC+WCZgPOg3vMi2Ih2n2YKAA4Cr8DuU36Ahno10Bv4q9kiae1hwWMn4Z/gvEBGRUt4HrnDtbgore4qIZLEl/mLk62gOM+lcU/EX4/tgCRzSmqLAXRgErLQ6cNZMwFLDgaOMvznu1hNbn2RobHl4kW1Da5XeRKT1VRMEvBObA+JNfJXgap0ArOHab6OU8kZbgFVknu2WdwZ2rdFzr4d/v1QbXFwDywQEG9r36yqfT0Q6z9n4g9yvAaOqfL6b8fu4Zs6iFulkB+E/p9+o8rnCY9QLKBw+J9JpzsMHZ/andxCmHMvgP6c3Vtkvye6r+IKLYSGNYnMCVlMYJOucgFEm4Gxs6qcoeJglCBgNHx5eZNvQ5VR/PCjS9KoJAo4FxrtbLTIKV6RwGPABFI77l8Z4FvhBsHwuvkBLNQbg3y8jq3yuo7E5AQEuBqZX+Xwi0nlexgJ3YENKDqry+cbg93FJQ2ZEJH+L4z+n1RzbLA/s4NqzgMuq7JdIq3saK+QIFoDZp4rn6ov/nMbndpP6GYIPnoXn9gsS1kWSMgFLzQlYbhAwCvbNjt0XCwJGCUqvx5aLmYhNS7Z/qQ1FWl0zDQc+Dx/p/z3wlxz70unOAf7j2ktTGBTM25LAnq69AFXiE5HKnR20D8d/B4mIFHMkPth/JYUVK0U6VTgv5pE013mmlKcnaGetDhwF/6KgXloQMFofBQ2zFgaZ4+6j+SeLzQkYBTNfD5ZLvR+/4O4PIb3vIm2hWXbOuwPbu/bHaH6mvC0ADsV/ARwBrJ1bbwqFJ+o34HfuIiLlegD4h2uPwYbCiIgUsziwn2v3oIuRIpG/AM+59krA53Psi1RnUdAuNidgFCzrIftw4Gh91sIg4ZyAUFkQsJvSQ9Sj7O6lgU1KbCvS0pohCDgIOCtYPhmbA0Ly9SA2LwLYTr8Zqt4NAg4Mls/NqR8i0j5+FrS/jZ9qQEQkyTfwQ8v+Cvwvx76INJMebH7MiOZ2b11JmYBJU51EgcFwTtSsw4GjYGHWOQGjIGCW4cBRwO+1YF2pIcFrB+0VSmwr0tKaIQh4NDYBLMBT6IpqMzkeq/gF8BngSzn2BezKezRZ69/xlYxFRCp1IzY/IMDqwOdy7IuINLduCiufNsMFUpFmcjl+ru5tgck59kUqV2514DAIOD/2WFw1hUGgvEzAt4J1Q4psP5DCgiDLF9lWpOXlHQRcEl+FCOAYVF2tmXwAnBYsn0l6qna9dWNDgSPn5NQPEWkvCym8+HR0Xh0Rkaa3I7Cyaz8P3JFjX0Sa0ScUFso5NK+OSFWyZgLG5/eD0pmA8eHApeYEjA8HjoKBafM498EH/LIGAZemcCTIxCLbirS8vIOAP8JXZ7sZHUw1owuwA12wHeLhRbatp+2BVV37BeCWnPohIu3nMnzmwtbAGjn2RUSa12FB+zwK580SEfMLfOBob2BEjn2RymSdE7CaTMD4cOC0RJN4JuCc2Pq4QfiA3kfB9sWCgMvElpcrsq1Iy8szCLgG8DXXXgAcl2NfJN18Cv83J5HPl7kOvEWkXj4GLg2WlbkgInErYxcJwC4aXJVjX0Sa2Sv4i/WD8Od70jrKzQQsZ07AtMIgWTMBo8zBtCDg4KA9E8tOhdKZgKExRbYVaXl5BgFPxu9MfgU8m2NfpLg/Afe69nAaX715JWAb1/4Y+E2DX19E2t/5+APdfVHmgogUOgifXXIl/sRSRHq7MGgfioputZpyMwHD4cDFMgG7E36m0sIgxTIBI1mDgKNKLKdRxqC0pLyCgGvgi0zMAX6cUz8ku5OC9pHAEg187cPw79XLgRkNfG0R6QyvAbe5tjIXRCQ0kMJ9wsV5dUSkRdyJn05oRWDLHPsi5QuDgLWsDhyuyzonYHw4cDmZgLPIFgSMzmujoiPDSc9kDF0IfDrDdiJNpZog4D3YPH434z/EWZ0WvPaFFE7aKc3pIfycjUOBb5f581Px75dHy/i5IfgD7x7gl2W+rohIVhcE7UMp7zvyPvw+bnaJbUUkHy/iP6cvl9g29BX8SeI9wNM17pdIu+mhMFh+cBk/OxP/OX2glp2Sol4F3nTtj4L1lRYGScoEDAN9UfxgfsJjoXgmYHSfNodgJcOBo/17FLjuAkYW2R5gPPA5bJqqvOssiJQlbcLOLI6o8OfWAXZy7ZnAGVX0QRrrJGxYbheWDXge8H7Gn30W+GIFr7kvMMy170TDxkWkfv4KPAesgmUubOPWZXFs6U1EJGd/dLdyhQEMXYwUyebXwKlYdv1O2Lxrbxb9CfMelZ0zSHXuAO7HLno8F6yvtDBIUiZd/4Ttosy+tKBe9DNR0DAKAg5M2BZ8EHAOFsCMgoCLJ28OFAYB13btUcCUIj+zLxYYXQfYAQtai7SEPKLWp+DnhbgQeDeHPkhlHgVud+0hlJ8NWImDgvYFqVuJiFSvB7goWFaBEBHZAFjftd8BbsqvKyIt5SPgBtfuC/xfjn2RbKJsv4XBulpWB04aDlzrOQGjIOBMd/9JbH2SKAj4Iv53LzUv4DZBe/3UrUSaUKODgOsBO7r2TOCsBr++VO97+IpRh1Pf6klbAmu69mvAX+r4WiIiAFfgDxx3ACbm1xURaQJhFuAlFA57E5HiwszZA8g2z5rkJykIWMvqwEmZgKWCgFGG4NzYfanCIPEgYJZMwA+Aaa49vMj2YIUrI6uV2FakqTQ6CHgaPgvw51i6t7SWf+ODcYOBY+r4WmEWzoUUfiGJiNTDNOC3rt1NYTayiHSW4cAerr0ACwKKSHb/xM4dAJZCw3ybXaWZgFmrA1eSCRgPAkYZgVkzAWfG1ieJgoBTgemuPbzI9oOwOQEjCgJKS2lkEHBjYDvX/gT4WQNfW2rrZHw24CHAknV4jWXwBwpzsarAIiKNcD5+H/d/+KvKItJZvoH//N9MtvnMRKRQOM1GOQVCpPHKzQQstzpwLTIBswYBZ8Xuix3LDXX308iWCbgCPrEJYPki24o0nWqCgIdik72eSvqHMHRS0C6noIQ0n8eAP7v2YLLNDbgM/v2ye4btD8J/sVyD3i8i0jhPY9V+AUZgk2SXcgB+H1dsyImI5GdT/Od0oxLbdgEHBssqCCJSmauxDCuArSidNTUc/zn9Zv26JTHrAZNcO8xya9ScgGlDxavNBCwnCDgDHwQclrwpYEVuQgOx48VS1gheSyQ31QQBD8Lmh/se6dV5IqsB27v2x8DZVbyuNIdTg/YBlD7pXQb/fvlyiW0HUPilf2HZvRMRqU5YiOiwDNt/E7+P0wGeSHPaBP85/VSJbbcGVnbtF4G/1bFfIu1sNvCbYHn/EtuPwH9O96tTn6S39YAVXXuZYH3emYCVFgaJZwKmDQceiO/rDPxw4GJBwLHufmawblyR7SM/QYUupQk0ajjwt/Eps5di1aKktT2GPyAeTm2v1O2BH2L8MFaVWESkkf6IH/q3NpZBJCKdIxy2eCGwKK+OiLSBC/DTbOxH8fnZJH/h/q6WmYBhoC8K/s1PeCxUbmGQKDkpChZGgbq0TMDwwm2YCTg8ZXvw56nP4Ps/PmXbyLpYUtTewM4lthWpq0YEAccAe7n2QuAXDXhNaYwwo/Moknf2lQgLguhqiYjkYQFwcbB8aNqGItJ2xgE7uvZs4Moc+yLSDl4A7nHt4cCe+XVFMugJ2rXMBIzWzQ9eo9ZzAsaDgKWGA8eDgFkKg0SZgO8AU1y7VCbgvvikqB2LbShSb40IAh6O/5D+Hni5Aa8pjXE78KRrLwfsWoPnXBfY0LXfx94zIiJ5+BX+oPPLlL7KKyLt4SD8yeq1+PnMRKRy4byaWabZkPzUOxMwrCZcryBg1uHAaZmAxYYDR5mAU7BAIFjiUzFrB+1S01GI1FW9g4CDsAOpiCoCt5cerMhL5JgaPOeRQfsi/A5fRKTR3gNudO1+WKVgEWlvfbGqwBEVBBGpjZuAt1x7bfxFf2k+9c4EnBesK1YYpIv0IOCA3psDPuOv3EzA+e5npsfWJ4kCfu8CH7j2qCLbdwFrBcurFemPSN3VOwi4H/4DcR/wzzq/njTe1fgrIOsDn6niuUbhKwfHh+KJiOQhnJLgQNKr14lIe9gFn/X7OJqXWKRWFmBzw0cOTttQclevTMB+sW3ABwH70DvQ2B8/hDY+J2DfhO2h8jkBZ7j7j919saKXUSXgD8gWBFyKwuHF3cDyRbaPjKSwSItITdQzCNiNzRMXUUXg9jSXwpPko6t4rgPwqd3hpPwiInl5CB8EWAr4Uo59EZH6CwMT5+fWC5H29Ct8AGgPLMghzacWmYBZg4BhOz4kOMz2izIAw1FiSUOC0+YEHIgPKIYqCQIu4e6n4oOAxd7LywXt6PedWGT7yFnAHcDoDNuKZFbPIOBOwEqu/TxwSx1fS/L1S+AT1/4CMKmC5+gD7B8sqyCIiDSLcDigCoSItK9JwOauPQ24PreeiLSnd4CbXXsgNmpMmk8tMgG7En6m2JyA4eORMAgYHw4cfzySFgTsIjkbMJorMDqXLScI+BHwoWsXywSMsvk+wtdHKBUE3Aj4GrAqcFqJbUXKUs8g4LeC9s8o3JlIe5mKr5zXBRxRwXPsCExw7f8Cf6++WyIiNXEt/iDv08DkHPsiIvVzID5T5HL8yaOI1E54YS38zEnzCM/bK80EDB+PFBsODNmCgJVmAkJyEHCIu48HAYeQHCvphw8QhpmAxYKAy7r71/FBwOVSto3siP9sfK7EtiJlqVcQcC3sRAnsg3FlkW2lPYSB3r0pXlEpyYFB+8Ka9EhEpDZmYwGByIFpG4pIyxoI7OvaPdiwRRGpvb9ho8TARo1tmWNfJFk4HLjSTEDoPY9yscIgSdtXkgkYBfri1YHBBwhDaUHALpIrCi8RtMMg4BIJ20aWdvevA2+49rgi2wNsFrQnACuU2F4ks3oFAcOy75fiI/HSvl4CbnPtIVj6clbLAtu69idYsRERkWZyEf5Cx14kHxiKSOvaHT/Z+73Ac/l1RaSt9QCXBMsH5NURSVXvTMBwm2oyAcsZDgzFMwGjAiIzgseSKgSPCNpT3Q2s8EdaVutYd/8W8LZrL5WyLe551outyzIKZTng2CL9EAGqCwJ+GpukcjQ2b0pkOPAV116EKrx2knAev8Mo3AE9gn+/fDP2cwfiv1iuoXDnKyLSDF4C7nHtYcCesce3xe/j3kFEmtEv8J/TeKZfGIjQsatIfV2Oz+jaGVgyeOw1/Od0x8Z2q6NdCbzi2n8J1kcBu256xw5qkQmYtTBIUiZgLYYDp80JCMnzAkYZfz3YHH8fueW+KduDf3+/iz9GLJYJOA4fnIx+39WLbA/2t7oJOAP4YYltpcNVEwSchi+LHaYMfx3/YboFvzOR9nc7hen9WwWPLcC/X8Kda18KJwXWgbeINKswaBAfEjwdv4/THLgizWk2/nMankhOAjZx7Q+APza4XyKd5kP856w/hecCi/Cf0+mN7VZHCzPswlF8C4N2PBuw3nMC9k/YrtxMwDn4WEU5w4EhOagXTXk1E/tdPwoeG9F7c6AwCJglE3AVd98D3OXapYKA2wBru/YB9A7Aivx/tR4O3AUcFCyrwmtn6aEwiJeliuZO+CshjwL/rnWnRERq5Cb8FdwNgHXz64qI1NDBQftyCk8yRaQ+wnOGA6hvwUrJJgryhYG/YkG9cjMBK60OvDB4jQVB/7IEARfhL/pkKQwyN+hTUhBwuLuPAtQfJTwWFw0HDjMBB6c8P1gyDcCb2PkxlJ4T8ItBexS+0n0pScO8pc3Veme7LbCya7+Ij1xL57gMn3b9BXzF3zRhNo0m4RaRZjafwkJXmsdIpPUNxAqagV3MvDTHvoh0knuBZ1x7IioQ0gySgoCNzARMCxrOi62Pgnrx4cD9gtcMsxmjdpZMwLA9hN6iTMDpwbbR75yUCdgfHxx8190iYxK2B185+EV8NeHlU7aNbBZb3rTE9mBD8V9Fw+47Tq2DgGHm1wVoSFQnmgZc59p9gP2LbDsRP2T4Y+D6+nVLRKQmLqawQEjSpNEi0jr2wJ+43YOf1kRE6u+yoB2fZkMaLymoV24QMMucgOE2C4PXSMsEjGdnz409HgmDgmEQMEpQyVIYBPyQ4GJBwGnuvidoJwUBR+PnyX8PeB8/PDktCLiMu38dP7XaKNIzBwfjE7GizMRPpWwbGY4V41wa+DXJAVJpU9UEAcdjEenl3fMsC2zvHpsFXFVd16SFnR+098d2yAPw75fR7rEw9f+3FM7BICLSjF4B7nTtIfgCIUvh93EaWiHSnIbhP6dRAF8FQUTycyWFBULGYYGl6HNabN40qa3F8UG4MJiWZThw1uHDSUHAcDltTsB4EDAtE3BgwjZQPBMwXhgkbGfJBAQfeEsKAo4M2h9gv0v0s2lBwGXd/RsU1leYkLL9ZPyx5xXufq2UbSNfxQdFR2MXxLJYpvQm0uyqCQLejqWnvox9GA7Bv/muxpfLls7zBPAP1x4N7IqVOY/eLxdgXwL7Bj9zSQP7JyJSjXDqgkPc/U34fdzY+A+ISFM4AP85/TqwJrCxe+wD7HMsIo3zIXCja/cFvoYFGaLP6e9y6lcn+go+iLVDsL4emYDx4b3RclomYHz7WmYC1iIIOM3dJwUBR7n7BcF277n7UpmAb2BzCEa/79Ip26/q7j8E/uTaS+ErGSfZIba8feJWhY7Ghg//IcO20sRqNRy4P3YwFflljZ5XWldYFOaQhMd3wV/d+yfweN17JCJSGzcDb7n2WsD6OfZFRCoXDj/8NSoIIpKHMAN3f1QgpBmEU3qVWxikkkzAebHHI+UOB07LBMwSBAyHA0dBwGKFQaYF66KA4DB6izIBp+KHARcLAnbhi2a+jv0vomPOtCBgNBT4Wfw8m2CV75P0wc8h+IG73xI/bDnJeOCH2OdzFyzJJwuNjmlCtdrJ7oB/Ez+MAjoCN2BzHgBsgs3/F1JBEBFpVQuwCqIRzWMk0nr6YcOhQAVBRPJ0H/C0ay9PtoIGUl89QbuawiBZ5gSE9EzAcguDpGUCZhkOPCtYV2km4PCE7aNMwA+DddE5clIQcAl8cDOqJBwFAccnbA8+CPiCe+73Y+vjVsNPifFjdz8SX5U4ydcp/PseVGTbyLZY32+i9/9WclSrIOB+QVsBHQG7QhNW0dwpaA8BtnDt6SjVX0Raz8X4A+OvoiudIq1mXfzQrbuxkycRyUcYhP9q6lbSKNVkAvYEy/Hto0BQ1jkB650J2BWsy1oYJAqezQjWTXP3wxO2jzIBPwjWRe1R9BbOgxkFAd9092mZgCu4+6iw1XPuPi2ot467n43FbqK/54Yp24Nl/4EPEG9Bcv8jw7AaEUticYDvFdk2tD7pwUupkVoFAaN00ulYBpgIwEX4HcXng/UT8enGv6Fwpysi0greAP7q2oMonPhZRJrfJkFbBUFE8nUFPlCzdY79EFNNJmC4XK9MwCholZYJOD/Wp7RMwIH4eEjScOBaBgHDTMAsQcD5wXZREDAtE3CCu3/J3b/o7ldM2X5td/8U9nv/xy2vm7L9SHzgMCoA2ofin9X9sABg5HDSqxtHjgf+hWUG71xi20ixIcySolZBwOiPfxWFqbTS2V4C7nHtcCc6IWhf1rDeiIjUVpj5nja5s4g0p+Xc/Qf4idRFJB/T8MUGlFmfv2oyAcEH+eo1J+Cc2OORgbHHI2mZgIODdtYgYBTICoOA0dDg4QnbR8U5yg0CTsH/H95290lBwBH4IcqvuvtSQcDV3P3jsfvJKdt/Gh/vOQM/7+BWKduDFfkB/zsMA75UZPtVsDkHwd43F5P89wytj2Xxv0NhwpGUUOuJV3UlVeKShodHO+yHsErCIiKt6C/Aa66dNNm0iDS/y1BBEJFmoCml8hVmVIVBwHpkAsYz+9KGA0fLWYcDR5mAs2Pr0zIBw2O3WmQCJhUGiYKAU4N1xYKAY939lGBdFARcit4mBO1X3X2UERifkz8SBQGjuTifdPdrpmy/kbt/AZvj7063nDZ/53L4bMPvYaNnAHZP2R7gBAqDxqOBg4tsPxIr1LcC9jf7PRZILGYyFvCchmUdZrE4sAZtdnGilkHAB4H/1vD5pD3cBLyb8pi+7EWklS2ksECIiLSWHjQiQaRZPIDOJfMUxgXC4cCNzARMmxMw63DgtEzAtCBgqUzApOGr5Q4HjuaezRoEjIbQhufPURBwJL0Dn8u6+1n4giCvuPth+CBkZHH83IJRRl8UBByDBd/iNnD3D7j7+939qgnPD5aV14UFk3+NrxOwJckXzYfiA4TP4v8+B5MefDsJHzAF+1+el7Itbts7sODkMOAnwLeLbA/wZSyA+ZS7rVpi+25gb+Bn7j5L4LAL+5/HA+Z1VcsgoLIAJck8bJ6PuGlYxF5EpJVdTO+DXxFpDX9DBUFEmomqdOcnjAvUKhOw2uHAaZmAacOB0zIBswwHzlIdeLGgT+VmAn4UrIuCdYMS+lQsCNhF72zAKAj4Oj54+0rw+PKx7VfGZ33+z92HwffVCjenC1jPtR9x9w8GjyUVE9km2H4KcKNbXgxfHDS0E/7v8C3gFNdehuQhxyOBA117Dn449uconO83dB6FcxSCVUaO/76RrYBr8f/TSdgxQ1pxlsHA7Vi9g6Pc/V0UH9K8C3YMMgULfJ5O78B2aCTwXWwk0rXArlQYz6smCBh+6KajgI6ku5jCLxOwKwKaP1JEWt07wK15d0JEKqIL2CLN5UoKAzjxzDCpn7RMwEqCgFGQL2thkFLVgdMyAbPOCVgqE3B+7DXSgoBDg/bHQXt68Hg8vpI0HDicHzCeDVgsCAgwLrb9Mu7+jWDde/jfISkICNb/qPrwtOA1Vo9tPxEfCHvU3U/BBxo3im3fBXzWte9w9//BFzfZlt52c/cvY0X3Lsf/ffdK2P4b+KDh4VgRkeg9e2zC9hsHrzEf/7v2xxc6CY3EgnjxIPY44Hp6v6/7YnOabhNbvzk2dHoovZ3ifiaq7DwU+A4WYF0mYfs9serPpwHbu+XfAw/jh16HihZM6erp6VlYbIMi4juKnrQNRbA3YtpcEyIirUz7N5HWoM+qSPPT5zQf4d89fm4fnffH/xfR+rTts66PXrtW6+N9rdX68Hcod329/qZpf4tS28dfu9K/aVIcKOl3y/I3DZ8ry/bhY8W2T3pvF4tlxfdBPQnLxbZPkva3SJP2tyi2fZZ19oQ9PT0K3omIiIiIiIiIiLSxvli6qIiIiIiIiIiIiLSp/wfkBTj/QwePSQAAAABJRU5ErkJggg=='
    return(diagram)


def exp_layout(exp_type):  #Generate a WCHORUS Layout
    if(exp_type=='WCHORUS'):  # define WCHORUS layout
        chorus_sequence = sequence_diagram()

        pulse1 = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                [sg.Input(key='tw1', size=(10,1))],
                [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                [sg.Input(key='rffactor1', size=(10,1))],
                [sg.Text('Shape parameter of P1:'), sg.Text(size=(10,1), key='N1_out')],
                [sg.Input(key='N1', size=(10,1))]]
        pulse2 = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                [sg.Input(key='tw2', size=(15,1))],
                [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                [sg.Input(key='rffactor2', size=(15,1))],
                [sg.Text('Shape parameter of P2:'), sg.Text(size=(10,1), key='N2_out')],
                [sg.Input(key='N2', size=(15,1))]]
        delay  = [[sg.Text('Duration of D3 (us):'), sg.Text(size=(10,1), key='tau3_out')],
                [sg.Input(key='tau3', size=(15,1))]]
        pulse3 = [[sg.Text('Duration of P3 (us):'), sg.Text(size=(10,1), key='tw3_out')],
                [sg.Input(key='tw3', size=(15,1))],
                [sg.Text('RF-Factor of P3:'), sg.Text(size=(10,1), key='rffactor3_out')],
                [sg.Input(key='rffactor3', size=(15,1))],
                [sg.Text('Shape parameter of P3:'), sg.Text(size=(10,1), key='N3_out')],
                [sg.Input(key='N3', size=(15,1))]]

        layout = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                [sg.Button('', image_data=chorus_sequence,
                button_color=(sg.theme_background_color(),sg.theme_background_color()),
                border_width=0, key='sequence')],
                [sg.Column(pulse1), sg.Column(pulse2), sg.Column(delay), sg.Column(pulse3)],
                [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta_out')],
                [sg.Input(key='delta')],
                [sg.Button('Start Simulation'), sg.Button('Exit')]
                ]
    elif(exp_type=='Double_Echo'):
        chorus_sequence = sequence_diagram()

        pulse1 = [[sg.Text('Duration of P1 (us):'), sg.Text(size=(10,1), key='tw1_out')],
                [sg.Input(key='tw1', size=(10,1))],
                [sg.Text('RF-Factor of P1:'), sg.Text(size=(10,1), key='rffactor1_out')],
                [sg.Input(key='rffactor1', size=(10,1))],
                [sg.Text('Shape parameter of P1:'), sg.Text(size=(10,1), key='N1_out')],
                [sg.Input(key='N1', size=(10,1))]]
        pulse2 = [[sg.Text('Duration of P2 (us):'), sg.Text(size=(10,1), key='tw2_out')],
                [sg.Input(key='tw2', size=(15,1))],
                [sg.Text('RF-Factor of P2:'), sg.Text(size=(10,1), key='rffactor2_out')],
                [sg.Input(key='rffactor2', size=(15,1))],
                [sg.Text('Shape parameter of P2:'), sg.Text(size=(10,1), key='N2_out')],
                [sg.Input(key='N2', size=(15,1))]]
        delay  = [[sg.Text('Duration of D3 (us):'), sg.Text(size=(10,1), key='tau3_out')],
                [sg.Input(key='tau3', size=(15,1))]]
        pulse3 = [[sg.Text('Duration of P3 (us):'), sg.Text(size=(10,1), key='tw3_out')],
                [sg.Input(key='tw3', size=(15,1))],
                [sg.Text('RF-Factor of P3:'), sg.Text(size=(10,1), key='rffactor3_out')],
                [sg.Input(key='rffactor3', size=(15,1))],
                [sg.Text('Shape parameter of P3:'), sg.Text(size=(10,1), key='N3_out')],
                [sg.Input(key='N3', size=(15,1))]]

        layout = [[sg.Text('Define Pulse Parameter here, then click start Simulation'), sg.Text(size=(15,1), key='text1_out')],
                [sg.Text('Shapes are saved in this folder and can be used by Bruker Spectrometer'), sg.Text(size=(15,1), key='text3_out')],
                [sg.Text('and the SIMPSON package.'), sg.Text(size=(15,1), key='text5_out')],
                [sg.Text('Simulations can also be restarted with changed Parameter from her.'), sg.Text(size=(15,1), key='text6_out')],
                [sg.Button('', image_data=chorus_sequence,
                button_color=(sg.theme_background_color(),sg.theme_background_color()),
                border_width=0, key='sequence')],
                [sg.Column(pulse1), sg.Column(pulse2), sg.Column(delay), sg.Column(pulse3)],
                [sg.Text('Sweepwidth (kHz):'), sg.Text(size=(15,1), key='delta_out')],
                [sg.Input(key='delta')],
                [sg.Button('Start Simulation'), sg.Button('Exit')]
                ]
    return layout


def selection_layout():  # Generates a selection Layout
    layout = [[sg.Text('Define Pulse Parameter here, then click start Simulation')],
              [sg.Listbox(values=('WCHORUS', 'Double_Echo'), size=(20, 12), key='exp_type', enable_events=True)],
              [sg.Button('Exit')]]
    return(layout)


# Create Sequence Selection Window
sequence_selection = sg.Window('Pulsesequence Phasecorrection Tool', selection_layout(), grab_anywhere=True)

# STEP3 - the event loop
while True:
    event, values = sequence_selection.read()   # Read the event that happened and the values dictionary
    print(event, values)
    if event == sg.WIN_CLOSED or event == 'Exit':     # If user closed window with X or if user clicked "Exit" button then exit
        break
    simulation_window = sg.Window('New Window Test', exp_layout(values['exp_type'][0]), grab_anywhere=True)
    while True:
        event, values = simulation_window.read()   # Read the event that happened and the values dictionary
        print(event, values)
        if event == sg.WIN_CLOSED or event == 'Exit':     # If user closed window with X or if user clicked "Exit" button then exit
            break
        if event == 'Start Simulation':
            # ADD SIMULATION PARAMETER SETUP HERE
            break
            # ADD SIMULATION HERE
    simulation_window.close()
        # # Update the "output" text element to be the value of "input" element
        # window['delta_out'].update(values['delta'])
        # window['tw1_out'].update(values['tw1'])
        # window['tw2_out'].update(values['tw2'])
        # window['tw3_out'].update(values['tw3'])
        # window['rffactor1_out'].update(values['rffactor1'])
        # window['rffactor2_out'].update(values['rffactor2'])
        # window['rffactor3_out'].update(values['rffactor3'])
        # window['tau3_out'].update(values['tau3'])
        # window['N1_out'].update(values['N1'])
        # window['N2_out'].update(values['N2'])
        # window['N3_out'].update(values['N3'])
        # # window['ss_offset_out'].update(values['ss_offset'])
        # # window['rms_limit_out'].update(values['rms_limit'])

        # delta       = values['delta']
        # tw1         = values['tw1']
        # tw2         = values['tw2']
        # tw3         = values['tw3']
        # rffactor1   = values['rffactor1']
        # rffactor2   = values['rffactor2']
        # rffactor3   = values['rffactor3']
        # tau1        = '0'
        # tau2        = '0'
        # tau3        = values['tau3']
        # N           = values['N1']

        # # ss_offset   = values['ss_offset']
        # # rms_limit   = int(values['rms_limit'])
        # ss_offset   = '50000'
        # rms_limit   = 500
        # poly_order = 42
        # start = 0.2
        # end = 0.8
        # outer_weight = 0.0
        # rmssteps = 0.01
        
        # create_simpson()
        # simulate_sequence()

        ########


sequence_selection.close()