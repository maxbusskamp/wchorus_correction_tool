#!/bin/bash

#                           filename  tw1  tw2   tw3  exp_type                 delta1 delta2 delta3 var11  var12  var13  rffactor1  rffactor2  rffactor3  tau1  tau2  tau3  'none'  ss_offset  var21  var22  var23  shape_type  var31  var32  var33  phaseoff1  phaseoff2  phaseoff3  compression
time simpson inputfile.tcl  test_seq  500  1250  1000 compressedCHORUS_cycled  300    300    300    40     0      0      0.265      1.0        1.0        0     250   0     'none'  25000      80     0      0      WURST       120    0      0      0          0          0          0
