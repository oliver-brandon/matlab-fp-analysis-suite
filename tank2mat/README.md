# tank2mat
## About
This collection of MATLAB scripts can be used to convert SINGLE FIBER TDT tank fiber photometry recordings into .mat files, split the .mat file into separate .mat files for each animal, and finally infuse PRL trial epocs (cue, correct_rewarded, correct_noreward, incorrect_rewarded, and incorrect_noreward) into the .mat files and save them.

## Instructions
Each step listed in "About" is its own script. Thus, to effectively use, the correct order is (1) tank2mat, (2) mat2mats, and (3) prl_epoc_infuser.    