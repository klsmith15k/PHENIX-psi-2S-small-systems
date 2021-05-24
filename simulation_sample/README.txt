These macros are run in a certain order to determine the second gaussian parameters for the data fits.

1. sanghoon simulations directory
    Place the text files to Sanghoon's embedded Run15pp simualtions for J/psi and psi(2S) here. The list can contain all files in the directory, whether or not they failed.  The list can be made using:
    ls /phenix/subsys/fvtx/shlim/simulation/Jpsi_PYTHIA_dimu_embedding/Jpsi_embed_pp_grp200/EMBED_C/PDST_Jpsi_embed_pp_grp200_*.root > filename.txt

2. list_sorter.C
    Then run the list_sorter.C , which prints a new text file to the above directory with only the simulations that completed.

3.  FillAnaTree directory
    Use Sanghoon's FillAnaTree module to run over the simulation files.  The mFillAnaTree.cc and mFillAnaTree.h files can be edited to include new cuts and histograms for your analysis

4. Run_fill_Ana_tree.C
   This will run the FillAnatree module

5. hadd
   The ROOT built-in hadd command can add together the group 200 and 201 files using 
   hadd new_filename.root jpsi_200.root jpsi_201.root

But the files are large, and this method only works if you load root6 first:  source /opt/sphenix/core/bin/sphenix_setup.csh -n
Also note that the FillAnaTree module should be run using root 5 (phenix environment), and not root6 (sphenix environment)

6. mass_sim_rebinning_macro.C 
    Run this macro to convert the 3D histograms from Sanghoon's AnaTree module into 1D histograms

7. CB_macro_sanghoon_psi2s_2nd_G_sigmafree_toy.C
    Run this macro after changing the input files.  You can comment out the tail parameter limits to test the correct values and also change the range of the fit. 

8.  Take the weighted average
    Since pp is a symmetric system, we take the wieghted average of the north and south arm results using the number of mass entries in the histograms ( GetEntries() )

9.  CB_macro_sanghoon_psi2s_2nd_G_sigmafree_ave.C
    Refit the mass histograms with the second gaussian values fixed to the weighted average to get the new ratio and fraction values.  Also use the same tail limits and fit range

10.  CB_macro_sanghoon_for_sim_check.C
     Change the input files and run the macro four times, changing the psi(2S) and the north_arm selections.  The output files are read into the toy model macro

11. random_sim_check_total.C
    Change the ratio and fraction entries at all instances in the macro as well as the tail limits and the fit range.  The north and south arm selections are in two places - double check you changed both

12. sim_check_ROOT/newlib_MC directory
    concatenate the text files into north and south arms using:
    cat S_sim*.txt > cat_all_S.txt

14. sim_check_ROOT directory
    Run plotting_macro_fvtx.C to see the toy model results for J/psi

15.  plotting_macro_fvtx_psi2s.C
     Run this macro to check the toy model results for the psi(2S) with the input gaussian parameters, fit range and tail parameter limits

16. data fits
    If the toy model has a systematic uncertainty you might expect, then the ratio and fraction values can be copied and pasted into the real data fitting macros and the same tail limits can also be used    
