# Main task:
  * Go over [Aim 1](/Notebook/P06.md)
  * Go over how to tackle Aim 2 objective. 

# Notes after dicussion with Avril

## Recap
* the “Cage ID” showed the largest R^2 value among all the tested variables in [Aim 1](/Notebook/P06.md); however, this category does not have biological implication. Hence we want to run the “Cage ID value again to check whether we need to control for the variable before re-running Aim 1 analysis
  * Test cage id to check if we need to control for it during the experiment because we want to look for biologically significant/useful results
  * Control for them = take note of “experimental error”
  * Run simple adonis dis~cageID (as justification of why we need to control for it)

## Re-subset the metadata input to the loop
* Variables
  * Age: there are multiple columns with ages in months, days, and in categorical characters
    * Only use “Age.New.Bin”: young, middle old
  * Experiment group and treatment type: both refer to cohouse or separate the mice (noticed that “Experiment group” also have succinate treatment - but was not explained in the original paper Cheney et al., 2022)
    * Subset metadata with “Separate” in the “Treatment.type” column and remove “Succinate_experiment” in the “Experiment.Group” column
    * We exclude the cohousing ones to prevent the microbiota of mice become similar over time
    * Expected number of mice in both the control and mutant genotype after subsetting:
      * control : 11 from “General” and 3 from “Cohousing_experiment” = 14 in total
        * Noted that those shown to be from “Cohousing_experiment” is the control for cohousing experiment, where the mice were housed separately
      * Mutant: 14 in total 
