# VLBI-FRB-injection

This script enables the injection of a given number of Fast Radio Bursts (FRBs) into radio interferometric visibility data. Additionally, if specific calibration tables are already present, it allows simulating systematic phase errors by corrupting the phases of one or more antennas, injecting a phase offset (only the parameter phi_0_c of the multi-band delay table for now) at one or more antennas.

## Usage:
If you have a UVFITS file available, running the main script from step "0" will convert it into a MEASUREMENT SET. You can skip this step by starting from step "1" if you already have a MEASUREMENT SET file.

Once inside the CASA environment

mysteps=[0,1,2,3,4,5]
execfile('VLBI_FRB_injection.py')

All parameters needed for the correct usage of the script are contained in the configuration file located in the directory "./config/conf.yaml". These parameters are extracted from the main script using the Python package yaml.

A sample calibration script looks like this:

```python
general:
    base: /path/to/base_path/
    listobs_path: /path/to/listobs/.listobs
    file_ToAs: /path/to/timeofarrivals.dat
    EXPERIMENT: test1.1

calibration:
    uvfile: /path/to/uvfits_file/.uvfits
    ms: /path/to/msfile/.ms
    base_calib: /path/to/calib_folder/
    targetFLAG: /path/to/flags/FLAGtarget.flagcmd
    bp_cal: J1829+4844
    ph_cal: J0529+3209
    target: R1_D
    sbdtab: ek051e.sbd
    mbdtab: ek051e.mbd
    bpasstab: ek051e.bpass
    Kselftab: J0529+3209_self.K
    pselftab: J0529+3209_self.p
    apselftab: J0529+3209_self.ap
    antID: EF,TR,WB,NT,O8,CM,DA,DU,KN,PI,DE,JM,MC

imaging:
    modeIM_single: DIRTY
    modeIM_concat: CLEAN
    imsize: 1280
    box_region: 626,626,654,654
    rms: 2.e-4

injection:
    nFRBs: 15
    fl_ch: 0.01
    fl_std: 0.025
    alpha_fl: -2
    modeFL: equal
    freq: 1.4GHz
    BW: 0.256GHz
    direction: J2000 05h31m58.700s +33d08m52.568s
    cell: 3.6mas
    t_start_inj: 2022/09/22/00:00:00
    t_stop_inj: 2022/09/22/05:48:39
    fld_inj: R1_D
    mode_inj: uniform

corruption:
    ant_not_used: All
    corr_ant: O8
    corruption: 1.0
    corrupt_table: mbd
```

t's important to set the paths under "general" and the "base_calib" path under "calibration," where the calibration tables should be located.
