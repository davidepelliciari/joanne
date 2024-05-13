# JOANNE: An injection pipeline of transient signals into real interferometric data

JOANNE,  inJection Of fAst traNsients iNto visibilitiEs, is a python-3 script which enables the injection of a given number of Fast Radio Bursts (FRBs) into radio interferometric visibility data.
Additionally, if specific calibration tables are already present, it allows simulating systematic phase errors by corrupting the phases of one or more antennas,
injecting a phase offset (only the parameter phi_0_c of the multi-band delay table for now) at one or more antennas. The script will split a CASA's MEASUREMENT SET (MS)
file into NFRBs MS files, of 1 second duration each around their time of arrivals (ToAs). Then, these splitted MS files will be combined via CASA's concat task.
A final dirty/cleaned image is created for the combined dataset, and the user can choose also to get single dirty images from each splitted visibilities.

## Usage:
If you have a UVFITS file available, running the main script from step "0" will convert it into a MEASUREMENT SET.
You can skip this step by starting from step "1" if you already have a MEASUREMENT SET file.

Once inside the CASA environment

```python
mysteps=[0,1,2,3,4,5]
execfile('go_joanne.py')
```

All parameters needed for the correct usage of the script are contained in the configuration file located in the directory "./config/conf.yaml".
These parameters are extracted from the main script using the Python package yaml.

A sample calibration script looks like this:

```python
general:
    listobs_path: /path/to/listobs/.listobs
    file_ToAs: /path/to/timeofarrivals.dat
    EXPERIMENT: test1.1
    output_path: /path/to/output_path/

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

It's important to set the paths under "general" and the "base_calib" path under "calibration," where the calibration tables should be located.
The configuration file must reside in the "./config/conf.yaml" directory. You can change this by modifying manually the main python script at 
row 750, just after the function definitions.

### Description of configuration file parameters

- ```file_ToAs```: this parameter permits to inject FRBs given an external file containing their time of arrivals (ToAs). The ToA format should be:

```python
2022/09/22/00:11:39
2022/09/22/00:28:51
...
```
in case the file does not exists, the script will generate randomly selected ToAs inside the observational scans (taken from listobs file).

- ```EXPERIMENT``` : the name of the given experiment / test you want to make (e.g. ```EXPERIMENT: 5FRBs_10Jy, TEST, ...```). All the output results will
be placed at ```base+EXPERIMENT``` directory.

- ```antID```: the ensamble of antennas used in your observation. You can find this from the listobs file, but you have to specify antID manually.

- ```modeIM_single, modeIM_concat```: select here if you want "DIRTY" or "CLEANED" images for individual splitted MS files. The same for the following parameters, but for the combined image.

- ```box_region, rms```: the pixel coordinates of a rectangular box inside which a gaussian fit will be performed after the cleaning of the combined image. ```rms``` defines
the flux threshold at which the global cleaning (CASA's ```tclean``` task) will stop.

- ```nFRBs, fl_ch, modeFL```: number of FRBs you want to inject. Each FRB will be a point-source lasting for a single integration time (the minimum possible duration).

The characteristic flux density of each FRB is determined by the ```fl_ch``` parameter. The user can choose different ways to generate FRBs fluxes, via the ```modeFL```
parameter. The implemented options are:

- equal: all FRBs will have same flux = ```fl_ch``` (units of Jy/beam);
- random: FRBs will have flux randomly distributed around a ```fl_ch``` flux with a ```fl_std``` standard deviation;
- powerlaw --> FRBs will have flux distributed as a power-law, with flux_min = ```fl_ch``` (i.e. the majority of bursts will have flux = ```fl_ch```), and a given
power-law index alpha. alpha must be > 0, otherwise a change of sign will take place.

The time of arrival of each FRB can be read from an external file, otherwise will be generated randomly inside a given scan from the main script.
You can choose the different ToA simulation modes via the ```mode_inj``` parameter. The implemented options are:
- from_file: read ToAs from external file. See above.
- uniform: the scans in which a given FRB will be injected is drawn randomly from the ensamble of scans of the given ```fld_inj``` field.
- heavenly: as "uniform", but the first and the last scan will contain an FRB.

- ```do_corrupt, ant_not_used```: determine the visibility corruption modality. The user can choose wheter
to corrupt or not with ```do_corrupt``` parameter. If this parameter is set to ```True```, then the other parameters will be effective. ```ant_not_used``` defines the antennas you want to
remove from the ```antID``` ensamble. NB: the same antenna must have been removed in the calibration steps used for obtaining the calibration
tables at ```base_calib``` path. If you specify ```ant_not_used: All```, all the antennas will be considered.

- ```corr_ant``` defines the antenna in which you want to inject a phase offset (i.e. the antenna of which calibration solution you want to corrupt).
  
- ```corruption``` defines the amount of corruption (in units of deg) to inject in a given ```corr_ant```antenna's phase solution.
  
- ```corrupt_table``` defines the type of calibration table in which the corruption will take place. For now only ```corrupt_table: mbd``` (mbd = multi-band delay) is implemented.
The corruption will take place in the ```phi_0_c``` parameter of the mbd table.

## Output format

The output of the code will be stored at ```output_path```+```EXPERIMENT``` directory, and includes combined and splitted MS files, along with images of the combined visibilities.
Additionaly, in the same directory the user will find a .txt file containing a 2D gaussian fit to the injected point source. The header of this output file is the following:

```python
nFRBs Time_interval_for_injections fitted_peak_flux rms_from_image S/N RA_fit +/- RA_err dec_fit +/- dec_err
```


## Dependencies

```python
numpy, matplotlib, scipy, math, random, re, yaml 6.0.1>, os, casatools, casavlbitools, datetime
```

For any questions regarding the pipeline please feel free to write to:
davide.pellciari@inaf.it
