### --  Script written by D. Pelliciari (Ph. D. XXXVII cycle, University of Bologna) -- ###
### -- Inject multi fake point sources in a MV visibility dataset in CASA

import numpy as np
import os
import os.path
from casatools import table
from scipy.constants import c
import matplotlib.pyplot as plt
from casatools import image
import math
from datetime import datetime, timedelta
import random
import re
import yaml


###############################################################################################
## just for better plotting (used only for plot_tab function)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=20)
plt.rc('axes', linewidth=2)

###############################################################################################
###############################################################################################


## read_config: extract parameters from a config.yaml file
def read_config(config="config/conf.yaml"):

    with open(config, 'r') as stream:
        data_loaded = yaml.load(stream, yaml.SafeLoader)

    general_params = []
    calibration_params = []
    imaging_params = []
    injection_params = []
    corruption_params = []

    # general paths
    general_params.append(data_loaded['general']['listobs_path'])
    general_params.append(data_loaded['general']['file_ToAs'])
    general_params.append(data_loaded['general']['EXPERIMENT'])
    general_params.append(data_loaded['general']['output_path'])

    # calibration parameters
    calibration_params.append(data_loaded['calibration']['uvfile'])
    calibration_params.append(data_loaded['calibration']['ms'])
    calibration_params.append(data_loaded['calibration']['base_calib'])
    calibration_params.append(data_loaded['calibration']['targetFLAG'])
    calibration_params.append(data_loaded['calibration']['bp_cal'])
    calibration_params.append(data_loaded['calibration']['ph_cal'])
    calibration_params.append(data_loaded['calibration']['target'])
    calibration_params.append(data_loaded['calibration']['sbdtab'])
    calibration_params.append(data_loaded['calibration']['mbdtab'])
    calibration_params.append(data_loaded['calibration']['bpasstab'])
    calibration_params.append(data_loaded['calibration']['Kselftab'])
    calibration_params.append(data_loaded['calibration']['pselftab'])
    calibration_params.append(data_loaded['calibration']['apselftab'])
    calibration_params.append(data_loaded['calibration']['antID'])

    # imaging parameters
    imaging_params.append(data_loaded['imaging']['modeIM_single'])
    imaging_params.append(data_loaded['imaging']['imsize'])
    imaging_params.append(data_loaded['imaging']['box_source'])
    imaging_params.append(data_loaded['imaging']['box_bkg'])
    imaging_params.append(data_loaded['imaging']['keep_single_ms'])

    # injection parameters
    injection_params.append(data_loaded['injection']['nFRBs'])
    injection_params.append(data_loaded['injection']['fl_ch'])
    injection_params.append(data_loaded['injection']['fl_std'])
    injection_params.append(data_loaded['injection']['alpha_fl'])
    injection_params.append(data_loaded['injection']['modeFL'])
    injection_params.append(data_loaded['injection']['freq'])
    injection_params.append(data_loaded['injection']['BW'])
    injection_params.append(data_loaded['injection']['direction'])
    injection_params.append(data_loaded['injection']['cell'])
    injection_params.append(data_loaded['injection']['t_start_inj'])
    injection_params.append(data_loaded['injection']['t_stop_inj'])
    injection_params.append(data_loaded['injection']['fld_inj'])
    injection_params.append(data_loaded['injection']['mode_inj'])

    # corruption parameters
    corruption_params.append(data_loaded['corruption']['do_corrupt'])
    corruption_params.append(data_loaded['corruption']['ant_not_used'])
    corruption_params.append(data_loaded['corruption']['corr_ant'])
    corruption_params.append(data_loaded['corruption']['corruption'])
    corruption_params.append(data_loaded['corruption']['corrupt_table'])
    corruption_params.append(data_loaded['corruption']['corr_mode'])
    corruption_params.append(data_loaded['corruption']['corr_min'])
    corruption_params.append(data_loaded['corruption']['corr_max'])
    corruption_params.append(data_loaded['corruption']['corr_parameter'])

    return general_params, calibration_params, imaging_params, injection_params, corruption_params


###############################################################################################


## get_sigmas: get a gaussian distribution around 0 mean and fact standard-deviation, with size N
def get_sigmas(fact, N):

    return np.random.normal(0., fact, N)


###############################################################################################


## inject_phase_offset: modify calibration table injecting an offset in phi0_c (to be implement: other parameters)
def inject_phase_offset(arr, ant, sigma, antID):

    N_ant = len(antID)
    newarr = np.zeros_like(arr)  # Crea una copia di arr

    if ant != "All" and ant != "None":

        for a in ant:

            id_ant = antID.index(a)
            sig = sigma[id_ant]
            array_ant = newarr[id_ant::N_ant]

            #noisy_array = [array_ant[i]+sig for i in range(0,len(array_ant))]
            noisy_array = [sig for i in range(0,len(array_ant))]

            #noisy_array = array_ant + np.random.normal(0, sig, len(array_ant))   #inject gaussian noise

            newarr[id_ant::N_ant] = noisy_array

    elif ant == "All":

        for ii in range(0, len(antID)):

            id_ant = ii
            array_ant = newarr[ii::N_ant]
            #noise_ant = [array_ant[i] + sigma[ii] for i in range(0,len(array_ant))]
            noise_ant = [sigma[ii] for i in range(0,len(array_ant))]
            #noise_ant = array_ant + np.random.normal(0, sigma[ii], len(array_ant))    #inject gaussian noise
            newarr[ii::N_ant] = noise_ant

    else:

        for ii in range(0, len(antID)):

            id_ant = ii
            array_ant = newarr[ii::N_ant]
            #noise_ant = [array_ant[i] + sigma[ii] for i in range(0,len(array_ant))]
            noise_ant = [0. for i in range(0,len(array_ant))]
            #noise_ant = array_ant + np.random.normal(0, sigma[ii], len(array_ant))    #inject gaussian noise
            newarr[ii::N_ant] = noise_ant
        
        print("## Doing nothing")

    return newarr


###############################################################################################


## plot_tab: plot calibration table after corruption
def plot_tab(phi, N_ant, antID):

    #print("phi che viene plottata: ")
    #print(phi)

    outfile = open("/home/davide/PHD/joanne/mytest/delay_test_distributions.dat", 'w')

    plt.figure(figsize=(12,8))

    for j in range(0,N_ant):
    
        phiNc = [phi[i] for i in range(j,len(phi), N_ant)]
    
        time = np.linspace(0,len(phiNc), len(phiNc))

        print("antID: ", antID[j], " media phiNC: ", np.mean(phiNc), " standard-dev: ", np.std(phiNc))
        outfile.write(str(antID[j])+" "+str(np.mean(phiNc))+" "+str(np.std(phiNc))+"\n")
        plt.plot(time, phiNc, marker='.', ls='--', markersize=10,label=antID[j])#, color=cl)


    #plt.ylim(-np.pi,np.pi)
    plt.xlabel("Scan number", fontsize=23)
    plt.ylabel(r"$\phi_{0,i}$", fontsize=23)
    #plt.text(0,45,'LL', fontsize=30)
    plt.legend(ncol=3, fontsize=16)
    plt.show()


###############################################################################################


## corruptVIS: corrupt a given calibration table (e.g. multi-band delay) and apply it to corrected visibilities
## mode: manual = choose manually the corruption factor by passing a "corruption" parameter in input
## mode: external = corrFACT is an array with all the corruptions.
## mode: TimeDep = as external, but the corrFACT array changes for each injection
def corruptVIS(ms, corrFACT, origtable, corr_table, corr_parameter, antID, mask, base_calib, mode='manual'):

    print("### corrupting visibilities!")
    new_corrected = []
    correct_data = []
    corr_bkg = []

    tb.open(ms, nomodify=False)     ######### open the MS with DATA - CORR. DATA - MODEL_DATA

    index_parameters_L = [0,1,2]
    index_parameters_R = [4,5,6]

    mapping = {"phase": 0, "delay": 1, "rate": 2}
    index = mapping.get(corr_parameter, 0)

    idx_L = index_parameters_L[index]
    idx_R = index_parameters_R[index]

    data_column = tb.getcol('DATA')
    corr_bkg = tb.getcol('CORRECTED_DATA')
    model_data = tb.getcol('MODEL_DATA')

    # SWAP the two

    tb.putcol('DATA', model_data)

    tb.close()

    if mode == 'manual':
        str_corrFACT = "_"+str(corrFACT)
        if corrFACT == 0.:
            str_corrFACT = ""
    elif mode == 'external' or mode == 'TimeDep':
        str_corrFACT = "_"+str(np.mean(corrFACT))

    N_ant = len(antID)
    antUSED = [antID[i] for i, m in enumerate(mask) if m == 1]

    allANT = all(ants == 1 for ants in mask[1:])
    noANTS = all(ants == 0 for ants in mask)

    if allANT:

        print("## Considering all antennas", allANT)
        antUSED = "All"
        str_corrFACT = 'All'
        #mask[0] = 0                ## UNCOMMENT THIS TO HAVE NULL EF CORRUPTIONS

    if noANTS:

        print("## No corruption applied", noANTS)
        antUSED = "None"

    print("### copying ", origtable, " into ", corr_table)
    os.system('cp -r '+origtable+" "+corr_table)
    print("### done with the copy")

    #tb.open(corr_table, nomodify=False)
    tb.open(origtable, nomodify=False)
    fpar = tb.getcol('FPARAM')

    phi0c_L = fpar[idx_L,0,:]
    phi0c_R = fpar[idx_R,0,:]

    #plot_tab(phi0c_L, N_ant, antID)


    if mode == 'manual':
        sigma = [float(corrFACT)*mask[ii] for ii in range(0,len(mask))]
    elif mode == 'external' or mode == 'TimeDep':
        sigma = [float(corrFACT[ii])*mask[ii] for ii in range(0,len(mask))]
    sigmaR = sigma

    print("### creating phi_noise, considering ", sigma)
    
    #print("ANTENNA DA USARE PER LA CORRUZIONE: ", antUSED)
    phi_noise = inject_phase_offset(phi0c_L, antUSED, sigma, antID)
    phi_noiseR = inject_phase_offset(phi0c_R, antUSED, sigmaR, antID)
    
    #plot_tab(phi_noise, N_ant, antID)

    newfpar = np.zeros_like(fpar)
    newfpar[idx_L,0,:] = phi_noise
    newfpar[idx_R,0,:] = phi_noiseR

    tb.putcol("FPARAM", newfpar)    # sostituisco la colonna 'FPARAM' con la nuova appena calcolata                                                                                                         
    tb.flush()   # flush the current content to disk

    # Apply corrupted calibration tables to the new DATA_COLUMN

    default(applycal)
    applycal(vis=ms,
        field='',
        gaintable=[corr_table],
        interp=['linear'],
        spwmap=[[0,0,0,0,0,0,0,0]])
    
    tb.open(ms, nomodify=False)
    correct_data = tb.getcol('CORRECTED_DATA')
    new_corrected = correct_data + corr_bkg

    tb.putcol('CORRECTED_DATA', new_corrected)
    tb.putcol('DATA', data_column)
    tb.flush()
    tb.close()


###############################################################################################


## getIMAGE: get the dirty image of a MS-file, and then exports a fits of that image
def getIMAGE(ms, antennae, imout, mode, imsize=1280, rms=None):

    if mode == 'DIRTY':
        print("## Making a dirty image")
        nitter = 0
        inter = False
        thresh = None
    else:
        print("## Making a cleaned image")
        nitter = 1000
        inter = False
        thresh = "{:.2f}".format(rms*1000.)
        thresh = thresh+'mJy'

    print("## Eliminating imout: ", imout)
    os.system('rm -r '+imout+'*')

    default(tclean)
    tclean(vis=ms, stokes='pseudoI',
           psfcutoff=0.5,
           cell=cell,
           #antenna='EF&*',
           #phasecenter=direction,
           datacolumn='corrected',
           imagename=imout,
           antenna=antennae,
           imsize=[imsize,imsize],
           deconvolver='clark',
           #weighting='standard',
           niter=nitter,
           #cycleniter=150,
           threshold=thresh,
           interactive=inter,
           #nsigma=3,
           parallel=False)

    print('## Done with tclean, exporting DIRTY image in FITS file..')

    default(exportfits)
    exportfits(imagename=imout+'.image', fitsimage=imout+'.fits', history=False)

    #rms=imstat(imagename=imout+'.image', box='53, 892, 1232, 1250')['rms'][0]
    #print("rms of your image: ", rms)


###############################################################################################


## modelPointSource: inject a point source from image plane into visibilities (MODEL_DATA column)
def modelPointSource(ms, flux, direction, freq, BW, cell, radir, decdir, im_name, out_name):

    os.system("rm -rf "+out_name)    
    cl.done()

    # add a continuum component at given direction, flux, freq. and with given "shape"

    #cl.addcomponent(dir=direction, flux=flux, fluxunit='Jy', freq=freq, shape="Gaussian"), majoraxis="10.0mas", minoraxis='10.0mas', positionangle='0.0rad')
    cl.addcomponent(dir=direction, flux=flux, fluxunit='Jy', freq=freq, shape="point")

    # define your image. It has to be the same size of the real one (from tclean)
    ia.fromshape(im_name+".im",[1280,1280,1,1],overwrite=True)
    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])

    # convert cell_size in "rad"
    cell_rad=qa.convert(qa.quantity(cell),"rad")['value']

    cs.setincrement([-cell_rad,cell_rad], 'direction')
    cs.setreferencevalue([qa.convert(radir,'rad')['value'],qa.convert(decdir,'rad')['value']],type="direction")
    cs.setreferencevalue(freq,'spectral')
    cs.setincrement(BW,'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(),subtract=False)
    #exportfits(imagename=im_name+".im",fitsimage=im_name+'.fits',overwrite=True, dropdeg=True)

    cl.rename(out_name)
    cl.done()

    print("## Fourier transforming back your component list..")

    # Convert the injected component into model visibility, inserting in MODEL_DATA column.

    default(ft)
    ft(vis=ms, complist=out_name, usescratch=True)
    print("## Done with ft task")


###############################################################################################


## getANT
def getANT(ANT):

    if len(ANT)==1:
        ANT = ANT[0]
        ant_toCorr = ANT
    else:
        ant_toCorr = ','.join(ANT)
        ANT = ''.join(ANT)

    return ANT, ant_toCorr


###############################################################################################


## parse_datetime: convert a string containing a date (YYYY/MM/DD/hh:mm:ss) into a datetime object
def parse_datetime(datetime_str):

    return datetime.strptime(datetime_str, "%Y/%m/%d/%H:%M:%S")


###############################################################################################


## get_ToAs_in_scan: returns time of arrival (ToA) of FRBs that has to be injected, from t_start and t_stop of a given scan
def get_ToAs_in_scan(t_start, t_stop):
    
    times = []
    start_time = datetime.strptime(t_start, '%Y/%m/%d/%H:%M:%S')
    stop_time = datetime.strptime(t_stop, '%Y/%m/%d/%H:%M:%S')
    start_time += timedelta(seconds=1)
    while start_time <= stop_time:
        times.append(start_time.strftime('%Y/%m/%d/%H:%M:%S'))
        start_time += timedelta(seconds=2)
    
    return times


###############################################################################################


## get_injected_times: returns num_events ToAs starting from the scans of the observations 
def get_injected_times(dirTOA, scan_data, num_events, mode='uniform', t_start=None, t_stop=None):
   
    num_scans = len(scan_data)
    event_counts = np.zeros(num_scans, dtype=int)

    if t_start is not None:
        t_start = parse_datetime(t_start)
    if t_stop is not None:
        t_stop = parse_datetime(t_stop)
    if t_start is None and t_stop is None:
        start_index = 0
        end_index = num_scans
    elif t_start is None:
        start_index = 0
        end_index = next(i for i, (start, stop) in enumerate(scan_data) if parse_datetime(stop) >= t_stop)
    elif t_stop is None:
        start_index = next(i for i, (start, stop) in enumerate(scan_data) if parse_datetime(start) >= t_start)
        end_index = num_scans
    else:
        start_index = next(i for i, (start, stop) in enumerate(scan_data) if parse_datetime(start) >= t_start)
        end_index = next(i for i, (start, stop) in enumerate(scan_data) if parse_datetime(stop) >= t_stop)

    scan_data = scan_data[start_index:end_index]
    print("## Number of scans considered for the injection: ", len(scan_data))
    event_counts = np.zeros(len(scan_data), dtype=int)
    num_scans_in_range = len(scan_data)

    if num_events > 3:
        if mode == 'uniform':
            total_events = num_events
            avg_events_per_scan = total_events / num_scans_in_range
            total_generated_events = 0
            while total_generated_events < total_events:
                for i in range(num_scans_in_range):
                    remaining_events = total_events - total_generated_events
                    events_for_scan = min(np.random.poisson(avg_events_per_scan), remaining_events)
                    event_counts[i] += events_for_scan
                    total_generated_events += events_for_scan

        if mode == 'heavenly':
            event_counts[0] = 1
            event_counts[-1] = 1
            num_scans_in_range = len(scan_data)-2
            total_events = num_events
            avg_events_per_scan = total_events / (num_scans_in_range)
            total_generated_events = 2
            while total_generated_events < total_events:
                for i in range(1,num_scans_in_range+1):
                    remaining_events = total_events - total_generated_events
                    events_for_scan = min(np.random.poisson(avg_events_per_scan), remaining_events)
                    event_counts[i] += events_for_scan
                    total_generated_events += events_for_scan

    elif num_events == 3 and mode == 'heavenly':
        event_counts[0] = 1
        event_counts[-1] = 1
        mid_index = int(num_scans_in_range / 2)
        event_counts[mid_index] = 1
    elif num_events == 2 and mode == 'heavenly':
        event_counts[0] = 1
        event_counts[-1] = 1
    else:
        mid_index = int(num_scans_in_range / 2)
        event_counts[mid_index] = 1

    print("## Number of instant point sources per scan: ")
    print(event_counts)
    sel_tstart = [sublist[0] for sublist in scan_data]
    sel_tstop = [sublist[1] for sublist in scan_data]
    sel_toas = []
    sel_scans = []
    output_file = open(dirTOA+"newTOAs.dat", 'w')

    for ii in range(0,len(scan_data)):
        toas_per_scan = get_ToAs_in_scan(sel_tstart[ii], sel_tstop[ii])
        num_events_per_scan = event_counts[ii]
        if num_events_per_scan > 0:
            selected_indices = random.sample(range(len(toas_per_scan)), min(num_events_per_scan, len(toas_per_scan)))
            for index in selected_indices:
                sel_toas.append(toas_per_scan[index])
                output_file.write(str(toas_per_scan[index])+"\n")
                
    print("## Number of time of arrivals: ", len(sel_toas))
    output_file.close()

    return sel_toas


###############################################################################################


## get_DeltaT: compute the difference in hh and mm between two dates
def get_DeltaT(t_start, t_stop):

    start = datetime.strptime(t_start, "%Y/%m/%d/%H:%M:%S")
    stop = datetime.strptime(t_stop, "%Y/%m/%d/%H:%M:%S")
    diff = stop - start
    hh = diff.seconds // 3600
    mm = (diff.seconds % 3600) // 60
    delta = str(hh)+"h"+str(mm)+"m"
    
    return delta


###############################################################################################


## is_non_empty_directory: check if the given path is an existing and non-empty directory
def is_non_empty_directory(path):

    if os.path.isdir(path):
        contents = os.listdir(path)
        if contents:
            return True

    return False



###############################################################################################


## find_date_to_add: find the date (format: dd-mmm(str)-yyyy) that has to be added to 2nd and 3rd column of listobs
def find_date_to_add(rows):

    for index, row in enumerate(reversed(rows)):
        match = re.search(r'\d{2}-[A-Za-z]{3}-\d{4}/', row)
        if match:
            return match.group(0), len(rows) - index - 1

    return None, None


###############################################################################################


## add_date: add a datetime to hours minutes and seconds of 2nd and 3rd cols of listobs
def add_date(row, date, mode):

    parts = row.split(" ", 1)
    if mode == 'multi':
        part1 = date + parts[0]
        part2 = date + parts[1]
    elif mode == 'single':
        part1 = parts[0]
        part2 = date + parts[1]

    return part1 + " " + part2


###############################################################################################


## skim_listobs: takes as input a CASA listob file and returns a path where a skimmed (csv-like) version of the listob is stored
def skim_listobs(input_file_path, output_file_path):

    outFile = open(output_file_path, "w")
    nrows_line_index = None
    with open(input_file_path, "r") as input_file:
        lines = input_file.readlines()
    for i, line in enumerate(lines):
        if "nRows = Total number of rows per scan" in line:
            nrows_line_index = i
            break
    outFile.writelines(lines[9:nrows_line_index])
    header = "Timerange_start_UTC Timerange_stop_UTC Scan FldId FieldName nRows SpwIds Average_Intervals"
    outFile.close()
    with open(output_file_path, 'r') as file:
        lines = file.readlines()
    updated_lines = []
    updated_lines.append(header + '\n')
    for ii in range(1, len(lines)):
        line = lines[ii]
        line = line.replace(" - ", " ")
        line = line.replace("[2, 2, 2, 2, 2, 2, 2, 2]", "[2,2,2,2,2,2,2,2]")
        updated_lines.append(line)
    date, idx_date = find_date_to_add(updated_lines)
    if date is None:
        print("## WARNING: no date time found to be added to listobs..")
        return
    unmod_rows = [add_date(updated_lines[ii].strip(), '', 'multi') for ii in range(1, idx_date+1)]
    first_dates = []
    for rows in unmod_rows:
        first_match = re.search(r'\d{2}-[A-Za-z]{3}-\d{4}/', rows)
        if first_match:
            first_dates.append(first_match.group(0))
            #print(first_match)
    unmod_rows = [add_date(updated_lines[ii].strip(), first_dates[ii-1], 'single') for ii in range(1, idx_date+1)]
    mod_rows = [add_date(updated_lines[ii].strip(), date, 'multi') for ii in range(idx_date+1, len(updated_lines))]
    with open(output_file_path, "w") as file_out:
        file_out.write(updated_lines[0])
        file_out.write("\n".join(unmod_rows))
        file_out.write("\n")
        file_out.write("\n".join(mod_rows))


###############################################################################################


## convert_time: convert a time string in another format (from 01-Sep-2022/20:02:00.0 to 2022/09/01/20:02:00)
def convert_time(time_str):

    date_str, hours_str = time_str.split('/')
    date = datetime.strptime(date_str, "%d-%b-%Y").strftime("%Y/%m/%d")
    converted_time = date+'/'+hours_str

    return converted_time[:-2]


###############################################################################################


## get_scans: reads a modified listobs file and returns an output file with a list of t_start and t_stop of each scan for a given field
def get_scans(input_file_path, FieldName = "R1_D"):

    with open(input_file_path, "r") as input_file:
        selected_scans = []
        for line in input_file:
            if "Scan" in line:
                continue
            columns = line.split()
            t_start, t_stop, Scan, FldId, FldName, nRows, SpwIds, Average_Intervals = columns
            if FldName == FieldName:
                selected_scans.append([convert_time(t_start), convert_time(t_stop)])

    return selected_scans


###############################################################################################


## getAntennae: from a string indicating which antenna do not use, returns the antennae string vector
def getAntennae(ANT_notUSED):

    ANT_notUSED_dir, ant_notUSED = getANT(ANT_notUSED)
    if ANT_notUSED == ['eMer']:
        ANT_notUSED = ["CM","DA","DU","KN","PI","DE","JM"]
    elif ANT_notUSED == ['None']:
        print("## No antenna selected! Changing ant to *&* (all)")
        antennae = '*&*'
    elif ANT_notUSED == ["CM","DA","DU","KN","PI","DE","JM"]:

        ANT_notUSED_dir = 'eMer'
        antennae='*;!CM,DA,DU,KN,PI,DE,JM'
    else:
        antennae='*;!'+ant_notUSED

    return antennae, ANT_notUSED_dir, ant_notUSED


###############################################################################################


## getAntennae_toCorrupt
def getAntennae_toCorrupt(corrFACT, corrANT, ANT_notUSED, doCorrupt, antID):

    CORRANT, ant_toCorr = getANT(corrANT)

    if ((corrANT == ANT_notUSED and corrANT != ['All']) or corrANT == ['None']):
        doCorrupt = False
        mask = [0 for ii in range(0,len(antID))]
        return mask, CORRANT, ant_toCorr, doCorrupt

    elif corrANT != ANT_notUSED:
        for cANT in corrANT:
            if cANT in ANT_notUSED:
                corrANT.remove(cANT)
        if len(corrANT) == 0 or corrFACT == 0:
            doCorrupt = False
            mask = [0 for ii in range(0,len(antID))]
            return mask, CORRANT, ant_toCorr, doCorrupt
        if corrANT == ['All']:
            mask = [1 for ii in range(0,len(antID))]
            return mask, CORRANT, ant_toCorr, doCorrupt
        else:
            mask = [1 if ant in corrANT else 0 for ant in antID]
            return mask, CORRANT, ant_toCorr, doCorrupt


###############################################################################################


## get_fluxes: get the fluxes for a given number of generated FRBs. fl_ch = characteristic flux used.
## modeFL can be:
## modeFL = equal --> All FRBs will have same flux = fl_ch
## modeFL = random --> FRBs will have flux randomly distributed around a fl_ch flux with a fl_std standard deviation
## modeFL = powerlaw --> FRBs will have flux distributed as a power-law, with flux_min = fl_ch (i.e. the majority of bursts will have flux = fl_ch), and a given
## power-law index alpha. alpha must be > 0, otherwise a change of sign will take place

def get_fluxes(nFRBs, modeFL, fl_ch, fl_std=None, alpha=None):

    if modeFL == 'equal':
        flux = [fl_ch for ii in range(0,nFRBs)]
        str_fl = 'eq'
        return flux, str_fl
    elif modeFL == 'random' and fl_std != None:
        flux = np.random.normal(fl_ch, fl_std, nFRBs)
        str_fl = 'rnd'
        return flux, str_fl
    elif modeFL == 'powerlaw' and alpha != None:
        if alpha > 0.:
            alpha = -alpha
        u = np.random.uniform(0, 1, size)
        flux = fl_ch*(1 - u) ** (1 / alpha)
        flux = PL_inverse_cdf(alpha_fl, fl_ch, nFRBs)
        str_fl = 'pl'
        return flux, str_fl
    else:
        print("## WARNING: bad flux mode. Considering same flux for all FRBs")
        if fl_ch!=None:
            print("## flux will be ", fl_ch, " Jy/beam")
            flux = [fl_ch for ii in range(0,nFRBs)]
            str_fl = 'eq'
            return flux, str_fl
        else:
            print("## considering 1 Jy/beam flux")
            flux = np.ones(nFRBs)
            str_fl = 'eq'
            return flux, str_fl


###############################################################################################


# extractFIT: extracts best-fit RA Dec (+- errors) from an image (.image) given a box_source. The information is extracted from a log file written by imfit task.
def extractFIT(im_name, box_source, logFILE):

    default(imfit)
    imfit(imagename=im_name, box=box_source, logfile=logFILE)

    ra_lines = ""
    dec_lines = ""
    pattern_ra = re.compile(r"--- ra:\s+(.*?)(?=\()")
    pattern_dec = re.compile(r"--- dec:\s+(.*?)(?=\()")

    with open(logFILE, 'r') as file:
        lines = file.readlines()
        ra_line = lines[20].strip()
        dec_line = lines[21].strip()
        ra_match = pattern_ra.search(ra_line)
        if ra_match:
            ra_info = ra_match.group(1).strip()
            ra_lines = "--- ra: " + ra_info
        else:
            print("## No correspondence for ra_line", ra_line)

    dec_info = dec_line.split(":")[1].strip()
    dec_lines = "--- dec: " + dec_info
    RESULTS_RA = re.sub(r'--- ra: ', '', ra_lines)
    RESULTS_DEC = re.sub(r'--- dec: ', '', dec_lines)

    return RESULTS_RA, RESULTS_DEC


###############################################################################################


## convert_dec: takes a declination string and returns a formatted string, ready to be read by astropy
def convert_dec(dec):
    match = re.search(r'([-+]?\d+\.\d+\.\d+\.\d+)', dec)
    if match:
        num = match.group(1)
        num_nozero = re.sub(r'\+0*', '+', num)
        num_nozero = re.sub(r'-0*', '-', num_nozero)
        formatted = re.sub(r'(\d+)\.(\d+)\.(\d+)\.(\d+)', r'\1:\2:\3.\4', num_nozero)
        return dec.replace(num, formatted)
    else:
        return None


# ---------------------------------------------------------------------------------------------------------------------------------------------- #
############################################################# Main code #########################################################################
# ---------------------------------------------------------------------------------------------------------------------------------------------- #


config_file = "./config/conf.yaml"
print(config_file)

general_params, calibration_params, imaging_params, injection_params, corruption_params = read_config(config_file)

# general parameters
listobs_path = general_params[0]
file_ToAs = general_params[1]
EXPERIMENT = general_params[2]
output_path = general_params[3]

# calibration parameters
uvfile = calibration_params[0]
ms = calibration_params[1]
base_calib = calibration_params[2]
targetFLAG = calibration_params[3]
bp_cal = calibration_params[4]
ph_cal = calibration_params[5]
target = calibration_params[6]
sbdtab = base_calib+calibration_params[7]
mbdtab = base_calib+calibration_params[8]
bpasstab = base_calib+calibration_params[9]
Kselftab = base_calib+calibration_params[10]
pselftab = base_calib+calibration_params[11]
apselftab = base_calib+calibration_params[12]
ants = calibration_params[13]
antID = ants.split(',')#


# imaging parameters
modeIM_single = imaging_params[0]
imsize = imaging_params[1]
box_source = imaging_params[2]
box_bkg = imaging_params[3]
keep_ms = imaging_params[4]

# injection parameters
nFRBs = injection_params[0]
fl_ch = injection_params[1]         ## uncomment for simulations purpouse
fl_std = injection_params[2]
alpha_fl = injection_params[3]
modeFL = injection_params[4]
freq = injection_params[5]
BW = injection_params[6]
direction = injection_params[7]
cell =  injection_params[8]
t_start_INJ = injection_params[9]
t_stop_INJ = injection_params[10]
fld_INJ = injection_params[11]
modeINJ = injection_params[12]

# visibility corruption parameters
doCorrupt = bool(corruption_params[0])
ant_not_used = corruption_params[1]
corr_ant = corruption_params[2]
corruption = corruption_params[3]           ## uncomment for simulations purpouse
corrupt_table = corruption_params[4]
mode_corr = corruption_params[5]
corr_min = corruption_params[6]
corr_max = corruption_params[7]
corr_parameter = corruption_params[8]

base_EXP = output_path+EXPERIMENT+"/"
outputFIT_path = base_EXP+EXPERIMENT+"_FITlog.txt"

newlistobs_path = output_path+'skimmedListobs.dat'

flux, str_fl = get_fluxes(nFRBs, modeFL, fl_ch, fl_std, alpha_fl)
str_fl = str_fl+str(fl_ch)
epoch, radir, decdir = direction.split(' ')
deltaT = get_DeltaT(t_start_INJ, t_stop_INJ)

if is_non_empty_directory(base_EXP) == False:

    print("## "+base_EXP+" is not an existing, empty directory!")
    print("## ..creating it")
    os.system('mkdir '+base_EXP)


if 'nsim' in locals() or 'nsim' in globals():
    if nsim is not None:
        try:
            outDIR_INJ = base_EXP + str(nFRBs) + "FRB_" + deltaT + "_" + str_fl + "Jy_" + str(nsim)+"_"
        except TypeErrofr:
            outDIR_INJ = base_EXP + str(nFRBs) + "FRB_" + deltaT + "_" + str_fl + "Jy_"
    else:
        outDIR_INJ = base_EXP + str(nFRBs) + "FRB_" + deltaT + "_" + str_fl + "Jy_"
else:
    outDIR_INJ = base_EXP + str(nFRBs) + "FRB_" + deltaT + "_" + str_fl + "Jy_"

ANT_notUSED = [ant_not_used]        ## specify which antenna has not been used for applycal
selANT = False
if ANT_notUSED != ['None']:
    selANT = True

corrANT = [item.strip() for item in corr_ant.split(',')]


print("## Not considering following antennas: ", ANT_notUSED)


antennae, ANT_notUSED_dir, ant_notUSED = getAntennae(ANT_notUSED)
mask, CORRANT, ant_toCorr, doCorrupt = getAntennae_toCorrupt(corruption, corrANT, ANT_notUSED, doCorrupt, antID)
print("mask, CORRANT, ant_toCorr, doCorrupt: ")
print(mask, CORRANT, ant_toCorr, doCorrupt)

if selANT:
    outDIR_INJ = outDIR_INJ+"notUSED_"+ANT_notUSED_dir
    if doCorrupt and mode_corr == 'manual':
        outDIR_INJ = outDIR_INJ+ANT_notUSED_dir+"notUsed_c"+CORRANT+"_corrupted_"+str(corruption)+"rad/"
    elif doCorrupt and mode_corr == 'external':
        corr_min_rad = corr_min*1.
        corr_max_rad = corr_max*1.
        if corr_parameter == 'phase':
            corr_min_rad = corr_min*0.0174533
            corr_max_rad = corr_max*0.0174533
        corruption = np.random.uniform(corr_min_rad, corr_max_rad, len(antID))
        str_corrupt = "{:.4f}".format(np.mean(corruption))
        outDIR_INJ = outDIR_INJ+ANT_notUSED_dir+"notUsed_c"+CORRANT+"_TimeIndep_"+str_corrupt+"rad/"
    elif doCorrupt and mode_corr == 'TimeDep':
        outDIR_INJ = outDIR_INJ+ANT_notUSED_dir+"notUsed_c"+CORRANT+"_TimeDep/"
    else:
        outDIR_INJ = outDIR_INJ+"/"
if doCorrupt and mode_corr == 'manual':
    outDIR_INJ = outDIR_INJ+CORRANT+"_corrupted_"+str(corruption)+"rad/"
elif doCorrupt and mode_corr == 'external':
    corr_min_rad = corr_min*1.
    corr_max_rad = corr_max*1.
    if corr_parameter == 'phase':
        corr_min_rad = corr_min*0.0174533
        corr_max_rad = corr_max*0.0174533
    corruption = np.random.uniform(corr_min_rad, corr_max_rad, len(antID))
    str_corrupt = "{:.4f}".format(np.mean(corruption))
    outDIR_INJ = outDIR_INJ+CORRANT+"_TimeIndep_"+str_corrupt+"rad/"
elif doCorrupt and mode_corr == 'TimeDep':
    outDIR_INJ = outDIR_INJ+CORRANT+"_TimeDep/"
else:
    outDIR_INJ = outDIR_INJ+"/"


print("final ANTENNAE: ", antennae)


# --------------------------------------------------------------------------------------------------------

if is_non_empty_directory(base_calib) == False:
    
    print("## "+base_calib+" is not an existing, non empty directory!")
    print("## Using all antennas instead, i.e. using calibration from All")
    antennae = '*&*'
    ANT_notUSED_dir = 'None'


thetable = base_calib+corrupt_table
expcode, corr_ant = corrupt_table.split(".")
corr_table = base_calib+expcode+"_inj."+corr_ant

print("## The corrupted table (copy of it) will be: ")
print(corr_table)

######################################################################

thesteps = []
step_title = {0: 'Import UVFITS in a CASA measurement set (MS) (importuvfits) + flag data (flagdata)',
                1: 'Divide the MS in multiple chunks at fake FRBs ToAs (split)',
                2: 'Apply calibration tables to the splitted MSs (applycal) + Clean other RFIs in CORRECTED data (rflag)',
                3: 'Create and inject a fake point source in splitted MSs',
                4: 'Concatenate the splitted MS files in a single MS (concat)',
                5: 'Get the dirty/cleaned image of the concatenated MS file (tclean)'
                }

thesteps = []

try:
    print('## List of steps to be executed: ', mysteps)
    thesteps = mysteps

except:
    print('## Global variable mystep not set')

if(thesteps==[]):

    thesteps = range(0,len(step_title))
    print('## Executing all steps: ', thesteps)


## Import UVFITS in MS

mystep = 0

if(mystep in thesteps):

    #print("rm -rf "+ms+"*")
    #os.system('rm -rf '+ms+"*")
    
    print("## Importing from .. "+uvfile)
    print("## MS file will be .. : "+ms)

    default(importuvfits)
    importuvfits(fitsfile=uvfile, vis=ms)

    default(flagdata)
    flagdata(mode='manual', vis=ms, autocorr=True)

    default(flagmanager)
    flagmanager(vis=ms, mode='save', versionname='originalMS')

    default(flagdata)
    flagdata(vis=ms, mode='manual', spw='0:0~6;58~63,1:0~6;58~63,2:0~6;58~63,3:0~6;58~63,4:0~6;58~63,5:0~6;58~63,6:0~6;58~63,7:0~6;58~63')

    default(flagmanager)
    flagmanager(vis=ms, mode='save', versionname='originalMS_v1')

    default(flagdata)
    flagdata(vis=ms, mode='list', inpfile=targetFLAG, antenna=antennae)

mystep = 1

if(mystep in thesteps):

    if listobs_path is not None and os.path.exists(listobs_path):
        print("## Considering listobs file.. ", listobs_path)

    else:
        print("## No LISTOBS path set. Creating it..")
        listobs_path = output_path+expcode+".listobs"
        listobs(vis=ms, listfile=listobs_path)
        print("listobs_path: ", listobs_path)
    
    if outDIR_INJ is not None and os.path.exists(outDIR_INJ):
        print("## Splitting MS files in directory: ", outDIR_INJ)
    
    else:
        print("## Creating directory for split.. ", str(outDIR_INJ))
        os.system("mkdir "+outDIR_INJ)

    if modeINJ != 'from_file':

        skim_listobs(listobs_path, newlistobs_path)
        scan_data = get_scans(newlistobs_path, fld_INJ)
        sel_toas = get_injected_times(outDIR_INJ, scan_data, nFRBs, modeINJ, t_start_INJ, t_stop_INJ)

    else:

        print("## Reading FRB ToAs from external file: ", file_ToAs)
        sel_toas = np.genfromtxt(file_ToAs, dtype=str)

    splitted_ms_names = []

    for ii in range(0,len(sel_toas)):

        outms = outDIR_INJ+expcode+"_INJ_"+str(ii)+".ms"
        toa = parse_datetime(sel_toas[ii])

        start_split_DT = toa - timedelta(seconds=0.5)
        stop_split_DT = toa + timedelta(seconds=0.5)
        start_split = start_split_DT.strftime('%Y/%m/%d/%H:%M:%S.%f')
        stop_split = stop_split_DT.strftime('%Y/%m/%d/%H:%M:%S.%f')
        timerange = start_split+"~"+stop_split

        print("## Timerange: ", timerange)
        default(split)
        split(vis=ms, outputvis=outms, field=fld_INJ, datacolumn='data', timerange=timerange)

        splitted_ms_names.append(outms)


mystep = 2

if(mystep in thesteps):

    for outms in splitted_ms_names:

        default(applycal)
        applycal(vis=outms,
            field='',
            antenna=antennae,
            gaintable=[sbdtab, bpasstab, mbdtab],
            interp=['nearest', 'nearest,nearest', 'linear'],
            spwmap=[[],[],[0,0,0,0,0,0,0,0]],
            gainfield=[bp_cal, bp_cal, ph_cal],
            applymode='calonly',
            parang=False)

        default(flagdata)
        flagdata(mode='rflag', vis=outms, field=fld_INJ, antenna=antennae, datacolumn='corrected', flagbackup=False, timedevscale=5.0, freqdevscale=5.0)

        #(outms, antennae, outms[:-3]+'_bkg', 'DIRTY')         ## if you want to look at the image pre-injection
        #rms = imstat(imagename=outms[:-3]+'_bkg.image', box='53, 892, 1232, 1250')['rms'][0]

        #print("## Injected burst ", outms, " rms of DIRTY image (pre-injection): ", rms)

mystep = 3

if(mystep in thesteps):

    print("## Fluxes will be: ")
    print(flux)

    if mode_corr == 'external':
        imout = outDIR_INJ+"injectedFRB_"+str_fl+"Jy_"+CORRANT+"_"+str(np.mean(corruption))+'rad_'+mode_corr     # final DIRTY image name

    for ii in range(0,len(splitted_ms_names)):

        fl = flux[ii]

        outms = splitted_ms_names[ii]

        im_name = outDIR_INJ+"TESTpoint_"+str_fl+"Jy"
        out_name = im_name+".cl"
        
        modelPointSource(outms, fl, direction, freq, BW, cell, radir, decdir, im_name, out_name)

        # if no corruption has to be applied then you want only to inject the point source into CORRECTED_DATA column
        if doCorrupt:

            if mode_corr == 'TimeDep':
                corr_min_rad = corr_min*1.
                corr_max_rad = corr_max*1.
                #corr_min_rad = corr_min*0.0174533
                #corr_max_rad = corr_max*0.0174533
                corruption = np.random.uniform(corr_min_rad, corr_max_rad, len(antID))
                imout = outDIR_INJ+"injectedFRB_"+str_fl+"Jy_"+CORRANT+"_"+str(np.mean(corruption))+'rad_'+mode_corr     # final DIRTY image name

            elif mode_corr == 'manual':
                imout = outDIR_INJ+"injectedFRB_"+str_fl+"Jy_"+CORRANT+"_"+str(corruption)+'rad'     # final DIRTY image name

            print("## Corrupting visibilities.. corruption: ", corruption, " rad")
            corruptVIS(outms, corruption, thetable, corr_table, corr_parameter, antID, mask, base_calib, mode_corr)
            #getIMAGE(outms, antennae, imout, modeIM_single, imsize)      ## get the DIRTY image of your injected-corrupted visibilities
            print("## Done!")

        else:

            imout = outms[:-3]+"_"+modeIM_single
            default(uvsub)
            uvsub(vis=outms, reverse=True)
            getIMAGE(outms, antennae, imout, modeIM_single, imsize)

mystep = 4

if(mystep in thesteps):

    concat_vis = outDIR_INJ+"CONCAT_"+str(nFRBs)+"FRB_"+deltaT+".ms"
    default(concat)
    concat(vis=splitted_ms_names, concatvis=concat_vis)

    if keep_ms != True:
        print("## Removing single MS files..")
        for ii in range(0,len(splitted_ms_names)):
            print("rm -rf "+splitted_ms_names[ii])
            #os.system("rm -rf "+splitted_ms_names[ii])
    else:
        print("## Keeping single MS files..")

mystep = 5

if(mystep in thesteps):

    print("## Ready for tclean..")
    imout = outDIR_INJ+"CONCAT_"+str(nFRBs)+"FRB_"+deltaT+"_"
    fitlog_file = outDIR_INJ+"imfit_logger.log"
    getIMAGE(concat_vis, antennae, imout+'DIRTY', 'DIRTY', imsize)      ## get the DIRTY image of your injected-corrupted visibilities
    peak = imstat(imagename=imout+"DIRTY.image", box=box_source)['max'][0]
    peak_dirty = "{:.4f}".format(peak)
    rms_used = peak/10.
    rmsUsed_str = "{:.4f}".format(rms_used)
    getIMAGE(concat_vis, antennae, imout+"CLEAN", "CLEAN", imsize, rms_used)      ## get the CLEANED image of your injected-corrupted visibilities
    RAfit, decfit = extractFIT(imout+"CLEAN.image", box_source, fitlog_file)
    decfit_form = convert_dec(decfit)
    peak_fromIm = imstat(imagename=imout+"CLEAN.image", box=box_source)['max'][0]
    peakstr_fromIm = "{:.4f}".format(peak_fromIm)
    rms_fromIm = imstat(imagename=imout+"CLEAN.image", box=box_bkg)['rms'][0]
    rms_str = "{:.4f}".format(rms_fromIm)
    print("## Writing fit to file: ", outputFIT_path)
    outFITfile = open(outputFIT_path, 'a')
    outFITfile.write(str(nFRBs)+" "+deltaT+" "+peak_dirty+" "+peakstr_fromIm+" "+rmsUsed_str+" "+rms_str+" "+RAfit+" "+decfit_form+"\n")
    outFITfile.close()
    print('rm -rf '+base_EXP+str(nFRBs)+"*")
    #os.system('rm -rf '+base_EXP+str(nFRBs)+"*")
