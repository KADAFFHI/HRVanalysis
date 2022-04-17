import os
import sys
import pandas as pd
import numpy as num
import matplotlib.pyplot as mplot
import math
import heartpy as hp
from scipy.signal import resample
import ntpath as pt
from scipy.stats.stats import pearsonr
from scipy import stats
from multiprocessing import Pool

# -----INSTRUCTION-----#

#   call program from terminal
#   use desired folder directory as argument
#   the program reads all csv files in the directory
#   program presents 2 choices: Process all or process specific. choose by typing in the displayed command
#   Process all: create csv files of RR values and txt files of additional measurements for all files in selected folder
#   Process specific: choose desired file by typing in the index,
#     program prints a poincare plot, EKG and other values over time

# ex: python file/location/RRmassAnalysis.py path/to/folder

def getDataMatrices(dirPath): # function read all files in directory, returns content of files as a matrix
    for root, dirs, files in os.walk(dirPath):
        return [(pd.read_csv(os.path.join(dirPath,file)),file) for file in files if os.path.splitext(file)[1] == ".csv"]

def histRR(RRldiff): # function to plot a histogram of the allotted data
    # ---REQUIRES FURTHER WORK --- #

    mplot.figure()
    size, scale = 1000, 10

    mplot.hist(RRldiff, color='#607c8e')
    mplot.title('Differance in RR-values')
    mplot.xlabel('Counts')

    mplot.ylabel('Differance')
    mplot.grid(axis='x', alpha=0.5)
    mplot.show()

    # --------------------------- #


def heartRate(wd, m): # function to plot EKG data over time and poincare plots using HeartPy
    hp.plotter(wd, m, figsize=(12, 4))
    hp.plot_poincare(wd, m)


def correlation(RRlist): # function to compute pearson correlation between RR values
    RR, RRp1 = [], []

    for idx, i in enumerate(RRlist):
        if idx % 2 == 0:
            RRp1.append(i)
        else:
            RR.append(i)

    if len(RRlist) % 2 != 0:
        RRp1.remove(RRp1[-1])

    #print('Pearson correlation: ', pearsonr(RR, RRp1))

    return pearsonr(RR, RRp1)


def HRplotRR(RRlist, RRindex, samplerate): # function to plot RR values over time
    import math
    """
    temp = [(100000 / r, idx) for idx, r in enumerate(RRlist) if r > 400]
    HRpBeat, skipIdx = [i[0] for i in temp], [i[1] for i in temp]
    idx = [RRindex[i][0] / 60000 for i in skipIdx]
    """

    idx = [index/130 for index in range(len(RRindex))]
    mplot.figure(figsize=(22, 11))

    max_Y, min_Y, max_X, min_X = math.ceil(max(RRlist)), math.floor(min(RRlist)), math.ceil(max(idx)), math.floor(
        min(idx))


    mplot.plot(idx, RRlist, 'b', label="RR values by heartbeat", linewidth=0.5)

    mplot.yticks(num.arange(min_Y, max_Y, 25))
    #mplot.xticks(num.arange(min_X, max_X))

    mplot.xlabel('Beats')
    mplot.ylabel('RR value')
    mplot.title(f'RR values over time')
    mplot.legend()
    mplot.grid()
    mplot.show()


def HRplotHP(dHR, file): # function to plot heartrate over time
    background_constant = ""
    idx = []
    for i in range(len(dHR)):
        idx.append(i * 5)


    mplot.figure(figsize=(22, 11))
    max_Y, min_Y, max_X, min_X = math.ceil(max(dHR)), math.floor(min(dHR)), math.ceil(max(idx)), math.floor(min(idx))

    mplot.plot(idx, dHR, 'g', label="HR/t", linewidth=1)
    mplot.yticks(num.arange(min_Y, max_Y, 5))
    mplot.xticks(num.arange(min_X, max_X, 10), rotation=45)
    mplot.xlabel('Time [s]')
    mplot.ylabel('HeartRate [BPM]')
    mplot.title(f'HeartRate over time: {file}')
    mplot.legend()
    mplot.grid()
    mplot.show()

def statistics(dataVector): # function returns [mean, stdiv]
    temp_mean = 0
    stdev  = stats.tstd(dataVector)
    for r in dataVector:
        temp_mean += r
    mean = temp_mean/len(dataVector)

    return [mean, stdev]

def RRanalysis(dataMatrix, sampleRate, target_file_path, mode):
    # depending on the "mode" variable the function behaves in 2 ways:
    # mode = 0: the function reads the allotted data and returns a .csv file of RR values and a .txt file
    # -with additional measurements.
    # mode = 1: the function reads the allotted data and plots it to the screen, allowing further inspection

    print("-- processing data")

    data = [dataMatrix[column] for column in dataMatrix.columns]
    df = data[1]
    

    #raw_data = [datapoint for datapoint in df[2]]

    ## filter and resample ECG data to fix peak-rejection issue
    filtered_data = hp.filter_signal(df, cutoff=0.05, sample_rate=sampleRate, filtertype='notch')
    # filtered_data = hp.filter_signal(filtered_data, cutoff=[1.11, 3.33], sample_rate=sampleRate, filtertype='bandpass')
    resampled_data = resample(filtered_data, len(filtered_data) * 2)
    ##---------------------------------------------------------

    savefileID = f'{(os.path.splitext(target_file_path)[0])}.RRdata.csv'  # filename for the finished file

    wdr, mr = hp.process(df, sampleRate)  # raw ECG data pre-filter
    wd, m = hp.process(resampled_data, sampleRate * 2)  # filtered and resampled ECG data
    wds, ms = hp.process_segmentwise(resampled_data, sampleRate * 2, segment_width=10, segment_overlap=0.5)
    RRlist = wd['RR_list']  # list of RR intervals
    RRdiff = wd['RR_diff']  # list of RR differences
    RRindex = wd['RR_indices']  # list of RR indices
    dHR = ms['bpm']

    # ------------------------------

    # aditional calculations ----------
    meanRR = statistics(RRlist)
    meanRR = f"{meanRR[0]} " + u"\u00B1" + f" {meanRR[1]}"
    meanHR = statistics(dHR)
    meanHR = f"{meanHR[0]} " + u"\u00B1" + f" {meanHR[1]}"
    pearsonCorrelation = correlation(RRlist)



    if mode == 0:
        print("-- storing data")

        # creating csv file with RR values
        with open(f"{savefileID}", "w") as f:
            f.write(f"\"RR values\",\"RR differances\",\"RR indices\"\n")
            for line in range(len(RRindex) - 1):
                try:
                    f.write(f"\"{RRlist[line]}\",\"{RRdiff[line]}\",\"{RRindex[line]}\"\n")

                except Exception as ex:
                    print(ex)

        # ----------------------------------

        f.close()

        # creating txt file with aditional data and measurements
        with open(f"{savefileID}_Measures.txt", "w") as f:
            f.write(f"Measures: ")
            for measure in m.keys():
                try:
                    f.write(f'\n{measure}: {m[measure]:f}')

                except Exception as ex:
                    print(ex)
            try:
                f.write(f"\nmax bpm: {max(dHR)}")
                f.write(f"\nmin bpm: {min(dHR)}")
                f.write(f"\nmean bpm: {meanHR}")
                f.write(f"\nmean RR: {meanRR}")
                f.write(f"\n\npearson correlation: {pearsonCorrelation}")
            except Exception as ex:
                print(ex)

            # ----------------------------------

        f.close()

    if mode == 1:
        # plotting values to screen
        heartRate(wd,m)
        HRplotRR(RRlist, RRindex, sampleRate)
        HRplotHP(dHR, savefileID)




def main(argv):
    sampleRate = 130

    print("\n\n-------RR_analysis v_2-------")


    if len(argv) > 0:
        file_directories = [arg for arg in argv if os.path.exists(arg)]
    else:
        print("Pass a directory as an argument")
        exit()

    if len(file_directories) == 0:
        print("No valid directories found.")
        exit()

    # multithreading to allow quick filescan

    if os.cpu_count() <= 2:
        threadCount = 1
    elif 2 < os.cpu_count() <= 4:
        threadCount = 2
    else:
        threadCount = os.cpu_count() - 2

    with Pool(processes=threadCount) as processPool:
        tskIoOps = processPool.map_async(getDataMatrices, file_directories)
        processPool.close()
        processPool.join()

    dataFileSets = tskIoOps.get()

    dataCollection = []


    for k in range(len(dataFileSets)):
        dataFileSet = dataFileSets[k]
        dataMatrices, fileNames = [dataFileSet[i][0] for i in range(len(dataFileSet))], [dataFileSet[i][1] for i in
                                                                                         range(len(dataFileSet))]


    plotRangeIndex = [i for i in range(len(dataMatrices))]



    # ----------- USER INTERFACE -------------------------------------------
    while True:
        try:
            ModeSelector = input(
                f"Process all[ALL], Process specific[SPEC] or exit[EXIT]: ")


            if ModeSelector == 'ALL':

                for matrix in range(len(dataMatrices)):

                    RRanalysis(dataMatrices[matrix], sampleRate, fileNames[matrix], 0)

            elif ModeSelector == "SPEC":

                print("choose desired file: \n")

                for file in range(len(fileNames)):
                    print("index: ", file, ", name: ", fileNames[file])

                try:
                    command = int(input("enter index of desired file: "))
                except:
                    command = -1

                if 0 <= command < len(fileNames):
                    RRanalysis(dataMatrices[command], sampleRate, fileNames[command],1)



            elif ModeSelector == 'EXIT':
                exit()

            else:
                print("Invalid input. Retry!")

        except Exception as ex:
            print(ex)
    #---------------------------------------------------------------------------------

if __name__ == '__main__':
    main([sys.argv[i] for i in range(1, len(sys.argv)) if len(sys.argv) > 1])
