import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import ascii
from matplotlib.patches import Rectangle
import numpy as np
import os
import seaborn as sns
import glob
import scipy.optimize as opt
import matplotlib.image as mpimg
import pandas as pd


#####    UNCOMMENT TO HAVE USER INPUT FILES AND PATH
#print('enter path to imm file')
#path = input()
#print('enter imm file name')
#filename = input()
#print('enter .batch file name')
#batchfile = input()


### HARD INPUT imm, batch and file path        JUST CHANGE PATH TO YOUR FOLDER WITH BATCH AND .IMM, DO NOT NEED .npy FILES
#path = "/home/avarga/Documents/RIT_HW/Data_Analy/multiple images in single file with batchinfo file/"
path = "/home/avarga/Documents/RIT_HW/Data_Analy/immtonpy/"
#filename = "au111-sunday_naf20_XCV_Lp25_3_00001-16010.imm"
##############
#filename = "ag001_NaF_p02M_MonLong1_00001-170020.imm"
#268
#150
#10
###########
#filename = "fe3al_long_04_00001-02070.imm"
#396
#477
#10
##########
#filename = "ErMnO3_S178_room_temp_00001-01210.imm"
#97
#107
#10
########
#filename = "au_001_batch2_00001-15003.imm"
#225
#76
#10

#########
filename = "au111-sunday_naf20_bp0_Lp25_00001-16010.imm"
#
#
#
#####
#filename = "au111-sunday_naf20_XCV_Lp25_3_00001-16010.imm"
#358
#54
#10
##########
#filename = "ErMnO3_S100_870C_150A_fast_00001-02570.imm"

#batchfile = "ErMnO3_S100_870C_150A_fast_0001-0860.batchinfo"
#batchfile = "au111-sunday_naf20_XCV_Lp25_3_0001-2198.batchinfo"
#batchfile = "au_001_batch2_0001-15003.batchinfo"
#batchfile = "ErMnO3_S178_room_temp_0001-1210.batchinfo"
#batchfile = "/home/avarga/Documents/RIT_HW/Data_Analy/immtonpy/ag001_NaF_p02M_MonLong1_0001-26158.batchinfo"
#batchfile = "/home/avarga/Documents/RIT_HW/Data_Analy/immtonpy/fe3al_long_04_0001-2070.batchinfo"
batchfile = "au111-sunday_naf20_bp0_Lp25_0001-3103.batchinfo"
fid = open(path+filename, mode = 'rb')
header_files = batchfile
##########################################################
#Attila - moved some of Sadie's teams code to handle reading batch file lines


with open(header_files, "r") as f:
    header_info = [line.strip() for line in f.readlines()]

# read all the lines in the .batchinfo file and pick out the "ndark0" and "ndarkend" parameters
for line in header_info:
    if line.startswith("ndark0"):
        ndark0 = int(line.split(" = ")[1])
    if line.startswith("ndataend"):
        nframes = int(line.split(" = ")[1])
    if line.startswith("cols"):
        ccols = int(line.split(" = ")[1])
    if line.startswith("rows"):
        rrows = int(line.split(" = ")[1])
    if line.startswith("ndata0"):
        ndata0 = int(line.split(" = ")[1])
    if line.startswith("ndataend"):
        ndataend = int(line.split(" = ")[1])
    elif line.startswith("ndarkend"):
        ndarkend = int(line.split(" = ")[1])
#print(ndata0)
#print(ndataend)

###########################################################################
#Attila- function immtotxt inputs total frame number and file as fid
#Takes .imm file and stores it as an array called imm, this is the return value


def immtotxt(fid,frame):
    # Finding first and last image index number from filename
    findUnderscore = filename.split('_')
    findDash = filename.split('-')
    firstImmIndex = int(findUnderscore[-1].split('-')[0])
    lastImmIndex = int(findDash[-1].split('.')[0])
    ImmIndexList = list(range(0,lastImmIndex+1))
    #immIndex = ImmIndexList[12]
    immIndex = frame
    fid.seek(0,0)
    modeflag = np.fromfile(fid, np.int32,count = 1)
    fid.seek(4,0)
    compressionflag = np.fromfile(fid, np.int32,count = 1)
    fid.seek(616,0)
    immversionflag  = np.fromfile(fid, np.int32,count = 1)

# Determine start position of the wanted image
    fid.seek(108,0)
    rows = np.fromfile(fid, np.int32, count=1)
    fid.seek(112,0)
    cols = np.fromfile(fid, np.int32, count=1)
    fid.seek(116,0)
    bytes = np.fromfile(fid, np.int32, count=1)

# Finding image size and image start
    compression = 0
    if(compression == 0):
        if(immversionflag>=11):
            immSize = bytes * rows * cols +1024
        else:
            immSize = 2 * rows * cols +1024
# Converting from array to integer numbers
    imageStart_array = (immIndex - 1)*immSize
    imageStart = int(imageStart_array)

# Load header of given image index
    fid.seek(imageStart,0)
    header = np.empty(shape=(53,2),dtype = object)
# LF - 11/16/2022 --------------------------------------------------------------
# MatLab ->  Python
# int    ->  np.int32
# uint32 ->  np.uint32
# uchar  ->  np.ubyte
# float  ->  np.single
# double ->  np.double
    header[0] = ('mode',        np.fromfile(fid, np.int32, count=1))
    header[1] = ('compression', np.fromfile(fid, np.int32, count=1))
    header[2] = ('date',        np.fromfile(fid, np.ubyte, count=32))
    header[3] = ('prefix',      np.fromfile(fid, np.ubyte, count=16))
    header[4] = ('number',      np.fromfile(fid, np.int32, count=1))
    header[5] = ('suffix',      np.fromfile(fid, np.ubyte, count=16))
    header[6] = ('monitor',     np.fromfile(fid, np.int32, count=1))
    header[7] = ('shutter',     np.fromfile(fid, np.int32, count=1))
    header[8] = ('row_beg',     np.fromfile(fid, np.int32, count=1))
    header[9] = ('row_end',     np.fromfile(fid, np.int32, count=1)-1)
    header[10] = ('col_beg',    np.fromfile(fid, np.int32, count=1))
    header[11] = ('col_end',    np.fromfile(fid, np.int32, count=1)-1)
    header[12] = ('row_bin',    np.fromfile(fid, np.int32, count=1))
    header[13] = ('col_bin',    np.fromfile(fid, np.int32, count=1))
    header[14] = ('rows',       np.fromfile(fid, np.int32, count=1))
    header[15] = ('cols',       np.fromfile(fid, np.int32, count=1))
    header[16] = ('bytes',      np.fromfile(fid, np.int32, count=1))
    header[17] = ('kinetics',   np.fromfile(fid, np.int32, count=1))
    header[18] = ('kinwinsize', np.fromfile(fid, np.int32, count=1))
    header[19] = ('elapsed',    np.fromfile(fid, np.double, count=1))
    header[20] = ('preset',     np.fromfile(fid, np.double, count=1))
    header[21] = ('topup',      np.fromfile(fid, np.int32, count=1))
    header[22] = ('inject',     np.fromfile(fid, np.int32, count=1))
    header[23] = ('dlen',       np.fromfile(fid, np.int32, count=1))
    header[24] = ('roi_number', np.fromfile(fid, np.int32, count=1))
    header[25] = ('buffer_number', np.fromfile(fid, np.uint32, count=1))
    header[26] = ('systick',       np.fromfile(fid, np.uint32, count=1))
    header[27] = ('pv1',       np.fromfile(fid, np.ubyte, count=40))
    header[28] = ('pv1VAL',    np.fromfile(fid, np.single, count=1))
    header[29] = ('pv2',       np.fromfile(fid, np.ubyte, count=40))
    header[30] = ('pv2VAL',    np.fromfile(fid, np.single, count=1))
    header[31] = ('pv3',       np.fromfile(fid, np.ubyte, count=40))
    header[32] = ('pv3VAL',    np.fromfile(fid, np.single, count=1))
    header[33] = ('pv4',       np.fromfile(fid, np.ubyte, count=40))
    header[34] = ('pv4VAL',    np.fromfile(fid, np.single, count=1))
    header[35] = ('pv5',       np.fromfile(fid, np.ubyte, count=40))
    header[36] = ('pv5VAL',    np.fromfile(fid, np.single, count=1))
    header[37] = ('pv6',       np.fromfile(fid, np.ubyte, count=40))
    header[38] = ('pv6VAL',    np.fromfile(fid, np.single, count=1))
    header[39] = ('pv7',       np.fromfile(fid, np.ubyte, count=40))
    header[40] = ('pv7VAL',    np.fromfile(fid, np.single, count=1))
    header[41] = ('pv8',       np.fromfile(fid, np.ubyte, count=40))
    header[42] = ('pv8VAL',    np.fromfile(fid, np.single, count=1))
    header[43] = ('pv9',       np.fromfile(fid, np.ubyte, count=40))
    header[44] = ('pv9VAL',    np.fromfile(fid, np.single, count=1))
    header[45] = ('pv10',      np.fromfile(fid, np.ubyte, count=40))
    header[46] = ('pv10VAL',   np.fromfile(fid, np.single, count=1))
    header[47] = ('imageserver', np.fromfile(fid, np.single, count=1))
    header[48] = ('CPUspeed',    np.fromfile(fid, np.single, count=1))
    header[49] = ('immversion',  np.fromfile(fid, np.int32, count=1))
    header[50] = ('corecotick',  np.fromfile(fid, np.uint32, count=1))
    header[51] = ('cameratype',  np.fromfile(fid, np.int32, count=1))
    header[52] = ('threshhold',  np.fromfile(fid, np.single, count=1))

# Building image
    if(compression == 0):
        fid.seek(imageStart+1024, 0)
        if(immversionflag >= 11):
            if (bytes == 2 ):
                imm = np.fromfile(fid, np.short, int((immSize-1024)/2))
            elif (bytes == 4):
                imm = np.fromfile(fid, np.short, int((immSize-1024)/2))
        else:
            imm = np.fromfile(fid, np.short, int((immSize-1024)/2))

    return imm
#####################################
#Attila - New method to take imm array and output to an array of only dark frames and of non dark frames

f_array = []
frame_array = []
dark_frames = []
#3000 for ag001
for n in range(ndark0, ndarkend+1):
    next_dframe = immtotxt(fid,n)
    dark_frames.append(next_dframe)
for j in range(ndark0, ndarkend+2000):
    nextf = immtotxt(fid,j)
    f_array.append(nextf)
#for k in range(ndata0,ndataend+1):
#    next_frame = immtotxt(fid,k)
#    frame_array.append(next_frame)
darks = np.reshape(dark_frames, (len(dark_frames),rrows,ccols))
lights = np.reshape(f_array, (len(f_array),rrows,ccols))
#print(np.shape(lights))
if (ndarkend>11):
   dnum = 21
elif (ndarkend<11):
    dnum = 11
print(ndark0,dnum,'here')
lights = lights[dnum::1]
#print(np.shape(lights))
averaged_dark = np.average(darks,axis=0)

clean = lights - averaged_dark
del(darks)
del(lights)
del(averaged_dark)
#######################################################################################################################
#Code from Sadie's team

#**************************************

#Attila - Commented out as I had to move stuff to make code work for arrays and not .npy files
#The method is still the same


#**************************************
'''
input_dir = "au111-sunday_naf20_bp0_Lp25_00001-16010/"
input_dir = "All_images/"
output_dir = "clean/"

# First I downloaded the .batchinfo file and and manually changed the extension to .txt - the location of the darks are indicated in the .batchinfo file

Open the batchinfo (header) file
header_files = glob.glob("*3103.txt")

with open(header_files[0], "r") as f:
    header_info = [line.strip() for line in f.readlines()]

# read all the lines in the .batchinfo file and pick out the "ndark0" and "ndarkend" parameters
for line in header_info:
    if line.startswith("ndark0"):
        ndark0 = int(line.split(" = ")[1])
    elif line.startswith("ndarkend"):
        ndarkend = int(line.split(" = ")[1])

# pull out only the dark frames using the parameters from the header
dark_files = [input_dir + '/au111-sunday_naf20_bp0_Lp25_00001-16010_'
              + str(n) + '.npy' for n in range(ndark0, ndarkend+1)]


# generate blank array in which to put the new averaged dark frame
sum_array = np.zeros((650, 100))
sum_array =+ dark_files

averaged_dark = np.divide(sum_array, len(dark_frames) * np.ones((650, 100)))
averaged_dark = np.average(dark_array)
# generate a heat map for the averaged dark
plt.figure(dpi=300)
sns.heatmap(averaged_dark, cmap='hot', square='True', xticklabels=False, yticklabels=False)
plt.show()


# subtract the averaged dark from every other frame and save cleaned datafiles to directory "clean" with new extension "_clean"

clean = frame_array - averaged_dark
'''

#####################################################################################################################################
#Attila - Removed load_data() as code just uses arrays
'''
def load_data():
    files = os.listdir("clean")    #Go to the clean folder
    file_count = len(files)-2 # CHANGE TO 2 LATER!!!!!!!!!!!!!!!
    count = 0
    start = 1
    x = 0
    y = 0
    r = 0

    DATA = [None] * file_count

    for file in os.listdir("clean"):    #go into clean directory
        filename = os.fsdecode(file)
        filename_split = filename.split('.')
        if filename_split[1] == 'npy':
            data = np.load("clean/"+f'{filename}')
            frame = filename_split[0].split('_')[-2]     #choose the number of data
            if count == 0:
                start = frame
            #print(data)
            #print(frame)

            DATA[int(frame)-int(start)] = data
            count += 1
        if filename == 'input.txt':
            region_data = np.loadtxt("clean/"+f'{filename}')
            x = int(region_data[0])
            y = int(region_data[1])
            r = int(region_data[2])
            region = (x,y,r)

    #print(f"Data loaded, {count} files")
    return region, DATA[10:]
'''
########################################################################################

#******************************************
#Attila - these functions are not changed except main() see below

#**********************************************
def load_data():
    files = os.listdir("clean")    #Go to the clean folder
    file_count = len(files)-2 # CHANGE TO 2 LATER!!!!!!!!!!!!!!!
    count = 0
    start = 1
    x = 0
    y = 0
    r = 0

    DATA = [None] * file_count

    for file in os.listdir("clean"):    #go into clean directory
        filename = os.fsdecode(file)
        filename_split = filename.split('.')
        if filename_split[1] == 'npy':
            data = np.load("clean/"+f'{filename}')
            frame = filename_split[0].split('_')[-2]     #choose the number of data
            if count == 0:
                start = frame
            #print(data)
            #print(frame)

            DATA[int(frame)-int(start)] = data
            count += 1
        if filename == 'input.txt':
            region_data = np.loadtxt("clean/"+f'{filename}')
            x = int(region_data[0])
            y = int(region_data[1])
            r = int(region_data[2])
            region = (x,y,r)

    print(f"Data loaded, {count} files")
    return region, DATA[10:]

def correlation(data, points):
    #print("points:", len(points))
    correl = []
    for point in points:
        data_points = []
        for frame in data:
            data_points.append(frame[point[1]][point[0]])
        data_points = np.array(data_points)
        numerator = cor_num(data_points)
        denominator = np.average(data_points) ** 2
        g2 = numerator / denominator
        correl.append(g2)
    correl = np.array(correl)
    correl_avg = np.average(correl, axis=0)
    print('hi')
    print(correl)
    print(len(correl[0]))
    print(len(correl))
    print(correl_avg)

    print(len(correl_avg)==len(correl[5]))

    #Add uncertainty here
    N = len(points)    #use the number of points in the region for now

    var_g2_array = []
    for a in range(len(correl[0])):   #loop through all data at one point
        diff_sq_sum = 0
        for b in range(len(correl)):   #loop through all the points
            diff_sq_sum += (correl[b][a]-correl_avg[a])**2

        var_g2_array.append(diff_sq_sum/(N*(N-1)))

    var_g2_array = np.asarray(var_g2_array)
    g2_err= np.sqrt(var_g2_array)
    print("done")
    print(var_g2_array[:12])
    print(g2_err[:12])

    return (correl_avg, g2_err)


def cor_num(data):
    out = []
    delta_t = range(1,len(data))
    for dt in delta_t:
        count = 0
        total = 0
        for t in range(len(data)):
            if t+dt < len(data):
                total += data[t]*data[t+dt]
                count += 1
        avg = total / count
        out.append(avg)
    return np.array(out)


def find_region(x,y,r,frame):
    print(frame)
    points = []

    val_sum = 0
    count = 1
    points.append([x,y])
    for i in range(1,r+1):
        for j in range(1,r+1):
            if np.sqrt( i**2 + j**2 ) <= r and x+i<len(frame[0]) and y+j<len(frame) and x-i>=0 and y-j>=0:
                val_sum += frame[y+j][x+i]
                count += 1
                points.append([x,y+j])
                points.append([x,y-j])
                points.append([x+i,y])
                points.append([x-i,y])
                points.append([x+i,y+j])
                points.append([x+i,y-j])
                points.append([x-i,y+j])
                points.append([x-i,y-j])

    val_avg = val_sum / count

    fig, ax = plt.subplots(figsize=(32,5))
    sns.heatmap(frame, cmap='hot', square='True', xticklabels=False)
    circle = plt.Circle((x,y), r, color='g', fill=False)
    ax.add_patch(circle)
    plt.title(f"Region: ({x}, {y}) r = {r}\nAverage Value : {np.round(val_avg,3)}")
    fig.savefig(f'{filename}'+'region.png')
    plt.show()

    print(np.array(points))
    return np.array(points)



def fitting(data):
    x = np.linspace(1, len(data[0]), len(data[0]))
    y = data[0][::-1]
    y_err = data[1][::-1]


    def func(x,A,gamma):
        return 1+A*np.exp(-2*gamma*x)

    parameters, covariance = opt.curve_fit(func, x, y, sigma = y_err)
    #print(covariance)
    perr = np.sqrt(np.diag(covariance))

    A = parameters[0]
    gamma = parameters[1]

    print("Parameters")
    print('A = '+str(A))
    print('Gamma = '+str(gamma))

    print("error")
    print(perr)
    plt.figure(figsize=(32,5))
    plt.errorbar(x,y, yerr= y_err, marker='.',markersize = 0.5, linestyle='None', capsize=3, label='Data',ecolor='black',elinewidth=0.1)
    plt.plot(x, func(x,*parameters), label='Fit', color='tab:orange')
    plt.title(f'A = {np.round(A,3)} $\pm$ {np.round(perr[0],3)}, $\Gamma$ = ({np.round(10000*gamma,1)} $\pm$ {np.round(10000*perr[1],1)}) x $10^{-2}$')
    plt.ylabel('$g_2$')
    plt.xlabel('$\Delta t$')
    plt.legend()
    plt.savefig(f'{filename}'+'FIT.png')
    plt.show()


def analysis(region,data):
    points = find_region(*region,data[0])
    corr = correlation(data, points)
    plot_correlations(*region,corr)
    fitting(corr)


def plot_correlations(x,y,r,corr):
    plt.figure(figsize=(32,5))
    t = np.linspace(1, len(corr[0]), len(corr[0]))
    #plt.scatter(t, corr[0][::-1], marker='.')
    plt.errorbar(t, corr[0][::-1], yerr= corr[1][::-1], marker='.',markersize = 0.5,ecolor='black',elinewidth=0.1)
    plt.title(f"Region: ({x}, {y}) r = {r}")
    plt.ylabel("Auto Correlation, $g_2$")
    plt.xlabel("$\Delta t$")
    plt.savefig(f'{filename}'+'g2.png')
    plt.show()


def main():
    data = clean
    fig, ax = plt.subplots(figsize=(32,5))
    ax.imshow(clean[0], cmap='hot')
    plt.title("First frame in file")
    #ax.xaxis.set_ticks(np.arange(0, len(data[0][0]), 50))
    #ax.yaxis.set_ticks(np.arange(0, len(data[0]), 50))
    plt.show()
    print('enter region x coordinate')
    x = int(input())
    print('enter region y coordinate')
    y = int(input())
    print('enter region r radius')
    r = int(input())
    region = (x,y,r)
    analysis(region, data)
    #fitting(data)
main()











