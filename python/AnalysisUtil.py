import tkinter.filedialog

import imageio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from os import listdir
from os.path import isfile, join
import tkinter
from tkinter import *
from scipy.ndimage import convolve, uniform_filter

class ButtonClickedEvent:
    @staticmethod
    def Open_Path(fd):

        filepath = tkinter.filedialog.askdirectory(initialdir=f"{fd}/")
        return filepath

    @staticmethod
    def Open_File(fd):

        filepath = tkinter.filedialog.askopenfilename(initialdir=f"{fd}/")
        return filepath

    @staticmethod
    def Read_Image(filepath, fileformat, filedtype, ImageSize):

        read_data = np.array([], dtype=np.float64)
        onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]

        for file_now in onlyfiles:
            if file_now[-3:] == fileformat:

                fid = open(f"{filepath}/{file_now}", "rb")
                read_data_now = np.fromfile(fid, dtype=filedtype, sep="")
                read_data_now = read_data_now.reshape(ImageSize)
                fid.close()

                if not read_data.any():
                    read_data = read_data_now[np.newaxis, :]
                    continue

                read_data = np.append(read_data, read_data_now[np.newaxis, :], axis=0)

        return np.array(read_data, dtype=np.float64)

    @staticmethod
    def Read_File(filepath, fileformat, filedtype, ImageSize):

        read_data = np.array([], dtype=np.float64)
        file_now = filepath

        if file_now[-3:] == fileformat:

            fid = open(f"{file_now}", "rb")
            read_data = np.fromfile(fid, dtype=filedtype, sep="")
            read_data = read_data.reshape(ImageSize)
            fid.close()

        return np.array(read_data, dtype=np.float64)

    @staticmethod
    def Save_File(filepath, filedtype, data):

        fn = filepath[filepath.rfind('/')+1:filepath.rfind('.')]

        filepath = tkinter.filedialog.asksaveasfilename(initialdir=f"{filepath}/",
                                                        initialfile = f'Line Calibrated {fn} W{data.shape[1]}H{data.shape[0]}',
                                                        title="Save as",
                                                        defaultextension=".raw",
                                                        filetypes=(("raw", ".raw"),
                                                                   ("tif", ".tiff"),
                                                                   ("all files", "*")))

        if filepath[-3:] == "raw":
            with open(filepath, 'wb') as f:
                f.write((data.astype(filedtype)).tobytes())
                f.close()
        elif filepath[-4:] == "tiff":
            imageio.imwrite(filepath, data.astype(filedtype))

    @staticmethod
    def Set_ROI(Data, ROI1, ROI2):
        X1, X2, Y1, Y2 = ROI1[0], ROI2[0], ROI1[1], ROI2[1]

        L, R = min(X1, X2), max(X1, X2)
        U, D = min(Y1, Y2), max(Y1, Y2)

        if Data.ndim == 2:
            return Data[U:D+1, L:R+1]

        if Data.ndim == 3:
            return Data[:, U:D + 1, L:R + 1]

    @staticmethod
    def Set_FOI(Data, FOI):
    #"""FOI: Frame of Interest"""

        return Data[FOI[0]-1:FOI[1]]

    @staticmethod
    def Division(ax, imageinfo, numRow, numCol):
        Plotting.DrawDivision(ax, imageinfo, numRow, numCol)

    @staticmethod
    def Calc(imageinfo, Differential = True):
        if Differential:
            imageinfo = DataProcessing.DifferentialImage(imageinfo)
        FrameNoise = DataProcessing.FrameNoise(data = imageinfo, Differential = Differential)
        TotalNoise = DataProcessing.TotalNoise(data = imageinfo, Differential = Differential)
        RowLineNoise = DataProcessing.LineNoise(data = imageinfo, Differential = Differential, Orientation = 'Row')
        ColLineNoise = DataProcessing.LineNoise(data = imageinfo, Differential = Differential, Orientation = 'Col')
        LineNoise = np.sqrt(np.square(RowLineNoise) + np.square(ColLineNoise))
        PixelNoise = DataProcessing.PixelNoise(TotalNoise, FrameNoise, LineNoise)
        return {"TotalNoise" : TotalNoise, "FrameNoise" : FrameNoise, "RowLineNoise" : RowLineNoise,
                "ColLineNoise" : ColLineNoise, "PixelNoise" : PixelNoise, "ImageInfo": imageinfo}

    @staticmethod
    def IQR_Method(imageinfo, NIQR, NIteration, Differential = True, ExcludingZero = True, HPF = True):
        if Differential:
            imageinfo = DataProcessing.DifferentialImage(imageinfo)

        if HPF:
            imageinfo = DataProcessing.Highpass_Filter(imageinfo)

        Mask = (imageinfo.copy()).astype(bool)
        Mask.fill(True)

        if ExcludingZero:
            Mask = (imageinfo!=0)

        for k in range(NIteration):
            MaskedArray = np.ma.masked_array(imageinfo, mask=~Mask)[np.ma.masked_array(imageinfo, mask=~Mask).mask==False].data
            Mask = Mask * DataProcessing.IQR_Mask(imageinfo = imageinfo, MaskedArray = MaskedArray, NIQR = NIQR)

        MaskedImage = np.ma.masked_array(data=imageinfo, mask=~Mask)

        FrameNoise = DataProcessing.FrameNoise(data = MaskedImage, Differential = Differential)
        TotalNoise = DataProcessing.TotalNoise(data = MaskedImage, Differential = Differential)
        RowLineNoise = DataProcessing.LineNoise(data = MaskedImage, Differential = Differential, Orientation = 'Row')
        ColLineNoise = DataProcessing.LineNoise(data = MaskedImage, Differential = Differential, Orientation = 'Col')
        LineNoise = np.sqrt(np.square(RowLineNoise) + np.square(ColLineNoise))
        PixelNoise = DataProcessing.PixelNoise(TotalNoise, FrameNoise, LineNoise)
        return {"TotalNoise" : TotalNoise, "FrameNoise" : FrameNoise, "RowLineNoise" : RowLineNoise,
                "ColLineNoise" : ColLineNoise, "PixelNoise" : PixelNoise, "ImageInfo" : imageinfo, "Mask": Mask}

    @staticmethod
    def DeNoise(imageInfo, nG, nD):
        r = int(imageInfo.shape[0] / nG)
        c = int(imageInfo.shape[1] / nD)

        for j in range(r):
            for k in range(c):
                tempData = DataProcessing.SelectBlock(imageInfo, nG*j, nG*(j+1), nD*k, nD*(k+1))
                imageInfo[nG*j:nG*(j+1), nD*k:nD*(k+1)] = DataProcessing.LineCalibration(tempData)

        return imageInfo


class DataProcessing:
    @staticmethod
    def TemporalAverage(data):
        # Return 2D Image averaged over Time
        return data.mean(axis=0) if data.ndim == 3 else data

    @staticmethod
    def SpatialAverage(data):
        # Return 1D Temporal Array Averaged over All Pixels
        return data.mean(axis=(1,2)) if data.ndim == 3 else data.mean()

    @staticmethod
    def DN2Coulomb(DN, IFS=2, LSBC = 0.125E-12, LSBV = 61.035E-6):
        return LSBC*(IFS+1)*LSBV*DN

    @staticmethod
    def Coulomb2Electron(Coulomb):
        return Coulomb / 1.602E-19

    @staticmethod
    def DN2Electron(DN, IFS=2, LSBC = 0.125E-12, LSBV = 61.035E-6):
        return DataProcessing.Coulomb2Electron(DataProcessing.DN2Coulomb(DN, IFS, LSBC, LSBV))

    @staticmethod
    def DifferentialImage(data):
        return np.diff(data, axis=0)

    @staticmethod
    def TotalNoise(data, Differential = True):
        X = DataProcessing.TemporalAverage(data)
        return np.std(X) / np.sqrt(2) if Differential else np.std(data)

    @staticmethod
    def FrameNoise(data, Differential = True):
        X = DataProcessing.SpatialAverage(data)
        return np.std(X) / np.sqrt(2) if Differential else np.std(X)

    @staticmethod
    def LineNoise(data, Differential = True, Orientation = 'Row'):
        X = DataProcessing.TemporalAverage(data)
        LineMeanX = DataProcessing.LineMean(X, Orientation)
        return np.std(LineMeanX) / np.sqrt(2) if Differential else np.std(LineMeanX)

    @staticmethod
    def PixelNoise(TotalNoise, FrameNoise, LineNoise):
        return np.sqrt(np.square(TotalNoise) - np.square(FrameNoise) - np.square(LineNoise))

    @staticmethod
    def IQR_Mask(imageinfo, MaskedArray, NIQR):
        Q1, Q2, Q3 = np.percentile(MaskedArray, [25, 50, 75], method='nearest')
        IQR = Q3 - Q1
        return (Q1 - NIQR * IQR < imageinfo) & (imageinfo < Q3 + NIQR * IQR)

    @staticmethod
    def Data2Histogram(data, Mask = None):
        if Mask is not None:
            MaskedArray = np.ma.masked_array(data, mask=~Mask)[np.ma.masked_array(data, mask=~Mask).mask == False].data
            hist = np.histogram(MaskedArray, bins=int(np.ma.masked_array(data, ~Mask).max() - np.ma.masked_array(data, ~Mask).min()))

        else:
            hist = np.histogram(data, bins=int(data.max() - data.min()))

        return np.transpose(np.array((hist[1][:-1], hist[0])))

    @staticmethod
    def Highpass_Filter(data):

        Filtered_Data = np.zeros_like(data)
        binomial_filter = 1/16 * np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]])

        if data.ndim == 3:
            for l in range(data.shape[0]):
                img = data[l, :, :].copy()
                img = uniform_filter(img, size=7, mode='reflect')
                img = uniform_filter(img, size=11, mode='reflect')
                img = convolve(img, weights = binomial_filter, mode='reflect')

                Filtered_Data[l, :, :] = data[l, :, :] - img
        if data.ndim == 2:
            img = data.copy()
            img = uniform_filter(img, size=7, mode='reflect')
            img = uniform_filter(img, size=11, mode='reflect')
            img = convolve(img, weights=binomial_filter, mode='reflect')

            Filtered_Data = data - img
        return Filtered_Data

    @staticmethod
    def SelectBlock(data, r1, r2, c1, c2):

        L, R, U, D = min(c1, c2), max(c1, c2), min(r1, r2), max(r1, r2)

        if data.ndim == 2:
            return data[U:D, L:R]

        if data.ndim == 3:
            return data[:, U:D, L:R]

    @staticmethod
    def LineMean(data, Orientation):
        return np.mean(data, axis=1) if Orientation == 'Row' else np.mean(data, axis=0)
    @staticmethod
    def LineCalibration(data, Orientation='Row'):
        return data + data.mean() - (DataProcessing.LineMean(data, Orientation))[:, None]


class Plotting:
    @staticmethod
    def MakeFigureWidget(frame, figuresize):
        fig, ax = plt.subplots(figsize=figuresize, tight_layout=True)
        tk_plt = FigureCanvasTkAgg(fig, frame)
        tk_plt.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
        plt.close(fig)
        return ax

    @staticmethod
    def ShowImage(imageinfo, ax):
        ax.cla()
        ax.imshow(imageinfo, cmap="gray", vmin=np.mean(imageinfo) - 3*np.std(imageinfo), vmax=np.mean(imageinfo) + 3*np.std(imageinfo), origin='lower')
        plt.pause(0.1)

    @staticmethod
    def DrawDivision(ax, imageInfo, numRow, numCol):

        totalRow, totalCol = imageInfo.shape[0], imageInfo.shape[1]

        plt.pause(0.1)
        plt.ion()

        for j in range(numRow - 1):
            ax.axhline(y=totalRow*(j+1)/numRow, xmin=0, xmax=totalCol-1, color='red')
        for k in range(numCol - 1):
            ax.axvline(x=totalCol*(k+1)/numCol, ymin=0, ymax=totalRow-1, color='red')


class UIConfiguration:

    @staticmethod
    def set_text(entry, text):
        entry.delete(0, END)
        entry.insert (0, text)
        return

    @staticmethod
    def Save2Clipboard(data):
        df = pd.DataFrame(data)
        df.to_clipboard(excel=True)

    @staticmethod
    def ButtonState(button, condition):
        if condition == True:
            button["state"] = 'normal'
        else:
            button["state"] = 'disable'
