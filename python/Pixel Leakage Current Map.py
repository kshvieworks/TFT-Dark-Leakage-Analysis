import numpy as np
import pandas as pd
import pathlib
import re
import xlsxwriter
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import ticker
import matplotlib.image as mImage
import multiprocessing as mp
import time
import tkinter
from tkinter import *
from os import listdir
from os.path import isfile, join

tkfont = {'Font': 'Calibri', 'FontSize': 10}
tickfontstyle = {'Font': 'Calibri', 'FontSize': 18}
fontstyle = {'Font': 'Calibri', 'FontSize': 24}
DefaultInfos = {'Title': 'Title', 'xAxisTitle': 'Dark Current Density, J$_{Dark}$ [A/cm$^2$]', 'yAxisTitle': 'Count',
                'xLim': (0, 1), 'yLim': (0, 1), 'MajorTickX': 1, 'MajorTickY': 1}

fd = pathlib.Path(__file__).parent.resolve()
Image_Size = [256, 256]
fw = int(1280/2)
fh = int(1280/2)
fs = (fw/200, fh/200)
fsPTC = (fw/200, fh/100)

class TFTDarkLeakageAnalysis:
    def __init__(self, window):
        self.window = window
        self.window.title("Pixel Leakage Current Map")
        # self.window.config(background='#FFFFFF')
        self.window.geometry(f"{fw+350}x{2*int(2*fh/3)}")
        self.window.resizable(False, False)

        # self.FOI = np.zeros(2, dtype=int)
        # self.ROI1 = np.zeros(2, dtype=int)
        # self.ROI2 = np.zeros(2, dtype=int)
        # self.Point1 = np.zeros(2, dtype=int)
        # self.Point2 = np.zeros(2, dtype=int)

        self.filepath1 = ""
        self.filepath2 = ""

        self.read_data1 = np.array([], dtype=np.float64)
        self.frame_average1 = np.array([], dtype=np.float64)
        self.Average1 = 0
        self.read_data2 = np.array([], dtype=np.float64)
        self.frame_average2 = np.array([], dtype=np.float64)
        self.Average2 = 0
        # self.ImageWidget = FALSE
        # self.ROIWidget = FALSE
        self.Output = np.array([], dtype=np.float64)

        self.__main__()

    def Open_Path(self, m):

        # self.__init__(window)
        if m == '1':
            self.filepath1 = tkinter.filedialog.askdirectory(initialdir=f"{fd}/")
            self.label1.configure(text=f"{self.filepath1[-100:]}")

        else:
            self.filepath2 = tkinter.filedialog.askdirectory(initialdir=f"{fd}/")
            self.label2.configure(text=f"{self.filepath2[-100:]}")

    def Read_Image(self, m):

        if m == '1':

            onlyfiles = [f for f in listdir(self.filepath1) if isfile(join(self.filepath1, f))]

            for file_now in onlyfiles:
                if file_now[-3:] == 'raw':

                    fid = open(f"{self.filepath1}/{file_now}", "rb")
                    read_data_now = np.fromfile(fid, dtype=np.uint16, sep="")
                    read_data_now = read_data_now.reshape(Image_Size)
                    fid.close()

                    if not self.read_data1.any():
                        self.read_data1 = read_data_now[np.newaxis, :]
                        continue

                    self.read_data1 = np.append(self.read_data1, read_data_now[np.newaxis, :], axis=0)

            self.frame_average1 = self.read_data1.sum(axis=0) / len(self.read_data1)

            if not hasattr(self, 'ImageWidget1'):
                self.ImageWidget1 = self.MakeFigureWidget(self.ImagePlotFrame1)
            self.Show_Image(self.frame_average1, self.ImageWidget1)

            self.set_text(self.FOI_EntryStart, '1')
            self.set_text(self.FOI_EntryEnd, f"{int(len(self.read_data1))}")
            self.set_text(self.ROI1_EntryX, '0')
            self.set_text(self.ROI1_EntryY, f"{int(self.frame_average1.shape[0] - 1)}")
            self.set_text(self.ROI2_EntryX, f"{int(self.frame_average1.shape[1] - 1)}")
            self.set_text(self.ROI2_EntryY, '0')

            self.label1_FOI.configure(text=f"{int(self.FOI_EntryStart.get()), int(self.FOI_EntryEnd.get())}")
            self.label1_ROI1.configure(text=f"{int(self.ROI1_EntryX.get()), int(self.ROI1_EntryY.get())}")
            self.label1_ROI2.configure(text=f"{int(self.ROI2_EntryX.get()), int(self.ROI2_EntryY.get())}")

        else:
            onlyfiles = [f for f in listdir(self.filepath2) if isfile(join(self.filepath2, f))]

            for file_now in onlyfiles:
                if file_now[-3:] == 'raw':

                    fid = open(f"{self.filepath2}/{file_now}", "rb")
                    read_data_now = np.fromfile(fid, dtype=np.uint16, sep="")
                    read_data_now = read_data_now.reshape(Image_Size)
                    fid.close()

                    if not self.read_data2.any():
                        self.read_data2 = read_data_now[np.newaxis, :]
                        continue

                    self.read_data2 = np.append(self.read_data2, read_data_now[np.newaxis, :], axis=0)

            self.frame_average2 = self.read_data2.sum(axis=0) / len(self.read_data2)

            if not hasattr(self, 'ImageWidget2'):
                self.ImageWidget2 = self.MakeFigureWidget(self.ImagePlotFrame2)
            self.Show_Image(self.frame_average2, self.ImageWidget2)

            self.set_text(self.FOI_EntryStart, '1')
            self.set_text(self.FOI_EntryEnd, f"{int(len(self.read_data2))}")
            self.set_text(self.ROI1_EntryX, '0')
            self.set_text(self.ROI1_EntryY, f"{int(self.frame_average2.shape[0] - 1)}")
            self.set_text(self.ROI2_EntryX, f"{int(self.frame_average2.shape[1] - 1)}")
            self.set_text(self.ROI2_EntryY, '0')

            self.label2_FOI.configure(text=f"{int(self.FOI_EntryStart.get()), int(self.FOI_EntryEnd.get())}")
            self.label2_ROI1.configure(text=f"{int(self.ROI1_EntryX.get()), int(self.ROI1_EntryY.get())}")
            self.label2_ROI2.configure(text=f"{int(self.ROI2_EntryX.get()), int(self.ROI2_EntryY.get())}")

    def Update_ROI(self, BTN, ROI_X, ROI_Y):

        if BTN == 'FOI1':
            self.label1_FOI.configure(text=f"({int(ROI_X)}, {int(ROI_Y)})")
        elif BTN == 'ROI1_1':
            self.label1_ROI1.configure(text=f"({int(ROI_X)}, {int(ROI_Y)})")
        elif BTN == 'ROI2_1':
            self.label1_ROI2.configure(text=f"({int(ROI_X)}, {int(ROI_Y)})")

        if BTN == 'FOI2':
            self.label2_FOI.configure(text=f"({int(ROI_X)}, {int(ROI_Y)})")
        elif BTN == 'ROI1_2':
            self.label2_ROI1.configure(text=f"({int(ROI_X)}, {int(ROI_Y)})")
        elif BTN == 'ROI2_2':
            self.label2_ROI2.configure(text=f"({int(ROI_X)}, {int(ROI_Y)})")

    def Show_ROI(self, BTN, FOI, ROI1, ROI2):

        FOI_Start = int(re.findall(r'\d+', FOI)[0])
        FOI_End = int(re.findall(r'\d+', FOI)[1])

        ROI1_X = int(re.findall(r'\d+', ROI1)[0])
        ROI1_Y = int(re.findall(r'\d+', ROI1)[1])

        ROI2_X = int(re.findall(r'\d+', ROI2)[0])
        ROI2_Y = int(re.findall(r'\d+', ROI2)[1])

        if ROI1_X <= ROI2_X:
            ROI_L = ROI1_X
            ROI_R = ROI2_X
        else:
            ROI_L = ROI2_X
            ROI_R = ROI1_X

        if ROI1_Y <= ROI2_Y:
            ROI_U = ROI1_Y
            ROI_D = ROI2_Y
        else:
            ROI_U = ROI2_Y
            ROI_D = ROI1_Y

        if BTN == '1':

            self.ROI_Data1 = self.read_data1[FOI_Start-1:FOI_End, ROI_U:ROI_D, ROI_L:ROI_R]
            self.ROI_Frame_Average1 = self.ROI_Data1.sum(axis=0) / len(self.ROI_Data1)

            if not hasattr(self, 'ROIWidget1'):
                self.ROIWidget1 = self.MakeFigureWidget(self.ROIPlotFrame1)
            self.Show_Image(self.ROI_Frame_Average1, self.ROIWidget1)

        else:
            self.ROI_Data2 = self.read_data2[FOI_Start-1:FOI_End, ROI_U:ROI_D, ROI_L:ROI_R]
            self.ROI_Frame_Average2 = self.ROI_Data2.sum(axis=0) / len(self.ROI_Data2)

            if not hasattr(self, 'ROIWidget2'):
                self.ROIWidget2 = self.MakeFigureWidget(self.ROIPlotFrame2)
            self.Show_Image(self.ROI_Frame_Average2, self.ROIWidget2)

            self.set_text(self.Config_EntryX, '0')
            self.set_text(self.Config_EntryY, '0')
            self.set_text(self.Config_EntrydX, self.ROI_Frame_Average2.shape[1])
            self.set_text(self.Config_EntrydY, self.ROI_Frame_Average2.shape[0])

            dT = np.abs(int(re.findall(r'\d+', self.filepath1)[-1]) - int(re.findall(r'\d+', self.filepath2)[-1]))
            self.set_text(self.Config_EntrydT, dT)

    def ConfigBTNEvent(self, X, Y, dX, dY, dT, IFS, W, H):

        '1. Difference of Pixel Value btw. Two Image'
        self.Output = np.abs(self.ROI_Frame_Average1 - self.ROI_Frame_Average2)

        '2. Apply Region of Interest'
        self.Output = self.Output[X:X+dX, Y:Y+dY]

        '3 Unit Conversion(DN to Current)'
        self.Output = self.DN2Coulomb(self.Output, IFS=IFS, LSBC=0.125E-12, LSBV=61.035E-6) / dT

        '4. Unit Conversion(Current to Current Density, A/cm$^2$)'
        self.Output = self.Current2CurrentDensity(self.Output, W, H)

    def CalculateBTNEvent(self, ax, Frame, row, col):

        nr = Frame.shape[0]
        nc = Frame.shape[1]
        mask = np.zeros((nr, nc))
        self.Average = np.zeros(row*col)
        Division = np.zeros((row*col, nr, nc))

        for j in range(row):
            for k in range(col):
                mask[int(nr*j/row):int(nr*(j+1)/row), int(nc*k/col):int(nc*(k+1)/col)] = 1
                Division[j*col + k] = Frame * mask
                self.Average[j*col + k] = self.Show_Division_Average(ax, Division[j*col + k], int(nr*(j+0.5)/row), int(nc*(k+0.5)/col))
                mask[:, :] = 0

        self.label_stddev.configure(text=f"{np.format_float_scientific(np.std(Frame), unique=False, precision=2)}")
        self.label_mean.configure(text=f"{np.format_float_scientific(np.mean(Frame), unique=False, precision=2)}")

    def Show_Division_Average(self, ax, Data, r, c):

        avg = Data.sum() / np.count_nonzero(Data)
        ax.text(c, r, int(avg), c='r', horizontalalignment='center', verticalalignment='center')
        return avg

    def DN2Coulomb(self, DN, IFS=2, LSBC = 0.125E-12, LSBV=61.035E-6):
        return LSBC*(IFS+1)*LSBV*DN

    def Current2CurrentDensity(self, Data, W, H):
        return Data/(W*1E-4*H*1E-4)
    def MakeFigureWidget(self, frame):
        fig, ax = plt.subplots(figsize=fs, tight_layout=True)
        tk_plt = FigureCanvasTkAgg(fig, frame)
        tk_plt.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
        plt.close(fig)
        return ax

    def Show_Image(self, imageinfo, ax):
        ax.cla()
        c = ax.imshow(imageinfo, cmap="gray", vmin=np.amin(imageinfo), vmax=np.amax(imageinfo), origin='lower')
        return c

    def Plot_Hist(self, Data):
        Data = Data.ravel()
        bins = 100

        # Data_Effective = Data[np.argwhere(np.abs(Data-np.mean(Data))<3*np.std(Data)).flatten()]
        params = f"$\mu$={np.format_float_scientific(np.mean(Data), unique=False, precision=2)} \n" \
                 f"$\sigma$={np.format_float_scientific(np.std(Data), unique=False, precision=2)}"

        fig, ax = plt.subplots(figsize=fs, tight_layout=True)
        ax.hist(Data, bins=bins, edgecolor='red', color='green', label = params, range = (0, 8.9E-9))
        ax.legend(loc='best', fontsize=fontstyle['FontSize'])
        self.FigureOptionSetting(ax)
        plt.pause(0.01)

        self.label_stddev.configure(text=f"{np.format_float_scientific(np.std(Data), unique=False, precision=3)}")
        self.label_mean.configure(text=f"{np.format_float_scientific(np.mean(Data), unique=False, precision=3)}")
# ax.set_xlim(0, 8.9E-9)
# ax.set_ylim(0, 49000)
# ax.xaxis.set_major_locator(MultipleLocator(2E-9))
# ax.yaxis.set_major_locator(MultipleLocator(10000))
# self.forceAspect(ax)

# Log scale
# fig, ax = plt.subplots(figsize=fs, tight_layout=True)
# ax.hist(Data, bins=bins, edgecolor='red', color='green', label=params, range = (0, 8.9E-9))
# ax.set_yscale('log')
# ax.set_ylim(0.1, 49000)
# ax.set_xlim(0, 8.9E-9)
# ax.xaxis.set_major_locator(MultipleLocator(2E-9))
# ax.legend(loc='best', fontsize=fontstyle['FontSize'])
# self.FigureOptionSetting(ax)
# plt.pause(0.01)

    def FigureOptionSetting(self, ax):

        ax.set_xlabel(DefaultInfos['xAxisTitle'], font=fontstyle['Font'], fontsize=fontstyle['FontSize'])
        ax.set_ylabel(DefaultInfos['yAxisTitle'], font=fontstyle['Font'], fontsize=fontstyle['FontSize'])
        # ax.set_xlim(DefaultInfos['xLim'][0], DefaultInfos['xLim'][1])
        # ax.set_ylim(DefaultInfos['yLim'][0], DefaultInfos['yLim'][1])
        # ax.xaxis.set_major_locator(MultipleLocator(DefaultInfos['MajorTickX']))
        # ax.yaxis.set_major_locator(MultipleLocator(DefaultInfos['MajorTickY']))

        ax.grid(True)
        ax.tick_params(axis='x', labelsize=tickfontstyle['FontSize'])
        ax.tick_params(axis='y', labelsize=tickfontstyle['FontSize'])
        plt.tight_layout()
        # self.forceAspect(ax)

    def forceAspect(self, ax, aspect=1):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_aspect(abs((xlim[1] - xlim[0]) / (ylim[1] - ylim[0])) / aspect)

    def DrawLeakage(self, Data):

        fig, ax = plt.subplots(figsize=fs, tight_layout=True)
        c = self.Show_Image(Data, ax)
        ax.cbar = fig.colorbar(c, ax=ax)
        ax.cbar.set_label(label='Leakage Current [nA/cm$^2$]', size=24)
        ax.set_xlabel('Horizontal Pixel Entry', font=fontstyle['Font'], fontsize=fontstyle['FontSize'])
        ax.set_ylabel('Vertical Pixel Entry', font=fontstyle['Font'], fontsize=fontstyle['FontSize'])
        ax.tick_params(axis='x', labelsize=tickfontstyle['FontSize'])
        ax.tick_params(axis='y', labelsize=tickfontstyle['FontSize'])
        ax.cbar.ax.tick_params(labelsize=16)
        ax.cbar.locator = ticker.MaxNLocator(nbins=2)
        ax.cbar.update_ticks()
        plt.pause(0.01)
# c.set_clim(0, 8.9E-9)

#Log Scale
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import LogFormatterExponent
# fig, ax = plt.subplots(figsize=fs, tight_layout=True)
# c = ax.imshow(Data, norm=LogNorm())
# ax.cbar = fig.colorbar(c, ax=ax)
# ax.cbar.set_label(label='Leakage Current [A/cm$^2$]', size=24)
# ax.set_xlabel('Horizontal Pixel Entry', font=fontstyle['Font'], fontsize=fontstyle['FontSize'])
# ax.set_ylabel('Vertical Pixel Entry', font=fontstyle['Font'], fontsize=fontstyle['FontSize'])
# ax.tick_params(axis='x', labelsize=tickfontstyle['FontSize'])
# ax.tick_params(axis='y', labelsize=tickfontstyle['FontSize'])
# ax.cbar.ax.tick_params(labelsize=16)
# ax.cbar.locator = ticker.MaxNLocator(nbins=2)
# ax.cbar.formatter = LogFormatterExponent(base=10)
# ax.cbar.update_ticks()
# plt.pause(0.01)

    def set_text(self, entry, text):
        entry.delete(0, END)
        entry.insert(0, text)
        return

    def SaveBTNEvent(self):

        filepath = tkinter.filedialog.asksaveasfilename(initialdir=f"{self.filepath}/",
                                                        title="Save as",
                                                        filetypes=(("Xlsx Files", ".xlsx"),
                                                                   ("all files", "*")))
        filepath = f"{filepath}.xlsx"

        mean = pd.DataFrame(self.Average)
        writer = pd.ExcelWriter(filepath, engine='xlsxwriter')
        mean.to_excel(writer, sheet_name='Mean')
        writer.close()

    def SaveClipboardBTNEvent(self):

        mean = pd.DataFrame(self.Average)
        mean.to_clipboard(excel=True)

    def __main__(self):

        self.InputFrame1 = tkinter.Frame(self.window, width=fw, height=fh+100)
        self.InputFrame1.grid(column=0, row=0)
        self.ImagePlotFrame1 = tkinter.Frame(self.InputFrame1, bg='white', width=fw/2, height=fh/2)
        self.ImagePlotFrame1.grid(column=0, row=0)
        self.ROIPlotFrame1 = tkinter.Frame(self.InputFrame1, bg='white', width=fw / 2, height=fh / 2)
        self.ROIPlotFrame1.grid(column=1, row=0)

        self.InputinfoFrame1 = tkinter.Frame(self.InputFrame1, width=fw, height=100)
        self.InputinfoFrame1.grid(column=0, row=1, columnspan = 2)

        self.InputFrame2 = tkinter.Frame(self.window, width=fw, height=fh+100)
        self.InputFrame2.grid(column=0, row=1)
        self.ImagePlotFrame2 = tkinter.Frame(self.InputFrame2, bg='white', width=fw/2, height=fh/2)
        self.ImagePlotFrame2.grid(column=0, row=1)
        self.ROIPlotFrame2 = tkinter.Frame(self.InputFrame2, bg='white', width=fw / 2, height=fh / 2)
        self.ROIPlotFrame2.grid(column=1, row=1)

        self.InputinfoFrame2 = tkinter.Frame(self.InputFrame2, width=fw, height=100)
        self.InputinfoFrame2.grid(column=0, row=0, columnspan = 2)

        labelspan = 1
        self.label1 = tkinter.Label(self.InputinfoFrame1)
        self.label1.grid(column=0, row=1, columnspan=10)
        self.OpenButton1 = tkinter.Button(self.InputinfoFrame1, text='Open Path', command=lambda m='1': self.Open_Path(m))
        self.OpenButton1.grid(column=0, row=2)

        self.label2 = tkinter.Label(self.InputinfoFrame2)
        self.label2.grid(column=0, row=0, columnspan=10)
        self.OpenButton2 = tkinter.Button(self.InputinfoFrame2, text='Open Path', command=lambda m='2': self.Open_Path(m))
        self.OpenButton2.grid(column=0, row=1)

        Readspan = 1
        self.ReadButton1 = tkinter.Button(self.InputinfoFrame1, text='Read Files', command=lambda m='1': self.Read_Image(m))
        self.ReadButton1.grid(column=labelspan, row=2)
        self.ReadButton2 = tkinter.Button(self.InputinfoFrame2, text='Read Files', command=lambda m='2': self.Read_Image(m))
        self.ReadButton2.grid(column=labelspan, row=1)

        FOIspan = 2
        self.FOI_EntryStart = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.FOI_EntryStart.grid(column=labelspan+Readspan, row=3)
        self.FOI_EntryStart.insert(0, '0')
        self.FOI_EntryEnd = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.FOI_EntryEnd.grid(column=labelspan+Readspan+1, row=3)
        self.FOI_EntryEnd.insert(0, '0')

        self.FOIBTN1 = tkinter.Button(self.InputinfoFrame1, text='Frame(Start, End)', command=lambda m='FOI1': self.Update_ROI(m, int(self.FOI_EntryStart.get()), int(self.FOI_EntryEnd.get())))
        self.FOIBTN1.grid(column=labelspan+Readspan, row=2, columnspan=FOIspan)
        self.FOIBTN2 = tkinter.Button(self.InputinfoFrame2, text='Frame(Start, End)', command=lambda m='FOI2': self.Update_ROI(m, int(self.FOI_EntryStart.get()), int(self.FOI_EntryEnd.get())))
        self.FOIBTN2.grid(column=labelspan+Readspan, row=1, columnspan=FOIspan)

        ROI1span = 2
        self.ROI1_EntryX = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.ROI1_EntryX.grid(column=labelspan+Readspan+FOIspan, row=3)
        self.ROI1_EntryX.insert(0, '0')
        self.ROI1_EntryY = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.ROI1_EntryY.grid(column=labelspan+Readspan+FOIspan+1, row=3)
        self.ROI1_EntryY.insert(0, '0')

        self.ROI1BTN1 = tkinter.Button(self.InputinfoFrame1, text='ROI1(Left, DN)', command=lambda m='ROI1_1': self.Update_ROI(m, int(self.ROI1_EntryX.get()), int(self.ROI1_EntryY.get())))
        self.ROI1BTN1.grid(column=labelspan+Readspan+FOIspan, row=2, columnspan=ROI1span)
        self.ROI1BTN2 = tkinter.Button(self.InputinfoFrame2, text='ROI1(Left, DN)', command=lambda m='ROI1_2': self.Update_ROI(m, int(self.ROI1_EntryX.get()), int(self.ROI1_EntryY.get())))
        self.ROI1BTN2.grid(column=labelspan+Readspan+FOIspan, row=1, columnspan=ROI1span)

        ROI2span = 2
        self.ROI2_EntryX = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.ROI2_EntryX.grid(column=labelspan+Readspan+FOIspan+ROI1span, row=3)
        self.ROI2_EntryX.insert(0, '0')
        self.ROI2_EntryY = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.ROI2_EntryY.grid(column=labelspan+Readspan+FOIspan+ROI1span+1, row=3)
        self.ROI2_EntryY.insert(0, '0')

        self.ROI2BTN1 = tkinter.Button(self.InputinfoFrame1, text='ROI2(Right, Up)', command=lambda m='ROI2_1': self.Update_ROI(m, int(self.ROI2_EntryX.get()), int(self.ROI2_EntryY.get())))
        self.ROI2BTN1.grid(column=labelspan+Readspan+FOIspan+ROI1span, row=2, columnspan=ROI2span)
        self.ROI2BTN2 = tkinter.Button(self.InputinfoFrame2, text='ROI2(Right, Up)', command=lambda m='ROI2_2': self.Update_ROI(m, int(self.ROI2_EntryX.get()), int(self.ROI2_EntryY.get())))
        self.ROI2BTN2.grid(column=labelspan+Readspan+FOIspan+ROI1span, row=1, columnspan=ROI2span)

        ROIShowspan = 2
        self.ROI_Show_BTN1 = tkinter.Button(self.InputinfoFrame1, text='Show ROI', command=lambda m='1': self.Show_ROI(m, self.label1_FOI.cget('text'), self.label1_ROI1.cget('text'), self.label1_ROI2.cget('text')))
        self.ROI_Show_BTN1.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=2, columnspan=ROIShowspan)
        self.label1_FOI_prompt = tkinter.Label(self.InputinfoFrame1, text='FOI')
        self.label1_FOI_prompt.grid(column=labelspan + Readspan + FOIspan + ROI1span + ROI2span, row=3)
        self.label1_ROI1_prompt = tkinter.Label(self.InputinfoFrame1, text='ROI1')
        self.label1_ROI1_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=4)
        self.label1_ROI2_prompt = tkinter.Label(self.InputinfoFrame1, text='ROI2')
        self.label1_ROI2_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=5)
        self.label1_FOI = tkinter.Label(self.InputinfoFrame1)
        self.label1_FOI.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=3)
        self.label1_ROI1 = tkinter.Label(self.InputinfoFrame1)
        self.label1_ROI1.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=4)
        self.label1_ROI2 = tkinter.Label(self.InputinfoFrame1)
        self.label1_ROI2.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=5)

        self.ROI_Show_BTN2 = tkinter.Button(self.InputinfoFrame2, text='Show ROI', command=lambda m='2': self.Show_ROI(m, self.label2_FOI.cget('text'), self.label2_ROI1.cget('text'), self.label2_ROI2.cget('text')))
        self.ROI_Show_BTN2.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=1, columnspan=ROIShowspan)
        self.label2_FOI_prompt = tkinter.Label(self.InputinfoFrame2, text='FOI')
        self.label2_FOI_prompt.grid(column=labelspan + Readspan + FOIspan + ROI1span + ROI2span, row=2)
        self.label2_ROI1_prompt = tkinter.Label(self.InputinfoFrame2, text='ROI1')
        self.label2_ROI1_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=3)
        self.label2_ROI2_prompt = tkinter.Label(self.InputinfoFrame2, text='ROI2')
        self.label2_ROI2_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=4)
        self.label2_FOI = tkinter.Label(self.InputinfoFrame2)
        self.label2_FOI.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=2)
        self.label2_ROI1 = tkinter.Label(self.InputinfoFrame2)
        self.label2_ROI1.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=3)
        self.label2_ROI2 = tkinter.Label(self.InputinfoFrame2)
        self.label2_ROI2.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=4)

        Divspan = 2
        self.Config_EntryX = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.Config_EntryX.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=3)
        self.Config_EntryX.insert(0, '1')
        self.Config_EntryY = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.Config_EntryY.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+1, row=3)
        self.Config_EntryY.insert(0, '1')

        self.Config_EntrydX = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.Config_EntrydX.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=4)
        self.Config_EntrydX.insert(0, '1')
        self.Config_EntrydY = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.Config_EntrydY.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+1, row=4)
        self.Config_EntrydY.insert(0, '1')

        self.Config_EntrydT = tkinter.Entry(self.InputinfoFrame1, width=4, textvariable="", relief="ridge")
        self.Config_EntrydT.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=5, columnspan=Divspan)
        self.Config_EntrydT.insert(0, '1')

        self.ConfigBTN = tkinter.Button(self.InputinfoFrame1, text='Config\nX, Y\ndX, dY\ndt', command=lambda: \
            self.ConfigBTNEvent(int(self.Config_EntryX.get()), int(self.Config_EntryY.get()), int(self.Config_EntrydX.get()), int(self.Config_EntrydY.get()), float(self.Config_EntrydT.get()), int(self.Config2_EntryIFS.get()), int(self.Config2_EntryWidth.get()), int(self.Config2_EntryHeight.get())))
        self.ConfigBTN.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=2, columnspan=Divspan)

        self.Config2BTN = tkinter.Button(self.InputinfoFrame2, text='Config\nIFS\nW, H [$\mu$m]')
        self.Config2BTN.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=1, columnspan=Divspan)
        self.Config2_EntryIFS = tkinter.Entry(self.InputinfoFrame2, width=4, textvariable="", relief="ridge")
        self.Config2_EntryIFS.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=2, columnspan=Divspan)
        self.Config2_EntryWidth = tkinter.Entry(self.InputinfoFrame2, width=4, textvariable="", relief="ridge")
        self.Config2_EntryWidth.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=3)
        self.Config2_EntryHeight = tkinter.Entry(self.InputinfoFrame2, width=4, textvariable="", relief="ridge")
        self.Config2_EntryHeight.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+1, row=3)

        Calcspan = 2
        self.HistogramBTN = tkinter.Button(self.InputinfoFrame1, text='Histogram', command=lambda: self.Plot_Hist(self.Output))
        self.HistogramBTN.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan, row=2, columnspan=Calcspan)
        self.mean_prompt = tkinter.Label(self.InputinfoFrame1, text='Mean')
        self.mean_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan, row=3)
        self.stddev_prompt = tkinter.Label(self.InputinfoFrame1, text='stddev')
        self.stddev_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan, row=4)
        self.label_mean = tkinter.Label(self.InputinfoFrame1)
        self.label_mean.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+1, row=3)
        self.label_stddev = tkinter.Label(self.InputinfoFrame1)
        self.label_stddev.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+1, row=4)

        Savespan = 1
        self.DrawLeakageButton = tkinter.Button(self.InputinfoFrame1, text='Draw Leakage Map', command=lambda: self.DrawLeakage(self.Output))
        self.DrawLeakageButton.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+Calcspan, row=2)
        # self.SavePath = tkinter.Label(self.InputinfoFrame1)
        # self.SavePath.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+Calcspan, row=1)

        self.SaveClipboardBoardBTN = tkinter.Button(self.InputinfoFrame1, text='Save Clipboard', command=self.SaveClipboardBTNEvent)
        self.SaveClipboardBoardBTN.grid(column=labelspan + Readspan + FOIspan + ROI1span + ROI2span + ROIShowspan + Divspan + Calcspan, row=3)

        # self.DrawPTC_BTN = tkinter.Button(self.OutputinfoFrame, text='Draw PTC', command=self.DrawPTC)
        # self.DrawPTC_BTN.grid(column=0, row=0)
        #
        # self.Point1_BTN = tkinter.Button(self.OutputinfoFrame, text='Point1(Left, DN)', command=self.Select_Point_1)
        # self.Point1_BTN.grid(column=0, row=1)
        # self.Point1_Label = tkinter.Label(self.OutputinfoFrame)
        # self.Point1_Label.grid(column=1, row=1)
        # self.Point2_BTN = tkinter.Button(self.OutputinfoFrame, text='Point2(Right, Up)', command=self.Select_Point_2)
        # self.Point2_BTN.grid(column=0, row=2)
        # self.Point2_Label = tkinter.Label(self.OutputinfoFrame)
        # self.Point2_Label.grid(column=1, row=2)
        # self.GainCalc_BTN = tkinter.Button(self.OutputinfoFrame, text='Calculate Camera Gain', command=self.GainCalc)
        # self.GainCalc_BTN.grid(column=0, row=3)
        # self.label_Line_prompt = tkinter.Label(self.OutputinfoFrame, text='Line Equation')
        # self.label_Line_prompt.grid(column=0, row=4)
        # self.label_Line = tkinter.Label(self.OutputinfoFrame)
        # self.label_Line.grid(column=1, row=4)
        # self.label_Gain_prompt = tkinter.Label(self.OutputinfoFrame, text='Camera Gain [DN/e]')
        # self.label_Gain_prompt.grid(column=0, row=5)
        # self.label_Gain = tkinter.Label(self.OutputinfoFrame)
        # self.label_Gain.grid(column=1, row=5)
        # self.label_Offset_prompt = tkinter.Label(self.OutputinfoFrame, text='Offset Noise [DN]')
        # self.label_Offset_prompt.grid(column=0, row=6)
        # self.label_Offset = tkinter.Label(self.OutputinfoFrame)
        # self.label_Offset.grid(column=1, row=6)


if __name__ == '__main__':
    window = tkinter.Tk()
    TFTDarkLeakageAnalysis(window)
    window.mainloop()