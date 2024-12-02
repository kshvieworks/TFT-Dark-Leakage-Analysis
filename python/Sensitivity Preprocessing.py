import numpy as np
import pandas as pd
import pathlib
import re
import xlsxwriter
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.image as mImage
import multiprocessing as mp
import time
import tkinter
from tkinter import *
from os import listdir
from os.path import isfile, join

fd = pathlib.Path(__file__).parent.resolve()
Image_Size = [1280, 1280]
fw = int(Image_Size[0]/2)
fh = int(Image_Size[1]/2)
fs = (fw/200, fh/200)
fsPTC = (fw/200, fh/100)

class TFTDarkLeakageAnalysis:
    def __init__(self, window):
        self.window = window
        self.window.title("Frame Averaging and Signal Mean")
        # self.window.config(background='#FFFFFF')
        self.window.geometry(f"{fw+350}x{int(2*fh/3)}")
        self.window.resizable(False, False)

        # self.FOI = np.zeros(2, dtype=int)
        self.ROI1 = np.zeros(2, dtype=int)
        self.ROI2 = np.zeros(2, dtype=int)
        self.Point1 = np.zeros(2, dtype=int)
        self.Point2 = np.zeros(2, dtype=int)
        self.filepath = ""
        self.read_data = np.array([], dtype=np.float64)
        self.dark_data = np.array([], dtype=np.float64)
        self.frame_average = np.array([], dtype=np.float64)
        self.variance_ij = np.array([], dtype=np.float64)
        self.Average = 0

        # self.ImageWidget = FALSE
        # self.ROIWidget = FALSE
        self.Output = FALSE

        self.__main__()

    def Open_Path(self):

        self.__init__(window)

        self.filepath = tkinter.filedialog.askopenfilename(initialdir=f"{fd}/")
        self.label.configure(text=f"{self.filepath[-20:]}")

    def Read_Image(self):

        file_now = self.filepath
        if file_now[-3:] == 'raw':

            fid = open(f"{file_now}", "rb")
            read_data_now = np.fromfile(fid, dtype=np.uint16, sep="")
            read_data_now = read_data_now.reshape(Image_Size)
            fid.close()

            self.read_data = read_data_now.copy()

        self.frame_average = self.read_data.copy()

        if not hasattr(self, 'ImageWidget'):
            self.ImageWidget = self.MakeFigureWidget(self.ImagePlotFrame)
        self.Show_Image(self.frame_average, self.ImageWidget)

        self.set_text(self.ROI1_EntryX, '0')
        self.set_text(self.ROI1_EntryY, f"{int(self.frame_average.shape[0] - 1)}")
        self.set_text(self.ROI2_EntryX, f"{int(self.frame_average.shape[1] - 1)}")
        self.set_text(self.ROI2_EntryY, '0')

        self.label_ROI1.configure(text=f"{int(self.ROI1_EntryX.get()), int(self.ROI1_EntryY.get())}")
        self.label_ROI2.configure(text=f"{int(self.ROI2_EntryX.get()), int(self.ROI2_EntryY.get())}")

        # self.Update_ROI(np.array([0, 0]), np.flip(self.frame_average.shape) - np.array([1, 1]))

    def DarkBTNEvent(self):
        fpath = tkinter.filedialog.askopenfilename(initialdir=f"{self.filepath}/")
        self.dklabel.configure(text=f"{fpath[-20:]}")

        file_now = fpath
        if file_now[-3:] == 'raw':

            fid = open(f"{file_now}", "rb")
            read_data_now = np.fromfile(fid, dtype=np.uint16, sep="")
            read_data_now = read_data_now.reshape(Image_Size)
            fid.close()

            self.dark_data = read_data_now.copy()

    def Update_ROI(self, BTN, ROI_X, ROI_Y):

        if BTN == 'ROI1':
            self.ROI1[0], self.ROI1[1] = int(ROI_X), int(ROI_Y)
            self.label_ROI1.configure(text=f"({int(self.ROI1[0])}, {int(self.ROI1[1])})")

        else:
            self.ROI2[0], self.ROI2[1] = int(ROI_X), int(ROI_Y)
            self.label_ROI2.configure(text=f"({int(self.ROI2[0])}, {int(self.ROI2[1])})")

    def Show_ROI(self):

        if self.ROI1[0] <= self.ROI2[0]:
            ROI_L = self.ROI1[0]
            ROI_R = self.ROI2[0]
        else:
            ROI_L = self.ROI2[0]
            ROI_R = self.ROI1[0]

        if self.ROI1[1] <= self.ROI2[1]:
            ROI_U = self.ROI1[1]
            ROI_D = self.ROI2[1]
        else:
            ROI_U = self.ROI2[1]
            ROI_D = self.ROI1[1]

        self.frame_average = self.frame_average - self.dark_data

        self.ROI_Data = self.frame_average[ROI_U:ROI_D, ROI_L:ROI_R]
        self.ROI_Frame_Average = self.ROI_Data.copy()

        if not hasattr(self, 'ROIWidget'):
            self.ROIWidget = self.MakeFigureWidget(self.ROIPlotFrame)
        self.Show_Image(self.ROI_Frame_Average, self.ROIWidget)

        asdf = 1

    def DivisionBTNEvent(self, ax, Frame, row, col):

        nr = Frame.shape[0]
        nc = Frame.shape[1]

        plt.pause(0.1)
        plt.ion()
        for j in range(row-1):
            ax.axhline(y=nr*(j+1)/row, xmin=0, xmax=nc-1, color='red')
        for k in range(col-1):
            ax.axvline(x=nc*(k+1)/col, ymin=0, ymax=nr-1, color='red')

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

    def MakeFigureWidget(self, frame):
        fig, ax = plt.subplots(figsize=fs, tight_layout=True)
        tk_plt = FigureCanvasTkAgg(fig, frame)
        tk_plt.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
        plt.close(fig)
        return ax

    def Show_Image(self, imageinfo, ax):
        ax.cla()
        ax.imshow(imageinfo, cmap="gray", vmin=np.amin(imageinfo), vmax=np.amax(imageinfo), origin='lower')

    def Plot_Hist(self, Data, frame):
        Data = Data.ravel()
        bins = 100
        fig, ax = plt.subplots(figsize=fs, tight_layout=True)
        ax.hist(Data, bins=bins, edgecolor='red', color='green')
        output_plt = FigureCanvasTkAgg(fig, frame)
        output_plt.get_tk_widget().pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=1)

        self.label_stddev.configure(text=f"{np.format_float_scientific(np.std(Data), unique=False, precision=3)}")
        self.label_mean.configure(text=f"{np.format_float_scientific(np.mean(Data), unique=False, precision=3)}")

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

        self.InputFrame = tkinter.Frame(self.window, width=fw, height=fh+100)
        self.InputFrame.grid(column=0, row=0)
        self.ImagePlotFrame = tkinter.Frame(self.InputFrame, bg='white', width=fw/2, height=fh/2)
        self.ImagePlotFrame.grid(column=0, row=0)
        self.ROIPlotFrame = tkinter.Frame(self.InputFrame, bg='white', width=fw / 2, height=fh / 2)
        self.ROIPlotFrame.grid(column=1, row=0)

        self.InputinfoFrame = tkinter.Frame(self.InputFrame, width=fw, height=100)
        self.InputinfoFrame.grid(column=0, row=1, columnspan = 2)

        # self.OutputFrame = tkinter.Frame(self.window, width=fw/2, height=fh/2)
        # self.OutputFrame.grid(column=1, row=0)
        # self.OutputPlotFrame = tkinter.Frame(self.OutputFrame, bg='white', width=fw/2, height=fh/2)
        # self.OutputPlotFrame.grid(column=0, row=0)
        # self.OutputinfoFrame = tkinter.Frame(self.OutputFrame, width=fh, height=fh)
        # self.OutputinfoFrame.grid(column=1, row=0)


        labelspan = 1
        self.label = tkinter.Label(self.InputinfoFrame)
        self.label.grid(column=0, row=1, columnspan=2)
        self.OpenButton = tkinter.Button(self.InputinfoFrame, text='Open Path', command=self.Open_Path)
        self.OpenButton.grid(column=0, row=2)

        Readspan = 1
        self.ReadButton = tkinter.Button(self.InputinfoFrame, text='Read Files', command=self.Read_Image)
        self.ReadButton.grid(column=labelspan, row=2)

        FOIspan = 2
        # self.FOI_EntryStart = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        # self.FOI_EntryStart.grid(column=labelspan+Readspan, row=3)
        # self.FOI_EntryStart.insert(0, '0')
        # self.FOI_EntryEnd = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        # self.FOI_EntryEnd.grid(column=labelspan+Readspan+1, row=3)
        # self.FOI_EntryEnd.insert(0, '0')
        self.dklabel = tkinter.Label(self.InputinfoFrame)
        self.dklabel.grid(column=2, row=1, columnspan=10)
        self.DarkBTN = tkinter.Button(self.InputinfoFrame, text='Dark File', command=self.DarkBTNEvent)
        self.DarkBTN.grid(column=labelspan+Readspan, row=2, columnspan=FOIspan)

        ROI1span = 2
        self.ROI1_EntryX = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.ROI1_EntryX.grid(column=labelspan+Readspan+FOIspan, row=3)
        self.ROI1_EntryX.insert(0, '0')
        self.ROI1_EntryY = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.ROI1_EntryY.grid(column=labelspan+Readspan+FOIspan+1, row=3)
        self.ROI1_EntryY.insert(0, '0')
        self.ROI1BTN = tkinter.Button(self.InputinfoFrame, text='ROI1(Left, DN)', command=lambda m='ROI1': self.Update_ROI(m, int(self.ROI1_EntryX.get()), int(self.ROI1_EntryY.get())))
        self.ROI1BTN.grid(column=labelspan+Readspan+FOIspan, row=2, columnspan=ROI1span)

        ROI2span = 2
        self.ROI2_EntryX = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.ROI2_EntryX.grid(column=labelspan+Readspan+FOIspan+ROI1span, row=3)
        self.ROI2_EntryX.insert(0, '0')
        self.ROI2_EntryY = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.ROI2_EntryY.grid(column=labelspan+Readspan+FOIspan+ROI1span+1, row=3)
        self.ROI2_EntryY.insert(0, '0')
        self.ROI2BTN = tkinter.Button(self.InputinfoFrame, text='ROI2(Right, Up)', command=lambda m='ROI2': self.Update_ROI(m, int(self.ROI2_EntryX.get()), int(self.ROI2_EntryY.get())))
        self.ROI2BTN.grid(column=labelspan+Readspan+FOIspan+ROI1span, row=2, columnspan=ROI2span)

        ROIShowspan = 2
        self.ROI_Show_BTN = tkinter.Button(self.InputinfoFrame, text='Show ROI', command=self.Show_ROI)
        self.ROI_Show_BTN.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=2, columnspan=ROIShowspan)
        # self.label_FOI_prompt = tkinter.Label(self.InputinfoFrame, text='FOI')
        # self.label_FOI_prompt.grid(column=labelspan + Readspan + FOIspan + ROI1span + ROI2span, row=3)
        self.label_ROI1_prompt = tkinter.Label(self.InputinfoFrame, text='ROI1')
        self.label_ROI1_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=3)
        self.label_ROI2_prompt = tkinter.Label(self.InputinfoFrame, text='ROI2')
        self.label_ROI2_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span, row=4)
        # self.label_FOI = tkinter.Label(self.InputinfoFrame)
        # self.label_FOI.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=3)
        self.label_ROI1 = tkinter.Label(self.InputinfoFrame)
        self.label_ROI1.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=3)
        self.label_ROI2 = tkinter.Label(self.InputinfoFrame)
        self.label_ROI2.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+1, row=4)

        Divspan = 2
        self.Division_EntryX = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Division_EntryX.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=3)
        self.Division_EntryX.insert(0, '1')
        self.Division_EntryY = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Division_EntryY.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+1, row=3)
        self.Division_EntryY.insert(0, '1')
        self.DivisionBTN = tkinter.Button(self.InputinfoFrame, text='Division(Column, Row)', command=lambda: self.DivisionBTNEvent(self.ROIWidget, self.ROI_Frame_Average, int(self.Division_EntryX.get()), int(self.Division_EntryY.get())))
        self.DivisionBTN.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan, row=2, columnspan=Divspan)

        Calcspan = 2
        self.CalculateBTN = tkinter.Button(self.InputinfoFrame, text='Calculate', command=lambda: self.CalculateBTNEvent(self.ROIWidget, self.ROI_Frame_Average, int(self.Division_EntryX.get()), int(self.Division_EntryY.get())))
        self.CalculateBTN.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan, row=2, columnspan=Calcspan)
        self.mean_prompt = tkinter.Label(self.InputinfoFrame, text='Mean')
        self.mean_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan, row=3)
        self.stddev_prompt = tkinter.Label(self.InputinfoFrame, text='stddev')
        self.stddev_prompt.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan, row=4)
        self.label_mean = tkinter.Label(self.InputinfoFrame)
        self.label_mean.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+1, row=3)
        self.label_stddev = tkinter.Label(self.InputinfoFrame)
        self.label_stddev.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+1, row=4)

        Savespan = 1
        self.SaveButton = tkinter.Button(self.InputinfoFrame, text='Save Files', command=self.SaveBTNEvent)
        self.SaveButton.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+Calcspan, row=2)
        self.SavePath = tkinter.Label(self.InputinfoFrame)
        self.SavePath.grid(column=labelspan+Readspan+FOIspan+ROI1span+ROI2span+ROIShowspan+Divspan+Calcspan, row=1)

        self.SaveClipboardBoardBTN = tkinter.Button(self.InputinfoFrame, text='Save Clipboard', command=self.SaveClipboardBTNEvent)
        self.SaveClipboardBoardBTN.grid(column=labelspan + Readspan + FOIspan + ROI1span + ROI2span + ROIShowspan + Divspan + Calcspan, row=3)


if __name__ == '__main__':
    window = tkinter.Tk()
    TFTDarkLeakageAnalysis(window)
    window.mainloop()