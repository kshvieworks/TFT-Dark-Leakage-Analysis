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

import AnalysisUtil as AU

fd = pathlib.Path(__file__).parent.resolve()
Image_Size = [1280, 1280]
fw = int(Image_Size[0]/2)
fh = int(Image_Size[1]/2)
fs = (fw/200, fh/200)
Image_Size = [207, 207]

class TFTDarkLeakageAnalysis:
    def __init__(self, window):
        self.window = window
        self.window.title("Frame Averaging and Signal Mean")
        # self.window.config(background='#FFFFFF')
        self.window.geometry(f"{fw+450}x{int(2*fh/3)+100}")
        self.window.resizable(True, True)

        self.FOI = np.zeros(2, dtype=int)
        self.ROI1 = np.zeros(2, dtype=int)
        self.ROI2 = np.zeros(2, dtype=int)
        self.Point1 = np.zeros(2, dtype=int)
        self.Point2 = np.zeros(2, dtype=int)
        self.filepath = ""
        self.read_data = np.array([], dtype=np.float64)
        self.frame_average = np.array([], dtype=np.float64)
        self.variance_ij = np.array([], dtype=np.float64)
        self.Average = 0

        self.IFS = IntVar()
        self.SystemGain = 0
        self.Differential = BooleanVar()
        self.NIQR = DoubleVar()
        self.NIteration = IntVar()
        self.ExcludingZero = BooleanVar()
        self.HPF = BooleanVar()
        self.Noise = None
        self.MaskedNoise = None

        # self.ImageWidget = FALSE
        # self.ROIWidget = FALSE
        self.Output = FALSE

        self.__main__()

    def Open_Path(self):

        self.__init__(window)

        self.filepath = AU.ButtonClickedEvent.Open_Path(fd)
        self.label1.configure(text=f"{self.filepath[-100:]}")

    def Read_Image(self):

        self.read_data = AU.ButtonClickedEvent.Read_Image(self.filepath, 'raw', np.uint16, Image_Size)
        self.frame_average = AU.DataProcessing.TemporalAverage(self.read_data)

        if not hasattr(self, 'ImageWidget'):
            self.ImageWidget = AU.Plotting.MakeFigureWidget(self.ImagePlotFrame, fs)

        AU.Plotting.ShowImage(self.frame_average, self.ImageWidget)

        AU.UIConfiguration.set_text(self.Entry3_1, '1')
        AU.UIConfiguration.set_text(self.Entry3_2, f"{int(len(self.read_data))}")
        AU.UIConfiguration.set_text(self.Entry4_1, '0')
        AU.UIConfiguration.set_text(self.Entry4_2, f"{int(self.frame_average.shape[0] - 1)}")
        AU.UIConfiguration.set_text(self.Entry5_1, f"{int(self.frame_average.shape[1] - 1)}")
        AU.UIConfiguration.set_text(self.Entry5_2, '0')
        self.Label6_1_2.configure(text=f"{int(self.Entry3_1.get()), int(self.Entry3_2.get())}")
        self.Label6_2_2.configure(text=f"{int(self.Entry4_1.get()), int(self.Entry4_2.get())}")
        self.Label6_3_2.configure(text=f"{int(self.Entry5_1.get()), int(self.Entry5_2.get())}")

        self.ROI_Data = self.read_data.copy()

    def Update_ROI(self, BTN, ROI_X, ROI_Y):

        if BTN == 'FOI':
            self.FOI[0], self.FOI[1] = int(ROI_X), int(ROI_Y)
            self.Label6_1_2.configure(text=f"({int(self.FOI[0])}, {int(self.FOI[1])})")
        elif BTN == 'ROI1':
            self.ROI1[0], self.ROI1[1] = int(ROI_X), int(ROI_Y)
            self.Label6_2_2.configure(text=f"({int(self.ROI1[0])}, {int(self.ROI1[1])})")
        else:
            self.ROI2[0], self.ROI2[1] = int(ROI_X), int(ROI_Y)
            self.Label6_3_2.configure(text=f"({int(self.ROI2[0])}, {int(self.ROI2[1])})")

    def Show_ROI(self, Frame):

        self.ROI_Data = AU.ButtonClickedEvent.Set_ROI(Frame, self.ROI1, self.ROI2)
        self.ROI_Data = AU.ButtonClickedEvent.Set_FOI(self.ROI_Data, self.FOI)

        ROI_Frame_Average = AU.DataProcessing.TemporalAverage(self.ROI_Data)

        if not hasattr(self, 'ROIWidget'):
            self.ROIWidget = AU.Plotting.MakeFigureWidget(self.ROIPlotFrame, fs)

        AU.Plotting.ShowImage(ROI_Frame_Average, self.ROIWidget)

    def Configurations(self, Frame):

        self.SystemGain = AU.DataProcessing.Coulomb2Electron(AU.DataProcessing.DN2Coulomb(1, self.IFS.get()))

        self.Label7_2_2.configure(
            text=f"{int(self.SystemGain)}")
        if self.Differential.get() == False:
            self.Label7_4_2.configure(text=f"{int(self.FOI[1] - self.FOI[0] + 1)}")
        else:
            self.Label7_4_2.configure(text=f"{int(self.FOI[1] - self.FOI[0])}")

    def Calculate(self, Frame, Widget):

        self.Noise = AU.ButtonClickedEvent.Calc(imageinfo = Frame, Differential = self.Differential.get())

        self.Label8_1_2.configure(text=f'{int(np.round(self.Noise["TotalNoise"] * self.SystemGain, 0))}')
        self.Label8_2_2.configure(text=f'{int(np.round(self.Noise["FrameNoise"] * self.SystemGain, 0))}')
        self.Label8_3_2.configure(text=f'{int(np.round(self.Noise["RowLineNoise"] * self.SystemGain, 0))}')
        self.Label8_4_2.configure(text=f'{int(np.round(self.Noise["ColLineNoise"] * self.SystemGain, 0))}')
        self.Label8_5_2.configure(text=f'{int(np.round(self.Noise["PixelNoise"] * self.SystemGain, 0))}')

        AU.Plotting.ShowImage(AU.DataProcessing.TemporalAverage(self.Noise["ImageInfo"]), Widget)
        AU.UIConfiguration.Save2Clipboard(AU.DataProcessing.Data2Histogram(self.Noise['ImageInfo']))
        plt.pause(0.1)
        plt.ion()

    def IQR(self, Frame, NIQR, NIteration, Widget):
        self.MaskedNoise = AU.ButtonClickedEvent.IQR_Method(Frame, NIQR, NIteration, self.Differential.get(), self.ExcludingZero.get(), self.HPF.get())

        self.Label10_1_2.configure(text=f'{int(np.round(self.MaskedNoise["TotalNoise"] * self.SystemGain, 0))}')
        self.Label10_2_2.configure(text=f'{int(np.round(self.MaskedNoise["FrameNoise"] * self.SystemGain, 0))}')
        self.Label10_3_2.configure(text=f'{int(np.round(self.MaskedNoise["RowLineNoise"] * self.SystemGain, 0))}')
        self.Label10_4_2.configure(text=f'{int(np.round(self.MaskedNoise["ColLineNoise"] * self.SystemGain, 0))}')
        self.Label10_5_2.configure(text=f'{int(np.round(self.MaskedNoise["PixelNoise"] * self.SystemGain, 0))}')

        AU.Plotting.ShowImage(AU.DataProcessing.TemporalAverage(self.MaskedNoise["Mask"]), Widget)
        AU.UIConfiguration.Save2Clipboard(AU.DataProcessing.Data2Histogram(self.MaskedNoise['ImageInfo'], self.MaskedNoise['Mask']))
        plt.pause(0.1)
        plt.ion()

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

    def SaveClipboard(self, data):

        AU.UIConfiguration.Save2Clipboard(data)

    def __main__(self):

        self.InputFrame = tkinter.Frame(self.window, width=fw, height=fh+100)
        self.InputFrame.grid(column=0, row=0)
        self.ImagePlotFrame = tkinter.Frame(self.InputFrame, bg='white', width=fw/2, height=fh/2)
        self.ImagePlotFrame.grid(column=0, row=0)
        self.ROIPlotFrame = tkinter.Frame(self.InputFrame, bg='white', width=fw / 2, height=fh / 2)
        self.ROIPlotFrame.grid(column=1, row=0)

        self.InputinfoFrame = tkinter.Frame(self.InputFrame, width=fw, height=100)
        self.InputinfoFrame.grid(column=0, row=1, columnspan = 2)

        col = 0

        Entry1Span = 1
        self.label1 = tkinter.Label(self.InputinfoFrame)
        self.label1.grid(column=col, row=1, columnspan=10)
        self.Button1 = tkinter.Button(self.InputinfoFrame, text='Open Path', command=self.Open_Path)
        self.Button1.grid(column=col, row=2)
        col = col + Entry1Span

        Entry2Span = 1
        self.Button2 = tkinter.Button(self.InputinfoFrame, text='Read Files', command=self.Read_Image)
        self.Button2.grid(column=col, row=2)
        col = col + Entry2Span

        Entry3span = 2
        self.Entry3_1 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Entry3_1.grid(column=col, row=3)
        self.Entry3_1.insert(0, '0')
        self.Entry3_2 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Entry3_2.grid(column=col+1, row=3)
        self.Entry3_2.insert(0, '0')
        self.Button3 = tkinter.Button(self.InputinfoFrame, text='Frame(Start, End)',
                                      command=lambda m='FOI': self.Update_ROI(m, int(self.Entry3_1.get()), int(self.Entry3_2.get())))
        self.Button3.grid(column=col, row=2, columnspan=Entry3span)
        col = col + Entry3span

        Entry4Span = 2
        self.Entry4_1 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Entry4_1.grid(column=col, row=3)
        self.Entry4_1.insert(0, '0')
        self.Entry4_2 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Entry4_2.grid(column=col+1, row=3)
        self.Entry4_2.insert(0, '0')
        self.Button4 = tkinter.Button(self.InputinfoFrame, text='ROI1(Left, DN)',
                                      command=lambda m='ROI1': self.Update_ROI(m, int(self.Entry4_1.get()), int(self.Entry4_2.get())))
        self.Button4.grid(column=col, row=2, columnspan=Entry4Span)
        col = col + Entry4Span

        Entry5Span = 2
        self.Entry5_1 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Entry5_1.grid(column=col, row=3)
        self.Entry5_1.insert(0, '0')
        self.Entry5_2 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable="", relief="ridge")
        self.Entry5_2.grid(column=col+1, row=3)
        self.Entry5_2.insert(0, '0')
        self.Button5 = tkinter.Button(self.InputinfoFrame, text='ROI2(Right, Up)',
                                      command=lambda m='ROI2': self.Update_ROI(m, int(self.Entry5_1.get()), int(self.Entry5_2.get())))
        self.Button5.grid(column=col, row=2, columnspan=Entry5Span)
        col = col + Entry5Span

        Entry6Span = 2
        self.Button6 = tkinter.Button(self.InputinfoFrame, text='Show ROI', command=lambda: self.Show_ROI(self.read_data))
        self.Button6.grid(column=col, row=2, columnspan=Entry6Span)
        self.Label6_1_1 = tkinter.Label(self.InputinfoFrame, text='FOI')
        self.Label6_1_1.grid(column=col, row=3)
        self.Label6_2_1 = tkinter.Label(self.InputinfoFrame, text='ROI1')
        self.Label6_2_1.grid(column=col, row=4)
        self.Label6_3_1 = tkinter.Label(self.InputinfoFrame, text='ROI2')
        self.Label6_3_1.grid(column=col, row=5)
        self.Label6_1_2 = tkinter.Label(self.InputinfoFrame)
        self.Label6_1_2.grid(column=col+1, row=3)
        self.Label6_2_2 = tkinter.Label(self.InputinfoFrame)
        self.Label6_2_2.grid(column=col+1, row=4)
        self.Label6_3_2 = tkinter.Label(self.InputinfoFrame)
        self.Label6_3_2.grid(column=col+1, row=5)
        col = col + Entry6Span

        Entry7Span = 2
        self.Label7_1_1 = tkinter.Label(self.InputinfoFrame, text='IFS')
        self.Label7_1_1.grid(column = col, row = 3)
        self.Label7_2_1 = tkinter.Label(self.InputinfoFrame, text='System Gain')
        self.Label7_2_1.grid(column = col, row = 4)
        self.Label7_3_1 = tkinter.Label(self.InputinfoFrame, text='Differential')
        self.Label7_3_1.grid(column = col, row = 5)
        self.Label7_4_1 = tkinter.Label(self.InputinfoFrame, text='Frames')
        self.Label7_4_1.grid(column = col, row = 6)

        self.Entry7_1_2 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable=self.IFS, relief="ridge")
        self.Entry7_1_2.grid(column=col + 1, row=3)
        self.Label7_2_2 = tkinter.Label(self.InputinfoFrame, text='')
        self.Label7_2_2.grid(column = col + 1, row = 4)

        self.CheckButton7_3_2 = tkinter.Checkbutton(self.InputinfoFrame, text="", variable=self.Differential)
        self.CheckButton7_3_2.select()
        self.CheckButton7_3_2.grid(column = col + 1, row = 5)
        self.Label7_4_2 = tkinter.Label(self.InputinfoFrame, text='')
        self.Label7_4_2.grid(column=col + 1, row=6)

        self.Button7 = tkinter.Button(self.InputinfoFrame, text='Configuration', command=lambda: self.Configurations(self.ROI_Data.copy()))
        self.Button7.grid(column=col, row=2, columnspan=Entry7Span)
        col = col + Entry7Span

        Entry8Span = 2
        self.Button8 = tkinter.Button(self.InputinfoFrame, text='Calculate',
                                      command=lambda: self.Calculate(self.ROI_Data.copy(), self.ROIWidget))
        self.Button8.grid(column=col, row=2, columnspan=Entry8Span)
        self.Label8_1_1 = tkinter.Label(self.InputinfoFrame, text='Total Noise')
        self.Label8_1_1.grid(column=col, row=3)
        self.Label8_2_1 = tkinter.Label(self.InputinfoFrame, text='Frame Noise')
        self.Label8_2_1.grid(column=col, row=4)
        self.Label8_3_1 = tkinter.Label(self.InputinfoFrame, text='Line Noise(R)')
        self.Label8_3_1.grid(column=col, row=5)
        self.Label8_4_1 = tkinter.Label(self.InputinfoFrame, text='Line Noise(C)')
        self.Label8_4_1.grid(column=col, row=6)
        self.Label8_5_1 = tkinter.Label(self.InputinfoFrame, text='Pixel Noise')
        self.Label8_5_1.grid(column=col, row=7)

        self.Label8_1_2 = tkinter.Label(self.InputinfoFrame)
        self.Label8_1_2.grid(column=col+1, row=3)
        self.Label8_2_2 = tkinter.Label(self.InputinfoFrame)
        self.Label8_2_2.grid(column=col+1, row=4)
        self.Label8_3_2 = tkinter.Label(self.InputinfoFrame)
        self.Label8_3_2.grid(column=col+1, row=5)
        self.Label8_4_2 = tkinter.Label(self.InputinfoFrame)
        self.Label8_4_2.grid(column=col+1, row=6)
        self.Label8_5_2 = tkinter.Label(self.InputinfoFrame)
        self.Label8_5_2.grid(column=col+1, row=7)

        col = col + Entry8Span

        Entry9Span = 2
        self.Button9 = tkinter.Button(self.InputinfoFrame, text='IQR Remove', command=lambda: self.IQR(self.ROI_Data.copy(),
                                                                                                       self.NIQR.get(),
                                                                                                       self.NIteration.get(),
                                                                                                       self.ROIWidget))
        self.Button9.grid(column=col, row=2, columnspan=Entry9Span)
        self.Label9_1_1 = tkinter.Label(self.InputinfoFrame, text='IQR')
        self.Label9_1_1.grid(column=col, row=3)
        self.Label9_2_1 = tkinter.Label(self.InputinfoFrame, text='Iterations')
        self.Label9_2_1.grid(column=col, row=4)
        self.Label9_3_1 = tkinter.Label(self.InputinfoFrame, text='Excluding 0')
        self.Label9_3_1.grid(column=col, row=5)
        self.Label9_4_1 = tkinter.Label(self.InputinfoFrame, text='HighPass Filter')
        self.Label9_4_1.grid(column=col, row=6)

        self.Entry9_1_2 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable=self.NIQR, relief="ridge")
        self.Entry9_1_2.grid(column=col + 1, row=3)
        self.Entry9_2_2 = tkinter.Entry(self.InputinfoFrame, width=4, textvariable=self.NIteration, relief="ridge")
        self.Entry9_2_2.grid(column=col + 1, row=4)
        self.CheckButton9_3_2 = tkinter.Checkbutton(self.InputinfoFrame, text="", variable=self.ExcludingZero)
        self.CheckButton9_3_2.grid(column = col + 1, row = 5)
        self.CheckButton9_4_2 = tkinter.Checkbutton(self.InputinfoFrame, text="", variable=self.HPF)
        self.CheckButton9_4_2.grid(column=col + 1, row=6)
        self.CheckButton9_4_2.select()

        col = col + Entry9Span

        Entry10Span = 2
        self.SaveButton = tkinter.Button(self.InputinfoFrame, text='Save Files', command=self.SaveBTNEvent)
        self.SaveButton.grid(column = col, row=2)
        self.SavePath = tkinter.Label(self.InputinfoFrame)
        self.SavePath.grid(column = col, row=1)

        self.SaveClipboardBoardBTN = tkinter.Button(self.InputinfoFrame, text='Save Clipboard')
        self.SaveClipboardBoardBTN.grid(column=col+1, row=2)

        self.Label10_1_1 = tkinter.Label(self.InputinfoFrame, text='Total Noise')
        self.Label10_1_1.grid(column=col, row=3)
        self.Label10_2_1 = tkinter.Label(self.InputinfoFrame, text='Frame Noise')
        self.Label10_2_1.grid(column=col, row=4)
        self.Label10_3_1 = tkinter.Label(self.InputinfoFrame, text='Line Noise(R)')
        self.Label10_3_1.grid(column=col, row=5)
        self.Label10_4_1 = tkinter.Label(self.InputinfoFrame, text='Line Noise(C)')
        self.Label10_4_1.grid(column=col, row=6)
        self.Label10_5_1 = tkinter.Label(self.InputinfoFrame, text='Pixel Noise')
        self.Label10_5_1.grid(column=col, row=7)

        self.Label10_1_2 = tkinter.Label(self.InputinfoFrame)
        self.Label10_1_2.grid(column=col+1, row=3)
        self.Label10_2_2 = tkinter.Label(self.InputinfoFrame)
        self.Label10_2_2.grid(column=col+1, row=4)
        self.Label10_3_2 = tkinter.Label(self.InputinfoFrame)
        self.Label10_3_2.grid(column=col+1, row=5)
        self.Label10_4_2 = tkinter.Label(self.InputinfoFrame)
        self.Label10_4_2.grid(column=col+1, row=6)
        self.Label10_5_2 = tkinter.Label(self.InputinfoFrame)
        self.Label10_5_2.grid(column=col+1, row=7)



if __name__ == '__main__':
    window = tkinter.Tk()
    TFTDarkLeakageAnalysis(window)
    window.mainloop()