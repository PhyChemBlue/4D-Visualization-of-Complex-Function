# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 04:08:43 2024

@author: PhyCh
"""

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.ttk as ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys

class GUI4D:
    def __init__(self):
        self.wd = tk.Tk()
        self.wd.title("4D Visualize")
        self.wd.geometry("1000x720")
        self.wd.resizable(False, False)
        self.canvas = tk.Canvas()
        self.inipic()
        self.btn()
        self.wd.protocol("WM_DELETE_WINDOW", self.wdcls)
        self.wd.mainloop()
    
    def inipic(self):
        #default
        self.fn = "z^2"
        self.co = "xy"
        self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
        self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        self.dtheta = 15
        self.f = lambda x: x**2
        self.normmu = (self.wmin + self.wmax) / 2
        self.normsigma = self.wmax - self.normmu
        self.x, self.y, self.z, self.w = self.xyzw()
        self.drawpic()
    
    def xyzw(self):
        #axis
        if self.co == 'xy':
            x = np.linspace(self.axmin, self.axmax, 100)
            y = np.linspace(self.aymin, self.aymax, 100)
            x, y = np.meshgrid(x, y)
        elif self.co == 'rt':
            r = np.linspace(self.axmin, self.axmax, 100)
            t = np.linspace(self.aymin, self.aymax, 100)
            r, t = np.meshgrid(r, t)
            x = r * np.cos(t)
            y = r * np.sin(t)
        a = x + 1j * y
        b = self.f(a)
        z = b.real
        w = b.imag
        return (x, y, z, w)
    
    def rtt(self, ax1, ax2):
        #rotate
        theta = self.dtheta * np.pi / 180
        axa = ax1 * np.cos(theta) + ax2 * np.sin(theta)
        axb = ax2 * np.cos(theta) - ax1 * np.sin(theta)
        return (axa, axb)
    
    def wclr(self, w1):
        #w -> facecolor
        w1 = (w1 - self.normmu) / self.normsigma
        w2 = np.zeros((len(w1), len(w1[0])), tuple)
        for k in range(len(w1)): #y
            for l in range(len(w1[0])): #x
                if w1[k][l] < -0.5 and w1[k][l] >= -1:
                    w2[k][l] = (0.0, 2 * w1[k][l] + 2, 1.0, 1.0)
                elif w1[k][l] < 0 and w1[k][l] >= -0.5:
                    w2[k][l] = (0.0, 1.0, -2 * w1[k][l], 1.0)
                elif w1[k][l] < 0.5 and w1[k][l] >= 0:
                    w2[k][l] = (2 * w1[k][l], 1.0, 0.0, 1.0)
                elif w1[k][l] <= 1 and w1[k][l] >= 0.5:
                    w2[k][l] = (1.0, 2 - 2 * w1[k][l], 0.0, 1.0)
                #<-1 or >1
                elif w1[k][l] < -1:
                    w2[k][l] = (0.0, 0.0, -1.0 / w1[k][l], 1.0)
                elif w1[k][l] > 1:
                    w2[k][l] = (1.0 / w1[k][l], 0.0, 0.0, 1.0)
                else:
                    w2[k][l] = (1.0, 1.0, 1.0, 1.0)
        return w2
    
    def pltcplxfun(self, wc):
        #plot1
        fig = plt.figure(figsize = (9.6, 4.8)) #960, 480
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.plot_surface(self.x, self.y, self.z, facecolors = wc)
        ax.view_init(elev = 30, azim = 300)
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)
        ax.set_zlim(self.zmin, self.zmax)
        #plot2
        ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax.plot_surface(self.x, self.y, self.z, facecolors = wc)
        ax.view_init(elev = 30, azim = 303)
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)
        ax.set_zlim(self.zmin, self.zmax)
        return fig
    
    def drawpic(self):
        plt.close('all')
        wc = self.wclr(self.w)
        fig = self.pltcplxfun(wc)
        self.canvas = FigureCanvasTkAgg(fig, self.wd)
        self.canvas.draw()
        self.canvas.get_tk_widget().place(x=20, y=20)
    
    def restoreDraw(self):
        #revert
        self.funsel()
        self.normmu = (self.wmin + self.wmax) / 2
        self.normsigma = self.wmax - self.normmu
        self.x, self.y, self.z, self.w = self.xyzw()
        self.drawpic()
    
    def funsel(self):
        #select function
        if self.fn == "z":
            self.f = lambda x: x
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -1, 1, -1, 1, -1, 1, -1, 1
        elif self.fn == "z^2":
            self.f = lambda x: x**2
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "z^3":
            self.f = lambda x: x**3
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "1/z":
            self.f = lambda x: 1/x
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -2, 2, -2, 2
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0.01, 2, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "1/(z^2+1)":
            self.f = lambda x: 1 / (x**2 + 1)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -2, 2, -2, 2
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 2, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "sqrt(z)":
            self.f = lambda x: x**0.5
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -1, 1, -1, 1, -1, 1, -1, 1
        elif self.fn == "z^i":
            self.f = lambda x: x**1j
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "exp(x)":
            self.f = lambda x: np.exp(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -8, 2, -2 * np.pi, 2 * np.pi
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 2, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -8, 8, -8, 8, -8, 8, -8, 8
        elif self.fn == "ln(z)":
            self.f = lambda x: np.log(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0.01, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "sin(z)":
            self.f = lambda x: np.sin(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -1, 1, -1, 1, -1, 1, -1, 1
        elif self.fn == "cos(z)":
            self.f = lambda x: np.cos(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "tan(z)":
            self.f = lambda x: np.tan(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -1, 1, -1, 1, -1, 1, -1, 1
        elif self.fn == "arcsin(z)":
            self.f = lambda x: np.arcsin(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -1, 1, -1, 1, -1, 1, -1, 1
        elif self.fn == "arccos(z)":
            self.f = lambda x: np.arccos(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -2, 2, -2, 2, -2, 2, -2, 2
        elif self.fn == "arctan(z)":
            self.f = lambda x: np.arctan(x)
            if self.co == 'xy':
                self.axmin, self.axmax, self.aymin, self.aymax = -1, 1, -1, 1
            elif self.co == 'rt':
                self.axmin, self.axmax, self.aymin, self.aymax = 0, 1, -np.pi, np.pi
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax, self.wmin, self.wmax = -1, 1, -1, 1, -1, 1, -1, 1
    
    def btn(self):
        #Rotate button
        self.bxy = tk.Button(self.wd, text="x->y", font=("Arial", 16), command=self.rtxy)
        self.bxy.place(x=160, y=560, width=120, height=40)
        self.byz = tk.Button(self.wd, text="y->z", font=("Arial", 16), command=self.rtyz)
        self.byz.place(x=300, y=560, width=120, height=40)
        self.bxz = tk.Button(self.wd, text="x->z", font=("Arial", 16), command=self.rtxz)
        self.bxz.place(x=440, y=560, width=120, height=40)
        self.bxw = tk.Button(self.wd, text="x->w", font=("Arial", 16), command=self.rtxw)
        self.bxw.place(x=580, y=560, width=120, height=40)
        self.byw = tk.Button(self.wd, text="y->w", font=("Arial", 16), command=self.rtyw)
        self.byw.place(x=720, y=560, width=120, height=40)
        self.bzw = tk.Button(self.wd, text="z->w", font=("Arial", 16), command=self.rtzw)
        self.bzw.place(x=860, y=560, width=120, height=40)
        self.bres = tk.Button(self.wd, text="Restore", font=("Arial", 16), command=self.rtres)
        self.bres.place(x=860, y=625, width=120, height=75)
        #Combobox
        self.fsv = tk.StringVar()
        self.cbfs = ttk.Combobox(self.wd, state="readonly", font=("Arial", 15), textvariable=self.fsv)
        self.cbfs['value'] = ("z", "z^2", "z^3", "1/z", "1/(z^2+1)", "sqrt(z)", "z^i", "exp(x)", "ln(z)", "sin(z)", "cos(z)", "tan(z)", "arcsin(z)", "arccos(z)", "arctan(z)")
        self.cbfs.current(1)
        self.cbfs.bind("<<ComboboxSelected>>", self.fscb)
        self.cbfs.place(x=20, y=660, width=120, height=40)
        self.asv = tk.StringVar()
        self.cbas = ttk.Combobox(self.wd, state="readonly", font=("Arial", 15), textvariable=self.asv)
        self.cbas['value'] = ("1°", "5°", "15°", "30°", "45°", "90°", "180°", "-90°", "-45°", "-30°", "-15°", "-5°", "-1°")
        self.cbas.current(2)
        self.cbas.bind("<<ComboboxSelected>>", self.ascb)
        self.cbas.place(x=20, y=560, width=120, height=40)
        self.cov = tk.StringVar()
        self.cbco = ttk.Combobox(self.wd, state="readonly", font=("Arial", 15), textvariable=self.cov)
        self.cbco['value'] = ("Cartesian", "Polar")
        self.cbco.current(0)
        self.cbco.bind("<<ComboboxSelected>>", self.cocb)
        self.cbco.place(x=160, y=660, width=120, height=40)
        #Zoom button
        self.bbg = tk.Button(self.wd, text="zoom in", font=("Arial", 16), command=self.zmbg)
        self.bbg.place(x=300, y=660, width=120, height=40)
        self.bsm = tk.Button(self.wd, text="zoom out", font=("Arial", 16), command=self.zmsm)
        self.bsm.place(x=440, y=660, width=120, height=40)
        self.bdb = tk.Button(self.wd, text="domain +", font=("Arial", 16), command=self.dmbg)
        self.bdb.place(x=580, y=660, width=120, height=40)
        self.bds = tk.Button(self.wd, text="domain -", font=("Arial", 16), command=self.dmsm)
        self.bds.place(x=720, y=660, width=120, height=40)
        #Label
        self.ldt = tk.Label(self.wd, text="Rotate Δθ", font=("Arial", 15))
        self.ldt.place(x=20, y=525, width=120, height=30)
        self.l3d = tk.Label(self.wd, text="3D Rotation", font=("Arial", 15), relief="groove")
        self.l3d.place(x=160, y=525, width=400, height=30)
        self.l4d = tk.Label(self.wd, text="4D Rotation", font=("Arial", 15), relief="groove")
        self.l4d.place(x=580, y=525, width=400, height=30)
        self.lfn = tk.Label(self.wd, text="f(z)=", font=("Arial", 15))
        self.lfn.place(x=20, y=625, width=120, height=30)
        self.lco = tk.Label(self.wd, text="Coordinate", font=("Arial", 15))
        self.lco.place(x=160, y=625, width=120, height=30)
        self.lzm = tk.Label(self.wd, text="Camera", font=("Arial", 15), relief="groove")
        self.lzm.place(x=300, y=625, width=260, height=30)
        self.ldm = tk.Label(self.wd, text="Plot Domain", font=("Arial", 15), relief="groove")
        self.ldm.place(x=580, y=625, width=260, height=30)
    
    def rtxy(self):
        self.x, self.y = self.rtt(self.x, self.y)
        self.drawpic()
    def rtyz(self):
        self.y, self.z = self.rtt(self.y, self.z)
        self.drawpic()
    def rtxz(self):
        self.x, self.z = self.rtt(self.x, self.z)
        self.drawpic()
    def rtxw(self):
        self.x, self.w = self.rtt(self.x, self.w)
        self.drawpic()
    def rtyw(self):
        self.y, self.w = self.rtt(self.y, self.w)
        self.drawpic()
    def rtzw(self):
        self.z, self.w = self.rtt(self.z, self.w)
        self.drawpic()
    def rtres(self):
        self.restoreDraw()
    
    def zmbg(self):
        self.xmin /= 2
        self.xmax /= 2
        self.ymin /= 2
        self.ymax /= 2
        self.zmin /= 2
        self.zmax /= 2
        self.wmin /= 2
        self.wmax /= 2
        self.normmu = (self.wmin + self.wmax) / 2
        self.normsigma = self.wmax - self.normmu
        self.drawpic()
    def zmsm(self):
        self.xmin *= 2
        self.xmax *= 2
        self.ymin *= 2
        self.ymax *= 2
        self.zmin *= 2
        self.zmax *= 2
        self.wmin *= 2
        self.wmax *= 2
        self.normmu = (self.wmin + self.wmax) / 2
        self.normsigma = self.wmax - self.normmu
        self.drawpic()
    def dmbg(self):
        self.axmin *= 2
        self.axmax *= 2
        self.aymin *= 2
        self.aymax *= 2
        self.x, self.y, self.z, self.w = self.xyzw()
        self.drawpic()
    def dmsm(self):
        self.axmin /= 2
        self.axmax /= 2
        self.aymin /= 2
        self.aymax /= 2
        self.x, self.y, self.z, self.w = self.xyzw()
        self.drawpic()
    
    def fscb(self, event):
        self.fn = self.fsv.get()
        self.restoreDraw()
    def ascb(self, event):
        self.an = self.asv.get()
        self.dtheta = int(self.an[:-1])
    def cocb(self, event):
        self.con = self.cov.get()
        if self.con == "Cartesian":
            self.co = "xy"
        elif self.con == "Polar":
            self.co = "rt"
        self.restoreDraw()
    
    def wdcls(self):
        #quit
        self.wd.destroy()
        sys.exit(0)

if __name__=="__main__":
    NDGUI = GUI4D()