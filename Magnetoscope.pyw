# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 03:11:55 2013

@author: Ali
"""

from ConfigFileReader import ConfigFileReader
import datetime
import sys
import matplotlib
import matplotlib.pylab as py
from matplotlib import pyplot as plt
from matplotlib.widgets import Button
from matplotlib.text import Text
from PyQt4 import QtGui
from PyQt4 import QtCore
from matplotlib.ticker import FuncFormatter
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
import myfilters
import math
import serial
import time
import os


global otimestamp
global tstamp
global fcoeffs
global N                       
global Area                   
global axfreq


otimestamp=0
tstamp=0

class LivePlot(FigureCanvas):
    """
    Matplotlib widget to plot the raw data
    """
    def __init__(self,parent,updatefunc):
        
        myf = "Probe.conf"
        mycfg = ConfigFileReader(myf)
        self.lockInnF = mycfg.getRotFreq()
        self.fsamp = mycfg.SamplingFreq()
        self.slope = mycfg.getCalibration()
        com = mycfg.getComPort()
        self.Dbug = mycfg.getDebugMode()
        self.newdataplotlimit = mycfg.NewDataPlotLimit()
        self.updatefunc = updatefunc
        self.xmax = mycfg.GetPlotWindowLength()
        self.MyFile = mycfg.getDebugFile()
        self.measurement = mycfg.GetNoOfMeas()
        
        self.newdata = 0
        self.oVphase=-999999
        self.wrapper=0
        self.wrapper1 = 0
        self.oCphase = -999999
        self.lastzeromark=0
        self.thiszeromark=0
        self.frequencyvalid=False
        self.frequencyvalid=True
        self.oldV=0
        self.linelength=44
        self.filtfreq = 0
        self.f = 0
        self.NB = 0
        
        self.a11=0
        self.a12=0
        self.a22=0
        self.b1=0
        self.b2=0
        self.a=0
        self.b=0
        self.d=0
        self.meanB = 0
        self.reset_mean = True

        

        
        self.cofs=myfilters.gbutter(0.5,self.fsamp,6)   
        self.myfilter0 = myfilters.gIIR(self.cofs)
        self.myfilter1 = myfilters.gIIR(self.cofs)
        self.myfilter2 = myfilters.gIIR(self.cofs)
        self.myfilter3 = myfilters.gIIR(self.cofs)
        
        # Phase determined cap filter
        self.myfilter4 = myfilters.gIIR(self.cofs)
        self.myfilter5 = myfilters.gIIR(self.cofs)
        self.cofs2=myfilters.gbutter(self.lockInnF/5,self.lockInnF,8)
        self.myfilter6 = myfilters.gIIR(self.cofs2)
        
        mylabels = ['V (mV)','Vc (V)','Cap-amp (mV)',\
            'Vamp (mV)','d (in)','f (Hz)',\
            'Vphi','Cphi','B (mT)']
        
        # Define the color scheme for the subplots
        mycolor = ['b-','r-','k-',
                 'b-','r-','k-',
                 'b-','r-','k-']
        mycolor2 =['b-','r-','k-',
                 'b-','r-','k-',
                 'b-','r-','k-']
                 
        # For annotation of mean and std deviation of plot
        self.myann=[]      
        self.ax=[]
        self.lines=[]
        self.lines2=[]
        
        self.fig = Figure()
        
        self.fig.subplots_adjust(left=0.085, right =0.95, wspace=0.5, hspace =0.3)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        
        self.canvas = FigureCanvas.__init__(self,self.fig)
        
        self.setParent(parent)
        self.main_frame = QtGui.QWidget()
        
        self.mainlayout = QtGui.QVBoxLayout()
        self.mainlayout.addWidget(self.canvas)
        
        self.button = QtGui.QPushButton('Update Frequency (f)')
        self.mainlayout.addWidget(self.button)
        self.button.setMaximumSize(200, 200)
        
        self.mainlayout.addStretch()
        self.button.setShortcut(0x46)
        self.button.clicked.connect(self.change_freq) 
        
        
        self.button1 = QtGui.QPushButton('Capture Data and Reset Mean (D)')
        self.mainlayout.addWidget(self.button1)
        self.button1.setMaximumSize(200, 200)
        self.mainlayout.addSpacing(1000)
        self.button1.setShortcut(0x44)
        self.button1.clicked.connect(self.BMean)
        
        self.buttonq = QtGui.QPushButton('Quit (q)')
        self.mainlayout.addWidget(self.buttonq)  
        self.buttonq.setMaximumSize(100, 100)
        self.buttonq.setShortcut(0x51)
        self.buttonq.clicked.connect(self.Quit)
        
        self.setLayout(self.mainlayout)
        
        # Define plot parameters.
        for i in range(9):
            # 3 rows 3 column ith plot
            self.ax.append(self.fig.add_subplot(3,3,i+1))   
            # show the last 20s(VoltageSamples_03-07-13_153318xmax defines the length of window for plot)
            self.ax[-1].set_xlim(0,self.xmax)   
            # Enable autoscale for both axis
            self.ax[-1].autoscale(enable=True, axis='both', tight=False)
            # Format the yaxis data such that the data appears on lhs
            self.ax[-1].yaxis.set_major_formatter(y_formatter)
            # Set the y-axis labels
            self.ax[-1].set_ylabel(mylabels[i])
            # Set the x-axis labels
            self.ax[-1].set_xlabel('time(s)')
            # Annotate top center, the mean and the dev of the plots
            an = self.ax[-1].annotate('Hi Fish',xy=(0.3,1.0),xycoords='axes fraction')
            self.myann.append(an)    
        
            newline, = self.ax[-1].plot([],[],mycolor[i],label=mylabels[i])
            self.lines.append(newline)
            
            newline, = self.ax[-1].plot([],[],mycolor2[i])
            self.lines2.append(newline)
       
        
        self.time=[]    # Collect the sampling time data from COM port
        self.Vdata=[]   # Collect EMF voltage data from COM port
        self.Cdata=[]   # Collect Cap-voltage data from COM port
        self.Rdata=[]   # Collect rotation data from zero crossings of Voltage data
        self.VAMP=[]    # EMF filtered amplitude
        self.VPHS=[]    # EMF filtered phase
        self.CAMP=[]    # Cap voltage filtered amplitude
        self.CPHS=[]    # Cap voltage filtered phase
        self.BAMP=[]    # B-field average amplitude
        self.CAMP1=[]   # Average Cap amplitude derived by using initial phase from voltage data
        self.freq=[]    # Get frequency of rotationfrom PyQt4 import QtCore
        self.Inphase=[] 
        self.Inquad=[]
        self.sampling=[]
        
        if (self.Dbug == True):
            print "Magnetoscope is in DEBUG Mode..."
           
            self.myfile = open (self.MyFile,'r')
            self.filedata=[]
            for lines in self.myfile:
                self.filedata.append(lines)
            self.myfile.close()
            self.recordNo=0
            self.timer = self.startTimer(1)               
        else:
            print "Starting Magnetoscope..."
            self.ser=serial.Serial(port=com,
                      baudrate=115200,
                      bytesize=serial.EIGHTBITS,
                      parity=serial.PARITY_NONE,
                      stopbits=serial.STOPBITS_ONE)

            self.ser.timeout = 3                      
            self.timer = self.startTimer(100)              
        self.fig.canvas.draw()  
        self.timerEvent(None)  
        
        
    def getdata(self):
        
        global otimestamp   # old timestamp 
        global tstamp       # new timestamp
      
        if (self.Dbug == True):
            #for fileline in self.myfile:
            i= self.filedata[self.recordNo].split()
            #print i
            self.recordNo+=1

                    
        else:
                samples = self.ser.readline()  
                i = samples.split()
        if (len(i)==3): 

            tstamp = otimestamp +float(i[0])*0.001
            otimestamp = tstamp
            # Error Check for Corrupted Data
            try:
                v_data =float(i[1]) # get voltage data
            except ValueError:
                print 'Corrupted Data'
                exit()
            try:
                c_data = float(i[2]) # get cap voltage data
            except ValueError:
                print 'Corrupted Data'
                exit()
           
            samp_data =float(i[0])  # sampling data
            
            return tstamp,v_data,c_data,samp_data
        # if the sample lengself.posbox.textC- Positionhanged.connect(self.CurrPos)th is not 3 then return 0's from PyQt4 import QtCore
        return 0,0,0,0
        
    def analysis(self,Volt,Time):
        ip = float(Volt)*math.cos(Time*self.lockInnF*2*np.pi)  
        iq = float(Volt)*math.sin(Time*self.lockInnF*2*np.pi)  
        return (ip, iq)
        
    def analysisC(self, CV,Time,vph):
        # Cap amplitude is calculated from the initial phase data of EMF plot
        # cap voltage leads emf by 90
        voffs=-0.087
        iq = float(CV)*math.sin(Time*self.lockInnF*2*np.pi + vph+voffs + math.pi/2)  
        return iq
    
    # Now update the frequency...        
    def update_freq(self,Time,Rot):
        
        if Rot==0:
            self.thiszeromark = Time
            # Calculate the new rotation frequency        
            frequen = 1.0/(-self.lastzeromark+self.thiszeromark)
            for count in range(1000):
               x= self.myfilter6.filtered(frequen)
            self.f = x
            self.lastzeromark=self.thiszeromark
        return self.f

    def change_freq(self):
        
        self.lockInnF = self.f   
        print self.lockInnF
         
    def nonzero(self,*arg):
        for i in arg:
            if i!=0: return True
        return False
    
    # DC average after filtering
    def Amp_Phase(self,fp,fq):
        Vaf =2.0*( np.sqrt(fp**2 +fq**2))
        # Calculate phase 
        Vpf = np.arctan2(fp, fq)
        
        return Vaf,Vpf
        
    def Mag_Amp(self,V,f):
        global N,Area
        N = 59 
        Area = 1.0*1.125*2.54e-2*2.54e-2   
        # V = N*A*B*2*pi*f
        BField1 = V/(N*Area*2*np.pi*f)
        
        return(BField1)
        
    def BMean(self):
        self.reset_mean = True
       
    def Quit(self):
        
        reply = QtGui.QMessageBox.question(self, 'Message',
            "Are you sure?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
             print "Quitting Application..."
             if (self.Dbug == True):
                 app.quit()
                 sys.exit(app.exec_())
             else:
                 self.ser.flush()
                 self.ser.close()
                 app.quit()
                 sys.exit(app.exec_())
        
       
    
    def timerEvent(self,evt):
        """
        This code will be called, whenever the timer triggers,
        don't rename this method
        """
        # Get the data every time the timer triggers
        #newdata=0
        if self.Dbug==True:
            haveData=True
        else:
            haveData=False
        while haveData==True or \
            (( self.Dbug==False) and (self.ser.inWaiting()>self.linelength)):
            t,V,C,s=self.getdata()
            haveData=False
            if np.sign(V)<np.sign(self.oldV):
                r=0
            else:
                r=1
                self.oldV=V 
            if self.nonzero(t,V,C):
                
                IP,IQ   = self.analysis(V,t)
                IPC,IQC = self.analysis(C,t)
                Vamp,Vphase = self.Amp_Phase(self.myfilter0.filtered(IP),\
                                            self.myfilter1.filtered(IQ))
                Bfield = self.Mag_Amp(Vamp,self.lockInnF)
                Camp,Cphase = self.Amp_Phase(self.myfilter2.filtered(IPC),
                                      self.myfilter3.filtered(IQC))
                Vphase = Vphase + self.wrapper*2*math.pi               
                Vphase_uw = Vphase
                if self.oVphase==-999999: self.oVphase=Vphase
                if Vphase-self.oVphase>math.pi:
                    Vphase=Vphase-2*math.pi
                    self.wrapper = self.wrapper -1
                elif self.oVphase-Vphase>math.pi:
                    Vphase=Vphase+2*math.pi
                    self.wrapper = self.wrapper +1                     
                self.oVphase=Vphase
                
                Cphase = Cphase + self.wrapper1*2*math.pi               

                if self.oCphase==-999999: self.oCphase=Cphase
                if Cphase-self.oCphase>math.pi:
                    Cphase=Cphase-2*math.pi
                    self.wrapper1 = self.wrapper1 -1
                elif self.oCphase-Cphase>math.pi:
                    Cphase=Cphase+2*math.pi
                    self.wrapper1 = self.wrapper1 +1                     
                self.oCphase=Cphase
                
                IA = self.analysisC(C,t,Vphase_uw)
                Camp1=2.0*self.myfilter4.filtered(IA)
                
                dist = Camp1/self.slope
                self.Vdata.append(V)
                self.Cdata.append(C)
                self.time.append(t)
                self.Inphase.append(IP)
                self.Inquad.append(IQ)
                self.Rdata.append(r)
                self.VAMP.append(Vamp)
                self.VPHS.append(Vphase)
                self.CAMP.append(dist)
                self.CPHS.append(Cphase)
                self.BAMP.append(Bfield)
                self.CAMP1.append(Camp1*1000)
                self.sampling.append(s)
                self.freq.append(self.f)
                self.newdata+=1
               
                # Clear the buffer                 
                if len(self.Vdata)>(self.xmax*100):
                    self.Vdata.pop(0)
                    self.time.pop(0)
                    self.Cdata.pop(0)
                    self.Rdata.pop(0)
                    self.VAMP.pop(0)
                    self.VPHS.pop(0)
                    self.CAMP.pop(0)
                    self.CPHS.pop(0)
                    self.freq.pop(0)
                    self.BAMP.pop(0)
                    self.CAMP1.pop(0)
                    self.sampling.pop(0)
                                    
                if self.reset_mean==True:
                    if self.NB>0 :
                        self.updatefunc(self.tB,self.t2B,self.NB\
                        ,self.tC,self.t2C,self.NC,True)
                    self.tB = 0
                    self.t2B = 0
                    self.NB = 0
                    self.tC = 0
                    self.t2C = 0
                    self.NC = 0
                    self.FB=[]
                    self.CA=[] 
                    self.reset_mean=False
                    
                    
                self.FB.append(Bfield)      
                self.CA.append(Camp1*1000)
                
                self.NB = self.NB + 1
                self.tB = self.tB + Bfield
                self.t2B = self.t2B + Bfield*Bfield
                
                self.NC = self.NC + 1
                self.tC =self.tC + Camp1*1000
                self.t2C =self.t2C + Camp1*Camp1*1e6
                
                while  self.NB>self.measurement:
                    self.tB = self.tB - self.FB[0]
                    self.t2B = self.t2B - self.FB[0]*self.FB[0]
                    self.NB = self.NB - 1
                    self.FB.pop(0)
                      
                #self.NB = len(self.FB)
                    
                while self.NC>1000:
                    self.tC = self.tC - self.CA[0]
                    self.t2C = self.t2C - self.CA[0]*self.CA[0]
                    self.NC = self.NC-1
                    self.CA.pop(0)
                    
        self.a11=0
        self.a12=0
        self.a22=0
        self.b1=0
        self.b2=0
        a=0
        b=0
        t0=self.time[0]
        for tn,Vphase in zip(self.time,self.VPHS):
            t=tn-t0
            self.a11+=t*t
            self.a12+=1*t
            self.a22+=1
            self.b1+=t*Vphase
            self.b2+=1*Vphase
        
        det = self.a11*self.a22-self.a12*self.a12
        if det!=0:
            b = 1.0/det * (self.a22*self.b1-self.a12*self.b2)
            a = 1.0/det * (-self.b1*self.a12+self.b2*self.a11)
        
        vphasefit=[]
        for t,phase in zip(self.time,self.VPHS):
            vphasefit.append(a+b*(t-t0))
        self.f=self.myfilter6.filtered(self.lockInnF+b/(2*math.pi))
     
        Data=zip(*zip(self.Vdata,self.Cdata,self.CAMP1,\
             self.VAMP,self.CAMP,self.freq,self.VPHS,self.CPHS,self.BAMP))
   
        mymeans=[]      
        mystds=[]      
        
        if self.newdata>self.newdataplotlimit:  
            self.updatefunc(self.tB,self.t2B,self.NB,self.tC,self.t2C,self.NC,False)
            self.newdata=0
            for i in range(len(Data)):
                self.ax[i].set_xlim(self.time[0],self.time[0]+self.xmax)
                if len(Data[i])>2:
                    mymin = min(Data[i])
                    mymax = max(Data[i])
                    mymeans.append(np.mean(Data[i]))
                    mystds.append(np.std(Data[i]))
                    
                    if mymax>mymin:
                        self.ax[i].set_ylim(mymin,mymax)
                        
                    #if (i==0):self.mlayout.setFocus() 
                    ll =len(self.time)
                    tets=range(10)
                    self.lines[i].set_data(tets,tets)
                    self.lines[i].set_data(self.time[0:ll],Data[i][0:ll])
                    if i==6:
                        self.lines2[i].set_data(self.time[0:ll],vphasefit[0:ll]) #this is not the problem
                        mystr="slope: {0:8.5f}".format(b)
                    else:
                       self.lines2[i].set_data([self.time[0],self.time[0]+self.xmax],[mymeans[i], mymeans[i]])
                       mystr="{0:8.5f} ({1:6.5f})".format(mymeans[i],mystds[i])
                    self.myann[i].set_text(mystr)
            self.fig.canvas.draw()
            
    
class StatusInfo(QtGui.QWidget):
        def __init__(self,parent):
            super(StatusInfo,self).__init__(parent)
            self.canvas = QtGui.QWidget()
           
            self.value=0   
          
            #self.setMaximumSize(400,400)
            self.setMaximumWidth(300)
            self.setMinimumWidth(300)
            
            self.sw =QtGui.QLabel("MAGNETOSCOPE V.1.0\nCreators:\nDr. Stephan Schlamminger\
                                    \nDr.Frank Seifert\nAlireza Panna\nLeon Chao")
           
            self.mainlayout = QtGui.QVBoxLayout()
           
            self.mainlayout.addWidget(self.sw)
           
            names = ["00/00/00","00:00:00","No of Measurements",\
                    "0","Mean B (mT)",\
                     "0.00","Std B (mT)", "0.00", \
                     "Cap-Amp (mV)", "0.00","Cap-Std (mV)",\
                    "0.00"]
            ss=20
            ls=48
            fontsizes=[ls,ls,ss,ls,ss,ls,ss,ls,ss,ls,ss,ls]
            self.labels=[]
            
            self.mainlayout.addStretch()
            for na,fs in zip(names,fontsizes):
                self.labels.append( QtGui.QLabel(na))
#                if fs==ls:
                self.labels[-1].setMinimumSize(280,int(fs*1.5))
                self.labels[-1].setMaximumSize(280,int(fs*1.5))
                myFont=QtGui.QFont("Times",fs)
                self.labels[-1].setFont(myFont)
                self.mainlayout.addWidget(self.labels[-1])

            self.lastr=""
            self.myedit=QtGui.QLabel("Press 'D' to capture")
            #self.myedit=QtGui.QText
            self.myedit.setBaseSize(280,100)
            #self.myedit.baseSize(QtCore.QSize(80,100))
            self.myedit.setMaximumSize(280,100)
            self.myedit.setMinimumSize(280,100)
            self.myedit.adjustSize()
            self.mainlayout.addWidget(self.myedit)
            
            self.setMaximumWidth(300)
            self.setMinimumWidth(300)

            #self.mainlayout.addStretch()
            self.setLayout(self.mainlayout)
            
            

        def makestrings(self,L,mymean,mystd):
            if L>3:
                le = L
#                mv = np.mean(myArray)
#                st = np.std(myArself.mlayout.setFocus() ray)
                
                di = int(math.ceil(-math.log10(mystd/math.sqrt(le)))+1)
                mele = math.ceil(math.log10(abs(mymean)))
                if mele<0:
                    mele=2
                if di<0:
                    di=1
                mele=int(mele+1+di)
                if mele<di:
                    mele=int(di+1)
                
                length_str = "{0:4d}".format(le)
                try:
                    mean_str= "{0:{le}.{si}f}".format(mymean,le=mele,si=di)
                except:
                    print "le=",mele," si=",di
                    mean_str="Error"
                try:
                    sigma_str = "{0:{le}.{si}f}".format(mystd,le=di+2,si=di)
                except:
                    print "le=",mele," si=",di
                    sigma_str="Error"
            else:
                length_str="No Data"
                mean_str="No Data"
                sigma_str="No Data"
                
            return length_str,mean_str,sigma_str

            

        def calc_mean_sig(self, S,S2,L):
                
                m = S/L
                V = 1.0/(L-1)*(S2-S*S/L)
                try:
                    Sig= math.sqrt(V)
                except:
                    print V,L,S2,S
                return L,m,Sig

        def updateMe(self,tB,t2B,NB,tC,t2C,NC,writeToFile):
            
                if NC<3 :return
                if NB<3 :return

                L,M,S = self.calc_mean_sig(tB,t2B,NB)            
                s2,s4,s6 = self.makestrings(L,M,S)
                L,M,S = self.calc_mean_sig(tC,t2C,NC)            
                du,s8,s10 = self.makestrings(L,M,S)
            
                now = datetime.datetime.now()
                dastr = now.strftime("%m/%d/%y")
                tistr = now.strftime("%H:%M:%S")
               
                self.labels[0].setText(dastr)   
                self.labels[1].setText(tistr)   
                self.labels[3].setText(s2)   
                self.labels[5].setText(s4)
                self.labels[7].setText(s6)
                self.labels[9].setText(s8)
                self.labels[11].setText(s10)
                
                if writeToFile==True:
                    dastr2 = now.strftime("%m%d%y")
                    fn = "MP_{0}.dat".format(dastr2)
                    if os.path.isfile(fn)==False:
                        myfile = open(fn,'w')
                        header="#Date Time Length B[mT] sigma[mT]  Cap[mV] sigma[mV]\n"
                        myfile.write(header)
                        
                    myfile = open(fn,'a')
                    mystr="{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\n".format(
                        dastr,tistr,s2,\
                        s4,s6,s8, \
                        s10)
                    myfile.write(mystr)
                    myfile.close()
                    mystr="{0} {1} mT {2} mV \n".format(
                        tistr,\
                        s4, s8)
                    self.lastr = self.lastr+mystr
                    mylines =self.lastr.split("\n")
                    test = [i  for i in mylines[-5:]]
                    testr = "\n".join(test)
                    self.myedit.setText(testr)
            
        def updateText(self,value):
            s=str(value)
            self.label1.setText(s)
            
        def add_one(self):
            self.value+=1
            self.updateText(self.value)
            
    
class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        
        self.main_widget = QtGui.QWidget(self) 
        self.main_widget.setAutoFillBackground(True)
        self.main_widget.setPalette(QtGui.QPalette(QtGui.QColor(255,255,255)))
        self.main_widget.setFont(QtGui.QFont("Helvetica", 9.5))
        hbl = QtGui.QHBoxLayout(self.main_widget)
        inf = StatusInfo(self.main_widget)
        qmc = LivePlot(self.main_widget,inf.updateMe)
        hbl.addWidget(qmc) 
        hbl.addWidget(inf) 
        # set the focus on the main widget
        self.main_widget.setFocus() 
        # set the central widget of MainWindow to main_widget 
        self.setCentralWidget(self.main_widget) 
        
app = QtGui.QApplication(sys.argv)   

aw = ApplicationWindow()
aw.show()
#aw.showFullScreen()
sys.exit(app.exec_())
#time.sleep(1)
#ser.flush()
#ser.close()