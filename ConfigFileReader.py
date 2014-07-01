

class ConfigFileReader:
    def __init__(self,fname):
        self.fname = fname
        self.mydict={}
        self.readfile()
        
        
    def cleanup(self,astring) :
        b = astring.replace('\n','')
        c= b.replace('\t','')
        d = c.rstrip(' ').lstrip(' ')
        e = d.upper()
        return e
        
        
    def readfile(self):
        #try:
        if 1==1:
            myf = file(self.fname,'r')
            for lines in myf:
                if lines[0]=='#':
                    continue
                validstr = lines.split('#')[0]
                if len(validstr)<3:
                    continue
                svstr = validstr.split('=')
                if len(svstr)!=2:
                    continue
                kk = svstr[0]
                vv =  svstr[1]
                kk =self.cleanup(kk)
                vv =self.cleanup(vv)
                self.mydict[kk] = vv
            myf.close()
        #except:
         #   print "I can't read ",self.fname

    def parseBool(self,instr):
        a = instr.upper()
        flist = ['F','FALSE','NONE']
        if a in flist:
            return False
        else:
            return True

    def getRotFreq(self):
        kk ='Rotation Frequency'.upper()
        kk = kk.lstrip(' ').rstrip(' ')
        def_value = 9.9
        if kk in self.mydict:
            return float(self.mydict[kk])
        else:
            return def_value
            
    def getDebugMode(self):
        kk ='Debug Mode'.upper()
        kk = kk.lstrip(' ').rstrip(' ')
        def_value = False
        if kk in self.mydict:
            return self.parseBool(self.mydict[kk])
        else:
            return def_value
            
    def getCalibration (self):
         cc = 'New Calibration Value'.upper()
         cc = cc.lstrip(' ').rstrip(' ')
         def_value = 4.075
         if cc in self.mydict:
             return float(self.mydict[cc])
         else:
             return def_value
             
    def SamplingFreq(self):
        ss = 'ADC Sampling Frequency'.upper()
        ss = ss.lstrip(' ').rstrip(' ')
        def_value = 100
        if ss in self.mydict:
            return float(self.mydict[ss])
        else: 
            return def_value
            
    def FilterOrder(self):
        fo = 'Butterworth Filter Order for EMF'.upper()
        fc = 'Butterworth Filter Order for Capacitor'.upper()
        fo = fo.lstrip(' ').rstrip(' ')
        fc = fc.lstrip(' ').rstrip(' ')
        def_value = 6
        if fo and fc in self.mydict:
            return int (self.mydict[fo]), int (self.mydict[fc])
        else:
            return def_value
            
    def getComPort(self):
        cp = 'New Serial Port'.upper()
        cp = cp.lstrip(' ').rstrip(' ')
        def_value = 'COM4'
        if cp in self.mydict:
            return str((self.mydict[cp]))
        else: 
            return (def_value)
            
    def getCalMode(self):
        cm = 'Calibration Mode'.upper()
        def_value = False
        if cm in self.mydict:
            return self.parseBool(self.mydict[cm])
        else :
            return def_value
            
    def getDebugFile(self):
        
        ff = "Filename".upper()
        def_value = 'VoltageSamples_03-07-13_153318.dat' 
        if ff in self.mydict:
            return str(self.mydict[ff])
        else: 
            return (def_value)
            
    def NewDataPlotLimit(self):
        
        nd = "New Data Plot Limit".upper()
        def_value = 10
        if nd in self.mydict:
            return int (self.mydict[nd])
        else:
            return def_value
    
    def GetPlotWindowLength(self):
        pl = "Window Length for sub-plots".upper()
        def_value = 10
        if pl in self.mydict:
            return int(self.mydict[pl])
        else:
            return def_value
            
    def GetNoOfMeas(self):
        nm = "No of Measurements".upper()
        def_value = 1000
        if nm in self.mydict:
            return int(self.mydict[nm])
        else:
            return def_value