import serial
import struct
import time

#Should Modify:
#Add method to set display mode
#Add method to set control mode (software vs. pot. analog in)
#Add method to class change serial timeout
#Use calibration line for output instead of PWR_BITS scaling
#Change readHex() to use method at bottom (for 0-padding on returns)
#Provide a build mode which does not print 
    #Imporved time performance and outputs as returned values
    #Can use buildFLG to supress prints and take it in as arg 
    #Can keep printDev() for printing if really needed
    #Throws actual errors instead of print lines

class TLS_Device(object):
    '''Class for controlling the Throlabs TLS001
    ***Device not setting Keys/Intr bits correctly so some items are omitted
        from this code to avoid confusion
        - The output of the device depends solely on the 'enable' bit
    '''
    DELAY       = .05       #Number of seconds to wait after writing a message
    PWR_BITS    = 32767     #Max power in bits (for power-write formatting)
    
    def __init__(self, devnm):
        '''Create serial object and instantiate instance variable
        *devnm should be a string like '/dev/ttyUSB0'
        '''
        #Create Serial Object for Communication (keep closed though)
        self.ser    = serial.Serial()     
        self.ser.baudrate    = 115200
        self.ser.port        = devnm
        self.ser.timeout     = 0.5
        
        #Other Instance Variables
        self.SN     = 'DevNotOpenedYet'     #Device serial number
                #self.SN also serves as flag to check if device has been opened
        self.TY     = 'DevNotOpenedYet'     #Device type
        self.FW     = 'DevNotOpenedYet'     #Device firmware version
        self.MCRT   = 'DevNotOpenedYet'     #Device maximum current (mA)
        self.MPWR   = 'DevNotOpenedYet'     #Device maximum power (mW)
        self.WL     = 'DevNotOpenedYet'     #Device operating wavelength (nm)
        
    def open(self):
        '''Opens connection to device
        Also queries the device to obtain basic information
            This serves to confirm communication
        *Does not reopen device if already open
        '''
        #Open port if not already open
        if self.ser.isOpen():
            print('(SN:%i) is already open' %self.SN)
        else:
            #print('Connecting to : %s...' %self.ser.name)
            self.ser.open()
            #Request Device Information
            self.reqInfo()
        
            #Request Maximum Limits
            self.reqLim()
        
            #print('Device is a   : %s \n'   %self.TY +
            #      'Serial Number : %i \n'   %self.SN +
            #      'Frameware vs. : %s \n'   %self.FW )
              
    def close(self):
        '''Closes the device connection
        '''
        self.ser.close()
        
    def identify(self):
        '''Makes device flash screen and LED for 3 seconds
        Useful for identifying connected device without checking SN
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Send identify command
        self.write('23 02 00 00 50 01')
        time.sleep(3)                   #Wait until identify is complete
        self.readAll()                  #Flush the buffer since identify 
                                            #returns unknown data chain
                
    def write(self, hexMSG, byteAppend = None):
        '''Sends a message to the device
        hexMSG should be the message as a hex string (ex. '05 00 00 00 50 01')
        Can append a byte (or byte array) to the end of hexMSG with byteAppend
        *Data requests using 'write' should be followed by a read
            Otherwise unread items in buffer may cause problems
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
            
        msg = bytearray.fromhex(hexMSG)     #convert to bytes
        
        #Append bytes to end
        if byteAppend != None:
            msg = b''.join((msg, byteAppend))
            
        #Send message using pySerial (serial)
        self.ser.write(msg)
        
    def readAll(self):
        '''Returns the full read buffer
        Also serves as a 'flush' function to clear buffer itself
            Useful for debugging reads to ensure read data is as expected
        Returns the read data as raw bytes in bytearray
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        time.sleep(self.DELAY)                  #Give device time to respond
        rd  = self.ser.readlines()
        
        return b''.join(rd)

    def readHex(self, nBytes = -1):
        '''Read from the device and fromat to hex values with spacing
            Makes output human legible and compatible with Thorlabs APT sheet
        Can specify number of bytes to read, otherwise will read full buffer
        *Useful for debugging.
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Read as specified by nBytes if provided
        if nBytes == -1:
            rd  = self.ser.readlines()
            rd  = b''.join(rd)
        else: 
            rd  = self.ser.read(nBytes)
            rd  = b''.join(rd)
              
        #Convert to readable hex with ' ' spacing
        #IF BYTEARRAY
        tmp = ''
        for i in rd:
            tmp += str(hex(i))
        return tmp.replace('0x', ' ')

    def printDev(self):
        '''Prints the device properties
        Very useful for debugging
        Works even if the port is closed
        *If port is closed, returns instance variables
        *If port is open, returns instance variables AND status bits
        '''
        #Exit if values have not been instantiated
        if self.SN == 'DevNotOpenedYet':
            print('ERROR::: Variables not instantiated yet:\n' +
                  '  open() must be called at least once before printDev()')
            return
        
        #Create string for basic print
        msgOut = ('Device is a    : %s \n' %self.TY   +
                  'Serial Number  : %i \n' %self.SN   +
                  'Frameware vs.  : %s \n' %self.FW   +
                  'Wavelength (nm): %i \n' %self.WL   + 
                  'Max Current(mA): %f \n' %self.MCRT +
                  'Max Power  (mW): %f \n' %self.MPWR +
                  'On port        : %s '   %self.ser.name)
        
        #Check if device is open and print accordingly
        if not self.ser.isOpen():
            print('WARNING::: Device is not open\n' +
                  '  Cannot return status bit values')
        else:
            stats  = self.reqStatBits(1)    #send 1 to request all data
            pwr    = self.reqPowerAct()
            #Key/Intr omitted since device is not setting the bits correctly
            #There is no way to operate in CL so omitted to avoid confusion
            msgOut = (msgOut + '\n' +
                      'Output State   : %s \n' %stats[1] +
                      #'Keyswitch State: %s \n' %stats[2] + 
                      #'Interlock State: %s \n' %stats[3] +
                      #'CL Mode        : %s \n' %stats[4] +
                      'Display Mode   : %s \n' %stats[2] +  #change stats[5]
                      'Power      (mW): %f '   %pwr)
        print(msgOut)
        
    def enableOut(self):
        '''Enables the device output
        *Throws a warning if interlock or keyswitch are disabled
            -NOT WORKING since the device is not setting the bits correctly
        *Does not prevent repeated 'enable' writes since that's not problematic
        Returns boolean TRUE if device reports that it is now ENABLED
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Send enable command
        self.write('11 08 00 00 50 01')
        
        #Check status bits to determine keyswitch and interlock states
        #DEVICE IS NOT SETTING THE KEY/INTR BITS CORRECTLY RIGHT NOW
        stats = self.reqStatBits(1)     #send 1 to request all data
        #isKeys = stats[2]
        #isIntr = stats[3]
        #if not isKeys:
        #    print('Warning::: Keyswitch is not enabled\n' +
        #          '  Laser will not output with key in \'OFF\' position')
        #if not isIntr:
        #    print('Warning::: Interlock is not enabled\n' +
        #          '  Laser will not output without interlock inserted')
            
        return stats[1]

    def disableOut(self):
        '''Disables the device output
        *Doesn't prevent repeated 'disable' writes since that's not problematic
        Returns boolean TRUE if device reports that it is now DISABLED
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Send disable command
        self.write('12 08 00 00 50 01')
        
        #Check status bits to determin if disabled    
        return not self.isEnableOut()     
    
    def isEnableOut(self):
        '''Checks the enable state of the device
        *ONLY CHECKS THE ENABLE BIT.
            If keyswitch or interlock are not enabled, there may be no output
            ^^Statement not true since device bits not setting correctly
        *To check if device should be outputting (enabled AND Keys + Intr)
            Use isOutput()
            ^^Not implemented since the bits aren't setting correctly
                Thus, this function will return the actual output of the device
        Returns: boolean TRUE if device is enabled
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Request and read enable bit (using reqStatBits)
        return self.reqStatBits(1)[1]
    
    #OMITTED: since bits not setting correctly. Thus, isEnableOut() accurately
    #           reflects the device output state
    #def isOutput(self):
    #    '''Checks if all output conditions are met (enabled AND Keys + Intr)
    #    Useful for confirming that device should be outputting
    #    *Use reqAmp() to see if power is below lasing threshold (~.3mA)
    #    '''
    #    pass
    
    def setPower(self, newPWR):
        '''Sets device power to newPWR in mW
        *Power will be set but output must be enabled separately
        Checks that newPWR is within device limits and pushes number if not
        Returns new power reported by device(float with 5 decimal precision)
            (Actual power, Set Power)
            *RETURNED values when device is disabled will not make sense
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Check that newPWR is within limits and push to correct if needed
        tmpfl = False
        if newPWR>self.MPWR:
            newPWR = self.MPWR
            tmpfl = True
        elif newPWR<0:
            newPWR = 0
            tmpfl = True
            
        if tmpfl:
            print('WARNING::: Setpoint beyond limits\n' +
                  '  Pushed setpoint to be within [0 and %f]' %self.MPWR)
        
        #Scale and convert to little-endian hex
        hxPWR = round(self.PWR_BITS*(newPWR/self.MPWR))
        hxPWR = struct.pack('<h', hxPWR)
        
        #Send new power
        self.write('00 08 04 00 D0 01 01 00', hxPWR) #MGMSG_LA_SET_PARAMS cmd
        
        #Check if device is enabled
        if not self.isEnableOut():
            print('WARNING::: Device is not enabled\n' +
                  '  Solution: call enableOut() to enable the device')
        
        #Check actual power
        return (self.reqPowerAct(), self.reqPowerSet())
        
    def setPowerPrct(self, newPWR):
        '''Sets device power to newPWR as a percentage
        *newPWR should be a fractional value (ex. .5 or .923)
        *Power will be set but output must be enabled separately
        Returns new power reported by device(float with 5 decimal precision)
            (Actual power, Set Power)
            *RETURNED values when device is disabled will not make sense
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
            
        #Check that input was a fractional percentage
        if newPWR > 1:
            print('ERROR::: setPowerPrct() takes fractional percentages \n' +
                  '  solution: for example, instead of 50, use .5')
            return
        elif newPWR < 0:
            print('WARNING::: Negative value provided\n' +
                  '  Pushed setpoint to 0')
            newPWR = 0
            
        #Convert to little-endian hex
        hxPWR = round(self.PWR_BITS*newPWR)
        hxPWR = struct.pack('<h', hxPWR)
        
        #Send new power
        self.write('00 08 04 00 D0 01 01 00', hxPWR) #MGMSG_LA_SET_PARAMS cmd
        
        #Check if device is enabled
        if not self.isEnableOut():
            print('WARNING::: Device is not enabled\n' +
                  '  Solution: call enableOut() to enable the device')
        
        #Check actual power
        return (self.reqPowerAct(), self.reqPowerSet())
        
    def reqPowerSet(self):
        '''Requests the power setpoint provided by the device
        *This value may differ from the actual power because the device does 
            not necessarily reach your desired power.
        *Empirical testing shows that the reported output does not match the 
            setpoint exactly even when within range.
        Returns Power Setpoint in mW (float with 5 decimal precision)
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Request and read information
        self.write('01 08 01 00 50 01')     #MGMSG_LA_REQ_PARAMS command
        rd = self.readAll()
        
        #Format output
        hxVal = struct.unpack('<8xh', rd)[0]
        return round(((hxVal)/self.PWR_BITS)*self.MPWR, 5) 
    
    def reqPowerAct(self):
        '''Requests the actual power as provided by the device
        *This value may differ from the power setpoint because the device
            does not necessarily reach your desired power.
        *Empirical testing shows that the reported output does not 
            match the setpoint exactly even when within range.
        Returns Actual Power in mW (float with 5 decimal precision)
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
            
        #Request and read information
        self.write('01 08 03 00 50 01')     #MGMSG_LA_REQ_PARAMS command
        rd  = self.readAll()
        
        #Format output
        hxVal = struct.unpack('<10xh', rd)[0]
        return round(((hxVal)/self.PWR_BITS)*self.MPWR, 5)

    def reqAmpAct(self):
        '''Requests the actual current output
        *There is no way to check the 'current setpoint'
        Returns the current in mA
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Request and read information
        self.write('01 08 03 00 50 01')     #MGMSG_LA_REQ_PARAMS command
        rd  = self.readAll()
        
        #Format output
        hxVal = struct.unpack('<8xh2x', rd)[0]
        return hxVal/100        
                        
    def reqInfo(self):
        '''Reads device information and updates variables
        *These values won't actually change so accessing them from the 
            instance variable is more efficient than repeating a call to 
            reqInfo()
        *To simply display values, use devPrint()
        Returns: Serial number, model number, firmware
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Request and read device information
        self.write('05 00 00 00 50 01')     #MGMSG_HW_REQ_INFO command
        rd  = self.readAll()
        
        #Format and set instance variables
        self.SN, FW, self.TY = struct.unpack('<L10xL19s', rd[6:43]) 
        self.TY    = str(self.TY)[2:-1]        
        FW         = str(FW) 
        self.FW    = FW[0:2] + '.' + FW[2:4] + '.' + FW[4:6]   #Format/Save FW
        
        return (self.SN, self.TY, self.FW)
        
    def reqLim(self):
        '''Reads device operational limits (and wavelength), updates variables
        *These values won't actually change so accessing them from the 
            instance variable is more efficient than repeating a call to 
            reqLim()
        *To simply display the values, use devPrint()
        Returns: Max Current, Max Power, Operating Wavelength
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Request and read device limits (and wavelength)
        self.write('01 08 09 00 50 01')     #MGMSG_LA_REQ_PARAMS command
        rd  = self.readAll()
        
        #Format and set instance variables
        self.MCRT, self.MPWR, self.WL =  struct.unpack('<3h', rd[8:]) 
        self.MPWR  = self.MPWR/10000        #Convert to mW
        self.MCRT  = self.MCRT/100          #Convert to mA

        return (self.MCRT, self.MPWR, self.WL)        
        
    def reqStatBits(self, retFull = None):
        '''Reads device status bits
        RETURNS: just Status bits if no argument provided (ie. reqStatBits())
        RETURNS: if an argument is provided (ie. reqStatBits(_something_)) 
          returns: status bits, output state, display units mode
        *Device does not seem to be setting the keys/intr bits correctly
            Thus, these are currently ommitted in the code to avoid confusion
        *Device does not have a way to go into CL mode...
            Thus, this is omitted from output to avoid confusion
        Status Bits: string of 8 bits where (starting with rightmost bit 1)
             bit 1 = Laser Output    (1 = enabled, 0 = disabled)
             bit 2 = Key Swith       (1 = enabled, 0 = disabled)**
             bit 3 = Control Mode    (1 = power [CL], 0 = current [OL])**
             bit 4 = Interlock       (1 = enabled, 0 = disabled)**
             bit 5 = Disp Units (mA) (1 = mA     , 0 = other)
             bit 6 = Disp Units (mW) (1 = mW     , 0 = other)
             bit 7 = Disp Units (dBm)(1 = dBm    , 0 = other)
             bit 8 = "Future Use"
         ex: 00101011 = display on mW, interlock enabled, keyswitch enabled, 
               output enabled.
        '''
        #Check if port is open
        if not self.ser.isOpen():
            print('ERROR::: Device must be open\n' + 
                  '  solution: call open()')
            return
        
        #Request and read status bits
        self.write('01 08 07 00 50 01')     #MGMSG_LA_REQ_PARAMS command
        rd  = self.readAll()
        
        #Format result
        stat = struct.unpack('<I', rd[-4:])[0]   #Translate Output
        stat = format(stat, '08b')     #Convert to binary string with 0 padding
        
        #Interperet result
        isEnab = bool(int(stat[7]))             #Output state
        isKeys = bool(int(stat[6]))             #Keyswitch stat
        isMdCL = bool(int(stat[5]))             #Control mode (Closed vs. Open)
        isIntr = bool(int(stat[4]))             #Interlock state
        dispMD = ''                         #Display mode (mW, mA, dBm, unkown)
        if bool(int(stat[3])):
            dispMD = 'mA'
        elif bool(int(stat[2])):
            dispMD = 'mW'
        elif bool(int(stat[1])):
            dispMD = 'dBM'
        else:
            dispMD = 'Unkown'
            
        if retFull != None:
            #return (stat, isEnab, isKeys, isIntr, isMdCL, dispMD)
            return (stat, isEnab, dispMD)
        else:
            return (stat)
           
'''
tmp8
Out[541]: bytearray(b'\x06\x00T\x00\x81"\x89S\x9a\x05ION001 \x00,\x00\x02\x019\x00Brushless DC Motor ION Drive\x00\x00\x11\x00\x01\x00\x00\x00\x01\x00')

tmp8 = str(tmp8)

tmp8
Out[543]: '\x06\x00T\x00\x81"\x89S\x9a\x05ION001 \x00,\x00\x02\x019\x00Brushless DC Motor ION Drive\x00\x00\x11\x00\x01\x00\x00\x00\x01\x00'

tmp8.encode('hex')
Out[544]: '06005400812289539a05494f4e30303120002c000201390042727573686c657373204443204d6f746f7220494f4e20447269766500001100010000000100'
'''
