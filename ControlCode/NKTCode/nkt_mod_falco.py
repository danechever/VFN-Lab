####################################################################
#  NKT Laser Python Driver
#
#  Author:  S. Felipe Fregoso / C. Masato Nakano
#  Date:    6/26/2018
#  Rev:     Initial
#           6/24/2018 Initial Working Python-only implm.
####################################################################

# 07/30/18 -- Barebones documentation complete, but more specific details are lacking
#             all get functions are working as expected, but some more testing on set
#             functions with different parameters should be done before initial release

####################################################################
#
#  Usage guide:
#
#    Initially, import this driver into python with:
#
#      from nkt_mod import *
#    
#    Connection to a device is opened by initializing an instance of
#    the Nkt() class to a label, with two parameters 1) IP address,
#    2) port number:
#      
#      ex. a = Nkt('192.168.100.41',2105)
#    
#    Methods are called upon an instance:
#      
#      ex. a.get_emission()
#    
#    Be sure to close the connection upon finishing a session:
#      
#      ex. a.close()
#
####################################################################

import serial
import socket
import time
import ctypes
import struct

# Receive-state macros
Hunting_SOT = 0
Hunting_EOT = 1
MessageReady = 2
Timeout_Error = 3
CRC_Error = 4
Garbage_Error = 5
Overrun_Error = 6
Content_Error = 7
Port_Lost = 8

# Verbosity levels: Generally, verbose is for reply analysis, and very verbose is for
# telegram construction analysis
verbose = False
vverbose = False

###################################################
#
#   module:  SerialDevice
#
#   description:  Handles read and write to any serial port or
#                 serial device using sockets.
#
###################################################

def comm_open (port, arg, parity, stopbits, bytesize, rtscts):
    try:
        fd = serial.Serial (port=port,
                            baudrate=arg,
                            parity=parity,
                            stopbits=stopbits,
                            bytesize=bytesize,
                            rtscts=rtscts,
                            timeout=0.1)
        fd.flushInput()
        fd.flushOutput()
        if fd.isOpen() == False:
            fd = None
    except (serial.SerialException, ValueError) as e:
        fd = None
        raise e

    return fd

def comm_close(fd):
    try:
        fd.close()
    except (serial.SerialException) as e:
        print 'Unable to close file'

def validate_ip(s):
    """ Validates whether input parameter s is a valid ip address
    """
    a = s.split('.')
    if len(a) != 4:
        return False
    for x in a:
        if not x.isdigit():
            return False
        i = int(x)
        if i < 0 or i > 255:
            return False
    return True

def get_bit(val, idx):
    return ((val&(1<<idx))!=0)

def print_statusbits(status_val):
    """ This function takes a 16 bit input and prints an interpretation of each status bit.
    """
    print "Status Bits:"
    print "    Emission LED on:            " + ("TRUE" if get_bit(status_val,0) else "FALSE")
    print "    Interlock off:              " + ("TRUE" if get_bit(status_val,1) else "FALSE")
    print "    Interlock power failure:    " + ("TRUE" if get_bit(status_val,2) else "FALSE")
    print "    Interlock loop off:         " + ("TRUE" if get_bit(status_val,3) else "FALSE")
    print "    External disable:           " + ("TRUE" if get_bit(status_val,4) else "FALSE")
    print "    Supply voltage low:         " + ("TRUE" if get_bit(status_val,5) else "FALSE")
    print "    Module temp range:          " + ("TRUE" if get_bit(status_val,6) else "FALSE")
    print "    USB log error code present: " + ("TRUE" if get_bit(status_val,14) else "FALSE")
    print "    Error code present:         " + ("TRUE" if get_bit(status_val,15) else "FALSE")
    return

def print_varia_statusbits(status_val):
    print "VARIA Status Bits:"
    print "    -                 :         " + ("TRUE" if get_bit(status_val,0) else "FALSE")
    print "    Interlock off     :         " + ("TRUE" if get_bit(status_val,1) else "FALSE")
    print "    Interlock loop in :         " + ("TRUE" if get_bit(status_val,2) else "FALSE")
    print "    Interlock loop out:         " + ("TRUE" if get_bit(status_val,3) else "FALSE")
    print "    -                 :         " + ("TRUE" if get_bit(status_val,4) else "FALSE")
    print "    Supply voltage low:         " + ("TRUE" if get_bit(status_val,5) else "FALSE")
    print "    Module temp range :         " + ("TRUE" if get_bit(status_val,6) else "FALSE")
    print "    -                 :         " + ("TRUE" if get_bit(status_val,7) else "FALSE")
    print "    Shutter sensor 1  :         " + ("TRUE" if get_bit(status_val,8) else "FALSE")
    print "    Shutter sensor 2  :         " + ("TRUE" if get_bit(status_val,9) else "FALSE")
    print "    -                 :         " + ("TRUE" if get_bit(status_val,10) else "FALSE")
    print "    -                 :         " + ("TRUE" if get_bit(status_val,11) else "FALSE")
    print "    Filter 1 moving   :         " + ("TRUE" if get_bit(status_val,12) else "FALSE")
    print "    Filter 2 moving   :         " + ("TRUE" if get_bit(status_val,13) else "FALSE")
    print "    Filter 3 moving   :         " + ("TRUE" if get_bit(status_val,14) else "FALSE")
    print "    Error code present:         " + ("TRUE" if get_bit(status_val,15) else "FALSE")
    return

def statusbitstolist(status_val,bitlen):
    """ This function returns a list of booleans representing the value of each bit in an input
        integer up until the bitlen-th bit.
    """
    ret = []
    for i in xrange(bitlen-1,-1,-1):
        ret.append(get_bit(status_val,i))
    return ret

class CommError(Exception):
    """Base class for communication exception in this module"""
    pass

class CommTimeoutException(CommError):
    """Exception raised for timeout errors"""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

class SerialDevice(object):

    def close (self):
        self.fd.close()
        self.fd = None

    def fd_read(self, numbytes):
        try:
            if type(self.fd) is socket._socketobject:
                return self.fd.recv(numbytes)
            else:
                return self.fd.read(numbytes)
        except Exception as e:
            print('caught this error: ' + repr(e))
            raise CommTimeoutException('Timeout on read', repr(e))

    def fd_write(self, buf):
        try:
            if type(self.fd) is socket._socketobject:
                return self.fd.send(buf)
            else:
                return self.fd.write(buf)
        except Exception as e:
            raise CommTimeoutException('Write timeout', str(e))

    def fd_bytestoread(self):
        try:
            if type(self.fd) is socket._socketobject:
                raise
            else:
                return self.fd.available()
        except Exception as e:
            raise

    ########################################################

    def write(self, cmd):
        """ Write command
        Write the command to the communication port but don't read the port"
        """
        self.send_cmd(self.build_cmd(cmd), response_expected=False)

    def query(self, cmd):
        """ Write command and return response.
        Write to the communication port and read the response
        """
        cmdstr = self.build_cmd(cmd)
        resp = self.send_cmd(cmdstr, response_expected=True)
        return resp

    def build_cmd (self, mnemonic, data=''):
        """Build a command.
        Append the data string and the terminating character to a mnemonic
        string.
        """
        cmdstr = mnemonic + ' ' + data # + self.terminator
        return cmdstr

    def send_cmd(self, cmd, response_expected=True):
        """ Send command to device
        Sends the command to the device and if 'response_expected' is True
        it will wait for a reply.  The reading ends if terminator characters
        are reached or if the read timesout.
        """
        response = ''

        if self.verbose:
            print ("Command Sent: %s"%(cmd))
        try:
            cmd += self.write_terminator
            self.fd_write(cmd)
        except (serial.SerialTimeoutException, socket.timeout) as e:
            print ("Error. Serial write timed out: %s"%(str(e)))
            raise
            return None
        if response_expected:
            while True:
                try:
                    response += self.fd_read(1)
                    if len(response) <= 0:
                        print "Error.  No data received"
                        return response
                    #print response

                    if response[-1*len(self.read_terminator):] == self.read_terminator:
                        response = response[:-1*len(self.read_terminator)]
                        if self.verbose:
                            print "Response: %s" % response
                        break
                except serial.SerialTimeoutException as e:
                    print "Error. Serial read timedout"
                    return None

        # sleep after last transmit 50ms per User's Manual
        # also sleep after completed read of response
        time.sleep(0.05)
        return response

    def close(self):
        comm_close(self.fd)
        self.fd = None

    def __init__(self, port, arg, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, bytesize=serial.EIGHTBITS, rtscts=True, timeout=0.1, verbose=False, write_terminator='', read_terminator=''):

        self.fd = None
        self.write_terminator = write_terminator
        self.read_terminator  = read_terminator

        ### Determine if port is an IP address
        if (validate_ip(port)):
            if verbose: print("IP address supplied: IP=%s, port=%d\n"%(port,arg));
            try:
                self.fd = socket.socket()
                self.fd.settimeout(timeout)
                self.fd.connect((port, arg))
            except socket.error as e:
                raise CommError('connection', 'Unable to connect to socket-server: %s'%(str(e)))
                self.fd = None
        else:
            ### Connect to serial port
            try:
                self.fd = comm_open (port, arg, parity, stopbits, bytesize, rtscts)
		print('Connected to serial port')
            except Exception as e:
                raise CommError('connection', 'Unable to connect to serial port: %s'%(str(e)))
                self.fd = None

class Nkt(SerialDevice):
    def __init__(self, port, arg, query_id=False, verbose=True):
        #Initialize variables
        self.verbose = verbose
        self.terminator = ''
        self.txBuffer = '' #char buffer
        self.txCRC = 0 #short int
        self.masterId = 66
        self.xtrmModuleAddr = 15
        super(Nkt, self).__init__(port, arg, rtscts=False, write_terminator='', read_terminator='')

    def addToTxMsgData2(self, data, escParse=False, updCRC=False):
        """ Adds a byte of data (parameter 'data') to the currently building telegram
        """
        if updCRC: 
            self.txCRC = self.calcCRC16(data, self.txCRC)
        if escParse:
            if data == 0x0D or data == 0x0A or data == 0x5E:
                if vverbose: print('escParse: data=%d'%(data))
                self.txBuffer = self.txBuffer + struct.pack('B', 0x5E)
                data += 0x40
        self.txBuffer = self.txBuffer + struct.pack('B', data)
        if vverbose:
            print("%s"%(':'.join("{:02x}".format(ord(c)) for c in self.txBuffer)))

    def addToTxMsgData(self, idx, cnt, data=None, escParse=False, updCRC=False):
        """ Adds a multi-byte list of data (parameter 'data') to the currently building
            telegram
        """
        for n in xrange(cnt):
            self.addToTxMsgData2(ord(data[n+idx]),escParse,updCRC) #likely not pythonic

    def calcCRC16(self, data, oldCRC):
        """ Adds CRC to the data encoding
        """
        data = ctypes.c_ubyte(data).value

        oldCRC = ctypes.c_ushort(oldCRC >> 8 | oldCRC << 8).value
        oldCRC = oldCRC ^ data
        oldCRC = oldCRC ^ ctypes.c_ubyte((oldCRC & 0xFF) >> 4).value
        oldCRC = oldCRC ^ ctypes.c_ushort((oldCRC << 8) << 4).value
        oldCRC = oldCRC ^ ctypes.c_ushort(((oldCRC & 0xFF) << 4) << 1).value
        return oldCRC

    def BytesToRead(self):
        """ Returns the number of bytes to be read in the current working port.
        """
        buf = array.array('h',[0])
        ret = fcntl.ioctl(self.fd,FIONREAD,buf)
        if(ret<0):
            return -1
        else:
            return buf[0]

    def flush(self): #broken
        if type(self.fd) is socket._socketobject:
            while True:
                numBytes = self.BytesToRead()
                if(numBytes>0):
                    self.fd_read(1)
                else:
                    break
        else:
            tcflush(self.fd,TCIFLUSH)
        
    def send_message(self, deviceId, registerId, msgType, data = None):
        """ Sends a telegram to the device and register specified in parameters deviceId
            and registerId respectively, consisting of the data specified in parameter data
        """
        self.txBuffer = '' #initialize empty buffer
        self.txCRC = 0 #initialize crc to 0

        if vverbose:
            print("+--------------------------------------------+\n"
            +"|             Building Telegram              |\n"
            +"+--------------------------------------------+")

        self.addToTxMsgData2(0x0D, False, False) #sot
        self.addToTxMsgData2(deviceId, True, True)
        self.addToTxMsgData2(self.masterId, True, True)
        self.addToTxMsgData2(msgType, True, True)
        self.addToTxMsgData2(registerId, True, True)
        if msgType == 0x05 and data is not None: # 0x05 = Write message - 0x04 = Read message
            self.addToTxMsgData(0,len(data),data,True,True)
        self.addToTxMsgData2(ctypes.c_ubyte(self.txCRC >> 8).value, True, False)
        self.addToTxMsgData2(ctypes.c_ubyte(self.txCRC).value, True, False)
        self.addToTxMsgData2(0x0A, False, False) #eot

        #if self.verbose: print(self.txBuffer)
        if vverbose: print("String Telegram: "+repr(self.txBuffer)+'\n')
        self.fd_write(self.txBuffer)

    def receive_message(self, deviceId, registerId, msgType, payload=None):
        """ Receives a telegram from the device and register specified in parameters deviceId
            and registerId respectively, and stores it in string parameter payload
        """
        RxBuffer = '' #data without beginning/ending markers
        #read 1024 bytes
        crc = 0
        reply = self.fd_read(1024)
        state = Hunting_SOT
        if vverbose:
            print("Received: %s\n"%(':'.join("{:02x}".format(ord(c)) for c in reply)))
            print("+--------------------------------------------+\n"
            +"|              Parsing Telegram              |\n"
            +"+--------------------------------------------+")
        for c in range(len(reply)):
            if vverbose: print("Byte %d = %d"%(c,ord(reply[c])))
            if (state==Hunting_SOT):
                if (ord(reply[c])==0x0d):
                    if vverbose: print("Found Start of Telegram")
                    escape = False
                    crc = 0
                    state = Hunting_EOT
            elif (state==Hunting_EOT):
                if (ord(reply[c])==0x0a):
                    if vverbose: print("Found End of Telegram\n")
                    if (len(RxBuffer)>=5): #Dest+Src+Type+msbCRC+lsbCRC==Minimum telegram length
                        if (crc==0): #We have collected a message with valid CRC - Check contents
                            if verbose:
                                print("+--------------------------------------------+\n"
                                +"|             Checking Contents              |\n"
                                +"+--------------------------------------------+")
                                print("RxBuffer[0] = 0x%X; masterId = 0x%X ; RxBuffer[1] = 0x%X; deviceId = 0x%X;"
                                " RxBuffer[2] = 0x%X ; msgType = 0x%X; RxBuffer[3] = 0x%X ; registerId = 0x%X"
                                %(ord(RxBuffer[0]),self.masterId,ord(RxBuffer[1]),deviceId,ord(RxBuffer[2]),msgType,
                                ord(RxBuffer[3]),registerId))
                            if (ord(RxBuffer[0])==self.masterId and ord(RxBuffer[1])==deviceId and ord(RxBuffer[2])==
                            msgType and ord(RxBuffer[3])==registerId):
                                RxBuffer = RxBuffer[:-2] + chr(0) + chr(0) #Remove CRC
                                payload.append(RxBuffer[4:-2]) #Modify payload
                                state = MessageReady
                            else:
                                state = Content_Error
                        else:
                            state = CRC_Error
                    else:
                        Garbage_Error
                else:
                    if (ord(reply[c])==0x5e):
                        escape = True
                    else:
                        if escape:
                            RxBuffer = RxBuffer + chr(ord(reply[c]) - 0x40)
                            escape = False
                        else:
                            RxBuffer = RxBuffer + reply[c]
                        crc = self.calcCRC16(ord(RxBuffer[len(RxBuffer)-1]),crc)
        if verbose: print("RxState = %d\n"%state)
        return state

    def read_UINTbyte(self, deviceId, registerId):
        """ Reads a telegram containing a one-byte message from the device specified by
            parameters deviceId and registerId
        """
        tempArray = []
        retVal = 0
        ### Flush port
        #self.flush()
        ### Send Message
        self.send_message(deviceId, registerId, 0x04)
        if (self.receive_message(deviceId,registerId,0x08,tempArray)==MessageReady):
            if verbose:
                print("+--------------------------------------------+\n"
                +"|                Parsing Reply               |\n"
                +"+--------------------------------------------+")
                print("tempArray size = %d"%len(tempArray[0]))
                print("tempArray[0] = %d\n"%ord(tempArray[0][0]))
            ### Interpret as Byte
            if (len(tempArray[0])==1):
                retVal = ord(tempArray[0][0])
        return retVal

    def write_UINTbyte(self, deviceId, registerId, data = None):
        """ Writes a one-byte message consisting of parameter data to the device specified by
            parameters deviceId and registerId
        """
        tempArray = []
        retVal = False
        outArray = b""
        outArray += chr(data)
        self.send_message(deviceId,registerId,0x05,outArray)
        if(self.receive_message(deviceId,registerId,0x03,tempArray)==MessageReady):
            retVal = True
        return retVal
        
    def read_UINT16(self, deviceId, registerId):
        """ Reads a telegram containing a one-byte message from the device specified by
            parameters deviceId and registerId
        """
        tempArray = []
        retVal = 0
        ### Flush port
        #self.flush()
        ### Send Message
        self.send_message(deviceId, registerId, 0x04)
 
        ### Get reply
	if (self.receive_message(deviceId,registerId,0x08,tempArray)==MessageReady):
            if verbose:
                print("+--------------------------------------------+\n"
                +"|                Parsing Reply               |\n"
                +"+--------------------------------------------+")
                print("tempArray size = %d"%len(tempArray[0]))
                print("tempArray[1] = %d"%ord(tempArray[0][1]))
                print("tempArray[0] = %d\n"%ord(tempArray[0][0]))
        ### Interpret as 16BIT
            if (len(tempArray[0])==2):
                retVal = ord(tempArray[0][1])<<8|ord(tempArray[0][0])
        return retVal

    def write_UINT16(self, deviceId, registerId, data = None):
        """ Writes a two-byte message consisting of parameter data to the device specified by
            parameters deviceId and registerId
        """
        tempArray = []
        retVal = False
        outArray = b""
        outArray += chr(data&0x00ff)
        outArray += chr((data&0xff00)>>8)
        self.send_message(deviceId,registerId,0x05,outArray)
        if(self.receive_message(deviceId,registerId,0x03,tempArray)==MessageReady):
            retVal = True
        return retVal

    def get_emission(self):
        """ This function gets the emission status of the queried device. A return of 1 means
            that the device is on, and a return of 0 means that the device is off.
        """

        #print("-------------------- Getting Emission --------------------")
        ret = self.read_UINTbyte(15, 0x30)
        if ret:
            #print("Emission = ON\n")
            return True
        else:
            #print("Emission = OFF\n")
            return False

    def set_emission(self,off_on):
        """ This function sets the emission status of the directed device. An input of boolean
            True turns on the device, while an input of boolean False turns off the device.
        """
        #print("-------------------- Setting Emission --------------------")
        if off_on:
            data = 3
            status = "ON"
        else:
            data = 0
            status = "OFF"
        ret = self.write_UINTbyte(15,0x30,data)
        #print("Emission set to : "+status+"\n")
        return int(ret)

    def get_interlock(self):
        """ This function gets the interlock status of the queried device. A return of >0
        """
        print("------------------- Getting Interlock --------------------")
        ret = self.read_UINT16(15,0x32)
        if ret:
            print("Interlock : ON\n")
            return 1
        else:
            print("Interlock : OFF\n")
            return 0

    def set_interlock(self,data):
        print("------------------- Setting Interlock --------------------")
        ret = self.write_UINT16(15,0x32,ctypes.c_ushort(data))
        if ret:
            if data:
                print("Interlock set to ON\n")
            else:
                print("Interlock set to OFF\n")
            return 0
        else:
            print("Interlock set error\n")
            return -1

    def get_powerlevel(self):
        """ This function gets the power level of the queried device and returns its value.
        """
        #print("------------------- Getting Power Level ------------------\n")
        ret = self.read_UINT16(15, 0x37)
        #print("Power Level = %0.1f%%\n"%(float(ret)*0.1))
        return float(ret) * 0.1

    def set_powerlevel(self,power_level):
        """ This function sets the power level of the directed device to its input value.
        """
        #print("------------------- Setting Power Level ------------------\n")
        data = int(power_level/0.1)
        ret = self.write_UINT16(15,0x37,data)
        #print("Power Level set to : %0.1f%%\n"%(power_level))
        return 0 if ret else -1

    def get_inlettemp(self):
        """ This function gets the inlet temperature of the queried device and returns its value.
        """
        print("------------------- Getting Inlet Temp -------------------\n")
        ret = self.read_UINT16(15,0x11)
        print("Inlet Temp = %0.1f\n"%(float(ret)*0.1))
        return float(ret) * 0.1

    def get_currentlevel(self):
        """ This function gets the current level of the queried device and returns its value.
        """
        print("----------------- Getting Current Level ------------------\n")
        ret = self.read_UINT16(15,0x38)
        print("Current Level = %0.1f\n"%(float(ret)*0.1))
        return float(ret) * 0.1

    def set_currentlevel(self,current_level):
        """ This function sets the current level of the directed device to its input value.
        """
        print("----------------- Setting Current Level ------------------\n")
        data = int(current_level/0.1)
        ret = self.write_UINT16(15,0x38,data)
        print("Current Level set to : %0.1f%%\n"%(power_level))
        return 0 if ret else -1

    def get_setupbits(self):
        """ This function gets the setup bits from the queried device. For the purpose of
            these devices, a return value of 1 means that the laser in question is in power
            mode, and a return value of 0 refers to the laser being in current mode.
        """
        print("------------------- Getting Setup Bits -------------------\n")
        ret = self.read_UINT16(15,0x31)
        print("Setup Bits : %d\n"%(ret))
        return ret

    def set_setupbits(self,mode):
        """ This function sets the setup bits of the directed device to the input. The meaning
            of the bits is described above.
        """
        print("------------------- Setting Setup Bits -------------------\n")
        if mode < 0 or mode > 1:
            print("Error : set_setupbits - mode must be 0 or 1\n")
            return -1
        data = int(mode)
        ret = self.write_UINT16(15,0x31,data)
        if ret:
            if data is 0:
                print("Setup Bits set to Current Mode\n")
            if data is 1:
                print("Setup Bits set to Power Mode\n")
            return 0
        else: return -1

    def get_statusbits(self):
        """ This function gets the status bits from the queried device.
            Status bits are defined on line 80.
        """
        print("------------------ Getting Status Bits -------------------\n")
        bitwise = self.read_UINT16(15,0x66)
        print("Status Bits : %d\n"%(bitwise))
        print_statusbits(bitwise)
        ret = (bitwise, statusbitstolist(bitwise,9))

    def set_statusbits(self,status_val):
        """ This function sets the status bits of the directed device to the input.
            Status bits are defined on line 87.
        """
        print("------------------ Setting Status Bits -------------------\n")
        data = int(status_val)
        ret = self.write_UINT16(15,0x66,data)
        print("Status Bits set\n")
        print_statusbits(data)
        return 0 if ret else -1

    def get_pulsepickerratio(self):
        """ This function gets the pulse picker ratio from the queried device.
        """
        print("-------------- Getting Pulse Picker Ratio ----------------\n")
        ret = self.read_UINT16(15,0x34)
        print("Pulse Picker Ratio : %d\n"%(ret))
        return ret

    def set_pulsepickerratio(self,val):
        """ This function sets the pulse picker ratio of the directed device.
        """
        print("-------------- Setting Pulse Picker Ratio ----------------\n")
        data = val
        ret = self.write_UINT16(15,0x34,data)
        print("Pulse Picker Ratio set to : %d\n"%(data))
        return 0 if ret else -1

    def get_pulsepickerdelay(self):
        """ This function gets the pulse picker delay from the queried device.
        """
        print("-------------- Getting Pulse Picker Delay ----------------\n")
        ret = self.read_UINTbyte(15,0x35)
        print("Pulse Picker Delay : %0.1f\n"%(float(ret)*0.25))
        return float(ret) * 0.25

    def set_pulsepickerdelay(self,ns):
        """ This function sets the pulse picker delay of the directed device.
        """
        print("-------------- Setting Pulse Picker Delay ----------------\n")
        data = int(ns/0.25)
        ret = self.write_UINTbyte(15,0x35,data)
        print("Pulse Picker Delay set to : %0.1f\n"%(ns))
        return 0 if ret else -1

    def get_watchdoginterval(self):
        """ This function gets the watchdog interval from the queried device.
        """
        print("-------------- Getting Watchdog Interval -----------------\n")
        ret = self.read_UINTbyte(15,0x36)
        print("Watchdog Interval : %d s"%(ret))
        return ret

    def set_watchdoginterval(self,seconds):
        """ This function sets the watchdog interval of the directed device to the
            input value seconds.
        """
        print("-------------- Setting Watchdog Interval -----------------\n")
        data = seconds
        ret = self.write_UINTbyte(15,0x36,data)
        print("Watchdog Interval set to : %d s"%(data))
        return 0 if ret else -1

    def get_varia_monitorinput(self):
        """ This function gets the monitor input from the queried Varia device.
        """
        print("----------------- Getting Monitor Input ------------------\n")
        ret = self.read_UINT16(16,0x13)
        print("Varia Monitor Input : %0.1f\n"%(float(ret)*0.1))
        return float(ret) * 0.1

    def get_varia_ndsetpoint(self):
        """ This function gets the ND setpoint from the queried Varia device and
            returns its value in nm. (?)
        """
        print("------------------ Getting ND Setpoint -------------------\n")
        ret = self.read_UINT16(16,0x32)
        print("ND Setpoint : %0.1f\n"%(float(ret)*0.1))
        return float(ret) * 0.1

    def set_varia_ndsetpoint(self,nd):
        """ This function takes an input parameter nd and sets the ND setpoint of
            the directed Varia device to that value.
        """
        print("------------------ Setting ND Setpoint -------------------\n")
        data = int(nd/0.1)
        ret = self.write_UINT16(16,0x32,data)
        print("ND Setpoint set to : %0.1f\n"%(nd))
        return 0 if ret else -1
        
    def get_varia_swpsetpoint(self):
        """ This function gets the SWP setpoint from the queried Varia device and
            returns its value in nm.
        """
        #print("------------------ Getting SWP Setpoint ------------------\n")
        ret = self.read_UINT16(16,0x33)
        #print("SWP Setpoint : %0.1f nm\n"%(ret*0.1))
        return ret*0.1

    def set_varia_swpsetpoint(self,nm):
        """ This function takes an input parameter nm and sets the SWP setpoint of
            the directed Varia device to that value.
        """
        #print("------------------ Setting SWP Setpoint ------------------\n")
        data = int(nm/0.1)
        ret = self.write_UINT16(16,0x33,data)
        #print("SWP Setpoint set to : %0.1f nm\n"%(nm))
        return 0 if ret else -1

    def get_varia_lwpsetpoint(self):
        """ This function gets the LWP setpoint from the queried Varia device and
            returns its value in nm.
        """
        #print("------------------ Getting LWP Setpoint ------------------\n")
        ret = self.read_UINT16(16,0x34)
        #print("LWP Setpoint : %0.1f nm\n"%(ret*0.1))
        return ret*0.1

    def set_varia_lwpsetpoint(self,nm):
        """ This function takes an input parameter nm and sets the LWP setpoint of
            the directed Varia device to that value.
        """
        #print("------------------ Setting LWP Setpoint ------------------\n")
        data = int(nm/0.1)
        ret = self.write_UINT16(16,0x34,data)
        #print("LWP Setpoint set to : %0.1f nm\n"%(nm))
        return 0 if ret else -1

    def get_varia_statusbits(self):
        """ This function gets the status bits from the queried Varia device.
            Status bits are defined on line 102.
        """
        print("--------------- Getting Varia Status Bits ----------------\n")
        bitwise = self.read_UINT16(15,0x66)
        print_varia_statusbits(bitwise)
        ret = (bitwise, statusbitstolist(bitwise,16))
        return ret

if __name__ == "__main__":
    """
    a = Nkt()
    a.get_powerlevel()
    a.get_emission()
    a.get_inlettemp()
    a.get_currentlevel()
    a.get_interlock()
    a.get_setupbits()
    a.get_statusbits()
    a.set_powerlevel(6)
    a.get_powerlevel()
    a.set_emission(False)
    a.get_emission()
    a.get_varia_swpsetpoint()
    a.get_varia_lwpsetpoint()
    """
