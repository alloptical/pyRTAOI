# communicate with prairie view via tcp socket
# (similar to prairieLink)
# see PraireView user manual for how to setup commands
#  zz 2018
from socket import *

import numpy as np
import time
import array
import os

# test1 = [[1,2],[1,2,3]]
# test2 = [(1,2),(1,2,3)]
# test3 = 'test3';

# print('test1 '+str(len(test1))+str(type(test1[0])))
# print('test2 '+str(len(test2))+str(type(test2[0])))
# print(type(test3))

def make_command(*args):
    # CMD = [int(ord(c)) for c in str_cmd]
    CMD = list()
    for arg in args:
        # print('arg = ' + str(arg)+'; length =' +str(len(arg)))
        if type(arg)==tuple:
            for unit_cmd in arg:
                # print('unit cmd = ' + str(unit_cmd))
                CMD.append(1)
                [CMD.append(int(ord(c))) for c in unit_cmd]
        else:        
            CMD.append(1)
            [CMD.append(int(ord(c))) for c in arg]
    CMD.append(13)
    CMD.append(10)

    return bytes(CMD)




class pvlink(socket):
    def __init__(self,TCP_IP,TCP_PORT):
        super().__init__(AF_INET, SOCK_STREAM)
        self.settimeout(3)

        # make commands
        self.EXIT = make_command('-x')
        self._ts = make_command('-ts')
        self._spp = make_command('-spp')
        self._ss = make_command('-ss')
        self._lv = make_command('-lv')
        self._abort = make_command('-stop')

        # parameters for aquisition
        self.SCAN_SETTINGS = [ ('-sam','Resonant Galvo'),
            ('-fa','1'), # ('-c','1','Off'), ('-c','2','On'), 
            ('-dw'), ('-lbs','true','0'), ('-srd','true','25')]


        try:
            self.connect((TCP_IP, TCP_PORT))
            print('pv socket created')
        except Exception as e:
            print('initiating pl error: '+str(e))

    def __del__(self):
        self.send_done(self.EXIT)
        self.close()
        print('pv socket closed')
            

    def _readline(self):
        data = ""
        raw = bytes(49)
        while raw!= b'\n':
            try:
                raw = self.recv(1)
                data += raw.decode('utf-8')
            except Exception as e:
                print('error reading line:' + str(e))
                break        
        return data        


    def send_done(self,CMD,*args):

        if type(CMD)!= bytes:
            CMD = make_command(CMD,*args)

        while True:    
            try:        
                self.sendall(CMD)
                # print(str(CMD) + 'sent')
            except Exception as e:
                print('error: '+str(e))
                return True
                break

            try:
                # print('checking ACK')
                data = self._readline()
                if 'ACK' in data: break
            except Exception as e:
                print('error: '+str(e))
                return  True
                break
        while True:
            try:
                data = self._readline()
                if 'DONE' in data: 
                    # print('CMD sent')
                    return False
                    break
            except Exception as e:
                print(str(CMD) + 'error: '+str(e))
                return True
                break

    def init_scan(self):
        for CMD in self.SCAN_SETTINGS:
            self.send_done(CMD)
            # print(str(CMD)+'sent')


    def get_frame_size(self):
        samplesPerPixel = np.int32(int(self.send_recv(self._spp)))
        pixelsPerLine = np.int32(self.send_recv('-gts','pixelsPerLine'))
        linesPerFrame = np.int32(self.send_recv('-gts','linesPerFrame'))

        return [samplesPerPixel, pixelsPerLine, linesPerFrame]

    def get_movie_name(self):
        savePath = str(self.send_recv(make_command('-gts', 'directory','1')))
        saveName = str(self.send_recv(make_command('-gts', 'directory','4')))
        saveIteration = str(self.send_recv(make_command('-gts', 'fileIteration','4')))[0:-2].zfill(4)
        full_name =  savePath[0:-2] + os.sep + saveName[0:-2]+'_'+saveIteration 

        return full_name

    def send_recv(self, CMD, *args):
        # NEED TO SYNCHRONISE HERE
        # print('This is send-recv function')
        if type(CMD)!= bytes:
            CMD = make_command(CMD,*args)

        while True:    
            try:
                # print('sending cmd: ' + str(CMD))
                self.sendall(CMD)
                
            except Exception as e:
                print('error: '+str(e))
                return False
                break
            try:
                # print('checking ACK')
                data = self._readline()
                # print('recv data:' + str(data))
                if 'ACK' in data: break
            except Exception as e:
                print('CMD' + str(CMD) + 'error: '+str(e))
                return False
                break
        try:
            rt_data = self._readline()
            # print('received rt_data = '+ rt_data)
        except Exception as e:
            print('error: '+str(e))
            return False
        while True:
            try:
                # print('checking Done')
                data = self._readline()
                # print('received data = '+ data)
                if 'DONE' in data: 
                    return rt_data
                    break
            except Exception as e:
                print('Recv CMD' + str(CMD) + 'error: '+str(e))
                return False
                break
'''
# example

# import prairie link
import pvlink
from pvlink import *
import ctypes
from ctypes import *

# initialise pvlink parameters in main
self.PV_IP = '128.40.156.161'  # bruker 2: TCP_IP = '128.40.202.220'
self.PV_PORT = 1236
self.flipEvenRows = np.int32(1)
self.pvBufferFrames = 1 
self.pl = []


# function to connect to prairie view
def connectPV(self):
    self.pl = pvlink(self.PV_IP,self.PV_PORT)
    print('pl created')
    ERROR = self.pl.send_done('-sam','Resonant Galvo')
    if(not ERROR):
        self.pl.init_scan()
        self.updateStatusBar('PV connected')
        self.FLAG_PV_CONNECTED = True

# start t-series
if not self.pl.send_done(self.pl._ts): # t-series
    self.FLAG_SCANNING = True
    print('started t series')


'''