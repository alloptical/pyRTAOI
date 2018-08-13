# work with HoloBlink
# send coordinates of holographic spots via TCP
# ZZ 2018

from socket import *
import numpy as np
import time



class bLink(socket):

	def __init__(self,TCP_IP,TCP_PORT):
		super().__init__(AF_INET, SOCK_STREAM)
#		self.setsockopt(socket.TCP_NODELAY,1) # ADDED 20180421
		self.settimeout(60)
		self.CONNECTED = False
		try:
			self.connect((TCP_IP, TCP_PORT))

			while True:
				try:
					data = self._readline()
					if 'Hello' in data: 
						print(data[0:-2])
						self.CONNECTED = True
						break
				except Exception as e:
					print('connection error: '+str(e))
					break

			print('bLink socket created')
		except Exception as e:
			print('initiating bl error: '+str(e))

	def add_prefix(self,prefix,xx):
		Xmsg = bytes(";".join([str(x) for x in xx]),'utf-8')
		Xmsg = bytes(prefix+str(len(Xmsg)).zfill(4),'utf-8')+ Xmsg
		return Xmsg

# x y should change together - TO DO

	def send_coords(self,xx,yy):
		# send as two commands - this is stupid
		# self.sendall(self.add_prefix("X",xx)+self.add_prefix("Y",yy))

		# send as one, added 20180423
		self.sendall(self.add_prefix("C", xx+yy))


		print("msg sent")
		while self.CONNECTED:
			try:
				data = self._readline()
				if 'Done' in data: 
					print("reply recvd")
					return False
					break
			except Exception as e:
				print('receiving error: '+str(e))
				self.CONNECTED = False
				return True
				break

		return False

	def _readline(self):
		data = ""
		raw = bytes(49)
		while raw!= b'\n':
			try:
				raw = self.recv(1)
				data += raw.decode('utf-8')
			except Exception as e:
				print('error reading line:' + str(e))
				self.CONNECTED = False
				break		
		return data		

	def abort(self):
		self.close()


