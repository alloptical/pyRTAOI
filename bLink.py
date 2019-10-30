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
		self.settimeout(10)
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


	def send_coords(self,xx,yy):
		# send coordinates, and wait until blink updated
		self.sendall(self.add_prefix("C", xx+yy))
		while self.CONNECTED:
			try:
				data = self._readline()
				if 'Done' in data:
					return False
					break
			except Exception as e:
				print('receiving error: '+str(e))
				self.CONNECTED = False
				return True
				break

		return False

	def send_duration(self,duration):
		print('sending duration')
		self.sendall(self.add_prefix("D", [duration]))
		while self.CONNECTED:
			try:
				data = self._readline()
				if 'Done' in data:
					return False
					break
			except Exception as e:
				print('receiving error: '+str(e))
				self.CONNECTED = False
				return True
				break

		return False


	def send_coords_power(self,xx,yy,volt):
		# send coords and control voltage for aom; use this when holoblink is set to trigger photostim
		self.sendall(self.add_prefix("P", xx+yy+[volt]))
		while self.CONNECTED:
			try:
				data = self._readline()
				if 'Done' in data:
					return False
					break
			except Exception as e:
				print('receiving error: '+str(e))
				self.CONNECTED = False
				return True
				break

		return False

	def send_trigger_power(self,volt):
		print('sending trigger')
		try:
			self.sendall(self.add_prefix("T", [volt]))
		except Exception as e:
				print('sending trigger error: '+str(e))
				self.CONNECTED = False
				return True

		while self.CONNECTED:
			try:
				data = self._readline()
				if 'Done' in data:
					return False
					break
			except Exception as e:
				print('receiving error: '+str(e))
				self.CONNECTED = False
				return True
				break

		return False

	def send_only(self,xx,yy):
		# send coordinates, do not wait
		self.sendall(self.add_prefix("C", xx+yy))
		print("msg sent")


	def rev_only(self):
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


