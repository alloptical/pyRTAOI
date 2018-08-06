def getStaFrameLut(self):
   start_frames = self.stim_frames -p['staPreFrame']
   
   if np.all(i>=0 for i in start_frames)==False: 
       start_frames = start_frames[start_frames>0]
       self.updateStatusBar('First stims not recorded, check stim settings')
       
   staLUT = np.zeros(p['stimStopFrame'])
   num_sta_frames = p['staPreFrame']+p['staPostFrame']

   for i in range(0, p['numberStims']):
       print(self.stim_frames[i])
       np.put(staLUT,np.arange(self.stim_frames[i]-1,self.stim_frames[i]+num_sta_frames),np.arange(1,num_sta_frames+1))
  
   self.staLUT = staLUT
   print(staLUT)