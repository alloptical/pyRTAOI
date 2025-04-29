# pyRTAOI
# Introduction
A python-based real-time all-optical interface (**pyRTAOI**) for closed-loop control of neural activity.<br/> 

The software integrates a calcium imaging analysis toolbox (CaImAn, Giovannucci et al., 2019, https://github.com/flatironinstitute/CaImAn), a custom hologram control software (HoloBlink in this repository) with a two-photon microscope control system (Prairie View, Bruker Corporation). Holographic photostimulation can be delivered to pre-defined or online-detected ensembles of neurons based on the calcium activity recorded from a single neuron or population of neurons. It also provides users with access to real-time motion-corrected and denoised imaging frames, manual curation and selection of ROIs for readout and stimulation, configuration of photostimulation and sensory stimulation protocols (e.g. laser power, frequency, amplitude, duration, etc.), automatic photostimulation and imaging laser alignment check, direct control of relevant parameters settings in CaImAn and the Prairie View microscope system.<br/>
 <img width="1287" alt="pyRTAOIscreenshot" src="https://github.com/user-attachments/assets/1342b6df-3d86-43d7-ac66-5bc129f87cf4">

 
<img width="745" alt="HoloblinkScreenshot" src="https://github.com/user-attachments/assets/449b3fd5-1d19-43b9-902f-446b5d3870f5">

# System requirements
* Software has been tested on a Windows 10 desktop with microscope control system (PrairieView, rev40, Bruker Corporation).<br/>
* Software platforms used in the package include: Python 3.6 (with PyQt5), Visual Studio 2013 (64 bit, with Qt 5.9) and MATLAB (2017b).<br/>
* Access to the raw image data stream depends on PrairieLink (Bruker Corporation).<br/>
* Phase masks are uploaded to the spatial light modulator (OverDrive Plus SLM, Meadowlark Optics) using the Blink_SDK dll (Meadowlark Optics).<br/>
* NVIDIA GPUs have been tested for raw data preprocessing and hologram calculation (CUDA8.0, NVDIA GeForce GTX 1080Ti).<br/>
* Analog voltage outputs are generated using NI-DAQmx (device used: PCI-6713, National Instruments).<br/>

# Installation notes
On microscope control PC:<br/>
*	Install Prairie View 5.4 (tested with rev40).<br/>
*	Install Python (Conda package) and CaImAn following instructions in https://github.com/flatironinstitute/CaImAn.<br/>
* Install Python dependencies: nidaqwx, pyqtgraph, scikit-image (skimage; tested with version 0.14.0), numpy (tested with version 1.16.0).<br/>
*	Install NI-DAQmx (tested with version 15.5)<br/>
*	Install CUDA toolkit 8.0.<br/>
*	Download pyRTAOI. Change the NI device ID and TCP/IP address (PV_IP, PV_PORT, BLINK_IP and BLINK_PORT).<br/>
*	(Optional, if using pyRTAOI to trigger photostimuli) connect an analog voltage output (tested with PCI-6713, National Instruments) to the photostimulation trigger input.<br/>

On SLM control PC:
*	Install Visual Studio 2013 (with Service Pack 1).<br/>
*	Install CUDA toolkit 8.0.<br/>
*	Install Measurement Studio (version 15) for Visual Studio 2013.<br/>
*	Install NI-DAQmx (tested with version 15.5).<br/>
*	Download folder 'HoloBlink'. Change the NI device ID. Rebuild the solution.<br/>
*	Connect an analog voltage output (tested with PCI-6713, National Instruments) to the photostimulation trigger input and AOM voltage control port.<br/>

# How to use
You can setup experiment parameters in the GUI and save into configuration files. Before running pyRTAOI, start the microscope control (PrairieView) and SLM control software (Holobink). Then start pyRTAOI, connect to PrairieView and Holoblink, load your configuration file, and start all-optical experiments.<br/>

**To start pyRTAOI**:

```
conda activate caiman # activate virtual environment
python pyRTAOI.py # open GUI
```
## Example protocols 

### 1.	Imaging-photostimulation calibration check

 This step aims to burn spirals in an fluorescent slide to check co-alignment of imaging and photostimulation beams before running an experiment.
  * Start PrairieView, start Markpoint
  *	Run Holobink (on the PC connected to SLM)
  *	Run pyRTAOI (in Spyder), 
  *	Click 'Connect Priarie'  and 'Connect Blink'.
  *	Load configuration file from Configs folder: 20190722_calcheck.cfg
  *	Click 'RunCalCheck' to burn spots.
 <img width="587" alt="image" src="https://github.com/user-attachments/assets/755a5151-7f85-46a0-a2b1-fe7889480874" />
 <img width="770" alt="image" src="https://github.com/user-attachments/assets/35b99380-3523-4cb3-9336-b161a4e71e07" />
 
  * Burnt spots should align with the orange targets if the imaging and photostimulation paths are well-aligned:
 <img width="510" alt="image" src="https://github.com/user-attachments/assets/008bfd86-92a2-4276-b0cb-0f51466fe341" />
 

### 2.	Cell detection

First, we need to take a short film to initiase CaImAn. 
 *	Set t-series frames (#Reps) to 500 in PrairieView
 *	Load configuration file (20190804_init_zoom2_16x.cfg).
<img width="458" alt="image" src="https://github.com/user-attachments/assets/0e306fac-8c03-4c31-8e4b-5d4f4e97e4ce" />

 *	Click ‘Take Ref Movie’ to grab a short film of FOV, then click ‘Initialise (auto save)’ when done.
   <img width="439" alt="image" src="https://github.com/user-attachments/assets/db8358c1-7948-4ed9-b9e2-1951355a43c3" />

Next, we take a long film to detect cells in FOV online. 
•	Load configuration file (20190804_tex_celldetection_zoom2_16x.cfg) then click 'update trials'
•	Click Run
You can initialse with the result file and repeat this step until enough cells are detected.
<img width="1054" alt="image" src="https://github.com/user-attachments/assets/437eb606-7988-4c24-b43c-ffafebefff14" />
Go to the 'plot' tab to check the cells:
<img width="916" alt="image" src="https://github.com/user-attachments/assets/8d25fff6-8b2f-404c-aee3-0589171273e7" />


### 3.  Cell identification and activity analysis (during behaviour)

This step will record activity from the detected cells as an animal is involved in a behavior task. Here we use PyBehaviour (by Lloyd Russell) to setup a behavioral task session and record task performance.
 *	Load configuration file (20200926_tex_cellidentification_zoom2_16x.config)
 *	Load the ‘trialOrder’ file (saved out from pyBehaviour) and click ‘Plot trial type’ to check the number of trials.
<img width="604" alt="image" src="https://github.com/user-attachments/assets/8c0375e6-dca8-4ccf-a004-b2fa025a8a19" />

 *	Change t-series frames in PrairieView according to the ‘Frames needed’ suggested by pyRTAOI.
 * Click Run.
Result files (in both .pkl and .mat formats) will be saved once the session ends (in ‘pyrtaoi_results’ folder), which contain the spatial footprint and inferred activity traces as extracted by CaImAn and the parameters used. You can use them to analyse the neuronal functional identities to decide when and which cells to stimulate in the next sessions. 
OnlineProcTexture.m is an example Matlab script to find trial-selective neurons in the barrel cortex of mice performing a texture discrimination task.
<img width="849" alt="image" src="https://github.com/user-attachments/assets/27fbce1a-065c-4430-b237-a82161bd7cfc" />


### 4. Automatic recruit newly-detected cells as photostimulation targets:

This allows cells to be added as photostimulation targets in the next trials as they are detected online.
 * Load opsin image (drag-drop to display window) to only include opsin-positive cells
<img width="467" alt="image" src="https://github.com/user-attachments/assets/6e6b482c-0f67-433f-ace5-5f10e29429ef" />

 * Check 'auto add' to photostimulate newly-detected cells
<img width="841" alt="image" src="https://github.com/user-attachments/assets/f729b790-777a-4f8d-842b-eb510de893de" />

<img width="829" alt="image" src="https://github.com/user-attachments/assets/a230d56a-c515-424b-9130-32955c457529" />

   
### 5.	Photo-excitabiliy test 

This is an optional step to test which cells-of-interest are responsive to photpstimulation by stimulating them one-by-one in a random sequence with several repeats.
 *	Load configuration file: 20190804_sequence_photo_zoom2_16x
 *	Select the cells to test by clicking on FOV, or (recommended) load from a configuration file by clicking 'Load trigger target config'. (e.g. select the ‘OutputParam’ file generated by the Matlab script in the texture example).
 *	Click 'Config Seq Stim'
 *	Check 'Stim from blink'
 *	Click 'Update Params' (make sure HoloBlink updates)
<img width="519" alt="image" src="https://github.com/user-attachments/assets/cf3c13de-9d9d-4cb4-adbe-e1e64404f462" />

 *	Change t-series frames accordingly in PrairieView and run mark points.
 *	Click 'Run' in pyRTAOI.

See Demo_photostim_sequence.mp4 

You can use the results from this session to find which cells are photo-excitable. In our example, run OnlinePhotoexcitability.m, then load the result file in the photoexcitability test section in OnlineProcTexture to filter the target cells by there photostimulation response.

<img width="924" alt="image" src="https://github.com/user-attachments/assets/bff1ed9f-0e42-4e44-9d82-0db72a047a7f" />


### 6.	Closed-loop all-optical experiment

You can load different experimental configurations or setup your own. Here is how to stimulate specific ensembles of neurons based on the population activity trajectory (i.e. the weighted sum of the activity of a group of neurons) during behavior. See OnlineProcTexture.m for an example of how to generate the experiment configuration files (file names that contain ‘_PyBehaviour_’, ‘_RTAOiPyBehaviour’ and ‘_OutputParams_’). Again, we use pyBehaviour to setup and record the behavioural task.
 * Setup behavioural session in pyBehaviour (load the _PyBehaviour_ file).
 * Load configuration file (2020926_tex_photo_zoom2_16x.cfg)
 * Click 'Load stimOrder (Matlab)' and select the _RTAOiPyBehaviour_ file.

   <img width="563" alt="image" src="https://github.com/user-attachments/assets/e231b572-9ea5-40a1-b3c4-bc1a5677321c" />

 * Click 'Load trigger target config' and select the _OutputParams_ file

   <img width="501" alt="image" src="https://github.com/user-attachments/assets/e11d4b43-63ad-4e10-8fe9-e3637c72bb8b" />

 * Check 'Stim from blink'
 * Click 'Update Params' (make sure HoloBlink updates)
  <img width="557" alt="image" src="https://github.com/user-attachments/assets/3aa96a49-1989-4adb-b2a4-debe0df84502" />
 * Change t-series frames accordingly in PrairieView and run mark points.
 * Click Run
Note: you can abort the session at anytime by aborting the acquisition in PrairieView. pyRTAOI will stop the processing threads and save out result files automatically when no more frames are available in buffer.

<img width="752" alt="image" src="https://github.com/user-attachments/assets/c56eb05a-beac-4d60-bb8c-5548e60e14e0" />

See 
Demo_trajectory_photo.mp4

## Developers
* Zihui Zhang (UCL, Stanford)
* Patrycja Dzialecka (UCL, Imperial College London)

## License

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

