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

See ‘Demo_ClosedLoopAllOpticalTexture’ for detailed steps for running a closed-loop all-optical experiment during beheviour.

**To start pyRTAOI**:

```
conda activate caiman # activate virtual environment
python pyRTAOI.py # open GUI
```


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

