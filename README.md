# Welcome to CUHK Seismology Group Shared Repository

## Guidelines to contribute the repository
**_NOTE:_** Please upload minimum amount of (waveform) data to keep the repository light. 

Demos:
  Demo use of builded packages/scripts
 
 ---
### Basic Usage of Github
Git(hub) is a code version-control platform to maintain the development and collaborations of code. The following shows the basic instruction to use and contribute the project.

#### First-time user
**_NOTE:_** This repository is for sharing of codes. Please do not work directly on the repository. It is suggested to make a copy of the useful code. Or else, your work will be shared in here.
Setup a "cloned" version of the repository to your computer

* git clone [url]

Synchronize the latest code from github to your local computer repository

* git pull

Contribute after amending codes in your local machine

* git add .
* git commint -m "msg of the commit"
* git push

### Common python package requirement
obspy
numpy
pandas
matploblib

#run setup.py to add the path to your python environment

---
### List of useful packages for seismology analysis

- Earthquake Detection
  - [PhaseNet](https://github.com/wayneweiqiang/PhaseNet.git)
  - [GeneralizedPhaseDetection](https://github.com/interseismic/generalized-phase-detection)
  - EQscancorr
  
- Earthquake Phase Association and Location (& Relocation)
  - [REAL: Rapid Earthquake Assocation](https://github.com/Dal-mzhang/REAL.git)
  - Velest
  - Hypoinverse
  - HypoDD
  - Glowclast
  
- Earthquake Source Study
  - [HASH - First-motion DC inversion](Earthquake_Source)
  
- Analysis
  - topodd
  - Abient Noise
    - [MSnoise](http://www.msnoise.org/)
    
- Plotting
  - GMT
  - PYGMT
  - Matplotlib

- Data Management
  - ASDF
  - Antelope

- Data Source
  - [IRIS](Data_Management/IRIS_fetch/)
 

---
## Contributors
Last but not least, thanks all for contribution the community
*   @zijinping
*   @jwjeremy
