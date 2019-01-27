# OPM processing Guide for SPM

The code in this toolbox can be used to create SPM MEEG objects from an arbitrary  data source or to simulate MEG data. The main function used in this toolbox is  `spm_opm_create`. Examples are given of how this function might be used in various contexts.  



## Table of contents
1. [Preliminaries and Warranty](#a)
2. [Data Import](#b)
	1. [Sensor level OPM data](#b1)
	2. [Source level OPM data](#b2)
3. [Preprocessing](#d)
	1. [Filtering](#d1) 
	2. [Synthetic Gradiometry](#d2)
	3. [Epoching](#d3)
4. [Sensor Space Analysis](#e)
	1. [Evoked Response](#e1)
	2. [Time Frequency Analysis](#e2)
5. [Source Space Analysis](#f)
	1. [Dipole Fitting](#f1)
	2. [Distributed Source Inversion on a Mesh](#f2)
	3. [Distributed Source Inversion in a Volume](#f3)
6. [A Whole Script](#h)
7. [Simulation](#c)
	1. [MNI space ](#c1)
	2. [Whole head MNI space](#c2)
	3. [Individual Subject](#c3)
	4. [Individual Subject with fixed Sensor Positions](#c4)
	5. [Individual Subject with Custom Cortical Mesh](#c5)	
8. [Contributing Code](#g)



<a name="a"></a>
## Preliminaries and Warranty

For  this code to run you must have [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) added to Matlab's path. To run these examples you will need to run this code snippet. This will add SPM12 and the simulation toolbox to the path. It will also change the directory to the [test data folder](https://github.com/tierneytim/OPM/tree/master/testData).


```matlab
%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM')
dir = 'OPM\testData';
cd(dir)
```

However, before you run this code you should know that it is licensed under a GNU license which means no warranty or liability. Please read the [License](LICENSE) before use.




<a name="b"></a>
## Data Import



<a name="b1"></a>
### Sensor level OPM data

While OPM data may come in many native formats currently the UCL native file format is a simple binary file that contains the magnetometer output. In order to read this file and assign appropriate labels, units and channel types to the dataset some metadata is required. This should be provided in the form of a channels.tsv and a meg.json file. These files should conform to the standards recommended by  [BIDS](https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/02-magnetoencephalography.html) specification for MEG. Code to create a dataset suitable for sensor level analysis is given below. The files can be found in the [test data folder](https://github.com/tierneytim/OPM/tree/master/testData).

```matlab
S =[];
S.data = 'meg.bin';
S.channels='channels.tsv';
S.meg='meg.json'
D = spm_opm_create(S);
```

### Source level OPM data
To go further than a sensor level analysis more metadata is required. Two files are necessary to achieve this aim. The first is a 'positions.tsv' file. This should give coordinates and orientations of the sensors. The coordinate system of these sensors should be defined in the coordsystem.json file. The format of which is described in [BIDS specification](https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/02-magnetoencephalography.html) and example file for both of these is available in the [test data folder](https://github.com/tierneytim/OPM/tree/master/testData). Lastly, a structural MRI is required to compute the MEG forward model. Code is given below as an example of how these files can be utilised.

```matlab
S =[];
S.data = 'meg.bin';
S.coordystem='coordsystem.json';
S.positions='positions.tsv';
S.channels='channels.tsv';
S.meg='meg.json';
S.sMRI='T1w.nii';
D = spm_opm_create(S);
```


<a name="d"></a>
## Preprocessing

A number of preprocessing steps can be optionally applied to OPM data in any order you want(similar to any MEG dataset). Here is an example of one pipeline.


<a name="d1"></a>
### Filtering

Temporal filtering can be applied to remove the effects of high or low frequency interference. In this case we filter between 1 and 80 Hz using a 2nd order butterworth filter. The following code snippet should filter the data we created in the previous step. 

```matlab
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [1 80];
S.dir = 'twopass';
S.order = 2;
D = spm_eeg_filter(S);
```
<a name="d2"></a>
### Synthetic Gradiometry


OPMs are magnetometers and not gradiometers which makes them somewhat more susceptible to environmental interference than gradiometers. To mitigate this effect we construct `Synthetic Gradiometers` by regressing the signal in a set of reference sensors from the signal in the scalp sensors. This can be done on a trial by trial basis or across the whole scanning session. If epoched data is supplied the it will take place on a trial by trial basis.  In this case we will use the whole session data An example of how you might do this is given below. 


``` matlab
S=[]
S.D=D;
S.confounds={'REF'};
D = spm_opm_synth_gradiometer(S);
```

Crucially this function works in an object oriented fashion. The S.confounds argument searches for channel types that with the corresponding label and models these channels as effects of no interest. This framework should therefore be easily extended to include motion estimate from optical tracking cameras simply by changing the S.confounds variable. 

<a name="d3"></a>
### Epoching

To epoch the data we just need to tell SPM the time windows(in ms) we want to analyze and SPM will extract the data around all triggers in order to epoch the data. In order for this to work there needs to be at least 1 channel in the MEEG object that has the type `TRIG`. 


``` Matlab
S =[];
S.D=D;
S.timewin=[-100 300];
S.condLabels= {'Median Nerve'};
D= spm_opm_epoch_trigger(S);
```


<a name="e"></a>
## Sensor Space Analysis

<a name="e1"></a>
### Evoked Responses

Once the data is preprocessed calculating an evoked response  is easy. It just requires a few lines of code which are given below resulting in a nice pretty picture.


```matlab
S =[]
S.D=D;
D =spm_eeg_average(S);
```
We can also apply a baseline correction to this dataset as well

```matlab
S=[];
S.D=D;
S.timewin=[-100 -20];
D = spm_eeg_bc(S);
```

<p align="center">
<img src="readme/evokedresponse.PNG" width="600"/>
</p>


<a name="e2"></a>
### Time Frequency Analysis
Coming soon!

<a name="f"></a>
## Source Space Analysis

<a name="f1"></a>
### Dipole Fitting
Coming soon!

<a name="f2"></a>
### Distributed Source Inversion on a Mesh
Coming soon!

<a name="f3"></a>
### Distributed Source Inversion in a Volume
Coming soon!


<a name="h"></a>
## An example script

Here we put all the steps together into 1 script for you to try out.

```matlab
%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM')
dir = '\OPM\testData';
cd(dir)
%% read data
S =[];
S.data = 'meg.bin';
S.sMRI='T1w.nii';
D = spm_opm_create(S);

%% filter the data
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [1 80];
S.dir = 'twopass';
S.order = 2;
D = spm_eeg_filter(S);

%% denoising
S=[];
S.D=D;
S.confounds={'REF'};
D = spm_opm_synth_gradiometer(S);

%% epoch the data
S =[];
S.D=D;
S.timewin=[-100 300];
D= spm_opm_epoch_trigger(S);

%% Average
S =[];
S.D=D;
D =spm_eeg_average(S);

%% baseline correct
S=[];
S.D=D;
S.timewin=[-100 -20];
D = spm_eeg_bc(S);

```

<a name="c"></a>
## Simulation


<a name="c1"></a>
### MNI Space
If you want to simulate MEG data you need to supply sensor positions and  orientations which should be in the same coordinate space as some brain image. If you don not have a brain image or sensor positions(and orientations) you can simulate data on an average template brain with automatically generated positions of fixed spacing using `spm_opm_create`. The following code snippet automatically generates sensors in this average space that are a fixed `space` apart. In this case the spacing is 15mm.

```matlab
S =[];  
S.space = 15;  
D = spm_opm_create(S);  
```
Once you run this code you will generate a figure like this. 

<p align="center">
<img src="readme/mni15.png" width="600"/>
</p>


<a name="c2"></a>
### Whole-Head MNI space
For some simulations you may want the entire scalp surface to be covered. In this case you just need to set the `wholehead` flag to 1.
```matlab
S =[];
S.space = 15;
S.wholehead=1;
D = spm_opm_create(S);
```
Running this code snippet should generate a figure like this.

<p align="center">
<img src="readme/mni15Whole.PNG" width="600"/>
</p>

<a name="c3"></a>
### Individual Subject
In some cases you may have already have an individual brain image that you want to simulate data in but no positions or orientations. In this case you just need  to give your desired sensor spacing with the `space` argument and the filepath to an MRI file with the `sMRI` argument. The slightly different orientation of this brain is due to the fact that this brain is not in the MNI space.

```matlab
S =[];
S.space = 15;
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);
```

<p align="center">
<img src="readme/native15.PNG" width="600" height="452"/>
</p>

<a name="c4"></a>
### Individual subject with fixed sensor positions
If you already have positions and orientations for a certain brain image this can be accounted for as well. This information can be supplied in the form of a tab-delimited text file where the first six columns give position and orientation(x,y,z,x,y,z) and the final column gives a label for the sensor. 
The code to incorporate this information just requires to specify the `pos argument 

```matlab
S =[];
S.pos = 'SEF_coarse';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);
```
<p align="center">
<img src="readme/coregSensors.PNG" width="600"/>
</p>

<a name="c5"></a>
### Individual Subject with Custom Cortical Mesh
Other times you may not wish to use the default meshes supplied by SPM and may wish to provide custom meshes. This is easily done with the `cortex` argument. Any mesh can be supplied but it is up to the user to ensure the meshes are in the same coordinate  space.
 ```matlab
S =[];
S.space = 15;
S.cortex='testCustom.gii';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

```

<p align="center">
<img src="readme/dragonNative.PNG" width="600" height="449" />
<p>

<a name="g"></a>
## Contributing Code
Contributions are always welcome. However, there are a few general rules that should be 
followed if you want to contribute code. 

- Main functions should begin with `spm_opm`. 
- They should accept one argument: `S`. 
- All arguments must be described in a help section within the file.  
- The default values should be clearly identified. 
- The use of `nargin` and  `varargin` are strictly forbidden. 
- The return value should usually be an MEEG object called `D` 
- A copy & paste example must be provided 
- No line of code should be longer than 75 columns(characters) in Matlab. 

For helper functions all these rules can be relaxed bar the use of `nargin`, `varargin` 
and the 
number of characters per line of code.
