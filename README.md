# FreeTSE

Calibration-free traffic state estimation method using detectors and connected vehicles data

## Executable version with manual

[Download](https://github.com/toruseo/FreeTSE/releases/download/v1.0.0/FreeTSE_executable.zip)

## How to use

- Executable version (for Windows, non-Python user): Run `run.bat` in the executable version.
	- Select an INI file that defines an estimation scenario.
- Code version (for Python user): Run or import `FreeTSE.py`.

### Executable version: Details

An estimation scenario should be defined by an INI file.
It should specify the input data files and the spatiotemporal resolution of the estimation.
The format can be easily understood by looking at the sample files.

An annotated example of INI file is below:
<details>
<summary>Click here to see the code</summary>

```ini
[Data]
name = ngsim_trajectory
;name of the scenario
speed = ./dat/ngsim_sampled_trajectories.csv
;path to speed data
speed_label_t = t
;column name for time
speed_label_x = x
;column name for position
speed_label_v = v
;column name for speed
flow = ./dat/ngsim_grid_flow_200m.csv
;path to flow data. if it does not exist, specify "None"
flow_label_t = t
flow_label_x = x
flow_label_q = q

[GroundTruth]
;This section specify ground truth data for validation. if it does not exist, delete it.
true_flow = ./dat/ngsim_grid_flow_400m.csv
;path to ground truth flow data for validation
flow_label_t = t
flow_label_x = x
flow_label_q = q

[Unit]
time = s
;not functioning
space = m
;not functioning

[Resolution]
time_min = 0
;initial time. unit: seconds
time_max = 800
;last time
space_min = 0
;upstream-end position. unit: meters
space_max = 500
;downstream-end position
dt = 4
;temporal resolution of estimation
dx = 100
;spatial resolution of estimation
;it is very preferable that dx/dt is almost equal  to the maximum speed of the traffic
number_of_lanes = 5
;not essential
```
</details>

### Code version: Details

FreeTSE module can be used with the following code:
<details>
<summary>Click here to see the code</summary>

```python
from FreeTSE import *

tse = FreeTSE()
tse.set_scenario(
	name = "ngsim_trajectory",	#name of the scenario
	dt = 4,		#temporal resolution of estimation. unit: seconds
	dx = 100,	#spatial resolution of estimation. unit: meters
	mint = 0,	#initial time
	maxt = 800,	#last time
	minx = 0,	#upstream-end position
	maxx = 500,	#downstream-end position
	number_of_lanes = 5,	#not essential
	speed_data_name = "./dat/ngsim_sampled_trajectories.csv",	#path to speed data
	speed_label_t = "t",	#column name for time
	speed_label_x = "x",	#column name for position
	speed_label_v = "v",	#column name for speed
	density_data_name = None,	#path to density data. if it does not exist, specify "None"
	density_label_t = "t",
	density_label_x = "x",
	density_label_k = "k",
	flow_data_name = "./dat/ngsim_grid_flow_200m.csv",	#path to flow data. if it does not exist, specify "None"
	flow_label_t = "t",
	flow_label_x = "x",
	flow_label_q = "q",
	density_dat_true_name = None,	#path to ground truth density data for validation. if it does not exist, specify "None"
	true_density_label_t = "t",
	true_density_label_x = "x",
	true_density_label_k = "k",
	flow_dat_true_name = "./dat/ngsim_grid_flow_400m.csv",	#path to ground truth flow data for validation. if it does not exist, specify "None"
	true_flow_label_t = "t",
	true_flow_label_x = "x",
	true_flow_label_q = "q"
)
tse.estimation()
tse.accuracy_evaluation()	#evaluate the accuracy if ground truth data exists
fname = "res_test"
tse.save_results(fname+".csv")	#save the results as CSV
tse.visualize(smooth=1, true=1, speed=1, timeseries=1, inputdata=1, save=1, fname=fname)	#visualize the results using matplotlib

q, k, v = tse.get_results()	#get the estimation results as 2 dimensional numpy.array. the 1st index of the array is time and the 2nd is position.
print("flow", q)
print("density", k)
print("speed", v)
```
</details>

Files: 
- `FreeTSE.py`: The main code
- `dat/*.ini`: Sample scenarios
- `dat/*.csv`: Sample datasets
- `gui.py`: GUI to run `FreeTSE.py`
- `util.py`: Miscellaneous utilities

## Data formats

The unit of variables must be based on seconds and meters (e.g., the speed must be in m/s and the flow must be in vehicles/s).
The numerical value of the position should be set so that the upstream position is represented by a smaller value.

### Speed data from connected vehicles (CSV)

The first row should be a header row that describes data type of each column.
The second row and below are the main body of data.
It should contain the traffic state variable (speed) and the time and position where the data was collected.

Speed data is assumed to be continuously collected over the entire area.
Missing data will be automatically interpolated.
You can either use trajectory data or aggregated average speed data.

A dataset may look like:
```
t,   x,     v
0.0, 379.6, 16.3
1.0, 396.3, 16.5
2.0, 413.0, 16.8
3.0, 429.7, 16.7
4.0, 446.5, 16.7
5.0, 462.9, 16.4
(cont.)
```

### Flow data from detectors (CSV)

The first row should be a header row that describes data type of each column.
The second row and below are the main body of data.
It should contain the traffic state variable (flow or density) and the time and position where the data was collected.

Flow data is assumed to be continuously collected at a limited number of locations.
Alternatively, you may use density data.

A dataset may look like:
```
t,  x,   q
0,  200, 2.41
4,  200, 2.42
8,  200, 2.88
12, 200, 2.70
16, 200, 2.33
20, 200, 2.59
(cont.)
```

## Mechanism

It estimates flow and density in the entire time-space domain based on speed data collected by connected vehicles and flow (or density) data collected by limited number of detectors.
It utilizes a parameter-free data-driven traffic flow model based on the conservation law and Kalman filtering and smoothing.
For the details, see
- Seo, T. Calibration-free traffic state estimation method using single detector and connected vehicles with Kalman filtering and RTS smoothing. In IEEE 23rd International Conference on Intelligent Transportation Systems, Web conference, 2020. [link](http://dx.doi.org/10.1109/ITSC45102.2020.9294229)

## Credits/Acknowledgements

- The sample files are based on the NGSIM dataset
- The executable version includes the embeddable version of Python
- Part of this work was financially supported by the Japan Society for the Promotion of Science (KAKENHI Grant-in-Aid for Young Scientists (B) 16K18164, KAKENHI Grant-in-Aid for Scientific Research (B) 20H02267)

## Further reading

- Seo, T. Calibration-free traffic state estimation method using single detector and connected vehicles with Kalman filtering and RTS smoothing. In IEEE 23rd International Conference on Intelligent Transportation Systems, Web conference, 2020.
- Seo, T., Bayen, A. M., Kusakabe, T., and Asakura, Y. Traffic state estimation on highway: A comprehensive survey. Annual Reviews in Control, Vol. 43, pp. 128-151, 2017.

## License

This software is released under the MIT License, see LICENSE.txt.

## Author

Toru Seo (Tokyo Institute of Technology)
- https://toruseo.github.io/
- http://seo.cv.ens.titech.ac.jp/
