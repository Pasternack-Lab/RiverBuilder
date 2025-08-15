# RiverBuilderv1.3 <img src="https://github.com/RiverBuilder/RiverBuilder/blob/master/logo.png" width="24">

## Dependencies

* [Python3](https://www.python.org/downloads/)
* [numpy](https://numpy.org/devdocs/user/quickstart.html)
* [matplotlib](https://matplotlib.org)

*These dependencies are already satisfied in GIS, so no worries about them if you are using GIS.*

## Quick Start
Please refer `riverbuilder_v1.3_Input_template.txt` for getting an idea on preparing your input file.

The main command to run the script is `python3 -m riverbuilder <inputfile.txt> <outputFolder> <your message>`. 

* `<inputfile.txt>` provides essential parameters for the program to run. Its format is explained in *riverbuilderInput_format.txt* file, and there are plenty of example input files in the `/examples_v1pt2` folder.

* `<outputFolder>` sets the path of script output files. It is optional. If this argument is omitted, program will create a `riverbuilder_output` folder in your current working folder, and store output files there.

* `<your message>` is a message will be added to the output. It is also optional. If `<outputFolder>` is not set, this can't be set.

More explanations can be found in [manual](https://github.com/RiverBuilder/RiverBuilder/blob/master/RiverBuilder_User_Manual_1.2.0_final.pdf).

### In MacOS or Linux

Type the following commands in a terminal.

```
# Download the repository
git clone https://github.com/RiverBuilder/RiverBuilder
cd RiverBuilder

# Install dependencies
pip3 install numpy
pip3 install matplotlib

# Run the script
python3 -m riverbuilder samples/S1/S1.txt samples/S1/S1 "A simple test sample."
```

And you should see the following screen output.

```
Start Parsing Inputs...
Start Creating River Channel...

It takes 0 seconds to build the channel.
Start Creating Valley...
It takes 0 seconds to build the valley.

Output files are saved to <Path to your folder>/RiverBuilder/samples/S1/S1
A simple test sample.
Given input file is: samples/S1/S1.txt
Given output folder is: samples/S1/S1.txt

Datum is set to: 1000.0
Length is set to: 420.0
Alert! X Resolution is not defined by user, use default value 1 instead.
X Resolution is set to: 1
Alert! Channel XS Points is not defined by user, use default value 21 instead.
Channel XS Points is set to: 21
Valley Slope (Sv) is set to: 0.002
Critical Shields Stress (t*50) is set to: 0.0
Inner Channel Lateral Offset Minimum is set to: 30.0
Inner Channel Depth Minimum is set to: 2.11
Median Sediment Size (D50) is set to: 0.0
Left Valley Boundary Lateral Offset Minimum is set to: 10.0
Left Valley Boundary Height Offset is set to: 20.0
Alert! Right Valley Boundary Lateral Offset Minimum is not defined by user, use default value 10 instead.
Right Valley Boundary Lateral Offset Minimum is set to: 10
Alert! Right Valley Boundary Height Offset is not defined by user, use default value 20 instead.
Right Valley Boundary Height Offset is set to: 20
Alert! Inner Channel Average Bankfull Width is not defined by user, use default value None instead.
Inner Channel Average Bankfull Width is set to: None
Alert! Inner Channel Average Bankfull Depth is not defined by user, use default value None instead.
Inner Channel Average Bankfull Depth is set to: None
Alert! River Slope is not defined by user, use default value None instead.
River Slope is set to: None
User defined funtions are:
MASK0 (ALL)
SIN1 (2.5, 1, 0, MASK0)

Alert! Can't find function Meandering Centerline Function in user-defined functions. Ignore function Meandering Centerline Function.
Reshape not needed for river centerline.
Creating Meandering Center line with Function:None
Use user defined Inner Channel Depth Minimum.
Alert! Can't find function Centerline Curvature Function in user-defined functions. Ignore function Centerline Curvature Function.
Creating Inner Channel Banks with left bank function: SIN1
                             with right bank function: SIN1
                             with thalweg elevation function: SIN1

Alert! Can't find function Valley Centerline Function in user-defined functions. Ignore function Valley Centerline Function.
Alert! Can't find function R1 Valley Breakline Function in user-defined functions. Ignore function R1 Valley Breakline Function.
Creating R1 Valley Breakline with constant width: 25.0
Alert! Can't find function L1 Valley Breakline Function in user-defined functions. Ignore function L1 Valley Breakline Function.
Creating L1 Valley Breakline with constant width: 25.0

Sinuosity:1.0
Channel Slope:0.002
Average Width of Inner Channel:35.0
Average Height of Inner Channel:4.823
Coefficient of Variation (W_ic):0.101
Coefficient of Variation (H_ic):0.406
```

### In Windows

1. Download the github folder.
2. Open a terminal and switch to RiverBuilder folder.
3. Run the following command:
```
<YOUR_PATH_TO_PYTHON3>\python.exe -m riverbuilder samples\S1\S1.txt samples\S1\S1 "A simple test sample."
```
4. A similar output should be printed.

## Documentation

The full manual is available [here](https://github.com/RiverBuilder/RiverBuilder/blob/master/RiverBuilder_User_Manual_1.2.0_FINAL_Rev1.pdf).

## License

MIT

## Disclaimer

No warranty is expressed or implied regarding the usefulness or completeness of the information contained in this manual or the associated software. References to commercial products do not imply endorsement by the authors. The concepts, materials, and methods presented in this manual are for informational purposes only. The authors have made substantial effort to ensure accuracy, but there is uncertainty and the authors shall not be held liable for calculations and/or decisions made on the basis of application of this software and manual. The information is provided "as is" and anyone who chooses to use the information is responsible for his or her own choices as to what to do with the data.

## Acknowledgment

The underlying research used to develop River Builder came from a variety of sources, including UC Davis, the USDA National Institute of Food and Agriculture [Hatch project number #CA-D-LAW-7034-H], the Ticho Foundation (unrestricted donation), US Army Corps of Engineers, and Yuba County Water Agency (Award #201016094). Dr. Rocko Brown and Dr. Greg Pasternack also contributed their own resources in developing different aspects of the underlying research. We thank postdoctoral scholars Dr. Colin Byrne and Dr. Sebastian Schwindt for conversations about various aspects of this version of the software as we were developing it. Other members of the Pasternack lab group gave input from time to time during 2019 and 2020.
