## RiverSTICH

_River-STICH_ (**S**urvey **T**ransect **I**nterpolation to reconstruct 3D **Ch**annels)

Updated on 6/10/2025

RiverSTICH converts traditional transect-based survey data into descriptive reach-scale attributes and variability functions and parameters that can then be used by [RiverBuilder](https://github.com/Pasternack-Lab/RiverBuilder) to construct a modular 3D synthetic river channel.


<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

* numpy
* pandas
* openpyxl
* simpledbf
* matplotlib
* scipy

<!-- USAGE EXAMPLES -->
## Usage Examples

Here, we present two examples using different types of XS survey data to demonstrate how RiverSTICH works.

#### Example 1. Auto level survey to RiverSTICH (main_SFE_Leggett.py)

- Input (/survey/SFE_Leggett)
    - Field survey data sheet (SFE_Leggett.xlsx, see [Survey_protocols.docx](/survey/SFE_Leggett/Survey_protocols.docx) for more information)
        - Equal-space transect survey
            <p align="center" width="100%">
            <img width="80%" src="/survey/SFE_Leggett/survey1.png" alt="input1">
            </p>
        - Longitudinal profile survey
            <!---![Figure 2.](/survey/SFE_Leggett/survey2.png)-->
            <p align="center" width="100%">
            <img width="60%" src="/survey/SFE_Leggett/survey2.png" alt="input2">
            </p>
        - Additional riffle crests and pool troughs transect survey
            <!---![Figure 3.](/survey/SFE_Leggett/survey3.png)-->
            <p align="center" width="100%">
            <img width="80%" src="/survey/SFE_Leggett/survey3.png" alt="input3">
            </p>
- Output (/output/SFE_Leggett)
    - X-Y contour plot (before and after transformation)
        - Note that the transformation was done to make the left and right bank margins symmetric (Black dots: thalweg, Blue dots: left and right bank margins)
            <!---![Figure 4. X-Y contour plot, before transformation]( /output/SFE_Leggett/XY_before_transformation.png)-->
            <p align="center" width="100%">
            <img width="60%" src="/output/SFE_Leggett/XY_before_transformation.png" alt="output1">
            </p>
            <!---![Figure 5. X-Y contour plot, after transformation](/output/SFE_Leggett/XY.png)-->
            <p align="center" width="100%">
            <img width="60%" src="/output/SFE_Leggett/XY.png" alt="output2">
            </p>
    - X-Y and X-Z interpolated contour plot 
            <!---![Figure 6. X-Y and X-Z interpolated contour plot](/output/SFE_Leggett/XYZ_contours.png)-->
            <p align="center" width="100%">
            <img width="60%" src="/output/SFE_Leggett/XYZ_contours.png" alt="output3">
            </p>
    - A channel attribute table of RiverSTICH channel, including reach-average bankfull depth, bed slope, and bankfull water surface elevation slope (channel_attributes.xlsx)
    - Interpolated contour series (SFE_Leggett_RB_metrics.xlsx)
        - These geomorphic variability functions (GVFs) will be used for RiverBuilder channel generation.

#### Example 2. X, Y, Z topographic survey to RiverSTICH (main_M1.py)

- Input (/survey/M1)
    - Field survey data 
        - Point shape file for topography (M1.shp)
            <!-- ![Figure 11.](/survey/M1/M1.png) -->
            <p align="center" width="100%">
            <img width="80%" src="/survey/M1/M1.png" alt="input2">
            </p>
        - Point shape file water surface elevation for baseflow (M1_base.shp)
            - This could be obtained from 2d hydraulic model output, through flume experiments or surveying
            <!-- ![Figure 12.](/survey/M1/M1_base.png) -->
            <p align="center" width="100%">
            <img width="80%" src="/survey/M1/M1_base.png" alt="input2">
            </p>
- Output (/output/M1)
    - Cross-section bed and water surface profile for width extraction
        <!--- ![Figure 13.](/output/M1/XS/x_0.png) -->
        <p align="center" width="100%">
        <img width="60%" src="/output/M1/XS/x_0.png" alt="output1">
        </p>
    - X-Y contour plot (XY.png)
        - Note that the transformation was not needed as the exact X, Y, Z coordinates were provided as input
        <!--- ![Figure 14.](/output/M1/XY.png) -->
        <p align="center" width="100%">
        <img width="60%" src="/output/M1/XY.png" alt="output2">
        </p>
    - X-Y and X-Z interpolated contour plot
        <!--- ![Figure 15.](/output/M1/XYZ_contours.png)-->
        <p align="center" width="100%">
        <img width="60%" src="/output/M1/XYZ_contours.png" alt="output3">
        </p>
    - A channel attribute table of RiverSTICH channel, including reach-average bankfull depth, bed slope, and bankfull water surface elevation slope (channel_attributes.xlsx)
    - Interpolated contour series (M1_RB_metrics.xlsx)
        - These geomorphic variability functions (GVFs) will be used for RiverBuilder channel generation.
        
<!---
<p align="center" width="100%">
<img width="50%" src="/SFE_Leggett_hand_param_calc/HAND_BM/SRCs_extended.png" alt="output3">
</p>
-->


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.



<!-- CONTACT -->
## Contact

Anzy Lee anzy.lee@usu.edu

GitHub repository: [https://github.com/USU-WET-Lab/RiverSTICH](https://github.com/USU-WET-Lab/RiverSTICH)


<!-- ACKNOWLEDGMENTS 
## Acknowledgments


<p align="right">(<a href="#readme-top">back to top</a>)</p>

