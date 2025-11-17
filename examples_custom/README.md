## River Builder with custom Geomorphic Variability Functions (GVFs)

Updated on 11/16/2025

In [River Builder](https://github.com/Pasternack-Lab/RiverBuilder) (RB), GVFs can be implemented in one of two ways: (1) by directly specifying mathematical functions for each GVF or (2) generating a straight River Builder channel based on reach-average attributes, and then applying GVFs externally by imposing GVF series from [RiverSTICH](https://github.com/USU-WET-Lab/RiverSTICH) to the River Builder output (e.g., table of x, y, z coordinates of topography). Here, we present the second approach which performs better for a channel with a highly sinuous thalweg function.


<!-- GETTING STARTED -->
## Getting Started

Here, we present an example developed for a V-shape River Builder terrain with feature-based XS survey for SFE Leggett to demonstrate how this approach works.

### Prerequisites

* numpy
* pandas
* matplotlib
* scipy
* arcpy (for RB_to_terrain_custom.py)

<!-- USAGE EXAMPLES -->
## Workflow

1. Write a RB txt input file to generate a straight RB channel based on reach-average attributes by running RB_gcs_custom.py (in /tools)
- Input (/tools/gcs)
    - SFE_Leggett_FB_V.xlsx
        - Interpolated contour series (e.g., SFE_Leggett_RB_metrics.xlsx)
        - This is the final output of the first example of [RiverSTICH](https://github.com/USU-WET-Lab/RiverSTICH).
- Output (/examples_custom/SFE_Leggett_FB_V)
    - SFE_Leggett_FB_V.txt
        - The main input txt file for River Builder with reach-average attributes
     
2. Generate a straight RB channel by running RB_run_custom.py (in /tools)
- Input (/examples_custom/SFE_Leggett_FB_V)
    - SFE_Leggett_FB_V.txt
        - The output of previous step
- Output (/examples_custom/SFE_Leggett_FB_V/SFE_Leggett_FB_V)
    - SRVtopo.csv
        - The main RB output (e.g., table of x, y, z coordinates of topography)


            <!-- ![Figure 1.](/SFE_Leggett_FB_V/SFE_Leggett_FB_V/SRVlevels_xy.png) -->
            <p align="center" width="100%">
            <img width="80%" src="/SFE_Leggett_FB_V/SFE_Leggett_FB_V/SRVlevels_xy.png" alt="input2">
            </p>

            <!-- ![Figure 2.](/SFE_Leggett_FB_V/SFE_Leggett_FB_V/SRVlevels_xz.png) -->
            <p align="center" width="100%">
            <img width="80%" src="/SFE_Leggett_FB_V/SFE_Leggett_FB_V/SRVlevels_xz.png" alt="input2">
            </p>
   

        
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

