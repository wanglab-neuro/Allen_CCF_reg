<p align="center">
  <img src="./misc/demo_header.png" alt="slider" width="800px"/>
</p>

[# Constructing An Adult Orofacial Premotor Atlas In Allen Mouse CCF](https://www.biorxiv.org/content/10.1101/2021.02.18.431923v1.abstract). 


## Contents
1. [Updates](#updates)
2. [Dependencies](#dependencies)
3. [Usage](#usage)
4. [References](#reference)
<!-- 5. [Questions](#questions) -->

## Dependencies
* Image Processing Toolbox
* Curve Fitting Toolbox
* [optional] Parallel Computing Toolbox
* The function [TVL1denoise](https://www.mathworks.com/matlabcentral/fileexchange/57604-tv-l1-image-denoising-algorithm) is provided in a subdirectory and needs to be in the path.

## Usage
0. Download the [Atlas.zip](https://drive.google.com/file/d/1-s8XfBQZxoolgtlJvwWK4daDHjxN6Jaa/view?usp=sharing). Unzip the archive and place `AllenAtlas.mat` in the `utils` folder.
1. Run `cell_count.m` (add code directory to Matlab's path if calling `cell_count` from another directory, e.g., your data folder)
2. In the atlas viewer, press keyboard arrows to switch between ML, DV and AP adjustment modes.
3. Use the mouse's scrolling wheel to make required adjustments in each axis.  

![image](https://user-images.githubusercontent.com/1872756/176929669-26424095-616f-4fc5-8177-1767ac6e7783.png)

4. Press "h" to draw a polygon around the outline you want the slice to be registered to in the atlas viewer (This is to improve the registration accuracy).

![image](https://user-images.githubusercontent.com/1872756/176953360-6d5f2d7c-ac58-4b84-8c18-c09823a4481d.png)

After completing the polygon by joining both ends, double-click to validate it. 

5. Press "p" to open and jump to the slice viewer.
6. 
  * Press "c" to enter manual mode and then click on the neurons. You will see an updated dot on the slice image if click succeeds.  
  or   
  * Press "y" to enter automatic cell detection mode.
7. Press "s" to save all the info.
  
## Reference
Please cite [the publication](https://doi.org/10.7554/eLife.67291) if it helps your research.
```
@article{10.7554/eLife.67291,
  title = {Constructing An Adult Orofacial Premotor Atlas In Allen Mouse CCF},
  author = {Takatoh, Jun and Park, Jae Hong and Lu, Jinghao and Li, Shun and Thompson, PM and Han, Bao-Xia and Zhao, Shengli and Kleinfeld, David and Friedman, Beth and Wang, Fan},
  doi = {10.7554/eLife.67291},
  url = {https://doi.org/10.7554/eLife.67291},
  journal = {eLife},
  year = 2021,
}
```
***Archive***
```
@article{takatoh2021constructing,
  title={Constructing An Adult Orofacial Premotor Atlas In Allen Mouse CCF},
  author={Takatoh, Jun and Park, Jae Hong and Lu, Jinghao and Li, Shun and Thompson, PM and Han, Bao-Xia and Zhao, Shengli and Kleinfeld, David and Friedman, Beth and Wang, Fan},
  journal={bioRxiv},
  year={2021},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Updates

7/1 Updates to procedure documentation

<!-- 
## Questions?
Please email to [`placeholder`](mailto:placeholder) for additional questions. -->
