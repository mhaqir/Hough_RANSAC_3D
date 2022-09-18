This repository contains the implementations for the 3D RANSAC and Hough-transform algorithms introduced in the paper: Hough-transform and Extended RANSAC Algorithms for Automatic Detection of 3D Building Roof Planes From Lidar Data.



I implemented a code to generate 3D points on a hyperplane. Given a desired number of 3D points and a hyperplane parameters, a set of 3D points that represent a hyperplane will be generated. There is a parameter epsilon that determines the boundaries of noise that will be added to the 3D points. The noise will be added independently for each dimension and for each data point.


##### RANSAC output for hyperplane detection
<!-- ![eps_0.05](https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac3d_outputs/noisy_0.05.jpg =250x250) ![eps_0.5](https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac3d_outputs/noisy_0.5.jpg =250x250) -->

Each row shows an original hyperplane with two different maximum noise value and the detected hyperplane.

<img src="https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac3d_outputs/noisy_0.05.jpg" width="350" height="350" /> <img src="https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac3d_outputs/noisy_0.5.jpg" width="350" height="350" />

<img src="https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac3d_outputs/ransac_a_2_b_1_c_-3_d_2_eps_0.05.jpg" width="350" height="350" /> <img src="https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac3d_outputs/ransac_a_2_b_1_c_-3_d_2_eps_0.5.jpg" width="350" height="350" />

For evaluating the output, I also calculated the angle between the normal vector of the original and detected hyperplanes. For the hyperplane 2x + y - 3z + 3 = 0, two sets of points that one contains no noise and the other with 0.05 maximum noise in each dimension were generated, see the output [here](https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/ransac_a_2_b_1_c_-3_d_3_eps_0.05.txt).

##### Hough-transform output for hyperplane detection
See the output [here](https://github.com/mhaqir/Hough_RANSAC_3D/blob/main/Hough_a_2_b_1_c_-3_d_3.txt). For this algorithm, I reported all of the detected hyperplanes with their confidence score which is the number of points that supports that hyperplane.
