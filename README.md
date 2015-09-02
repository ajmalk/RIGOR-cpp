RIGOR: Reusing Inference in Graph Cuts for generating Object Regions
====================================================================
Ahmad Humayun, Fuxin Li, James M. Rehg  
Georgia Institute of Technology  
[http://cpl.cc.gatech.edu/projects/RIGOR/](http://cpl.cc.gatech.edu/projects/RIGOR/)
--------------------------------------------------------------------
IEEE Computer Vision and Pattern Recognition (CVPR) 2014

This code is accompanied with a [LICENSE](rigor/src/master/LICENSE)

--------------------------------------------------------------------

This fork is maintained by Ajmal Kunnummal

RIGOR is a cutting edge algorithm that generates overlapping segment proposals in images that can then be fed into image recognition algorithms to speed up object recognition in images. RIGOR was proposed and implemented in 2014 and at the time was the best performing such algorithm. The long-term goal is to design a GPU driven highly parallelized min-cut algorithm that can produce high quality object proposals for real time object recognition. For this to be feasible, the ‘setup’ stages of the algorithm where the graphs are generated from the images have to be extremely fast. Initially the setup stages were done in Matlab and took over a second to run. We have improved the running time of the algorithm drastically by implementing the ‘setup’ part in C++.

This fork contains the current C++ version of Rigor. More information is contained in the paper. 