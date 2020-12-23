# TIEGCM 2.1
High-resolution version of TIEGCM with horizontal resolution up to 0.625 degrees.

If you use this code, please cite the paper [Dang T. et al. GMD, (2020)](https://doi.org/10.5194/gmd-2020-243).

Please cantact the code owner Tong DANG(dangt@ustc.edu.cn) for further questions.

### BUILD & RUN
The compilation and run of TIEGCM 2.0r are based on the typical [NCAR TIEGCM2.0](https://www.hao.ucar.edu/modeling/tgcm/tie.php). 
Please refer to the [TIEGCM userguide](https://www.hao.ucar.edu/modeling/tgcm/tiegcm2.0/userguide/userguide.pdf) for the basic settings and test runs.

In TIEGCM2.1, one could choose the resolution of 1.25 or 0.625 degrees in the job script by
```
   set modelres = 0.625
```
or
```
   set modelres = 1.25
```

