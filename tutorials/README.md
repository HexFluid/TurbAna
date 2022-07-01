## Tutorials

We present four case studies in fluid mechanics using *TurbAna*:
- 2D backward-facing step [[1](#dns-bstep),[2](#ddes-bstep)]: **01_bstep.py**
- 2D transonic bump [[3](#les-bump)]: **02_bump.py**
- 2D turbine film cooling [[4](#ddes-cool)]: **03_cooling.py**
- 2D shock/boundary layer interaction [[5](#dns-sbli)]: **04_SBLI.py**

This page briefly introduces the flow physics of each case and the major results of *TurbAna*. The tutorial scripts are well-commented, self-explanatory, and easily transferrable to other cases. If you are using any data from the tutorial cases, please explicitly mention the corresponding publications.

---
### Case 1: 2D backward-facing step
The backward-facing step flow is composed of an attached boundary layer, a separated region (due to a sudden change of geometry) and a re-attached boundary layer. It is one of the standard test cases for DES-type methods. The data analyzed in this case is time-averaged and spanwise-averaged DNS data from Le et al. [[1](#dns-bstep)] and DDES data from He et al. [[2](#ddes-bstep)]. The Reynolds numbers for the DNS and the DDES data are *Re<sub>H</sub>=5100* and *Re<sub>H</sub>=36000*, respectively. The figure below illustrates the flow schematics.

<figure>
<p align="left">
    <img alt="bstep_schematics" src="../docs/figs/bstep_schematic.png" width="300" />
    <figcaption align="left"> Flow schematic of backward-facing step </figcaption>
</p>
</figure>

After reading the turbulence statistics data and defining it as a native *ReynoldsStressTensor* type, the turbulence anisotropy tensor and its eigenvalues are automatically calculated. The turbulence anisotropy can be visualized in Lumley Triangle, Turbulence Triangle and Barycentric Map by using methods of *.LumleyTriCoor()*, *.TurbTriCoor()* and *.BaryTriCoor()*, respectively. The figures below present an example visualization of the profile at x/H=4 and 0<y/H<1.

<figure>
<p align="left">
    <img alt="bstep_LumleyTri" src="bstep_data/results/x=4.0H_Lumley_tri.png" height="200" />
    <img alt="bstep_TurbTri" src="bstep_data/results/x=4.0H_turb_tri.png" height="200" />
    <img alt="bstep_BaryMap" src="bstep_data/results/x=4.0H_bary_map.png" height="200" />
    <figcaption align="left"> (left): Lumley triangle; (medium): turbulence triangle; (right): barycentric map </figcaption>
</p>
</figure>

It is also possible to visualize the turbulence anisotropy in the flow domain. By using the method of *.AniRGB(c_off,c_exp)*, the turbulence anisotropy state is represented by a color, which is indicated by a barycentric color map. With *c<sub>off</sub>=0.65* and *c<sub>exp</sub>=5*, the Reynolds stress anisotropy componentality contour and the barycentric colormap are obtained as the following figures.

<figure>
<p align="left">
    <img alt="bstep_AniCont" src="bstep_data/results/anisotropy_contour.png" height="200" />
    <img alt="bstep_AniContMap" src="bstep_data/results/anisotropy_colormap.png" height="120" />
    <figcaption align="left"> (left): Reynolds stress anisotropy componentality contour; (right): barycentric colormap </figcaption>
</p>
</figure>

With the mean flow field and the gradient field data read-in as a native *MeanFlowField* type variable and a *MeanGradField* type variable, the turbulent viscosity field can be derived. This is achieved by using the build-in function of *.calc_EddyVisc()* with additional input including a reference strain rate value and a constitutive relation (e.g., Bousinessq, QCR). The left figure below shows the derived turbulent-to-laminar viscosity ratio using the QCR and *S<sub>ref</sub>=30s<sup>-1</sup>*; the right figure shows the limiter *f<sub>lim</sub>* applied to the derivation, where *f<sub>lim</sub>=0,1,2,3* represent no limiter applied, *S<sub>ref</sub>* applied, *ν<sub>t</sub> ≥ 0* applied or both limiters applied, respectively.

<figure>
<p align="left">
    <img alt="bstep_NutrCont" src="bstep_data/results/sref=30_QCR2013V_ViscRatio_contour.png" height="200" align="center"/>
    <img alt="bstep_NutrLimitery" src="bstep_data/results/sref=30_QCR2013V_ViscRatio_limiter.png" height="200" align="center"/>
    <figcaption align="left"> (left): viscosity ratio contour; (right): limiter contour </figcaption>
</p>
</figure>

---
### Case 2: 2D transonic bump
The axisymmetric transonic bump flow is featured by streamline curvature and shock-induced flow separation. The data analyzed in this case is time-averaged and circumferentially-averaged LES data from Uzun and Malik [[3](#les-bump)]. The Reynolds numbers of the problem is *Re<sub>c</sub>=2763000*. The figure below illustrates the flow schematics.

<figure>
<p align="left">
    <img alt="bump_schematics" src="../docs/figs/bump_schematic.png" width="240" />
    <figcaption align="left"> Flow schematic of transonic bump </figcaption>
</p>
</figure>

With the method of *.AniRGB(c_off=0.65,c_exp=5)*, the Reynolds stress anisotropy componentality contour and the barycentric colormap are obtained as the following figures.

<figure>
<p align="left">
    <img alt="bump_AniCont" src="bump_data/results/anisotropy_contour.png" height="200" align="center"/>
    <img alt="bump_AniContMap" src="bump_data/results/anisotropy_colormap.png" height="120" align="center"/>
    <figcaption align="left"> (left): Reynolds stress anisotropy componentality contour; (right): barycentric colormap </figcaption>
</p>
</figure>

With the function of *.calc_EddyVisc()* using the QCR and *S<sub>ref</sub>=1000s<sup>-1</sup>*, the turbulent-to-laminar viscosity ratio and the limiter *f<sub>lim</sub>* contours are obtained as follows.

<figure>
<p align="left">
    <img alt="bump_NutrCont" src="bump_data/results/sref=1000_QCR2013V_ViscRatio_contour.png" height="200" align="center"/>
    <img alt="bump_NutrLimitery" src="bump_data/results/sref=1000_QCR2013V_ViscRatio_limiter.png" height="200" align="center"/>
    <figcaption align="left"> (left): viscosity ratio contour; (right): limiter contour </figcaption>
</p>
</figure>


---
### Case 3: 2D turbine film cooling
The 2D trailing edge cutback film cooling flow is feature by the mixing between the hot mainstream and the coolant in the trailing edge wake region. The data analyzed in this case is a time-averaged 2D slice of the DDES data from Wang and Yan [[4](#ddes-cool)], and the working condition is at a Reynolds number of *Re<sub>H</sub>=6200*, a Strouhal number of 0.20 and a blowing ratio of 0.5. The figure below illustrates the flow schematics.

<figure>
<p align="left">
    <img alt="cooling_schematics" src="../docs/figs/cooling_schematic.png" width="300" />
    <figcaption align="left"> Flow schematic of turbine film cooling </figcaption>
</p>
</figure>

With the method of *.AniRGB(c_off=0.65,c_exp=5)*, the Reynolds stress anisotropy componentality contour and the barycentric colormap are obtained as the following figures.

<figure>
<p align="left">
    <img alt="cooling_AniCont" src="cooling_data/results/anisotropy_contour.png" height="200" align="center"/>
    <img alt="cooling_AniContMap" src="cooling_data/results/anisotropy_colormap.png" height="120" align="center"/>
    <figcaption align="left"> (left): Reynolds stress anisotropy componentality contour; (right): barycentric colormap </figcaption>
</p>
</figure>

With the function of *.calc_EddyVisc()* using the QCR and *S<sub>ref</sub>=1000s<sup>-1</sup>*, the turbulent-to-laminar viscosity ratio and the limiter *f<sub>lim</sub>* contours are obtained as follows.

<figure>
<p align="left">
    <img alt="cooling_NutrCont" src="cooling_data/results/sref=1000_QCR2013V_ViscRatio_contour.png" height="200" align="center"/>
    <img alt="cooling_NutrLimitery" src="cooling_data/results/sref=1000_QCR2013V_ViscRatio_limiter.png" height="200" align="center"/>
    <figcaption align="left"> (left): viscosity ratio contour; (right): limiter contour </figcaption>
</p>
</figure>

---
### Case 4: 2D shock/boundary layer interaction
The 2D shock/boundary layer interaction flow is featured by shock-induced flow separation. The data analyzed in this case is time-averaged and spanwise-averaged DNS data from Pirozzoli and Bernardini [[5](#dns-sbli)]. The Reynolds numbers of the problem is *Re<sub>θ</sub>=2300*. The figure below illustrates the flow schematics.

<figure>
<p align="left">
    <img alt="SBLI_schematics" src="../docs/figs/SBLI_schematic.png" width="300" />
    <figcaption align="left"> Flow schematic of shock/boundary layer interaction </figcaption>
</p>
</figure>

With the method of *.AniRGB(c_off=0.65,c_exp=5)*, the Reynolds stress anisotropy componentality contour and the barycentric colormap are obtained as the following figures.

<figure>
<p align="left">
    <img alt="SBLI_AniCont" src="SBLI_data/results/anisotropy_contour.png" height="200" align="center"/>
    <img alt="SBLI_AniContMap" src="SBLI_data/results/anisotropy_colormap.png" height="120" align="center"/>
    <figcaption align="left"> (left): Reynolds stress anisotropy componentality contour; (right): barycentric colormap </figcaption>
</p>
</figure>

With the function of *.calc_EddyVisc()* using the QCR and *S<sub>ref</sub>=100s<sup>-1</sup>*, the turbulent-to-laminar viscosity ratio and the limiter *f<sub>lim</sub>* contours are obtained as follows.

<figure>
<p align="left">
    <img alt="SBLI_NutrCont" src="SBLI_data/results/sref=100_QCR2013V_ViscRatio_contour.png" height="200" align="center"/>
    <img alt="SBLI_NutrLimitery" src="SBLI_data/results/sref=100_QCR2013V_ViscRatio_limiter.png" height="200" align="center"/>
    <figcaption align="left"> (left): viscosity ratio contour; (right): limiter contour </figcaption>
</p>
</figure>


---
## References
[<a id="dns-bstep">1</a>] Le, H., Moin, P., & Kim, J. (1997). Direct numerical simulation of turbulent flow over a backward-facing step. Journal of Fluid Mechanics, 330, 349-374. [[DOI](https://doi.org/10.1017/S0022112096003941)]

[<a id="ddes-bstep">2</a>] He, X., Zhao, F., & Vahdati, M. (2022). Detached eddy simulation: recent development and application to compressor tip leakage flow. ASME Journal of Turbomachinery, 144(1), 011009. [[DOI](https://doi.org/10.1115/1.4052019)][[preprint](https://www.researchgate.net/publication/347355348_Detached_Eddy_Simulation_Recent_Development_and_Application_to_Compressor_Tip_Leakage_Flow)]

[<a id="les-bump">3</a>] Uzun, A., & Malik, M. (2019). Wall-resolved large-eddy simulations of transonic shock-induced flow separation. AIAA Journal, 57(5), 1955-1972. [[DOI](https://doi.org/10.2514/1.J057850.)]

[<a id="ddes-cool">4</a>] Wang, R., & Yan, X. (2021). Delayed-detached eddy simulations of film cooling effect on trailing edge cutback with land extensions. ASME Journal of Engineering for Gas Turbines and Power, 143(11), 111004. [[DOI](https://doi.org/10.1115/1.4051865)]

[<a id="dns-sbli">5</a>] Pirozzoli, S., & Bernardini, M. (2011). Direct numerical simulation database for impinging shock wave/turbulent boundary-layer interaction. AIAA Journal, 49(6), 1307-1312. [[DOI](https://doi.org/10.2514/1.J050901)]
