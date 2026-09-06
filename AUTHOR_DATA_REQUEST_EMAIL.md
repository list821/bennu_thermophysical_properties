# 给论文作者的数据与代码索取邮件

**收件人（To）：** danidg@lpl.arizona.edu; jemery2@utk.edu  
**主题（Subject）：** Request for processed data and model materials for reproducing the 2019 Bennu thermal analysis

---

Dear Dr. DellaGiustina and Dr. Emery,

My name is **[Your name]**, and I am **[your position/program]** at **[your institution]**. I am currently working on a research project concerning asteroid thermophysical modelling and am attempting to reproduce the analyses reported in:

DellaGiustina, D. N., Emery, J. P., et al. (2019), “Properties of rubble-pile asteroid (101955) Bennu from OSIRIS-REx imaging and thermal analysis,” *Nature Astronomy*, 3, 341-351. https://doi.org/10.1038/s41550-019-0731-1

I have downloaded the publicly available calibrated OVIRS and OTES products, observation-geometry products, SPICE kernels, and a public Bennu shape model from the NASA PDS/NAIF archives. I understand from the paper’s Data availability statement that the data supporting the plots and other findings may be obtained from the corresponding authors upon reasonable request.

If they can be shared, I would be very grateful for the following processed data and modelling materials. Native files are entirely acceptable; no reformatting or additional documentation is necessary.

### 1. Single-component thermophysical analysis

1. The exact OVIRS and OTES product IDs used in the analysis, together with any quality masks, excluded observations, or other selection criteria.
2. The processed OVIRS thermal-radiance data used for fitting, particularly the 3.5-4.0 μm data after reflected-sunlight subtraction and 1° rotational-phase binning, including the associated point and absolute-calibration uncertainties.
3. The processed OTES data used for fitting, particularly the approximately 14 μm radiance/flux series, rotational phases or time stamps, and associated uncertainties.
4. The exact observation-geometry table and rotational-phase convention used to associate the observations with the shape model.
5. The SPC v13 shape model used in the paper, preferably including facet vertices, facet connectivity, units, coordinate-frame definition, and any mapping or resampling used to obtain the approximately 12 m thermal facets.
6. The custom ATPM-based source code used for the analysis, together with configuration/input files and basic build or runtime instructions. If the source code cannot be distributed, a description of the numerical implementation and all non-default settings would still be extremely helpful.
7. The model grids and numerical results supporting Supplementary Figures 1 and 2, including model fluxes, χ² values as functions of thermal inertia and roughness, the fitted rotational-phase offsets, and the best-fit results for each observing date/instrument.
8. The 100 Monte Carlo/bootstrap realizations or their fitted-parameter table, together with the adopted individual-point and absolute-calibration error model.

### 2. Two-component thermal model

9. The boulder spatial-density or areal-fraction data used to constrain the two-component model.
10. The definition and implementation of the two-component mixing calculation, including the assumed fine-regolith thermal inertia, the tested boulder thermal-inertia grid, and the model results used to derive the 1,400 J m⁻² K⁻¹ s⁻½ upper limit for the average boulder thermal inertia.

### 3. Photometric analysis

11. The Image Photometric Data Information Files (IPDIFs), or equivalent tables containing image/product ID, filter, mean I/F, phase angle, incidence/emission geometry, and uncertainty.
12. The 1° phase-binned reflectance data underlying Figure 1, and the IDL/MPFIT photometric-modelling script or its configuration and complete fitted parameter tables for all five MapCam filters.

### 4. Boulder and particle-size analyses

13. The boulder catalogue used for Figure 3 and the supplementary size-frequency analyses, including locations, measured dimensions/areas, image IDs, completeness flags, and uncertainty information, if available.
14. The numerical cumulative size-frequency tables, fitted samples, statistical-completeness/KS-test results, and any scripts used to obtain the reported power-law slope.
15. The input parameter table and calculation code or spreadsheet used to convert thermal inertia to particle diameter and diurnal skin depth, including the adopted density, heat capacity, porosity, temperature, conductivity-model constants, and sampled parameter ranges.

If the complete set is not readily available, the highest-priority materials for my work are items 1-10: the processed OVIRS/OTES fitting arrays, their uncertainties, the v13 shape model, the custom ATPM implementation/configuration, the χ² grids, and the bootstrap outputs. A repository link, DOI, archive identifier, or any later/revised public data product would also be very helpful.

The requested materials would be used solely for academic reproduction and methodological comparison. I will cite the paper and any associated datasets or software as requested, and I would be glad to acknowledge the OSIRIS-REx team’s assistance. Please let me know if a data-use agreement or any additional information about my project is required.

Thank you very much for your time and consideration.

Best regards,

**[Your full name]**  
**[Degree/program or position]**  
**[Department and institution]**  
**[City, country]**  
**[Institutional email]**  
**[Optional ORCID / personal webpage]**

---

## 发送前填写提示

请替换所有方括号字段。建议使用学校邮箱发送，并在正文中保留“已经下载PDS公开数据”这句话，以说明你索取的是论文级处理中间数据而不是让作者代为下载公开产品。如果担心一次索取范围过大，可先只保留第1-10项；光度和巨石统计数据可在作者回复后第二次询问。
