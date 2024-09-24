# Retrospective Motion Correction for Cardiac Multi-Parametric Mapping with Dictionary Matching-Based Image Synthesis and a Low-Rank Constraint


# Abstract

**Purpose**

To develop a model-based motion correction (MoCo) method that does not need an analytical signal model to improve the quality of cardiac multi-parametric mapping.

**Methods**

The proposed method constructs a hybrid loss that includes a dictionary-matching loss and a signal low-rankness loss, where the former registers the multi-contrast original images to a set of motion-free synthetic images and the latter forces the deformed images to be spatiotemporally coherent. We compared the proposed method with non-MoCo, a pairwise registration method (Pairwise-MI), and a groupwise registration method (pTVreg) via a free-breathing Multimapping dataset of 15 healthy subjects, both quantitatively and qualitatively.

**Results**

The proposed method achieved the lowest contour tracking errors (epicardium: 2.00±0.39mm vs 4.93±2.29mm, 3.50±1.26mm, and 2.61±1.00mm, and endocardium: 1.84±0.34mm vs 4.93±2.40mm, 3.43±1.27mm, and 2.55±1.09mm for the proposed method, non-MoCo, Pairwise-MI, and pTVreg, respectively; all p<0.01) and the lowest dictionary matching errors among all methods. The proposed method also achieved the highest scores on the visual quality of mapping (T1: 4.74±0.33 vs 2.91±0.82, 3.58±0.87, and 3.97±1.05, and T2: 4.48±0.56 vs 2.59±0.81, 3.56±0.93, and 4.14±0.80 for the proposed method, non-MoCo, Pairwise-MI, and pTVreg, respectively; all p<0.01). Finally, the proposed method had similar T1 and T2 mean values and standard deviations relative to the breath-hold reference in nearly all myocardial segments, whereas all other methods led to significantly different T1 and T2 measures and increases of standard deviations in multiple segments.

**Conclusions**

The proposed method significantly improves the motion correction accuracy and mapping quality compared with non-MoCo and alternative image-based methods.


# How to use

Environment: We developed and tested the method on MATLAB R2023a.

Demo: [moco_main.m](moco_main.m) shows how to apply the proposed MoCo method to cardiac Multimapping (described in our paper). By changing the signal dictionary, you may conveniently extend the method to your own qMRI applications.  


# How to cite
If you decide to use our code, please consider referring to the following paper:

Chen H, Emu Y, Gao J, Chen Z, Aburas A, Hu C. Retrospective motion correction for cardiac multi-parametric mapping with dictionary matching-based image synthesis and a low-rank constraint. Magnetic Resonance in Medicine. 2024. doi:10.1002/mrm.30291
