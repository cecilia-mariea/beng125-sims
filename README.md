# Python Simulation Module for UCSD BENG125

Written for **Spring 2025** quarter session with Professor Jeff Hasty. 

## Network Description

The module was written to model the simplified dynamics of p53 autoregulation and ubiquitation by MDM2. The equation and parameters are explictly defined in `set_params_eqns.py` and changes to those statements should be reflected throughout the entire module. 

***NOTE***: The model was not fit to experimental data and is limited in many aspects ; some assumptions may not be biologically sound. The goal was to replicate the qualitative behaviors seen in literature to practice *dynamics* not biology. There is *no* physiological relevance or innovative features.

Although there is module functionality, each Python file has a main function for direct script calls that generate the figures that were included in our report. These may be useful as example cases for future use.

## Features By File

### 1D Functionality

- `solve_1D.py`
- `diffeq_1D.py`
- `bifurcations_1D.py`
- `param_sweep_1D.py`

### 2D Functionality

- `diffeq_2D.py`
- `phase_plane_2D.py`

No type annotations and my naming conventions get shaky throughout the module - sorry.

## Requirements

The scripts were developed in an environment with the follwing packages and their appropriate dependancies. For full intended functionality, have the following installed or equivalent versions. 

Most up-to-date, stable versions of these packages should work fine.

0. python 3.12.6
1. scipy 1.15.2
2. numpy 2.2.5
3. matplotlib 3.9.2
4. tqdm 6.4.1 - only used for parameter sweeping scripts that run for > 1 min.

## References

1. Nag S, Qin J, Srivenugopal KS, Wang M, Zhang R. The MDM2-p53 pathway revisited. J Biomed Res. 2013 Jul;27(4):254-71. doi: 10.7555/JBR.27.20130030. Epub 2013 Jun 6. PMID: 23885265; PMCID: PMC3721034.
2. Levine AJ. p53, the cellular gatekeeper for growth and division. Cell. 1997 Feb 7;88(3):323-31. doi: 10.1016/s0092-8674(00)81871-1. PMID: 9039259.
3. Wang H, Guo M, Wei H, Chen Y. Targeting p53 pathways: mechanisms, structures, and advances in therapy. Signal Transduct Target Ther. 2023 Mar 1;8(1):92. doi: 10.1038/s41392-023-01347-1. PMID: 36859359; PMCID: PMC9977964.
4. Lev Bar-Or R, Maya R, Segel LA, Alon U, Levine AJ, Oren M. Generation of oscillations by the p53-Mdm2 feedback loop: a theoretical and experimental study. Proc Natl Acad Sci U S A. 2000 Oct 10;97(21):11250-5. doi: 10.1073/pnas.210171597. PMID: 11016968; PMCID: PMC17186.
5. Ciliberto A, Novak B, Tyson JJ. Steady states and oscillations in the p53/Mdm2 network. Cell Cycle. 2005 Mar;4(3):488-93. doi: 10.4161/cc.4.3.1548. Epub 2005 Mar 18. PMID: 15725723.
6. Geva-Zatorsky N, Rosenfeld N, Itzkovitz S, Milo R, Sigal A, Dekel E, Yarnitzky T, Liron Y, Polak P, Lahav G, Alon U. Oscillations and variability in the p53 system. Mol Syst Biol. 2006;2:2006.0033. doi: 10.1038/msb4100068. Epub 2006 Jun 13. PMID: 16773083; PMCID: PMC1681500.
7. Eliaš J, Macnamara CK. Mathematical Modelling of p53 Signalling during DNA Damage Response: A Survey. Int J Mol Sci. 2021 Sep 30;22(19):10590. doi: 10.3390/ijms221910590. PMID: 34638930; PMCID: PMC8508851.
