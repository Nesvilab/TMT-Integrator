<div align="center">
<img src="https://raw.githubusercontent.com/Nesvilab/FragPipe/gh-pages/images/TMT-I_logo_transparent.png" width="200px"/>
</div>

TMT-Integrator extracts and combines channel abundances from multiple TMT or iTRAQ-labeled samples. It takes PSM tables generated by [Philosopher](https://philosopher.nesvilab.org/) as input, and generates quantification reports at the user-specified levels. TMT-Integrator currently provides four quantification options: gene, protein, peptide, and modified site levels. TMT-Integrator is included in [FragPipe](https://fragpipe.nesvilab.org/) with the Quant (Isobaric) options for complete analyses of isobaric labeling experiments.

There are six steps in TMT-Integrator, (1) PSM selection, (2) normalization and log transformation, (3) grouping, (4) outlier removal, (5) PSM aggregation to multiple levels, and (6) protein-level normalization.

<div align="center">
<img src="https://raw.githubusercontent.com/Nesvilab/FragPipe/gh-pages/images/TMT-I_overview_transparent.png" width="500px"/>
</div>

#### Best PSM selection

For improved quantification, we have the default option to include only PSMs that pass the following criteria:

* TMT/iTRAQ-labeled
* Reference channel intensity > 0
* Precursor ion purity >= 50%
* Summed (across all channels) MS2 intensity >= 5% (2.5% for phospho)
* PSM with the highest summed MS2 intensity among all PSMs identifying the same peptide ion in the sample/fraction (i.e. same LC-MS/MS run)
* Must not map (as unique or razor peptides) to contaminant proteins
* (Phospho: Must contain phosphorylation)


####	Normalization and log transformation
For each selected PSM, all channel intensities are divided by the reference channel, then log2 transformed, yielding normalized intensities.


####	PSM-level normalization (optional)
Within each multiplexed sample, selected PSMs are sorted and divided into groups by retention time. Within each retention time group, normalized intensities are subtracted from the median channel intensities.


####	Outlier removal (optional)
PSMs are grouped (by gene, protein, peptide, or modified site). Within these groups, the first quartile (Q1), third quartile (Q3) and interquartile range (IQR; Q3-Q1) are computed, and PSMs with intensities outside the boundaries of Q1-1.5xIQR and Q3+1.5xIQR are removed.


####	Protein-level normalization
There are two protein level normalization options: median centering or global normalization. With median centering, median normalized log2 intensity of the ith sample (m<sub>0</sub>) is subtracted from each channel so that the median ratio of each becomes 0. If global normalization is used, median centering is first performed, then the median deviation and median absolute deviation of centered values (m<sub>1</sub> and m<sub>2</sub>) is calculated and used to scale all values: y<sub>ij</sub>'=(y<sub>ij</sub>/m<sub>2</sub>)xm<sub>1</sub>+m<sub>0</sub>, where y<sub>ij</sub> is the log2 ratio in the ith group of the jth sample.                                                                           

## Download
TMT-Integrator is included in [FragPipe](https://fragpipe.nesvilab.org/). The latest standalone version can be downloaded [here](https://github.com/Nesvilab/TMT-Integrator/releases/latest).


## Use
For most users, we recommend running TMT-Integrator through [FragPipe](https://fragpipe.nesvilab.org/) graphical user interface. See the TMT tutorials below, or other tutorials listed on the [FragPipe homepage](https://fragpipe.nesvilab.org/).
- [Basic tutorial for TMT data](https://fragpipe.nesvilab.org/docs/tutorial_tmt.html)
- [Tutorial for analyzing TMT data with multiple plexes](https://fragpipe.nesvilab.org/docs/tutorial_tmt-2plexes.html)


For command line uses, see the Philosopher tutorials for running TMT analyses [step-by-step](https://github.com/Nesvilab/philosopher/wiki/Step-by-step-TMT-analysis) or with the [pipeline](https://github.com/Nesvilab/philosopher/wiki/Pipeline-mode-for-TMT-analysis) command. These examples can also be used to run the tool separately (note that both label-free/MS1 and isobaric quantification must be run prior to TMT-Integrator): 

`java -jar TMTIntegrator.jar TMTIntegrator.yaml PSM_Tables`

`java -jar -Xmx16g TMTIntegrator.jar TMTIntegrator.yaml /*_psm.tsv`
