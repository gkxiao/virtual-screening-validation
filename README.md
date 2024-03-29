<h2>Performance</h2>
<ol>
<li>Chaput2016.csv<sup>[1]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: BEDROC (alpha=80.5)</p>
<p>Software: GOLD，Glide, Surflex and FlexX</p>  
<li>Cleves2020.csv<sup>[2]</sup></li>
<p>Data set: DUDE+</p>
<p>Metric: ROC AUC and ER 1%</p>
<p>Software: Dock, Glide and Surflex</p>
<li>Eberhardt2021.csv<sup>[3]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC, BEDROC (alpha=20), EF at 1%, 5% and 10%</p>
<p>Software: AutoDock 1.2</p>
<li>Mysinger2012.csv<sup>[4]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC, logAUC and EF at 1%</p>
<p>Software: DOCK</p>
<li>Wang2019.csv<sup>[6]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC and BEDROC (alpha=80.5)</p>
<p>Software: GLIDE</p>

<li>Jiang2020.xlsx<sup>[10]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC, BEDROC (alpha=20.0,80.5,321.0) and EF at 0.5%, 1%, 2%, 8%, 20%</p>
<p>Software: AutoPH4</p>

<li>Cleves2019.csv<sup>[12]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC</p>
<p>Software: Surflex eSim(-pscreen), maximum AUC over the alternate methods</p>  

<li>Koes2014.csv<sup>[13]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC and BEDROC</p>
<p>Software: USR, ROCS and VAMS</p>

<li>Puertas-Martín2019.csv<sup>[14]</sup></li>
<p>Data set: DUDE</p>
<p>Metric: ROC AUC</p>
<p>Software: OptiPharm and WEGA</p>

<li>Shen2020.xlsx<sup>[15]</sup></li>
<p>Data set: DUDE, DEKOIS2.0, dataset III</p>
<p>Metric: ROC AUC, logAUC, BEDROC(alpha=80.5), EF at 0.1%，0.5%, 1%, 5%</p>
<p>Software: GLIDE, GOLD, LeDock</p>

<li>Jiang2021.xlsx<sup>[16]</sup></li>
<p>Data set: DUD-E, LIT-PCBA</p>
<p>Metric: ROC AUC, EF at 1%, 5%, 10%</p>
<p>Software: ROCS、Phase Shape、SHAFTS、WEGA、ShaEP、Shape-it、Align-it、LIGSIFT、LS-align</p> 

<li>Jocelyn2021.xlsx<sup>[17]</sup></li>
<p>Data set: DUD-E, LIT-PCBA</p>
<p>Metric: ROC AUC, EF at 1%</p>
<p>Software: GNINA 1.0 with scoring function: Affinity,Pose,Affinity-dense,Pose-dense,Affinity-General,Pose-General,Vina,Vinardo,RFScore-VS,RFScore-4</p>
</p> 
</ol>

<h2>Tools</h2>
<ol>
  <li>metrics.py</li>
  <p>metrics: ROC AUC, BEDAUC, enrichment_factor(EF) and logAUC</p>
  <p>metrics.py can be available from <a href="https://github.com/oddt/oddt">oddt</a>.</p>

  <li>ROCKER<sup>[9]</sup></li>
     <p>ROCKER is a visualization tool for ROC and semi-log ROC curve</p>
     <p>ROCKER can be available from: http://www.medchem.fi/rocker</p>
  
  <li>bootstrap_tldr.py<sup>[11]</sup></li>
      <p>bootstrap_tldr.py is a visualization tool for ROC and semi-log ROC curve</p>
      <p>bootstrap_tldr.py can be available from: https://dudez.docking.org</p>
</ol>

<h2>Reference</h2>
<ol>
<li>Chaput, L.; Martinez-Sanz, J.; Saettel, N.; Mouawad, L. Benchmark of Four Popular Virtual Screening Programs: Construction of the Active/Decoy Dataset Remains a Major Determinant of Measured Performance. J. Cheminform. 2016, 8 (1), 56. https://doi.org/10.1186/s13321-016-0167-x.</li>
  
<li>Cleves, A. E.; Jain, A. N. Structure- and Ligand-Based Virtual Screening on DUD-E + : Performance Dependence on Approximations to the Binding Pocket. J. Chem. Inf. Model. 2020, 60 (9), 4296–4310. https://doi.org/10.1021/acs.jcim.0c00115.</li>
<p>Download: https://www.jainlab.org/downloads/</p>
  
<li>Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. J. Chem. Inf. Model. 2021, acs.jcim.1c00203. https://doi.org/10.1021/acs.jcim.1c00203.</li>
  
<li>Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. J. Med. Chem. 2012, 55 (14), 6582–6594. https://doi.org/10.1021/jm300687e.</li>
  
<li>Giangreco, I.; Mukhopadhyay, A.; C. Cole, J. Validation of a Field-Based Ligand Screener Using a Novel Benchmarking Data Set for Assessing 3D-Based Virtual Screening Methods. J. Chem. Inf. Model. 2021, 61 (12), 5841–5852. https://doi.org/10.1021/acs.jcim.1c00866.</li>
  
<li>Wang, D.; Cui, C.; Ding, X.; Xiong, Z.; Zheng, M.; Luo, X.; Jiang, H.; Chen, K. Improving the Virtual Screening Ability of Target-Specific Scoring Functions Using Deep Learning Methods. 2019, 10 (August), 1–11. https://doi.org/10.3389/fphar.2019.00924.</li>
  
<li>Imrie, F.; Bradley, A. R.; Van Der Schaar, M.; Deane, C. M. Protein Family-Specific Models Using Deep Neural Networks and Transfer Learning Improve Virtual Screening and Highlight the Need for More Data. J. Chem. Inf. Model. 2018, 58 (11), 2319–2330. https://doi.org/10.1021/acs.jcim.8b00350.</li>
  
<li>Scoring - Calculate rank statistics. http://www.rdkit.org/docs/source/rdkit.ML.Scoring.Scoring.html</li>
  
<li>Lätti, S.; Niinivehmas, S.; Pentikäinen, O. T. Rocker: Open Source, Easy-to-Use Tool for AUC and Enrichment Calculations and ROC Visualization. J. Cheminform. 2016, 8 (1), 45. https://doi.org/10.1186/s13321-016-0158-y.</li>

<li>Jiang, S.; Feher, M.; Williams, C.; Cole, B.; Shaw, D. E. AutoPH4: An Automated Method for Generating Pharmacophore Models from Protein Binding Pockets. J. Chem. Inf. Model. 2020, 60 (9), 4326–4338. https://doi.org/10.1021/acs.jcim.0c00121.</li>

<li>Stein, R. M.; Yang, Y.; Balius, T. E.; O’Meara, M. J.; Lyu, J.; Young, J.; Tang, K.; Shoichet, B. K.; Irwin, J. J. Property-Unmatched Decoys in Docking Benchmarks. J. Chem. Inf. Model. 2021, 61 (2), 699–714. https://doi.org/10.1021/acs.jcim.0c00598.</li>
 
<li>Cleves, A. E.; Johnson, S. R.; Jain, A. N. Electrostatic-Field and Surface-Shape Similarity for Virtual Screening and Pose Prediction. J. Comput. Aided. Mol. Des. 2019, 33 (10), 865–886. https://doi.org/10.1007/s10822-019-00236-6.</li>

<li>Koes, D. R.; Camacho, C. J. Shape-Based Virtual Screening with Volumetric Aligned Molecular Shapes. J. Comput. Chem. 2014, 35 (25), 1824–1834. https://doi.org/10.1002/jcc.23690.</li>

<li>Puertas-Martín, S.; Redondo, J. L.; Ortigosa, P. M.; Pérez-Sánchez, H. OptiPharm: An Evolutionary Algorithm to Compare Shape Similarity. Sci. Rep. 2019, 9 (1), 1–24. https://doi.org/10.1038/s41598-018-37908-6.</li>

<li>Shen, C.; Hu, Y.; Wang, Z.; Zhang, X.; Pang, J.; Wang, G.; Zhong, H.; Xu, L.; Cao, D.; Hou, T. Beware of the Generic Machine Learning-Based Scoring Functions in Structure-Based Virtual Screening. 2020, 00 (April), 1–22. https://doi.org/10.1093/bib/bbaa070.</li>

<li>Jiang, Z.; Xu, J.; Yan, A.; Wang, L. A Comprehensive Comparative Assessment of 3D Molecular Similarity Tools in Ligand-Based Virtual Screening. Brief. Bioinform. 2021, 22 (6), 1–17. https://doi.org/10.1093/bib/bbab231.</li>

<li>Sunseri, J.; Koes, D.R. Virtual Screening with Gnina 1.0. Molecules 2021, 26, 7369. https://doi.org/10.3390/molecules26237369</li>
</ol>
