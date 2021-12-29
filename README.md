<h2>DUDE+</h2>
<p>Source: Cleves, A. E.; Jain, A. N. Structure- and Ligand-Based Virtual Screening on DUD-E + : Performance Dependence on Approximations to the Binding Pocket. J. Chem. Inf. Model. 2020, 0 (0), acs.jcim.0c00115. https://doi.org/10.1021/acs.jcim.0c00115.</p>
<p>Download: https://www.jainlab.org/downloads/</p>
<h2>Tools</h2>
<ol>
<li>Remove the brackets in the title tag</li>
<p>The duplicated molecules has brackets in the SDF title tag. We need to remove it before grouping molecules with aggregation funtion.</p>
<pre line="1" lang="python">
CHEMBL162 (1)
CHEMBL471526 (1)
CHEMBL452341 (1)
CHEMBL383300 (1)
CHEMBL383300 (2)
CHEMBL383300 (3)
</pre>
<p>sed command is good to do it:</p>
<pre line="1" lang="python">
>sed s/\ \([0-9]*\)//g foo.sdf >> oof.sdf
</pre>
<li>Add Title and label tags to the SDF file.</li>
<pre line="1" lang="python">
>add_label_to_sdf.py actives_final_dock.sdf actives_final_dock_label.sdf active
</pre>
<li>Group together of molecules (differnet protonation states,tautomers) with the same title stored in SDF, based on a aggregation functions such as max and min.</li>
<pre line="1" lang="python">
>extract_sdfscore2csv.py actives_final_dock_label.sdf  actives_score.csv Chemgauss4 min
</pre>
</ol>
