<h2>Tools</h2>
<ol>
<li>Add Title and label tags to the SDF file.</li>
<pre line="1" lang="python">
>add_label_to_sdf.py actives_final_dock.sdf actives_final_dock_label.sdf active
</pre>
<li>Group together of molecules with the same title stored in SDF, based on a aggregation functions such as max and min.</li>
<pre line="1" lang="python">
>extract_sdfscore2csv.py actives_final_dock_label.sdf  actives_score.csv Chemgauss4 min
</pre>
<li>Remove the brackets in the title tag</li>
<p>Some duplicated molecules has brackets in the title tag, remove it before grouping molecules with aggregation funtion.</p>
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
</ol>
