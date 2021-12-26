<h2>Tools</h2>
<ol>
<li>add_label_to_sdf.py</li>
<p>Add Title and label tags to the SDF file.</p>
<pre line="1" lang="python">
>add_label_to_sdf.py actives_final_dock.sdf actives_final_dock_label.sdf active
</pre>
<li>extract_sdfscore2csv.py</li>
<p>Group together of a molecule with the same title stored in SDF, based on a aggregate functions such as max and min.</p>
<pre line="1" lang="python">
>extract_sdfscore2csv.py actives_final_dock_label.sdf  actives_score.csv Chemgauss4 min
</pre>
</ol>
