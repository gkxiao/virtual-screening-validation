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
</ol>
