# SNE_CP2_PDE_Matlab
Solution of task C3 of the SNE benchmark CP2. Definition: https://www.sne-journal.org/fileadmin/user_upload_sne/benchmarks/CP/CP2-definition.pdf

* Usage: run(location, type, n)
  * location = {"local", "cluster"} (for location: cluster editing the cluster profile name -> cluster = parcluster("Seneca")
  * type = {"loop", "matrix", "dist", "distSlow"}
  * n = number of workers
