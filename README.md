# localsancestrybreeding

A pilot model for modelling breeding using local ancestry, as described in Lawson, Howard-McCombe, Beaumont and Senn 2023 "To unscramble an egg: admixed captive breeding populations can be rescued using local ancestry information" (under submission).

All functions required for the simulation are in [localancestryfns.R](localancestryfns.R).

Figures 2-3 are produced by [sim.R](sim.R) which performs a direct simulation on the same captive population using a variety of different selection measures.

Figure 4 is produced by [multisim.R](multisim.R) which simulates a range of different captive populations and evaluates the performance of the best few measures.

Figure 5 is run on real [MOSAIC](https://maths.ucd.ie/~mst/MOSAIC/) output, here based on data from the Scottish Wildcat described in [Jamieson et al. 2023](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4384594). [convert_mosaic.R](convertmosaic/convert_mosaic.R) converts the output to the same format as the simulations, whilst [realdata.R](realdata.R) performs simulations very similarly to those undertaken for Figures 2-3.

Procedure:

```{sh}
Rscript sim.R # Produces the simulation plot
cd convertmosaic
Rscript convert_mosaic.R # Convert from Mosaic and VCF format to a simpler panel-based format
cd ../
Rscript realdata.R # Produces the wildcat breeding simlulation plot
```
