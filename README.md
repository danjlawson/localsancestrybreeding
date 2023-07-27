# localsancestrybreeding

A pilot model for modelling breeding using local ancestry, as described in Lawson, Howard-McCombe, Beaumont and Senn 2023 "To unscramble an egg: admixed captive breeding populations can be rescued using local ancestry information" (under submission).

All functions required for the simulation are in [localancestryfns.R](localancestryfns.R).

Procedure:

```{sh}
Rscript sim.R # Produces the simulation plot
cd convertmosaic
Rscript convert_mosaic.R # Convert from Mosaic and VCF format to a simpler panel-based format
cd ../
Rscript realdata.R # Produces the wildcat breeding simlulation plot
```

You get two outputs, SimulatedBreeding.pdf and SimulatedBreedingWildcat_L5000.pdf.
