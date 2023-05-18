# localsancestrybreeding
Pilot model for breeding using local ancestry 

Procedure:

```{sh}
Rscript sim.R # Produces the simulation plot
cd convertmosaic
Rscript convert_mosaic.R # Convert from Mosaic and VCF format to a simpler panel-based format
cd ../
Rscript realdata.R # Produces the wildcat breeding simlulation plot
```

You get two outputs, SimulatedBreeding.pdf and SimulatedBreedingWildcat_L5000.pdf.
