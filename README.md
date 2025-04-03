# ComingUpShort
 Code for my paper "Coming up short", coming soon to a journal near you

This repository contains the code to do all the analysis and all the figure generation.

Requirements:
[MATLAB R2023b](https://au.mathworks.com/products/new_products/release2023b.html)
[ImageMagick](https://imagemagick.org/index.php)(not required per se, but to make all the combined plots it is very much needed) 
[plotSurfaceROIBoundary](https://github.com/StuartJO/plotSurfaceROIBoundary) (only for making the plots of degrees spatial embedding)

## Rerunning analysis

### Generative models

To rerun the analysis for the generative models, this is how you would run the code
```
for form = 1:2
    for expo = 1:2
        for mdl = 1:9
            run_fitGNM(mdl,form,expo)
        end
    end
 end

 for form = 1:2
    for expo = 1:2
        for mdl = 1:12
            run_fitGNM(mdl,form,expo)
        end
    end
 end

for mdl = 1:10
    run_fitGNM_FLaG(mdl,2,2)
end

for mdl = 1:10
    run_fitGNM_WB(mdl,2,2)
end
```
Note that running it this way will take an exceptionally long time as it takes 2.5 hours minimum for one model to fit. And there are 104 models being run.

### Empirical network similarity

The following code wil compute the similarity between individual empirical networks constructed using 10 parcellations and two tractography algorithms (takes about 16 hours to run):

```
GetEmpComp
```

To then extract the needed data from the raw model outputs run:
```
GetBestResults.m
```


### Network rewiring

The code to do this is run as follows:
```
RunIterativeRewire
```
Will take ages, much like everything else

## Making figures

Figures can all* be remade by running the following script:
```
MakeFigures
```

*<sup><sub>It won't make Figures 1,3 and S6 but will make some of the elements you need to make those plots</sub></sup>

## Preprocessing the data

To get the original data for the biophysiological measures and grou[ consensus connectome] and configure it, please see this paper:

[Hansen JY, Shafiei G, Voigt K, Liang EX, Cox SML, Leyton M, et al. (2023) Integrating multimodal and multiscale connectivity blueprints of the human cerebral cortex in health and disease. PLoS Biol 21(9): e3002314. https://doi.org/10.1371/journal.pbio.3002314](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002314)

You can download the data from [this directory](https://github.com/netneurolab/hansen_many_networks/tree/v1.0.0/data/Schaefer400)

Put the data in the folder ./data/Schaefer_Hansen

Then run the following:
```
GetHansenData
```

The individual networks were processed as part of the following two papers:

[Arnatkeviciute, A., Fulcher, B. D., Oldham, S., Tiego, J., Paquola, C., Gerring, Z., ... & Fornito, A. (2021). Genetic influences on hub connectivity of the human connectome. Nature communications, 12(1), 4237.](https://www.nature.com/articles/s41467-021-24306-2)

[Oldham, S., Fulcher, B. D., Aquino, K., Arnatkevičiūtė, A., Paquola, C., Shishegar, R., & Fornito, A. (2022). Modeling spatial, developmental, physiological, and topological constraints on human brain connectivity. Science advances, 8(22), eabm6127.](https://www.science.org/doi/full/10.1126/sciadv.abm6127)

It is all provided now in the data download. 

To get the Euclidean distances between parcels of the Schaefer parcellations run
```
GetScha7Dist
```

--
TO DO

Some variables names are inconsistent across files

Sort code into appropriate subfolders

Add headers to function

Put all the raw processed data somewhere (many GB)