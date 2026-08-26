# `convertR`: A climate4R package for unit conversion and variable derivation

`convertR` is an R package performing unit conversion and variable derivation using the [**climate4R**](http://meteo.unican.es/climate4r) Common Data Model.

## Installation and requirements

Physical quantity handling is done via the Unidata's [UDUNITS-2](https://www.unidata.ucar.edu/software/udunits/) software libraries, tailoring the functionalities of the R package `udunits2` to the **climate4R** framework. Thus, the C-based UDUNITS-2 package **must be installed in your system** prior to package installation ([link to installation instructions](https://www.unidata.ucar.edu/software/udunits/udunits-current/doc/udunits/udunits2.html#Installation)).

Then, `install_github` can be used:

```R
devtools::install_github("SantanderMetGroup/convertR")
```

## Unit conversion

The function `udConvertGrid` is used for unit conversion. An example follows:

```R
library(convertR)
data("ps.iberia")
transformeR::getGridUnits(ps.iberia)
## [1] "Pascals"
range(ps.iberia$Data)
## [1]  86657.5 105690.0
ps.mmHg <- udConvertGrid(ps.iberia, new.units = "mmHg")
transformeR::getGridUnits(ps.mmHg)
## [1] "mmHg"
range(ps.mmHg$Data)
## [1] 649.9846 792.7401
```

## Derivation of physical quantities

In the following table the available variable derivations are summarized, including the input variables. The nomenclature of the input variables to calculate the derivations is according to the [climate4R vocabulary](https://github.com/SantanderMetGroup/loadeR/blob/devel/inst/vocabulary.txt) (see `loadeR::c4R.vocabulary()`)

The (*) symbol indicates that the function is planned but not yet available.


| function  	| Definition                                      	| psl 	| ps 	| tas 	| dpds 	| tdps 	| zgs 	| huss 	| hurs 	| rsds 	| rlds 	| ua(s) | va(s) |
|-----------	|-------------------------------------------------	|-----	|----	|-----	|------	|------	|----	|------	|------	|------	|------	| ------	|------	|
| `psl2ps`    	| Sea-level pressure to surface pressure          	| X   	|    	| X   	|      	|      	| X  	|      	|      	|      	|      	|      |      	|
| `ps2psl`    	| Surface pressure to sea-level pressure          	|     	| X  	| X   	|      	|      	| X  	|      	|      	|      	|      	|      |      	|
| `rad2cc`    	| Radiation to cloud cover                        	|     	|    	|     	|      	|      	|    	|      	|      	| X    	| X    	|      |      	|
| `huss2hurs` 	| Specific humidity to relative humidity          	|     	| X  	| X   	|      	|      	|    	| X    	|      	|      	|      	|      |      	|
| `hurs2huss` 	| Relative humidity from specific humidity        	|     	| X  	| X   	|      	|      	|    	|      	| X    	|      	|      	|      |      	|
| `tas2ws`    	| Saturation vapor pressure from temperature      	|     	| X  	| X   	|      	|      	|    	|      	|      	|      	|      	|      |      	|
| `hurs2w`    	| Water vapor mixing ratio from relative humidity 	|     	| X  	| X   	|      	|      	|    	|      	| X    	|      	|      	|      |      	|
| `dpds2hurs`* | Relative humidity from dew-point depression     	|     	| X  	| X   	| X    	|      	|    	|      	|      	|      	|      	|      |      	|
| `tdps2hurs` 	| Relative humidity from dew-point temperature    	|     	|    	| X   	|      	| X    	|    	|      	|      	|      	|      	|      |      	|
| `hurs2tdps` 	| Dew-point temperature from relative humidity    	|     	|    	| X   	|      	|     	 |    	|      	| X    	|      	|      	|      |      	|
| `huss2pvp`  	| Partial vapor pressure from specific humidity   	|     	| X  	|     	|      	|      	|    	| X    	|      	|      	|      	|      |      	|
| `uava2ws`  	 | Wind speed from wind components   	              |     	|    	|     	|      	|      	|    	|      	|      	|      	|      	|  X   |   X  	|


The input variables do not require to be in specific units, as all the necessary unit conversions are internally undertaken by the unit converter utility.


 <a href="https://unidata.ucar.edu/software/udunits" title="UDUnits"><img src="https://unidata.ucar.edu/images/logos/badges/badge_udunits_100.jpg" width="100px"></a>



