

![banner](https://github.com/gfalbery/SpRanger/blob/master/front_logo.png)

# SpRanger

## Very in development, mostly useless at this point

Tools for deriving and manipulating IUCN species ranges and creating overlap matrices.

Get the ranges from here: https://www.iucnredlist.org/resources/spatial-data-download

Change them to rasters with this script:

`RangePlot` will plot a species' ranges on the earth.

Get proportional species overlaps with `PairsWisely`, which takes a rasterstack with one layer per species and returns a square overlap matrix where each row and column is a species and each cell is a proportion of raster grid cells that both species inhabit.

`GetIntersects` will take the rasterstack and iteratively produce a list of range overlap rasters, one range per species pair. I ###highly recommend### using `PairsWisely`first to get the range overlap matrix, then feeding that into the `Predicate` argument. It will save bags of time.
