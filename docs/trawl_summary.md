---
title: 'Trawl Summary'
author: "Jim Junker"
date: "11 May, 2020"
knit: (function(inputFile, encoding) { 
      out_dir <- '../docs';
      rmarkdown::render(inputFile,
                        encoding=encoding, 
                        output_file=file.path(dirname(inputFile), out_dir, 'trawl_summary.html')) })
output: 
  html_document:
      keep_md: true
      toc: true
      toc_float: true
      toc_depth: 3
      numbered_sections: true
---



#TERREBONNE BAY TRAWL ANALYSIS SUMMARY {.unlisted .unnumbered}

*insert terrebonne image here*

This document summarizes trawl catch data within Terrebonne Bay, Louisiana (*lat. long.*) conducted by the Education Department of the Louisiana Universities Marine Consortium (LUMCON) from 2007 - 2020.

# Metadata summary

## Sampling locations map
*insert sampling locations, i.e. trawls start and stop, here*

# Environmental Data

A suite of environmental data have been collected through <a href="https://lumcon.edu/environmental-monitoring/" target="_blank">LUMCON's environmental monitoring program</a> over the course of the study period. These variables include:

- Water Temperature
- Salinity
- Specific conductivity
- Dissolved Oxygen



These data vary in coverage and are being updated:

<img src="../figures/env_facet_plot.png" width="2000" />

More data will be necessary to fully complete analyses, but this shows the emerging environmental picture.

# Biological Data

The education department has been conducting trawl surveys pretty consistently since 2007. These sampling dates are pretty evenly distributed throughout the year

## Sampling Events

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/sampling events-2.png)<!-- -->

There are a lot of sampling trawls over the past decade (n = 359) and they are spread pretty consistently over the year, but with some periods of heavy sampling: 

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/month sampling events-1.png)<!-- -->

## Species sampling

Over the full dataset we can estimate the species accumulation by # of trawls

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/full richness curve-1.png)<!-- -->

And by the trawls within each year to get at the differences in sampling effort (assuming each trawl is the same). This suggests there are certain years (or trawl locations):

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

Squinting this suggests that, generally, there might be increasing # of species more broadly. 
