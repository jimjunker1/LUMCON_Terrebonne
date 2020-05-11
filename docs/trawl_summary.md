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

There are a number of reasons for data holes and more complete data will be necessary for any final analyses, but this shows the emerging environmental picture.

# Biological Data

The education department has been conducting trawl surveys pretty consistently since 2007.

## Sampling Events

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/sampling events-2.png)<!-- -->

There are a lot of sampling trawls over the past decade (n = 359) and they are spread over the year, with some periods of heavy sampling: 

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/month sampling events-1.png)<!-- -->

## Species sampling

Over the full dataset we can estimate the species accumulation by # of trawls from the `specaccum()` function in the <a href = "https://cran.r-project.org/web/packages/vegan/index.html" target = "_blank">vegan package</a>

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/full richness curve-1.png)<!-- -->

And by the trawls within each year to get at the differences in sampling effort (assuming each trawl is the same). This suggests there are certain years (or trawl locations):

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/SAC curves-1.png)<!-- -->

Squinting this suggests that, generally, there might be increasing # of species more broadly based on the species-effort curves. The 'brighter' colors appear to be higher. I think we still need to dig in further to see if this is due to changing taxonomic resolution across the sampling period, or an actual biological signal.

## Beta diversity

We can start to break this down further, and the inference that species richness might be increasing slightly over time has slight support (again, not sure yet if this is real or methodological). 

![](C:/Users/james/Documents/Projects/LUMCON_Terrebonne/text/../docs/trawl_summary_files/figure-html/beta partition-1.png)<!-- -->
 
The bottom panel can use a bit more explanation. This panel shows jaccard distance of each individual trawl from the initial trawl based on presence/absence data. Understanding the highlights and limitations of the analysis can give us a better picture of what is going on. First, a community in flux responding to directional change, we would expect the distance to increase over time. There is maybe slight evidence of this, but a lot of noise around that signal. Second, the relative effects of species turnover is quite variable over time and seems to show a larger seasonal signal relative to 
