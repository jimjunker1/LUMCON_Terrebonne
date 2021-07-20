#' @title diversity_clean_spp
#' @description This function identifies species names that are not found within
#'  the RESTORE saltmarsh species list. For missing species names, the function attempts
#'  to find correct names through four steps: fuzzy matching for minor transcription errors,
#'  fuzzy matching by common name, matching by known synonyms, and finally conversion from common name
#'  to scientific name.
#'  @param x data.frame or coercible to data.frame; ts the unique species names from data frame
#'  @param spp_col character; string of column name where species names are located. Default is 'species'
#'  @param common_col character; string of column name where common names are located. This is only
#'   necessary if species names are not found during initial fuzzy matching from species list. Default is 'common_name'
#'  @param collapse logical; if FALSE (the default), returns a tibble with output returned from all matching steps in
#'   individual columns. If TRUE, collapses the matched names into named vector and recodes spp_col with a warning.
#'  @param check logical; If FALSE (the default) the function will return a list or manipulated object depending on 'collapse'
#'  if TRUE, returns a logical for if all spp_col names are found in taxonomy_list
#'  @param taxonomy_list string or object of taxonomy list from which to check species names.
#'  @param ... other arguments to pass to sub functions.
#'  @returns list of matched names, manipulated dataframe (if collapse == TRUE), or logical (TRUE == all species found)

clean_trawl_spp <- function(x, spp_col = "species", common_col = "common_name", collapse = FALSE, 
                          check = FALSE, taxonomy_list = "./raw-data/taxonomy_valid_common.rds",...) {
  
  ## Load required functions for matching
  "%ni%" <- Negate("%in%")
  require(fuzzySim)
  # require(taxize)
  require(tibble)
  
  ## Parameter checks
  if (spp_col %ni% colnames(x)) stop("Error: spp_col is not found within the data.")
  
  # call in the taxonomy dataframe to reference good spp names
  if(class(taxonomy_list) %ni% c("character","tbl_df","tbl","data.frame", "matrix")){
    stop("Error: taxonomy_list must be a dataframe or coercible or character vector of file path.")
  }
  if(is.null(taxonomy_list)){
    stop("Error: taxonomy_list is empty. Must be dataframe or coercible or file path.")
  }
  if(is.character(taxonomy_list)){
    if(grepl(".rds", taxonomy_list, ignore.case = TRUE)){
      taxonomy_df = readRDS(file = taxonomy_list)
    } else if(grepl(".csv", taxonomy_list, ignore.case = TRUE)){
      taxonomy_df = read.csv(file = taxonomy_list, header = TRUE)
    }
  } else{
    taxonomy_df = as.data.frame(taxonomy_list)
  }
  
  ## Guts of function
  #### Level 1. direct matching of the scientific names to taxonomy list
  good_sci_names <- unlist(x[unlist(x[, spp_col]) %in% taxonomy_df$Species, spp_col]) %>%
    setNames(., nm = .)
  
  bad_sci_names <- unlist(x[unlist(x[, spp_col]) %ni% taxonomy_df$Species, spp_col]) %>% unique()
  
  ##### If just checking, return check results
  if(check){
    if(is.null(bad_sci_names)){
      return(cat("TRUE\n\nAll scientific names are present in taxonomy list. No actions done"))
    } else{
      return(cat("FALSE\n\nAll scientific names not present. set check == FALSE to attempt to fix names.\nSpecies not found:\n\n",
                 unique(bad_sci_names)))
    }
  }
  
  ##### If check == FALSE, continue function first doing a check for all names matched
  if (is.null(bad_sci_names) & !collapse) {
    return("All scientific names are present in taxonomy list. No actions done")
  } else if (is.null(bad_sci_names) & collapse) {
    message("All scientific names are present in taxonomy list. Returning data unchanged.")
    return(x)
  } else if (!is.null(bad_sci_names) & length(bad_sci_names) >= 1) {
    message(cat("Message: Some species names not found. ", "Species names inferred through fuzzy Matching:\n\n", unique(bad_sci_names)))
  }
  
  #### Level 2. If all names not present, do fuzzy matching of scientific names
  ##### fuzzy matching of the scientific names. Checking for small transcription errors. Return if all fixed.
  sci_name_fuzzy <- lapply(bad_sci_names, agrepl, x = taxonomy_df$Species, max.distance = 0.15) %>%
    lapply(., function(x) taxonomy_df$Species[x]) %>%
    setNames(., nm = bad_sci_names) %>%
    .[lengths(.) >= 1]
  ##### If names are genus names do not attempt to find.
  if(any(grepl(" sp.", names(sci_name_fuzzy)))){
    sci_name_fuzzy[grepl(" sp.", names(sci_name_fuzzy))] <- NULL
  }
  
  if (any(lengths(sci_name_fuzzy) == 0 | sapply(names(sci_name_fuzzy), is.na) | sapply(sci_name_fuzzy, is.null))) {
    sci_name_fuzzy[lengths(sci_name_fuzzy) == 0 | unlist(sapply(names(sci_name_fuzzy), is.na)) | unlist(sapply(sci_name_fuzzy, is.null)) ] <- NULL  }
  if(any(lengths(sci_name_fuzzy) >1)){
    sci_name_fuzzy[lengths(sci_name_fuzzy) > 1] = sci_name_fuzzy[lengths(sci_name_fuzzy) >1] %>%
      map(~.x %>% flatten %>% .[1])
  }
  
  sci_name_fuzzy = unlist(sci_name_fuzzy)
  
  ##### If all the names were found with fuzzy matching in agrepl with max.distance = 0.15m return 
  if (identical(length(sci_name_fuzzy), length(bad_sci_names))) {
    message(cat("Fuzzy matching found matches for all missing names."))
    
    out <- tibble(old_sci_names = unlist(x[, spp_col]), direct_match = NA_character_, sci_name_fuzzmatches = NA_character_) %>%
      dplyr::mutate(
        direct_match = case_when(
          old_sci_names %in% good_sci_names ~ recode(old_sci_names, !!!good_sci_names),
          TRUE ~ NA_character_
        ),
        sci_name_fuzzmatches = case_when(
          old_sci_names %in% names(sci_name_fuzzy) ~ recode(old_sci_names, !!!sci_name_fuzzy),
          TRUE ~ NA_character_
        )
      )
    return(out)
  }
  #### Level 3. Check for direct matches to common names and fill in species from there.
  ##### direct matching by common names
  if (common_col %ni% colnames(x)) stop("Error: common_col is not found within the data.")
  
  common_names_extracted <- x %>%
    dplyr::select(all_of(common_col)) %>%
    unlist() %>%
    trimws() %>%
    purrr::set_names() %>% 
    .[!duplicated(.)] %>%
    purrr::imap(~ purrr::map2(., taxonomy_df$ComName, ~ grepl(.x, .y) %>%
                                sum() %>%
                                as.logical()) %>%
                  unlist()) %>%
    purrr::map(~ taxonomy_df$Species[.x]) %>%
    unlist()
  
  if (!is.null(common_names_extracted)) {
    # identify duplicates and grab only the first when duplicated
    common_names_filtered <- common_names_extracted %>% .[grepl("[[:alpha:]]$", names(.), perl = TRUE)]
    common_dups_filtered <- common_names_extracted %>%
      .[grepl("\\w(?=1)", names(.), perl = TRUE)] %>%
      sapply(., function(x) gsub("\\s\\w.*", " sp.", x, ignore.case = TRUE), USE.NAMES = TRUE)
    new_names <- sapply(names(common_dups_filtered), function(x) gsub("1", "", x))
    
    common_dups_filtered <- setNames(common_dups_filtered, nm = new_names)
    
    # grab just the unique to change the key in future.
    common_names_extracted<- do.call(c, list(common_names_filtered, common_dups_filtered)) %>%
      .[!duplicated(.)]
  }
  
  
  if (is.null(common_names_extracted) & !collapse) {
    message("All common names were matched in taxonomy list.")
    ## Create a 
    out <- tibble(old_sci_names = unlist(x[, spp_col]), direct_match = NA_character_, sci_name_fuzzmatches = NA_character_,
                  old_common_names = unlist(x[,common_col]), common_match = NA_character_) %>%
      dplyr::mutate(
        direct_match = case_when(
          old_sci_names %in% good_sci_names ~ recode(old_sci_names, !!!good_sci_names),
          TRUE ~ NA_character_
        ),
        sci_name_fuzzmatches = case_when(
          old_sci_names %in% names(sci_name_fuzzy) ~ recode(old_sci_names, !!!sci_name_fuzzy),
          TRUE ~ NA_character_
        )
      )
    # return(out)
  } else if (is.null(common_names_extracted) & collapse) {
    message("All common names were matched in taxonomy list.")
    
    out <- tibble(old_sci_names = unlist(x[, spp_col]), direct_match = NA_character_, sci_name_fuzzmatches = NA_character_,
                  old_common_names = unlist(x[,common_col]), common_match = NA_character_) %>%
      dplyr::mutate(
        direct_match = case_when(
          old_sci_names %in% good_sci_names ~ recode(old_sci_names, !!!good_sci_names),
          TRUE ~ NA_character_
        ),
        sci_name_fuzzmatches = case_when(
          old_sci_names %in% names(sci_name_fuzzy) ~ recode(old_sci_names, !!!sci_name_fuzzy),
          TRUE ~ NA_character_
        ),
        common_match = case_when(
          old_sci_names %in% names(common_names_extracted) ~ recode(old_sci_names, !!!common_names_extracted),
          TRUE ~ NA_character_
        )
      )
  } else if (!is.null(common_names_extracted) & length(common_names_extracted) >= 1) {
    message(cat("Message: Some species names not found. ", "Species names inferred through fuzzy Matching of common_names."))
    out <- tibble(old_sci_names = unlist(x[, spp_col]), direct_match = NA_character_, sci_name_fuzzmatches = NA_character_,
                  old_common_names = unlist(x[,common_col]), common_match = NA_character_) %>%
      dplyr::mutate(
        direct_match = case_when(
          old_sci_names %in% good_sci_names ~ recode(old_sci_names, !!!good_sci_names),
          TRUE ~ NA_character_
        ),
        sci_name_fuzzmatches = case_when(
          old_sci_names %in% names(sci_name_fuzzy) ~ recode(old_sci_names, !!!sci_name_fuzzy),
          TRUE ~ NA_character_
        ),
        common_match = case_when(
          old_common_names %in% names(common_names_extracted) ~ recode(old_common_names, !!!common_names_extracted),
          TRUE ~ NA_character_
        )
      )
  }
  if(collapse){
    out = out %>%
      dplyr::mutate(names_out = case_when(!is.na(direct_match) ~ direct_match,
                                          is.na(direct_match) & !is.na(sci_name_fuzzmatches) ~ sci_name_fuzzmatches,
                                          is.na(direct_match) & is.na(sci_name_fuzzmatches) & !is.na(common_match) ~ common_match,
                                          TRUE ~ NA_character_))
    warning(paste0("Warning: original dataframe ",spp_col," replaced by direct & fuzzy matches,"))
    out_collapse = x %>% dplyr::mutate(across(all_of(spp_col), .x = out %>% dplyr::select(names_out))) 
    return(out_collapse)
  }
  return(out)
}
