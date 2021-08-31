#!/usr/bin/env Rscript

library(jsonlite, quietly = T)
library(stringr, quietly = T)
library(optparse, quietly = T)

#-------------------------------------------------------------------------------
#' Utility to evaluate phenotyping rules (i.e. phenotype prediction) 
#' on a given genotyping file (e.g. genetic report card or GRC).
#-------------------------------------------------------------------------------

option_list <- list(make_option(c("--datafile"), 
                                help="Path to the data/genotypes file, a tabular file with the required columns (by default we expect the first column to be sample names).", 
                                type="character"),
                    make_option(c("--rulesfile"), 
                                help="Path to the JSON file that encodes the rules to be used for evaluating each drug phenotype.", 
                                type="character"), 
                    make_option(c("--ruleout"), 
                                help="The field to be extracted from a rule when it evaluates to true.", 
                                default="phenotype", 
                                type="character"),
                    make_option(c("--samplecolumn"), 
                                help="The name of the column in the data file that designates the sample id (optional). If not provided, we assume the first column provides sample names.", 
                                default=NULL, 
                                type="character"),
                    make_option(c("--verbose"), 
                                help="A flag to obtain verbose detailed output in the standard output (warning: can generate very big logs, useful for debugging).", 
                                default=F, 
                                action="store_true"),
                    make_option(c("--output"), 
                                help="Path for the output file. This file will have a row per sample and a column per drug/phenotype evaluated.", 
                                default="output.tab"))

program_description <- "Evaluates a set of phenotyping rules on all the samples of the given data/genotypes file."
arg <- parse_args(OptionParser(description=program_description, option_list=option_list))

#-------------------------------------------------------------------------------
#' Evaluates the rules associated with each drug for all sample instances.
#' @param metadata A data.frame with all the relevant columns/fields, by default 
#'  we expect sample names as first column.
#' @param drugs_rules A list of drugs, each one being itself a list that 
#'  contains a list of rules to be applied sequentially. 
#' @param rule_return_field The field to return from the rule list when a rule 
#'  is evaluated to true (e.g. phenotype).
#' @param verbose Activates or deactivates detailed output.
#' @param protection_mapping Replacement mapping to protect rules and column 
#'  names from conflicts with R code (better to protect column names with 
#'  backticks in the JSON rules file).
#' @return A data.frame of values as returned by the rules evaluated as true 
#'  for each instance. A row per sample and a column per each drug evaluated.
#-------------------------------------------------------------------------------

evaluateDrugsOnAllInstances <- function(metadata, drugs_rules, 
                                        rule_return_field="phenotype", 
                                        verbose=T, protection_mapping=list(),
                                        sample_field=NULL) {
  if(is.null(sample_field))
    sample_field <- 1
  checkRuleOutCoherence(drugs_rules, rule_return_field)
  pt <- ProgressTracker(nrow(metadata), progress_notification_threshold=0.05)
  results <- t(sapply(1:nrow(metadata), function(i) {
    instance <- metadata[i,]
    if(verbose) {
      out("\nProcessing SAMPLE: " %+% instance[sample_field] %+% "\n")
    }
    pt(i)
    return(evaluateDrugs(drugs_rules, instance, rule_return_field, 
                         verbose, protection_mapping))
  }, simplify=T))
  return(as.data.frame(results, stringsAsFactor=F, row.names=metadata[,1]))
}

#-------------------------------------------------------------------------------
#' Evaluates the rules associated with each drug for a sample instance 
#'  (i.e. a data row).
#' @param drugs_rules A list of drugs, each one being itself a list that 
#'  contains a list of rules to be applied sequentially. 
#' @param instance A named vector or 1-row data.frame with the instance 
#'  to evaluate (i.e. sample).
#' @param rule_return_field The field to return from the rule list when a rule 
#'  is evaluated to true (e.g. phenotype).
#' @param verbose Activates or deactivates detailed output.
#' @param protection_mapping Replacement mapping to protect rules and column 
#'  names from conflicts with R code (better to protect column names with 
#'  backticks in the JSON rules file).
#' @return A vector of values as returned by the rules evaluated as true.
#-------------------------------------------------------------------------------

evaluateDrugs <- function(drugs_rules, instance, 
                          rule_return_field="phenotype", 
                          verbose=T, protection_mapping) {
  
  return(sapply(names(drugs_rules), function(drug) {
    if(verbose) {
      out("Evaluating DRUG : " %+% drug %+% "...")
    }
    dr <- (drugs_rules[[drug]])$rules
    for(r in dr) {
      if(executeRule(r, instance, drugs_rules, rule_return_field, verbose, 
                     protection_mapping))
        return(r[[rule_return_field]])
    }
    stop("Error: all rules evaluated to false; there was no default-type rule.")
    }, simplify=T))
}

#-------------------------------------------------------------------------------
#' A preprocessing utility to avoid conflicts with R special characters.
#' It is better to protect column names with backticks in the JSON rules file). 
#' @param s The string(s) to be protected from conflict.
#' @param replace_mapping A list relating patters to replacement strings.
#' @return The protected string (i.e. with replaced patterns).
#-------------------------------------------------------------------------------

protectText <- function(s, replace_mapping) {
  for(p in names(replace_mapping))
    s <- str_replace_all(s,p,replace_mapping[[p]])
  return(s)
}

#-------------------------------------------------------------------------------
#' Protects metadata columns/fields to avoid conflicts with R special characters.
#' @param fields The set of columns/fields to be protected.
#' @param replace_mapping A list relating patters to replacement strings.
#' @return The protected string (i.e. with replaced patterns).
#-------------------------------------------------------------------------------

protectMetadataFields <- function(fields, replace_mapping) {
  return(sapply(fields, function(f) { protectText(f,replace_mapping) },
                simplify=T))
}

#-------------------------------------------------------------------------------
#' A preprocessing utility to allow rules to reference the output obtained
#' after evaluating another drug. These are refereced by @drug_name@
#' in the JSON file.
#' @param rule A list representing the rule. It must have (at least) the 
#'  fields name and evaluation, and include a specific return field
#'  (e.g. phenotype).
#' @param instance A named vector or 1-row data.frame with the sample/instance 
#'  to evaluate.
#' @param drugs_rules A list of drugs, each one being itself a list that 
#'  contains a list of rules to be applied sequentially.
#' @param rule_return_field The field to return from the rule list when a rule 
#'  is evaluated to true (e.g. phenotype).
#' @param verbose Activates or deactivates detailed output.
#' @param protection_mapping Replacement mapping to protect rules and column 
#'  names from conflicts with R code (better to protect column names with 
#'  backticks in the JSON rules file).
#' @return The original rule in which drug references have been replaced by 
#'  the result of evaluating those drugs for the given instance.
#-------------------------------------------------------------------------------

preprocessRuleReferences <- function(rule, instance, drugs_rules, 
                                     rule_return_field="phenotype", 
                                     verbose=T, protection_mapping) {
  # Extract references @drug@.
  references <- str_extract_all(rule,"@([^@]*)@")[[1]]
  
  # No references.
  if(length(references) == 0) {
    if(verbose)
      out("No references to other drugs found....")
    return(rule)
  }
  
  # References.
  else {
    
    if(verbose)
      out("Preprocessing rule with references: " %+% rule %+% ": (" %+% 
            paste(references, collapse=",") %+% ")")
    
    for(ref in references) {
      drug <- str_replace_all(ref,"@","")
      # Check we are referencing a valid drug.
      if(!(drug %in% names(drugs_rules)))
        stop("Unknown drug reference to: " %+% drug)
      # Execute rules for the drug.
      dr <- (drugs_rules[[drug]])$rules
      value <- NULL
      for(r in dr) {
        if(executeRule(r, instance, drugs_rules, rule_return_field, verbose, 
                       protection_mapping)) {
          value <- (r[[rule_return_field]])
          break
        }
      }
      if(is.null(value))
        stop("Error: all rules evaluated to false and no default-type rule.")
      # Replace reference in rule.
      if(is.character(value))
        value <- "'" %+% value %+% "'"
      rule <- str_replace(rule, ref, value)
    }
    
    if(verbose)
      out("Pre-processed rule: " %+% rule)
  }
  
  return(rule)
}

#-------------------------------------------------------------------------------
#' Evaluates a rule on the given instance (e.g. sample). To prevent injection 
#' attacks (e.g. access to the shell) we overwrite the system functions from 
#' the base package but bear in mind this could be an entry point for a 
#' skilled attacker. Only use JSON configuration files from sources you trust.
#' @param r A list representing the rule (having name and evaluation fields).
#' @param instance A named vector or 1-row data.frame with the instance/sample 
#'  to evaluate.
#' @param drugs_rules A list of drugs, each one being itself a list that 
#'  contains a list of rules to be applied sequentially.
#' @param rule_return_field The field to return from the rule list when a rule 
#'  is evaluated to true (e.g. phenotype).
#' @param verbose Activate, deactivate detailed output.
#' @param protection_mapping Replacement mapping to protect rules and column 
#'  names from conflicts with R code (better to protect column names with 
#'  backticks in the JSON rules file).
#' @return The truth value of the rule evaluated over the instance.
#-------------------------------------------------------------------------------

executeRule <- function(r, instance, drugs_rules, rule_return_field, 
                        verbose=T, protection_mapping) {
  
  # Print the whole rule.
  if(verbose) {
    out("Executing rule: " %+% r$name %+% " {" %+% r$evaluation %+% "}...")
  }
    
  # Process rule references to other drugs.
  r <- preprocessRuleReferences(r$evaluation, instance, drugs_rules, 
                                rule_return_field, verbose, protection_mapping)
  
  # This could be an entry point for injection attacks.
  protected_rule <- protectText(r, replace_mapping=protection_mapping)
  if(verbose)
    out("Protected rule: " %+% protected_rule)
  result <- eval(parse(text=protected_rule), envir=as.list(instance))
  
  if(verbose)
    out("Result: " %+% result)
  
  return(result)
}

#-------------------------------------------------------------------------------
#' Checks that the given rule return field is present in all rules. 
#'  An error is raised if the field is missing.
#' @param drug_rules The list of rules to be applied (by drug).
#' @param rule_return_field The field to return when a rule is true.
#-------------------------------------------------------------------------------

checkRuleOutCoherence <- function(drugs_rules, rule_return_field) {
  invisible(sapply(names(drugs_rules), function(drug) {
    dr <- (drugs_rules[[drug]])$rules
    sapply(dr, function(r){
      if(!rule_return_field %in% names(r)) {
        stop("Error: rule output field (--ruleout) not found for rule " %+% 
               r$name %+% " in drug: "  %+% drug %+% ".")
      }
    })
  }))
}

#-------------------------------------------------------------------------------
#' Utilities and syntactic sugar.
#' A set of functions and operators to ease things.
#-------------------------------------------------------------------------------

# Loose _contains_ operator that works with strings.
"%contains%" <- function(a, b) { any(str_detect(a, b)) }

# Symmetric _any_ operator (a in b or b in a).
"%==%" <- function(a, b) { any(a %in% b || b %in% a) }

# Concatenation.
"%+%" <- function(a, b) { paste(a,b,sep="") } 
tl <- tolower

# Formatting strings (like Python).
"%format%" <- function(fmt, list) {
  pat <- "%\\{([^}]*)\\}"
  fmt2 <- gsub(pat, "%s", fmt)
  list2 <- list[strapplyc(fmt, pat)[[1]]]
  do.call("sprintf", c(fmt2, list2))
}

# Formatting output for strings.
out <- function(text, use_time=F) {
  if(use_time)
    writeLines("[" %+% format(Sys.time(), "%H:%M") %+% "] " %+% text)
  else
    writeLines(text)  
}

# Formatting output for data/objects.
printd <- function(data, desc="",...) {
  if(desc != "")
    writeLines(paste0("[",format(Sys.time(), "%H:%M"),"] ",desc))
  print(data,quote=F,...)
  cat("\n")
}

# Conversion between time units (for printing remaining time).
convertToLargestTimeUnit <- function(duration, units) {
  unit_udpate_map <- c("secs"="mins", "mins"="hours")
  while(duration > 60 && units != "hours") {
    duration <- duration/60
    units <- unit_udpate_map[[units]]
  }
  if(units == "hours" && duration > 24) {
    duration <- duration/24
    units <- "days"
  }
  return(list(duration=duration, units=c(units)))
}

#-------------------------------------------------------------------------------
#' Builds a progress tracker.
#' @param iterations The number of iterations to be run or a vector with 
#'  the values to iterate through. NULL if we want the tracker to automatically 
#'  track the number of iterations (assuming a call per iteration).
#' @param progress_notification_threshold The threshold for printing progress 
#'  (0.01 corresponds to 1%).
#' @param progress_by_iteration A flag to show progress by iteration step.
#' @return A function that can be called with the current iteration index or 
#'  the current iteration value.
#-------------------------------------------------------------------------------

ProgressTracker <- function(iterations, 
                            progress_notification_threshold=0.1, 
                            progress_by_iteration=F) {
  
  .last_notification_progress <- 0
  .current_iteration <- 0
  .initial_tic <- Sys.time()
  
  tracker <- function(i=NULL) {
    if(is.null(i)) {
      .current_iteration <<- .current_iteration + 1
      current_progress <- .current_iteration/iterations
    }
    else if(length(iterations) == 1)
      current_progress <- i/iterations
    else
      current_progress <- which(i==iterations)/length(iterations)
    delta <- current_progress - .last_notification_progress
    epsilon <- -0.000000000001
    if(progress_by_iteration || delta - progress_notification_threshold >= epsilon) {
      current_toc <- Sys.time()
      out(as.integer(current_progress*100) %+% "% done...")
      elapsed_progress <- current_progress - .last_notification_progress
      .last_notification_progress <<- current_progress
      elapsed_time <- current_toc - .initial_tic
      et <- convertToLargestTimeUnit(elapsed_time, attr(elapsed_time, "units"))
      out("Elapsed time: " %+% round(et$duration,3) %+% " "  %+% et$units)
      rt <- convertToLargestTimeUnit(round(((1-current_progress)*elapsed_time) /current_progress, 2), 
                                     attr(elapsed_time, "units"))
      out("Estimated remaining time: " %+% round(rt$duration,3) %+% " " %+% rt$units)
    }
  }
  return(tracker)
}

#-------------------------------------------------------------------------------
# Protect from potential code injections by replacing system call functions.
# We evalute code from the rules to allow for discrepancies in the encoding of
# different data files (can use any operator/function in R).
#-------------------------------------------------------------------------------

psystem <- function(...) {stop("Preventing code injection from rules")}
assignInNamespace("system", psystem, "base")
assignInNamespace("system2", psystem, "base")
assignInNamespace("source", psystem, "base")
assignInNamespace("library", psystem, "base")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Execution.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Parsing the JSON rules.
out("Parsing rules file: " %+% arg$rulesfile)
rules <- fromJSON(arg$rulesfile, simplifyVector=F)

# Loading the data/metadata file (avoids quotes/comments).
out("Loading data file: " %+% arg$datafile)
metadata <- read.table(arg$datafile, header=T, sep="\t", comment.char="", 
                       check.names=F, stringsAsFactors=F, quote="")

# Execution of the logic.
out("Applying rules, this may take a while...")
results <- evaluateDrugsOnAllInstances(metadata, rules$drugs, 
                                       rule_return_field=arg$ruleout, 
                                       verbose=arg$verbose, 
                                       protection_mapping=list(),
                                       sample_field=arg$samplecolumn)

# Final results.
out("Writing output at " %+% arg$output)
if(!is.null(arg$samplecolumn)) {
  sample_ids <- metadata[,arg$samplecolumn]
  results <- cbind(sample=sample_ids, results)
}
write.table(results, file=arg$output, quote=F, 
            sep="\t", row.names=F, col.names=T)
