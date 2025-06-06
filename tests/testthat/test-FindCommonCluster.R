#Setup
if(!exists("eu_ptms_list")){      #Check if global variables are already made as to not run MakeClusterList multiple times
  set.seed(1)                     #Set the seed (very important)
  load("../../data/ptmtable.rda") #Load Sample data

  sink("noprint")                 #Suppress print statements from function
  MakeClusterList(ptmtable)       #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)
  sink()
}

#Unit Testing - Test cluster sizes MAKE THIS A BETTER TEST (if any, this doesn't seem to do anything unique, I'd like to test to see if it writing plots or not)
test_that("Testing eu_ptms_list", {expect_equal(length(sapply(eu_ptms_list, function(x) dim(x)[1])), 6)})
test_that("Testing sed_ptms_list", {expect_equal(length(sapply(sed_ptms_list, function(x) dim(x)[1])), 6)})
test_that("Testing sp_ptms_list", {expect_equal(length(sapply(sp_ptms_list, function(x) dim(x)[1])), 88)})


#Clean Up
if(file.exists("noprint")) file.remove("noprint") #Clean up file created by sink()
