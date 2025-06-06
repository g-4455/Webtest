#Set the seed (very important)
set.seed(1)
load("../../data/ptmtable.rda")

#Create sample data
sink("noprint") #Supressess the print statements from spearman
testptmt <- SpearmanDissimilarity(ptmtable) #find some way to suppress print statements
sink()

#n is the number of digits that testptmt will round to. Note it should match the expected values on the right
n <- 1e-5

#Unit tests for correct output, ONLY TESTS FIRST ROW
test_that("Testing SpearmanDissimilarity Row 1, Col 1", {expect_equal(testptmt[1, 1], -2.59383, tolerance = n)})
test_that("Testing SpearmanDissimilarity Row 2, Col 2", {expect_equal(testptmt[1, 2], 25.35892, tolerance = n)})
test_that("Testing SpearmanDissimilarity Row 3, Col 3", {expect_equal(testptmt[1, 3], 6.842180, tolerance = n)})

if(file.exists("noprint")) file.remove("noprint") #Clean Up Sink file
