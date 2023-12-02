# source("gaclust.R")

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


# test_that("integration test", {
#   
#   # Get test data "BeatAML_dataset_input.RDS" and "known_L2fc.RDS"
#   dataset <- readRDS(test_path("fixtures", "BeatAML_dataset_input.RDS"))
#   known_l2fc <- readRDS(test_path("fixtures", "known_L2fc.RDS"))
#   
#   # Get a static population (currently cannot do this)
#   # population <- readRDS(test_path("fixtures", "population.RDS"))
#   
#   # call ga.clust with test data (300 genes)
#   obj_gene <- ga.clust(dataset = dataset, k = 2,
#     plot.internals = FALSE, seed.ga = 42,
#     pop.size = nrow(dataset),
#     known_lfc = known_l2fc)
#   
#   # expect something
#   expect_equal(3 * 3, 9) #remove
# })


# test_that("ga.clust retuns correct solution", {
#   # Ensure that the solution returned by ga.clust() is the same as the one returned by GA::ga()
#   
#   
#   
#   expect_equal(2 * 2, 4) #remove
# })
