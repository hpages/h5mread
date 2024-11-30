test_that("check_uaselection()", {
    check_uaselection <- h5mread:::check_uaselection

    ## specifying 'starts' only (no 'counts')

    expect_identical(check_uaselection(6:4), 6:4)
    expect_identical(check_uaselection(6:4, starts=list(NULL, 5, 3:4)),
                     c(6L, 1L, 2L))

    expect_identical(check_uaselection(integer(0)), integer(0))
    expect_identical(check_uaselection(integer(0), starts=NULL), integer(0))
    expect_identical(check_uaselection(integer(0), starts=list()), integer(0))

    expect_identical(check_uaselection(15), 15L)
    expect_identical(check_uaselection(15, starts=NULL), 15L)
    expect_identical(check_uaselection(15, starts=list(NULL)), 15L)

    expect_identical(check_uaselection(15, starts=list(c(6, 5, 2, 2, 10))), 5L)
    expect_identical(check_uaselection(15, starts=list(c(15:1, 1:15))), 30L)
    expect_identical(check_uaselection(3e9, starts=list(c(6, 5, 2, 2, 10))), 5L)
    expect_identical(check_uaselection(1e18, starts=list(1e18)), 1L)
    expect_identical(check_uaselection(3e9), 3e9)

    expect_error(check_uaselection(9, starts=list()))
    expect_error(check_uaselection(9, starts=list(0)))
    expect_error(check_uaselection(9, starts=list(NA)))
    expect_error(check_uaselection(9, starts=list(NA_integer_)))
    expect_error(check_uaselection(9, starts=list(NA_real_)))
    expect_error(check_uaselection(9, starts=list(NaN)))
    expect_error(check_uaselection(9, starts=list(Inf)))
    expect_error(check_uaselection(9, starts=list(-Inf)))
    expect_error(check_uaselection(9, starts=list(1e19)))
    expect_error(check_uaselection(9, starts=list(10)))

    ## specifying 'starts' and 'counts'

    expect_identical(check_uaselection(integer(0), list(), list()), integer(0))
    expect_identical(check_uaselection(15, starts=list(NULL),
                                           counts=list(NULL)), 15L)
    expect_identical(check_uaselection(15, starts=list(c(6, 5, 2, 2, 10)),
                                           counts=list(NULL)), 5L)
    expect_identical(check_uaselection(15, starts=list(c(14, 5, 8)),
                                           counts=list(c( 2, 0, 3))), 5L)

    expect_error(check_uaselection(-1, starts=list(NULL), counts=list()))
    expect_error(check_uaselection(-1, starts=list(), counts=list(NULL)))
    expect_error(check_uaselection(-1, starts=list(NULL), counts=list(3)))
    expect_error(check_uaselection(-1, starts=list(6), counts=list(-1)))
    expect_error(check_uaselection(15, starts=list(11), counts=list(6)))
    expect_error(check_uaselection(-1, starts=list(1), counts=list(3e9)))
    expect_error(check_uaselection(-1, starts=list(c(3, 5, 2)),
                                       counts=list(1e9, 1e9, 1e9)))
})

test_that("check_ordered_uaselection()", {
    check_ordered_uaselection <- h5mread:::check_ordered_uaselection

    # TODO!
})

test_that("reduce_uaselection()", {
    reduce_uaselection <- h5mread:::reduce_uaselection

    expect_identical(reduce_uaselection(c(99, 99)), NULL)  # no reduction

    ## specifying 'starts' only (no 'counts')

    starts <- list(NULL, c(2, 6))
    expect_identical(reduce_uaselection(c(99, 99), starts), NULL)  # no reduction
    starts <- list(NULL, integer(0))
    expect_identical(reduce_uaselection(c(99, 99), starts), NULL)  # no reduction
    starts <- list(NULL, c(2:5, 7:10))
    current <- reduce_uaselection(c(99, 99), starts)
    expect_identical(current[[1L]], list(NULL, c(2L, 7L)))
    expect_identical(current[[2L]], list(NULL, c(4L, 4L)))

    current <- reduce_uaselection(c(99, 99), starts=list(7:10, c(1:2, 5)))
    expect_identical(current[[1L]], list(7L, c(1L, 5L)))
    expect_identical(current[[2L]], list(4L, c(2L, 1L)))

    starts <- list(NULL, c(4, 6:7), c(2, 5), integer(0), numeric(0))
    current <- reduce_uaselection(rep(99, 5), starts)
    expect_identical(current[[1L]],
                     list(NULL, c(4L, 6L), c(2L, 5L), integer(0), integer(0)))
    expect_identical(current[[2L]],
                     list(NULL, c(1L, 2L), NULL, NULL, NULL))

    ## specifying 'starts' and 'counts'

    dim <- c(99, 99)
    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, NULL)
    expect_identical(reduce_uaselection(dim, starts, counts),  # no reduction
                     NULL)

    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, c(3, 4))
    expect_identical(reduce_uaselection(dim, starts, counts),  # no reduction
                     NULL)

    starts <- list(NULL, 5)
    counts <- list(NULL, 0)
    current <- reduce_uaselection(dim, starts, counts)
    expect_identical(current[[1L]], list(NULL, integer(0)))
    expect_identical(current[[2L]], list(NULL, integer(0)))

    starts <- list(NULL, c(6, 9))
    counts <- list(NULL, c(2, 0))
    current <- reduce_uaselection(dim, starts, counts)
    expect_identical(current[[1L]], list(NULL, 6L))
    expect_identical(current[[2L]], list(NULL, 2L))

    starts <- list(NULL, c(2, 5, 5, 6, 11, 11))
    counts <- list(NULL, c(3, 0, 1, 4,  0,  5))
    current <- reduce_uaselection(dim, starts, counts)
    expect_identical(current[[1L]], list(NULL, c(2L, 11L)))
    expect_identical(current[[2L]], list(NULL, c(8L,  5L)))

    starts <- list(NULL, c(2, 6, 10), c(5:10, 15, 3e9-10:1), c( 99, 2e9))
    counts <- list(NULL, c(3, 4,  3),                  NULL, c(6e8, 5e8))
    current <- reduce_uaselection(c(99, 99, 3e9, 3e9), starts, counts)
    expect_identical(current[[1L]],
                     list(NULL, c(2L, 6L), c(5, 15, 3e9-10), c(99L, 2e9L)))
    expect_identical(current[[2L]],
                     list(NULL, c(3L, 7L), c(6L, 1L, 10L), c(6e8L, 5e8L)))
})

test_that("test_map_starts_to_chunks()", {
    map_starts_to_chunks <- h5mread:::map_starts_to_chunks

    ## 1 dimension

    current <- map_starts_to_chunks(list(13:22), 85, 1)
    expected <- list(list(as.double(1:10)), list(as.double(12:21)))
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(13:22), 85, 10)
    expected <- list(list(c(8, 10)), list(c(1, 2)))
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(c(1:15, 49:51)), 85, 10)
    expected <- list(list(c(10, 15, 17, 18)), list(c(0, 1, 4, 5)))
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(2*(10:35)), 85, 10)
    expected <- list(list(1 + 5*(0:5)), list(as.double(1:6)))
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(c(6e9, 6e9 + 1)), 9e9, 3)
    expected <- list(list(c(1, 2)), list(c(1999999999, 2000000000)))
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(8e9), 9e9, 3)
    expected <- list(list(1), list(2666666666))
    expect_identical(current, expected)

    ## more dimensions

    current <- map_starts_to_chunks(list(NULL, 13:22, NULL),
                                    c(0, 85, 999), c(0, 10, 1))
    expected <- list(list(NULL, c(8, 10), NULL), list(NULL, c(1, 2), NULL))
    expect_identical(current, expected)

    ## edge cases

    current <- map_starts_to_chunks(list(), integer(0), integer(0))
    expected <- list(list(), list())
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(NULL), 0, 0)
    expected <- list(list(NULL), list(NULL))
    expect_identical(current, expected)

    expect_error(map_starts_to_chunks(list(NULL), 1, 0))

    current <- map_starts_to_chunks(list(NULL), 0, 1)
    expected <- list(list(NULL), list(NULL))
    expect_identical(current, expected)

    current <- map_starts_to_chunks(list(integer(0)), 0, 1)
    expected <- list(list(numeric(0)), list(numeric(0)))
    expect_identical(current, expected)
})

