test_that("h5mread() on 2D dataset", {
    do_2D_tests <- function(m, M, noreduce=FALSE, as.integer=FALSE, method=0L) {
        read <- function(starts=NULL, counts=NULL, as.vector=NA)
            h5mread(M@seed@filepath, M@seed@name,
                    starts=starts, counts=counts, noreduce=noreduce,
                    as.vector=as.vector, as.integer=as.integer, method=method)

        current <- read()
        expect_identical(current, m)
        current <- read(as.vector=TRUE)
        expect_identical(current, as.vector(m))

        current <- read(list(NULL, NULL))
        expect_identical(current, m)

        current <- read(list(integer(0), integer(0)))
        expected  <- m[NULL, NULL, drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(integer(0), integer(0)), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        ## With indices strictly sorted

        current <- read(list(c(2:5, 7:10), NULL))
        expected <- m[c(2:5, 7:10), , drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(c(2:5, 7:10), NULL), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        current <- read(list(NULL, 1:2))
        expected <- m[ , 1:2, drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(NULL, 1:2), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        current <- read(list(7:10, c(1:2, 5)))
        expected <- m[7:10, c(1:2, 5), drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(7:10, c(1:2, 5)), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        ## With indices in any order and with duplicates

        i <- c(2:6, 6:3, 1, 1, 1, 9:8)

        current <- read(list(i, NULL))
        expect_identical(current, m[i, , drop=FALSE])
        expect_error(read(list(i, NULL), as.vector=TRUE))

        current <- read(list(i, c(6:5, 5)))
        expect_identical(current, m[i, c(6:5, 5), drop=FALSE])
        expect_error(read(list(i, c(6:5, 5)), as.vector=TRUE))

        ## Only methods 1 and 3 support 'counts'.
        if (!(method %in% c(1L, 3L)))
            return()

        starts <- list(integer(0), 4L)
        counts <- list(integer(0), 2L)
        current <- read(starts, counts)
        expected <- m[NULL, 4:5, drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        starts <- list(5L, integer(0))
        counts <- list(4L, integer(0))
        current <- read(starts, counts)
        expected <- m[5:8, NULL, drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        starts <- list(c(2L, 5L), 4L)
        counts <- list(c(3L, 4L), 2L)
        current <- read(starts, counts)
        expected <- m[2:8, 4:5, drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))
    }

    do_2D_sparse_tests <- function(M, as.integer=FALSE) {
        test_with <- function(starts=NULL) {
            coo <- h5mread(M@seed@filepath, M@seed@name,
                           starts=starts,
                           as.integer=as.integer, as.sparse=TRUE)
            expect_true(is(coo, "COO_SparseArray"))
            current <- as.array(coo)
            expected <- h5mread(M@seed@filepath, M@seed@name,
                                starts=starts, as.integer=as.integer)
            expect_identical(current, expected)
        }
        test_with()
        test_with(list(NULL, NULL))
	test_with(list(integer(0), NULL))
	test_with(list(NULL, integer(0)))
	test_with(list(integer(0), integer(0)))
        test_with(list(c(2:5, 7:10), NULL))
        test_with(list(NULL, 1:2))
        test_with(list(7:10, c(1:2, 5)))
    }

    chunkdims <- list(0,         # no chunking (i.e. contiguous data)
                      c(10, 6),
                      c(10, 1),
                      c( 1, 6),
                      c( 4, 5),
                      c( 4, 3),
                      c( 2, 2))

    ## with an integer matrix

    m0 <- matrix(1:60, ncol=6)
    for (chunkdim in chunkdims) {
        M0 <- HDF5Array::writeHDF5Array(m0, chunkdim=chunkdim)
        do_2D_tests(m0, M0, method=1L)
        do_2D_tests(m0, M0, noreduce=TRUE, method=1L)
        do_2D_tests(m0, M0, method=2L)
        do_2D_tests(m0, M0, method=3L)
        do_2D_tests(m0, M0, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m0, M0, method=4L)
            do_2D_tests(m0, M0, method=5L)
            do_2D_tests(m0, M0, method=6L)
            do_2D_sparse_tests(M0)
        }
        do_2D_tests(m0, M0)
    }

    ## with a logical matrix

    m1 <- m0 %% 3L == 0L
    for (chunkdim in chunkdims) {
        M1 <- HDF5Array::writeHDF5Array(m1, chunkdim=chunkdim)
        do_2D_tests(m1, M1, method=1L)
        do_2D_tests(m1, M1, noreduce=TRUE, method=1L)
        do_2D_tests(m1, M1, method=2L)
        do_2D_tests(m1, M1, method=3L)
        do_2D_tests(m1, M1, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m1, M1, method=4L)
            do_2D_tests(m1, M1, method=5L)
            do_2D_tests(m1, M1, method=6L)
            do_2D_sparse_tests(M1)
        }
        do_2D_tests(m1, M1)
        storage.mode(m1) <- "integer"
        do_2D_tests(m1, M1, as.integer=TRUE)
        do_2D_tests(m1, M1, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m1, M1, as.integer=TRUE, method=4L)
            do_2D_tests(m1, M1, as.integer=TRUE, method=5L)
            do_2D_tests(m1, M1, as.integer=TRUE, method=6L)
            do_2D_sparse_tests(M1, as.integer=TRUE)
        }
    }

    ## with a numeric matrix

    m2 <- matrix(10 * runif(60), ncol=6)
    for (chunkdim in chunkdims) {
        M2 <- HDF5Array::writeHDF5Array(m2, chunkdim=chunkdim)
        do_2D_tests(m2, M2, method=1L)
        do_2D_tests(m2, M2, noreduce=TRUE, method=1L)
        do_2D_tests(m2, M2, method=2L)
        do_2D_tests(m2, M2, method=3L)
        do_2D_tests(m2, M2, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m2, M2, method=4L)
            do_2D_tests(m2, M2, method=5L)
            do_2D_tests(m2, M2, method=6L)
            do_2D_sparse_tests(M2)
        }
        do_2D_tests(m2, M2)
        storage.mode(m2) <- "integer"
        do_2D_tests(m2, M2, as.integer=TRUE)
        do_2D_tests(m2, M2, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m2, M2, as.integer=TRUE, method=4L)
            do_2D_tests(m2, M2, as.integer=TRUE, method=5L)
            do_2D_tests(m2, M2, as.integer=TRUE, method=6L)
            do_2D_sparse_tests(M2, as.integer=TRUE)
        }
    }

    ## with a character matrix

    m3 <- matrix(as.character(1:60), ncol=6)
    for (chunkdim in chunkdims) {
        M3 <- HDF5Array::writeHDF5Array(m3, chunkdim=chunkdim)
        do_2D_tests(m3, M3, method=4L)
        if (!identical(chunkdim, 0)) {
            do_2D_sparse_tests(M3)
        }
    }

    m3[cbind(5:10, 6:1)] <- NA_character_
    for (chunkdim in chunkdims) {
        M3 <- HDF5Array::writeHDF5Array(m3, chunkdim=chunkdim)
        do_2D_tests(m3, M3, method=4L)
        if (!identical(chunkdim, 0)) {
            do_2D_sparse_tests(M3)
        }
    }

    ## with a raw matrix

    m4 <- m0
    storage.mode(m4) <- "raw"
    for (chunkdim in chunkdims) {
        M4 <- HDF5Array::writeHDF5Array(m4, chunkdim=chunkdim)
        do_2D_tests(m4, M4, method=1L)
        do_2D_tests(m4, M4, noreduce=TRUE, method=1L)
        do_2D_tests(m4, M4, method=2L)
        do_2D_tests(m4, M4, method=3L)
        do_2D_tests(m4, M4, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m4, M4, method=4L)
            do_2D_tests(m4, M4, method=5L)
            do_2D_tests(m4, M4, method=6L)
            do_2D_sparse_tests(M4)
        }
        do_2D_tests(m4, M4)
        do_2D_tests(m0, M4, as.integer=TRUE)
        do_2D_tests(m0, M4, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m0, M4, as.integer=TRUE, method=4L)
            do_2D_tests(m0, M4, as.integer=TRUE, method=5L)
            do_2D_tests(m0, M4, as.integer=TRUE, method=6L)
            do_2D_sparse_tests(M4, as.integer=TRUE)
        }
    }
})

test_that("h5mread() on 3D dataset", {
    DIM <- c(10, 15, 6)

    do_3D_tests <- function(a, A, noreduce=FALSE, as.integer=FALSE, method=0L) {
        read <- function(starts=NULL, counts=NULL, as.vector=NA)
            h5mread(A@seed@filepath, A@seed@name,
                    starts=starts, counts=counts, noreduce=noreduce,
                    as.vector=as.vector, as.integer=as.integer, method=method)

        current <- read()
        expect_identical(current, a)
        current <- read(as.vector=TRUE)
        expect_identical(current, as.vector(a))

        current <- read(list(NULL, NULL, NULL))
        expect_identical(current, a)

        current <- read(list(integer(0), integer(0), integer(0)))
        expect_identical(current, a[NULL, NULL, NULL, drop=FALSE])

        ## With indices strictly sorted

        current <- read(list(c(2:5, 7:10), NULL, NULL))
        expected <- a[c(2:5, 7:10), , , drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(c(2:5, 7:10), NULL, NULL), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        current <- read(list(NULL, c(11:12, 14), NULL))
        expected <- a[ , c(11:12, 14), , drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(NULL, c(11:12, 14), NULL), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        current <- read(list(NULL, NULL, 1:2))
        expected <- a[ , , 1:2, drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(NULL, NULL, 1:2), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        current <- read(list(7:10, c(11:12, 14), c(1:2, 5)))
        expected <- a[7:10, c(11:12, 14), c(1:2, 5), drop=FALSE]
        expect_identical(current, expected)
        current <- read(list(7:10, c(11:12, 14), c(1:2, 5)), as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        ## With indices in any order and with duplicates

        i <- c(2:6, 6:3, 1, 1, 1, 9:8)

        current <- read(list(i, NULL, NULL))
        expect_identical(current, a[i, , , drop=FALSE])
        expect_error(read(list(i, NULL, NULL), as.vector=TRUE))

        current <- read(list(i, NULL, c(6:5, 5)))
        expect_identical(current, a[i, , c(6:5, 5), drop=FALSE])
        expect_error(read(list(i, NULL, c(6:5, 5)), as.vector=TRUE))

        ## Only methods 1 and 3 support 'counts'.
        if (!(method %in% c(1L, 3L)))
            return()

        starts <- list(integer(0), NULL, 4L)
        counts <- list(integer(0), NULL, 2L)
        current <- read(starts, counts)
        expected <- a[NULL, , 4:5, drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        starts <- list(5L, integer(0), NULL)
        counts <- list(4L, integer(0), NULL)
        current <- read(starts, counts)
        expected <- a[5:8, NULL, , drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        starts <- list(c(2L, 5L), NULL, 4L)
        counts <- list(c(3L, 4L), NULL, 2L)
        current <- read(starts, counts)
        expected <- a[2:8, , 4:5, drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))

        starts <- list(c(2L, 5L), 11L, 4L)
        counts <- list(c(3L, 4L),  5L, 2L)
        current <- read(starts, counts)
        expected <- a[2:8, 11:15, 4:5, drop=FALSE]
        expect_identical(current, expected)
        current <- read(starts, counts, as.vector=TRUE)
        expect_identical(current, as.vector(expected))
    }

    do_3D_sparse_tests <- function(A, as.integer=FALSE) {
        test_with <- function(starts=NULL) {
            coo <- h5mread(A@seed@filepath, A@seed@name,
                           starts=starts,
                           as.integer=as.integer, as.sparse=TRUE)
            expect_true(is(coo, "COO_SparseArray"))
            current <- as.array(coo)
            expected <- h5mread(A@seed@filepath, A@seed@name,
                              starts=starts, as.integer=as.integer)
            expect_identical(current, expected)
        }
        test_with()
        test_with(list(NULL, NULL, NULL))
	test_with(list(integer(0), NULL, NULL))
	test_with(list(NULL, integer(0), NULL))
	test_with(list(NULL, NULL, integer(0)))
	test_with(list(NULL, integer(0), integer(0)))
	test_with(list(integer(0), integer(0), integer(0)))
        test_with(list(c(2:5, 7:10), NULL, NULL))
        test_with(list(NULL, c(11:12, 14), NULL))
        test_with(list(NULL, NULL, 1:2))
        test_with(list(7:10, c(11:12, 14), c(1:2, 5)))
    }

    chunkdims <- list(0,             # no chunking (i.e. contiguous data)
                      c(10, 15, 6),
                      c(10, 15, 1),
                      c(10,  1, 6),
                      c( 1, 15, 6),
                      c(10,  1, 1),
                      c( 1, 15, 1),
                      c( 1,  1, 6),
                      c( 4,  5, 5),
                      c( 2,  2, 2))

    ## with an integer array

    a0 <- array(1:900, DIM)
    for (chunkdim in chunkdims) {
        A0 <- HDF5Array::writeHDF5Array(a0, chunkdim=chunkdim)
        do_3D_tests(a0, A0, method=1L)
        do_3D_tests(a0, A0, noreduce=TRUE, method=1L)
        do_3D_tests(a0, A0, method=2L)
        do_3D_tests(a0, A0, method=3L)
        do_3D_tests(a0, A0, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a0, A0, method=4L)
            do_3D_tests(a0, A0, method=5L)
            do_3D_tests(a0, A0, method=6L)
            do_3D_sparse_tests(A0)
        }
        do_3D_tests(a0, A0)
    }

    ## with a logical array

    a1 <- a0 %% 3L == 0L
    for (chunkdim in chunkdims) {
        A1 <- HDF5Array::writeHDF5Array(a1, chunkdim=chunkdim)
        do_3D_tests(a1, A1, method=1L)
        do_3D_tests(a1, A1, noreduce=TRUE, method=1L)
        do_3D_tests(a1, A1, method=2L)
        do_3D_tests(a1, A1, method=3L)
        do_3D_tests(a1, A1, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a1, A1, method=4L)
            do_3D_tests(a1, A1, method=5L)
            do_3D_tests(a1, A1, method=6L)
            do_3D_sparse_tests(A1)
        }
        do_3D_tests(a1, A1)
        storage.mode(a1) <- "integer"
        do_3D_tests(a1, A1, as.integer=TRUE)
        do_3D_tests(a1, A1, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a1, A1, as.integer=TRUE, method=4L)
            do_3D_tests(a1, A1, as.integer=TRUE, method=5L)
            do_3D_tests(a1, A1, as.integer=TRUE, method=6L)
            do_3D_sparse_tests(A1, as.integer=TRUE)
        }
    }

    ## with a numeric array

    a2 <- array(10 * runif(900), DIM)
    for (chunkdim in chunkdims) {
        A2 <- HDF5Array::writeHDF5Array(a2, chunkdim=chunkdim)
        do_3D_tests(a2, A2, method=1L)
        do_3D_tests(a2, A2, noreduce=TRUE, method=1L)
        do_3D_tests(a2, A2, method=2L)
        do_3D_tests(a2, A2, method=3L)
        do_3D_tests(a2, A2, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a2, A2, method=4L)
            do_3D_tests(a2, A2, method=5L)
            do_3D_tests(a2, A2, method=6L)
            do_3D_sparse_tests(A2)
        }
        do_3D_tests(a2, A2)
        storage.mode(a2) <- "integer"
        do_3D_tests(a2, A2, as.integer=TRUE)
        do_3D_tests(a2, A2, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a2, A2, as.integer=TRUE, method=4L)
            do_3D_tests(a2, A2, as.integer=TRUE, method=5L)
            do_3D_tests(a2, A2, as.integer=TRUE, method=6L)
            do_3D_sparse_tests(A2, as.integer=TRUE)
        }
    }

    ## with a character array

    a3 <- array(as.character(1:900), DIM)
    for (chunkdim in chunkdims) {
        A3 <- HDF5Array::writeHDF5Array(a3, chunkdim=chunkdim)
        do_3D_tests(a3, A3, method=4L)
        if (!identical(chunkdim, 0)) {
            do_3D_sparse_tests(A3)
        }
    }

    ## with a raw array

    a0 <- a0 %% 256L
    a4 <- a0
    storage.mode(a4) <- "raw"
    for (chunkdim in chunkdims) {
        A4 <- HDF5Array::writeHDF5Array(a4, chunkdim=chunkdim)
        do_3D_tests(a4, A4, method=1L)
        do_3D_tests(a4, A4, noreduce=TRUE, method=1L)
        do_3D_tests(a4, A4, method=2L)
        do_3D_tests(a4, A4, method=3L)
        do_3D_tests(a4, A4, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a4, A4, method=4L)
            do_3D_tests(a4, A4, method=5L)
            do_3D_tests(a4, A4, method=6L)
            do_3D_sparse_tests(A4)
        }
        do_3D_tests(a4, A4)
        do_3D_tests(a0, A4, as.integer=TRUE)
        do_3D_tests(a0, A4, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a0, A4, as.integer=TRUE, method=4L)
            do_3D_tests(a0, A4, as.integer=TRUE, method=5L)
            do_3D_tests(a0, A4, as.integer=TRUE, method=6L)
            do_3D_sparse_tests(A4)
        }
    }
})

test_that("h5mread() on 1D dataset", {
    a1 <- array(101:120)
    A1 <- HDF5Array::writeHDF5Array(a1, name="A1")

    expect_identical(h5mread(path(A1), "A1"), as.vector(a1))
    expect_identical(h5mread(path(A1), "A1", as.vector=FALSE), a1)

    starts <- list(c(8:1, 7:11))
    expected <- a1[c(8:1, 7:11)]
    current <- h5mread(path(A1), "A1", starts=starts)
    expect_identical(current, as.vector(expected))
    current <- h5mread(path(A1), "A1", starts=starts, as.vector=FALSE)
    expect_identical(current, expected)
})

