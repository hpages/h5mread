test_that("h5getdimscales() and h5setdimscales()", {
    h5getdimscales <- h5mread:::h5getdimscales
    h5setdimscales <- h5mread:::h5setdimscales

    h5file <- tempfile(fileext=".h5")
    h5createFile(h5file)
    h5createGroup(h5file, "stuff")
    h5createGroup(h5file, "more stuff")
    HDF5Array::writeHDF5Array(array(1:24, 4:1), h5file, "A")
    HDF5Array::writeHDF5Array(matrix(11:16, ncol=3), h5file, "B")
    HDF5Array::writeHDF5Array(matrix(letters, ncol=2), h5file, "stuff/C")
    HDF5Array::writeHDF5Array(matrix(runif(10), ncol=5), h5file, "more stuff/D")
    HDF5Array::writeHDF5Array(matrix(1:4, ncol=2), h5file, "more stuff/E")
    HDF5Array::writeHDF5Array(matrix(1:8, ncol=2), h5file, "E")

    expected <- rep(NA_character_, 4)
    current <- h5getdimscales(h5file, "A")
    expect_identical(current, expected)

    dimscales <- c("B", "stuff/C", NA, "more stuff/D")
    expected <- c(TRUE, TRUE, FALSE, TRUE)
    current <- h5setdimscales(h5file, "A", dimscales, dry.run=TRUE)
    expect_identical(current, expected)
    current <- h5setdimscales(h5file, "A", dimscales)
    expect_identical(current, expected)
    expect_identical(logical(4), h5setdimscales(h5file, "A", dimscales))

    expected <- paste0("/", dimscales); expected[is.na(dimscales)] <- NA
    current <- h5getdimscales(h5file, "A")
    expect_identical(current, expected)

    expected <- rep(NA_character_, 4)
    current <- h5getdimscales(h5file, "A", scalename="foo")
    expect_identical(current, expected)

    dimscales <- c(NA, "E", "more stuff/E", "E")
    expected <- c(FALSE, TRUE, TRUE, TRUE)
    current <- h5setdimscales(h5file, "A", dimscales, scalename="foo",
                              dry.run=TRUE)
    expect_identical(current, expected)
    current <- h5setdimscales(h5file, "A", dimscales, scalename="foo")
    expect_identical(current, expected)
    expect_identical(logical(4), h5setdimscales(h5file, "A", dimscales,
                                                scalename="foo"))

    expected <- paste0("/", dimscales); expected[is.na(dimscales)] <- NA
    current <- h5getdimscales(h5file, "A", scalename="foo")
    expect_identical(current, expected)

    for (name in paste0("Adim", 1:4))
        HDF5Array::writeHDF5Array(matrix(1:2), h5file, name)
    for (scalename in c(NA, "foo", "bar"))
        dimscales0 <- h5getdimscales(h5file, "A", scalename=scalename)
        for (scale1 in c(NA, "Adim1", "B"))
          for (scale2 in c(NA, "Adim2", "E"))
            for (scale3 in c(NA, "Adim3", "bogus"))
              for (scale4 in c(NA, "Adim4", "bogus")) {
                dimscales <- c(scale1, scale2, scale3, scale4)
                dimscales2 <- paste0("/", dimscales)
                ok <- is.na(dimscales) |
                      is.na(dimscales0) |
                      dimscales2 == dimscales0
                if (scale1 %in% "B" && !(scalename %in% NA))
                    ok[[1]] <- FALSE
                if (scale2 %in% "E" && !(scalename %in% "foo"))
                    ok[[2]] <- FALSE
                if (!all(ok) || "bogus" %in% dimscales) {
                    expect_error(h5setdimscales(h5file, "A", dimscales,
                                                scalename=scalename))
                    next
                }
                expected <- !is.na(dimscales)
                current <- h5setdimscales(h5file, "A", dimscales,
                                          scalename=scalename, dry.run=TRUE)
                expect_identical(current, expected)
              }
})

test_that("h5getdimlabels() and h5setdimlabels()", {
    h5getdimlabels <- h5mread:::h5getdimlabels
    h5setdimlabels <- h5mread:::h5setdimlabels

    h5file <- tempfile(fileext=".h5")
    h5createFile(h5file)
    h5createGroup(h5file, "stuff")
    h5createGroup(h5file, "more stuff")
    HDF5Array::writeHDF5Array(array(1:24, 4:1), h5file, "stuff/A")

    expect_identical(NULL, h5getdimlabels(h5file, "stuff/A"))

    dimlabels <- c("dim1", NA, NA, "dim4")
    h5setdimlabels(h5file, "stuff/A", dimlabels)
    expected <- c("dim1", "", "", "dim4")
    current <- h5getdimlabels(h5file, "stuff/A")
    expect_identical(current, expected)

    dimlabels <- c(NA, "dim2", NA, "XXXX")
    h5setdimlabels(h5file, "stuff/A", dimlabels)
    expected <- c("dim1", "dim2", "", "XXXX")
    current <- h5getdimlabels(h5file, "stuff/A")
    expect_identical(current, expected)

    dimlabels <- c("", "", NA, "")
    h5setdimlabels(h5file, "stuff/A", dimlabels)
    expect_identical(NULL, h5getdimlabels(h5file, "stuff/A"))
})

