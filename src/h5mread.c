/****************************************************************************
 *            Exploring alternate rhdf5::h5read() implementations           *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread.h"

#include "H5File.h"
#include "global_errmsg_buf.h"
#include "H5DSetDescriptor.h"
#include "uaselection.h"
#include "h5mread_startscounts.h"
#include "h5mread_index.h"
#include "h5mread_sparse.h"

#include "hdf5.h"


/****************************************************************************
 * C_get_h5mread_returned_type()
 *
 * The R type returned by h5mread() is determined by arguments 'filepath',
 * 'name', and 'as_integer'.
 */

/* --- .Call ENTRY POINT --- */
SEXP C_get_h5mread_returned_type(SEXP filepath, SEXP name, SEXP as_integer)
{
	int as_int, ret;
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;
	const H5TypeDescriptor *h5type;
	SEXPTYPE Rtype;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	ret = _init_H5DSetDescriptor(&h5dset, dset_id, as_int, 1);
	/* It's ok to close 'dset_id' **before** destroying its descriptor. */
	H5Dclose(dset_id);
	if (!isObject(filepath))
		H5Fclose(file_id);
	if (ret < 0)
		error("%s", _h5mread_global_errmsg_buf());

	h5type = h5dset.h5type;
	if (!h5type->Rtype_is_set) {
		_destroy_H5DSetDescriptor(&h5dset);
		PRINT_TO_ERRMSG_BUF(
			"h5mread() does not support this type "
			"of dataset yet, sorry. You can\n  "
			"use 'H5DSetDescriptor(filepath, name)' "
			"to see details about the dataset.");
		error("%s", _h5mread_global_errmsg_buf());
	}

	Rtype = h5type->Rtype;
	_destroy_H5DSetDescriptor(&h5dset);
	return ScalarString(type2str(Rtype));
}


/****************************************************************************
 * C_h5mread()
 */

/* Return -1 on error. */
static int select_method(const H5DSetDescriptor *h5dset,
			 SEXP starts, SEXP counts, int as_sparse, int method)
{
	int along;

	if (as_sparse) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "'as.sparse' is TRUE");
			return -1;
		}
		if (h5dset->h5chunkdim == NULL) {
			PRINT_TO_ERRMSG_BUF("'as.sparse=TRUE' is not supported "
					    "on a contiguous dataset");
			return -1;
		}
		if (method == 0) {
			method = 7;
		} else if (method != 7) {
			PRINT_TO_ERRMSG_BUF("only method 7 is supported "
					    "when 'as.sparse' is TRUE");
			return -1;
		}
		return method;
	}
	if (method < 0 || method > 6) {
		PRINT_TO_ERRMSG_BUF("'method' must be >= 0 and <= 6 "
				    "when 'as.sparse' is FALSE");
		return -1;
	}
	if (h5dset->h5type->Rtype == STRSXP) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "reading string data");
			return -1;
		}
		if (method == 0) {
			method = 4;
		} else if (method != 4 && method != 5) {
			PRINT_TO_ERRMSG_BUF("only methods 4 and 5 are "
					    "supported when reading "
					    "string data");
			return -1;
		}
		return method;
	}
	if (method == 0) {
		method = 1;
		/* March 27, 2019: My early testing (from Nov 2018) seemed
		   to indicate that method 6 was a better choice over method 4
		   when the layout is chunked and 'counts' is NULL. Turns out
		   that doing more testing today seems to indicate the opposite
		   i.e. method 4 now seems to perform better than method 6 on
		   all the datasets I've tested so far, including those used
		   by Pete Hickey here:
		     https://github.com/Bioconductor/DelayedArray/issues/13
		   and those used in the examples in man/h5mread.Rd.
		   Note sure what happened between Nov 2018 and today. Did I
		   do something stupid in my early testing? Did something
		   change in Rhdf5lib?
		   Anyway thanks to Pete for providing such a useful report.

		   Nov 26, 2019: I added method 5. Is like method 4 but
		   bypasses the intermediate buffer if a chunk is fully
		   selected. This is now preferred over methods 4 or 6. */
		if (h5dset->h5chunkdim != NULL &&
		    counts == R_NilValue &&
		    starts != R_NilValue)
		{
			for (along = 0; along < h5dset->ndim; along++) {
				if (VECTOR_ELT(starts, along) != R_NilValue) {
					//method = 6;
					//method = 4;
					method = 5;
					break;
				}
			}
		}
		return method;
	}
	if (method >= 4) {
		/* Make sure the data is chunked and 'counts' is NULL. */
		if (h5dset->h5chunkdim == NULL) {
			PRINT_TO_ERRMSG_BUF("methods 4, 5, and 6 cannot "
				"be used on a contiguous dataset (unless\n  "
				"it contains string data in which case "
				"methods 4 and 5 can be used)");
			return -1;
		}
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("methods 4, 5, and 6 can "
				"only be used when 'counts' is NULL");
			return -1;
		}
	}
	return method;
}

static int set_ans_dim(SEXP ans_dim, const size_t *ans_dim_buf,
		       int suggest_as_vec)
{
	int ndim = LENGTH(ans_dim);
	for (int along = 0; along < ndim; along++) {
		size_t d = ans_dim_buf[along];
		if (d <= INT_MAX) {
			INTEGER(ans_dim)[along] = (int) d;
			continue;
		}
		const char *s1 = "Too many elements (>= 2^31) "
				 "selected along dimension";
		const char *s2 = "of dataset. The\n  selection is "
				 "too large to fit in an R array.";
		if (suggest_as_vec) {
			PRINT_TO_ERRMSG_BUF(
				"%s %d %s Please reduce the size\n  "
				"of the selection or use 'as.vector=TRUE' "
				"to return it as an ordinary\n  vector.",
				s1, along + 1, s2);
		} else {
			PRINT_TO_ERRMSG_BUF("%s %d %s", s1, along + 1, s2);
		}
		return -1;
	}
	return 0;
}

/* If the H5 datatype that was used to store the logical data is an 8-bit
   or 16-bit signed integer type (e.g. H5T_STD_I8LE or H5T_STD_I16BE) then
   NA values got loaded as negative values that are not equal to NA_LOGICAL
   (e.g. as -128 for H5T_STD_I8LE and -2^16 for H5T_STD_I16BE).
   These values must be replaced with NA_LOGICAL. */
static void fix_logical_NAs(SEXP x)
{
	R_xlen_t x_len, i;
	int *x_p;

	x_len = XLENGTH(x);
	for (i = 0, x_p = LOGICAL(x); i < x_len; i++, x_p++) {
		if (*x_p < 0)
			*x_p = NA_LOGICAL;
	}
	return;
}

/* Replace "NA" strings in 'x' with character NAs (NA_character_).
   'x' is assumed to be NA-free. */
static void set_character_NAs(SEXP x)
{
	R_xlen_t x_len, i;
	int x_elt_len;
	SEXP x_elt;

	x_len = XLENGTH(x);
	for (i = 0; i < x_len; i++) {
		x_elt = STRING_ELT(x, i);
		x_elt_len = LENGTH(x_elt);
		if (x_elt_len == 2) {
			const char *s = CHAR(x_elt);
			if (s[0] == 'N' && s[1] == 'A')
				SET_STRING_ELT(x, i, NA_STRING);
		}
	}
	return;
}

/* Return R_NilValue on error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce,
		    int as_vec, int as_int, int as_sparse,
		    int method, int use_H5Dread_chunk)
{
	H5DSetDescriptor h5dset;
	if (_init_H5DSetDescriptor(&h5dset, dset_id, as_int, 0) < 0)
		return R_NilValue;

	SEXP ans = R_NilValue;
	int nprotected = 0;

	const H5TypeDescriptor *h5type = h5dset.h5type;
	if (!h5type->Rtype_is_set) {
		PRINT_TO_ERRMSG_BUF(
			"h5mread() does not support this type "
			"of dataset yet, sorry. You can\n  "
			"use 'H5DSetDescriptor(filepath, name)' "
			"to see details about the dataset.");
		goto done;
	}

	int ret = _shallow_check_uaselection(h5dset.ndim, starts, counts);
	if (ret < 0)
		goto done;

	if (as_vec == NA_INTEGER) {
		as_vec = h5dset.ndim == 1;
	} else if (as_vec && as_sparse) {
		warning("'as.vector' is ignored when 'as.sparse' is TRUE");
	}

	method = select_method(&h5dset, starts, counts, as_sparse, method);
	if (method < 0)
		goto done;
	if (use_H5Dread_chunk && method != 4 && method != 5) {
		PRINT_TO_ERRMSG_BUF("invalid use of 'use.H5Dread_chunk'");
		goto done;
	}

	size_t *ans_dim_buf = R_alloc0_size_t_array(h5dset.ndim);
	SEXP ans_dim;
	if (as_sparse || !as_vec) {
		ans_dim = PROTECT(NEW_INTEGER(h5dset.ndim));
		nprotected++;
	}

	if (method <= 3) {

		/* --- Methods 1, 2, 3 (as.sparse=FALSE) --- */

		SEXP startscounts = _compute_startscounts_ans_dim(&h5dset,
						starts, counts, noreduce,
						method, ans_dim_buf);
		if (startscounts == R_NilValue)
			goto done;
		if (!as_vec) {
			ret = set_ans_dim(ans_dim, ans_dim_buf, 1);
			if (ret < 0)
				goto done;
		}
		PROTECT(startscounts);
		ans = _h5mread_startscounts(&h5dset, startscounts, noreduce,
					    method, ans_dim_buf);
		UNPROTECT(1);

	} else if (method <= 6) {

		/* --- Methods 4, 5, 6 (counts=NULL, as.sparse=FALSE) --- */

		AllTChunks all_tchunks;
		ret = _init_AllTChunks(&all_tchunks, &h5dset, starts,
					     ans_dim_buf);
		if (ret < 0)
			goto done;
		if (!as_vec) {
			ret = set_ans_dim(ans_dim, ans_dim_buf, 1);
			if (ret < 0)
				goto done;
		}
		ans = _h5mread_index(&all_tchunks, method,
				     use_H5Dread_chunk, ans_dim_buf);

	} else {

		/* --- Method 7 (counts=NULL, as.sparse=TRUE) --- */

		AllTChunks all_tchunks;
		ret = _init_AllTChunks(&all_tchunks, &h5dset, starts,
				       ans_dim_buf);
		if (ret < 0)
			goto done;
		/* 'as_vec' ignored. */
		ret = set_ans_dim(ans_dim, ans_dim_buf, 0);
		if (ret < 0)
			goto done;
		/* Return 'list(NULL, nzcoo, nzdata)' or R_NilValue if
		   an error occured. */
		ans = _h5mread_sparse(&all_tchunks, ans_dim_buf);

	}

	if (ans != R_NilValue) {
		PROTECT(ans);
		nprotected++;
		if (as_sparse) {
			if (h5type->Rtype == LGLSXP) {
				fix_logical_NAs(VECTOR_ELT(ans, 2));
			} else if (h5type->Rtype == STRSXP &&
				   h5dset.as_na_attr)
			{
				set_character_NAs(VECTOR_ELT(ans, 2));
			}
			/* Final 'ans' is 'list(ans_dim, nzcoo, nzdata)'. */
			SET_VECTOR_ELT(ans, 0, ans_dim);
		} else {
			if (h5type->Rtype == LGLSXP)
				fix_logical_NAs(ans);
			if (!as_vec)
				SET_DIM(ans, ans_dim);
		}
	}

    done:
	_destroy_H5DSetDescriptor(&h5dset);
	UNPROTECT(nprotected);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce,
	       SEXP as_vector, SEXP as_integer, SEXP as_sparse,
	       SEXP method, SEXP use_H5Dread_chunk)
{
	/* Check 'noreduce'. */
	if (!(IS_LOGICAL(noreduce) && LENGTH(noreduce) == 1))
		error("'noreduce' must be TRUE or FALSE");
	int noreduce0 = LOGICAL(noreduce)[0];

	/* Check 'as_vector'. */
	if (!(IS_LOGICAL(as_vector) && LENGTH(as_vector) == 1))
		error("'as.vector' must be TRUE or FALSE");
	int as_vec = LOGICAL(as_vector)[0];

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as.integer' must be TRUE or FALSE");
	int as_int = LOGICAL(as_integer)[0];

	/* Check 'as_sparse'. */
	if (!(IS_LOGICAL(as_sparse) && LENGTH(as_sparse) == 1))
		error("'as.sparse' must be TRUE or FALSE");
	int as_sparse0 = LOGICAL(as_sparse)[0];

	/* Check 'method'. */
	if (!(IS_INTEGER(method) && LENGTH(method) == 1))
		error("'method' must be a single integer");
	int method0 = INTEGER(method)[0];

	/* Check 'use_H5Dread_chunk'. */
	if (!(IS_LOGICAL(use_H5Dread_chunk) && LENGTH(use_H5Dread_chunk) == 1))
		error("'use.H5Dread_chunk' must be TRUE or FALSE");
	int use_H5Dread_chunk0 = LOGICAL(use_H5Dread_chunk)[0];

	hid_t file_id = _get_file_id(filepath, 1);  /* read-only */
	hid_t dset_id = _get_dset_id(file_id, name, filepath);
	SEXP ans = PROTECT(h5mread(dset_id, starts, counts, noreduce0,
				   as_vec, as_int, as_sparse0,
				   method0, use_H5Dread_chunk0));
	H5Dclose(dset_id);
	if (!isObject(filepath))
		H5Fclose(file_id);
	UNPROTECT(1);
	if (ans == R_NilValue)
		error("%s", _h5mread_global_errmsg_buf());
	return ans;
}

