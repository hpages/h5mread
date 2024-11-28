/****************************************************************************
 *                    Workhorse behind h5mread method 7                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_sparse.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "h5mread_helpers.h"
#include "ChunkIterator.h"

#include <string.h>  /* for memcpy */
#include <limits.h>  /* for INT_MAX */
//#include <time.h>


/* Note that we cap both the length of the 'nzdata' buffer and the number of
   rows in the 'nzcoo' buffer to INT_MAX. This is to prevent 'nzcoo' from
   growing into a matrix with more than INT_MAX rows, which R does not support
   yet. */
#define	NZDATA_MAXLENGTH INT_MAX


/****************************************************************************
 * Fast append an element to an auto-extending buffer
 */

static inline void IntAE_fast_append(IntAE *ae, int val)
{
	/* We don't use IntAE_get_nelt() for maximum speed. */
	if (ae->_nelt == ae->_buflength)
		IntAE_extend(ae, increase_buflength(ae->_buflength));
	ae->elts[ae->_nelt++] = val;
	return;
}

static inline void DoubleAE_fast_append(DoubleAE *ae, double val)
{
	/* We don't use DoubleAE_get_nelt() for maximum speed. */
	if (ae->_nelt == ae->_buflength)
		DoubleAE_extend(ae, increase_buflength(ae->_buflength));
	ae->elts[ae->_nelt++] = val;
	return;
}

static inline void CharAE_fast_append(CharAE *ae, char c)
{
	/* We don't use CharAE_get_nelt() for maximum speed. */
	if (ae->_nelt == ae->_buflength)
		CharAE_extend(ae, increase_buflength(ae->_buflength));
	ae->elts[ae->_nelt++] = c;
	return;
}

static inline void CharAEAE_fast_append(CharAEAE *aeae, CharAE *ae)
{
	/* We don't use CharAEAE_get_nelt() for maximum speed. */
	if (aeae->_nelt == aeae->_buflength)
		CharAEAE_extend(aeae, increase_buflength(aeae->_buflength));
	aeae->elts[aeae->_nelt++] = ae;
	return;
}


/****************************************************************************
 * Manipulation of the 'nzcoo' and 'nzdata' buffers
 */

static void *new_nzdata_buf(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return new_IntAE(0, 0, 0);
	    case REALSXP:             return new_DoubleAE(0, 0, 0.0);
	    case STRSXP:              return new_CharAEAE(0, 0);
	    case RAWSXP:              return new_CharAE(0);
	}
	/* Should never happen. The early call to _init_H5DSetDescriptor()
	   in h5mread() already made sure that Rtype is supported. */
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return NULL;
}

static SEXP make_nzcoo_from_bufs(const IntAEAE *nzcoo_bufs)
{
	int ndim = IntAEAE_get_nelt(nzcoo_bufs);
	size_t nzcoo_nrow = IntAE_get_nelt(nzcoo_bufs->elts[0]);
	/* 'nzcoo_nrow' is guaranteed to be <= INT_MAX (see
	   NZDATA_MAXLENGTH above) otherwise earlier calls to
	   copy_selected_chunk_data_to_nzbufs() (see below) would
	   have raised an error. */
	SEXP nzcoo = PROTECT(allocMatrix(INTSXP, (int) nzcoo_nrow, ndim));
	int *out_p = INTEGER(nzcoo);
	for (int along = 0; along < ndim; along++) {
		memcpy(out_p, nzcoo_bufs->elts[along]->elts,
		       sizeof(int) * nzcoo_nrow);
		out_p += nzcoo_nrow;
	}
	UNPROTECT(1);
	return nzcoo;
}

static SEXP make_nzdata_from_buf(const void *nzdata_buf, SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP:  return new_LOGICAL_from_IntAE(nzdata_buf);
	    case INTSXP:  return new_INTEGER_from_IntAE(nzdata_buf);
	    case REALSXP: return new_NUMERIC_from_DoubleAE(nzdata_buf);
	    case STRSXP:  return new_CHARACTER_from_CharAEAE(nzdata_buf);
	    case RAWSXP:  return new_RAW_from_CharAE(nzdata_buf);
	}
	/* Should never happen. The early call to _init_H5DSetDescriptor()
	   in h5mread() already made sure that Rtype is supported. */
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return R_NilValue;
}

static int copy_nzcoo_and_nzdata_to_ans(SEXPTYPE Rtype,
		const IntAEAE *nzcoo_bufs, const void *nzdata_buf, SEXP ans)
{
	SEXP ans_elt;

	/* Move the data in 'nzcoo_bufs' to an ordinary matrix. */
	ans_elt = PROTECT(make_nzcoo_from_bufs(nzcoo_bufs));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(1);
	if (ans_elt == R_NilValue)  /* should never happen */
		return -1;
	/* Move the data in 'nzdata_buf' to an atomic vector. */
	ans_elt = PROTECT(make_nzdata_from_buf(nzdata_buf, Rtype));
	SET_VECTOR_ELT(ans, 2, ans_elt);
	UNPROTECT(1);
	if (ans_elt == R_NilValue)  /* should never happen */
		return -1;
	return 0;
}


/****************************************************************************
 * copy_selected_chunk_data_to_nzbufs()
 */

static inline void append_array_index_to_nzcoo_bufs(
		const H5Viewport *mem_vp,
		const size_t *inner_midx_buf,
		IntAEAE *nzcoo_bufs)
{
	/* We don't use IntAEAE_get_nelt() for maximum speed. */
	int ndim = nzcoo_bufs->_nelt;
	for (int along = 0; along < ndim; along++) {
		IntAE *nzcoo_buf = nzcoo_bufs->elts[along];
		size_t midx = mem_vp->off[along] + inner_midx_buf[along] + 1;
		/* In the context of _h5mread_sparse(), 'midx' is guaranteed
		   to be <= INT_MAX so cast to int is safe. */
		IntAE_fast_append(nzcoo_buf, (int) midx);
	}
	return;
}

static inline int CharAEAE_append_if_nonzero(CharAEAE *aeae, const char *s,
					     size_t n)
{
	size_t s_len;
	for (s_len = 0; s_len < n; s_len++)
		if (s[s_len] == 0)
			break;
	if (s_len == 0)
		return 0;
	/* We don't use CharAEAE_get_nelt() for maximum speed. */
	if (aeae->_nelt >= NZDATA_MAXLENGTH)
		return -1;
	CharAE *ae = new_CharAE(s_len);
	memcpy(ae->elts, s, s_len);
	/* We don't use CharAE_set_nelt() for maximum speed. */
	ae->_nelt = s_len;
	CharAEAE_fast_append(aeae, ae);
	return 1;
}

static long long int copy_selected_string_chunk_data_to_CharAEAE_buf(
		const ChunkIterator *chunk_iter, size_t *inner_midx_buf,
		const char *in,
		IntAEAE *nzcoo_bufs, CharAEAE *nzdata_buf)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	size_t h5type_size = h5dset->h5type->h5type_size;
	size_t in_offset;
	_init_in_offset(ndim,
			chunk_iter->index,
			h5dset->h5chunkdim,
			&chunk_iter->mem_vp,
			&chunk_iter->h5dset_vp,
			&in_offset);
	while (1) {
		const char *s = in + in_offset * h5type_size;
		int ret = CharAEAE_append_if_nonzero(nzdata_buf, s,
						     h5type_size);
		if (ret < 0) {
			PRINT_TO_ERRMSG_BUF("too many non-zero "
					    "values to load");
			return -1;
		}
		if (ret > 0)
			append_array_index_to_nzcoo_bufs(
					&chunk_iter->mem_vp,
					inner_midx_buf,
					nzcoo_bufs);
		int inner_moved_along = next_midx(ndim,
						  chunk_iter->mem_vp.dim,
						  inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset(ndim,
				 chunk_iter->index,
				 h5dset->h5chunkdim,
				 &chunk_iter->mem_vp,
				 inner_midx_buf,
				 inner_moved_along,
				 &in_offset);
	};
	return (long long int) CharAEAE_get_nelt(nzdata_buf);
}

/* The "fast walk" method cannot be used on a truncated chunk! It works
   properly only if the chunk data spans the full chunk data buffer
   (this is 'chunk_data_buf->data' and it gets passed to the 'in' pointer),
   that is, if the current chunk is a full-size chunk and not a "truncated"
   chunk (a.k.a. "partial edge chunk" in HDF5's terminology). */
#define	ARGS_AND_BODY_OF_COPY_FUNCTION(in_type, nzdatabuf_type)(	 \
		const ChunkIterator *chunk_iter, size_t *inner_midx_buf, \
		const in_type *in,					 \
		IntAEAE *nzcoo_bufs, nzdatabuf_type *nzdata_buf)	 \
{									 \
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;		 \
	int ndim = h5dset->ndim;					 \
	int do_fast_walk = _tchunk_is_fully_selected(ndim,		 \
					&chunk_iter->h5dset_vp,		 \
					&chunk_iter->mem_vp)		 \
			   && ! _tchunk_is_truncated(h5dset,		 \
					&chunk_iter->h5dset_vp);	 \
	if (do_fast_walk) {						 \
		while (1) {						 \
			in_type val = *in;				 \
			if (val != (in_type) 0) {			 \
				if (nzdata_buf->_nelt >=		 \
				    NZDATA_MAXLENGTH)			 \
				{					 \
					PRINT_TO_ERRMSG_BUF(		 \
						"too many non-zero "	 \
						"values to load");	 \
					return -1;			 \
				}					 \
				nzdatabuf_type ## _fast_append(		 \
						nzdata_buf, val);	 \
				append_array_index_to_nzcoo_bufs(	 \
						&chunk_iter->mem_vp,	 \
						inner_midx_buf,		 \
						nzcoo_bufs);		 \
			}						 \
			int inner_moved_along = next_midx(ndim,		 \
						chunk_iter->mem_vp.dim,	 \
						inner_midx_buf);	 \
			if (inner_moved_along == ndim)			 \
				break;					 \
			in++;						 \
		};							 \
	} else {							 \
		size_t in_offset;					 \
		_init_in_offset(ndim,					 \
				chunk_iter->index,			 \
				h5dset->h5chunkdim,			 \
				&chunk_iter->mem_vp,			 \
				&chunk_iter->h5dset_vp,			 \
				&in_offset);				 \
		while (1) {						 \
			in_type val = in[in_offset];			 \
			if (val != (in_type) 0) {			 \
				if (nzdata_buf->_nelt >=		 \
				    NZDATA_MAXLENGTH)			 \
				{					 \
					PRINT_TO_ERRMSG_BUF(		 \
						"too many non-zero "	 \
						"values to load");	 \
					return -1;			 \
				}					 \
				nzdatabuf_type ## _fast_append(		 \
						nzdata_buf, val);	 \
				append_array_index_to_nzcoo_bufs(	 \
						&chunk_iter->mem_vp,	 \
						inner_midx_buf,		 \
						nzcoo_bufs);		 \
			}						 \
			int inner_moved_along = next_midx(ndim,		 \
						chunk_iter->mem_vp.dim,	 \
						inner_midx_buf);	 \
			if (inner_moved_along == ndim)			 \
				break;					 \
			update_in_offset(ndim,				 \
					 chunk_iter->index,		 \
					 h5dset->h5chunkdim,		 \
					 &chunk_iter->mem_vp,		 \
					 inner_midx_buf,		 \
					 inner_moved_along,		 \
					 &in_offset);			 \
		};							 \
	}								 \
	return (long long int) nzdatabuf_type ## _get_nelt(nzdata_buf);	 \
}

/* copy_selected_XXX_chunk_data_to_IntAE_buf() functions: copy ints and
   any smaller standard native type to an IntAE buf. */
static long long int copy_selected_int_chunk_data_to_IntAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(int, IntAE)
static long long int copy_selected_char_chunk_data_to_IntAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(char, IntAE)
static long long int copy_selected_schar_chunk_data_to_IntAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(signed char, IntAE)
static long long int copy_selected_uchar_chunk_data_to_IntAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, IntAE)
static long long int copy_selected_short_chunk_data_to_IntAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(short, IntAE)
static long long int copy_selected_ushort_chunk_data_to_IntAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned short, IntAE)

/* copy_selected_XXX_chunk_data_to_DoubleAE_buf() functions: copy doubles
   and any smaller standard native type to a DoubleAE buf. */
static long long int copy_selected_double_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(double, DoubleAE)
static long long int copy_selected_char_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(char, DoubleAE)
static long long int copy_selected_schar_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(signed char, DoubleAE)
static long long int copy_selected_uchar_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, DoubleAE)
static long long int copy_selected_short_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(short, DoubleAE)
static long long int copy_selected_ushort_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned short, DoubleAE)
static long long int copy_selected_int_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(int, DoubleAE)
static long long int copy_selected_uint_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned int, DoubleAE)
static long long int copy_selected_long_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(long, DoubleAE)
static long long int copy_selected_ulong_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned long, DoubleAE)
static long long int copy_selected_llong_chunk_data_to_DoubleAE_buf // be safe
	ARGS_AND_BODY_OF_COPY_FUNCTION(long long, DoubleAE)
static long long int copy_selected_ullong_chunk_data_to_DoubleAE_buf // be safe
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned long long, DoubleAE)
static long long int copy_selected_float_chunk_data_to_DoubleAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(float, DoubleAE)

/* Copy unsigned chars to a CharAE buf. */
static long long int copy_selected_uchar_chunk_data_to_CharAE_buf
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, CharAE)

static int copy_selected_chunk_data_to_nzbufs(
		const ChunkIterator *chunk_iter,
		ChunkDataBuffer *chunk_data_buf,
		size_t *inner_midx_buf,
		IntAEAE *nzcoo_bufs, void *nzdata_buf)
{

	//clock_t t0;
	//double dt;

	//t0 = clock();
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	if (h5dset->h5type->Rtype == STRSXP) {
		//printf("- copying selected chunk character data ... ");
		long long int nvals =
			copy_selected_string_chunk_data_to_CharAEAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
		if (nvals < 0)
			return -1;
		//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
		//printf("ok (%lld value%s copied in %3.3f ms)\n",
		//       nvals, nvals == 1 ? "" : "s", dt);
		return 0;
	}
	int copy_without_type_casting =
				chunk_data_buf->data_type_id ==
				h5dset->h5type->native_type_id_for_Rtype;
	//printf("- copying selected chunk data %s type casting ... ",
	//       copy_without_type_casting ? "WITHOUT" : "WITH");
	long long int nvals;
	switch (h5dset->h5type->Rtype) {
	    case INTSXP: case LGLSXP:
		if (copy_without_type_casting) {
			nvals = copy_selected_int_chunk_data_to_IntAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_CHAR) {
			nvals = copy_selected_char_chunk_data_to_IntAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SCHAR) {
			nvals = copy_selected_schar_chunk_data_to_IntAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UCHAR) {
			nvals = copy_selected_uchar_chunk_data_to_IntAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SHORT) {
			nvals = copy_selected_short_chunk_data_to_IntAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_USHORT) {
			nvals = copy_selected_ushort_chunk_data_to_IntAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	    case REALSXP:
		if (copy_without_type_casting) {
			nvals = copy_selected_double_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_CHAR) {
			nvals = copy_selected_char_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SCHAR) {
			nvals = copy_selected_schar_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UCHAR) {
			nvals = copy_selected_uchar_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SHORT) {
			nvals = copy_selected_short_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_USHORT) {
			nvals = copy_selected_ushort_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_INT) {
			nvals = copy_selected_int_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UINT) {
			nvals = copy_selected_uint_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_LONG) {
			nvals = copy_selected_long_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_ULONG) {
			nvals = copy_selected_ulong_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_LLONG) {
			nvals = copy_selected_llong_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_ULLONG) {
			nvals = copy_selected_ullong_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_FLOAT) {
			nvals = copy_selected_float_chunk_data_to_DoubleAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	    case RAWSXP:
		if (copy_without_type_casting) {
			nvals = copy_selected_uchar_chunk_data_to_CharAE_buf(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data,
					nzcoo_bufs, nzdata_buf);
			break;
		}
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	    default:
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	}
	if (nvals < 0)
		return -1;
	//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
	//printf("ok (%lld value%s copied in %3.3f ms)\n",
	//       nvals, nvals == 1 ? "" : "s", dt);
	return 0;
}


/****************************************************************************
 * read_data_7()
 *
 * One call to _read_H5Viewport() (wrapper for H5Dread()) per chunk touched
 * by the user-supplied array selection.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Make one call to _read_H5Viewport() to load the **entire** chunk data
 *     to an intermediate buffer.
 *   - Gather the non-zero user-selected data found in the chunk into
 *     'nzcoo_bufs' and 'nzdata_buf'.
 */

static int read_data_7(ChunkIterator *chunk_iter,
		IntAEAE *nzcoo_bufs, void *nzdata_buf)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	size_t *inner_midx_buf = R_alloc0_size_t_array(ndim);

	ChunkDataBuffer chunk_data_buf;
	int ret = _init_ChunkDataBuffer(&chunk_data_buf, h5dset, 0);
	if (ret < 0)
		return ret;
	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		//_print_tchunk_info(chunk_iter);

		//clock_t t0 = clock();
		ret = _load_chunk(chunk_iter, &chunk_data_buf, 0);
		if (ret < 0)
			break;
		//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
		//printf("- load chunk: %3.3f ms\n", dt);

		ret = copy_selected_chunk_data_to_nzbufs(
				chunk_iter,
				&chunk_data_buf,
				inner_midx_buf,
				nzcoo_bufs, nzdata_buf);
		if (ret < 0)
			break;
	}
	_destroy_ChunkDataBuffer(&chunk_data_buf);
	return ret;
}


/****************************************************************************
 * _h5mread_sparse()
 *
 * Implements method 7.
 * Return 'list(NULL, nzcoo, nzdata)' or R_NilValue if an error occured.
 */

SEXP _h5mread_sparse(ChunkIterator *chunk_iter, const size_t *ans_dim)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	IntAEAE *nzcoo_bufs = new_IntAEAE(ndim, ndim);
	void *nzdata_buf = new_nzdata_buf(h5dset->h5type->Rtype);
	if (nzdata_buf == NULL)  /* should never happen */
		return R_NilValue;

	int ret = read_data_7(chunk_iter, nzcoo_bufs, nzdata_buf);
	if (ret < 0)
		return R_NilValue;

	SEXP ans = PROTECT(NEW_LIST(3));
	//clock_t t0 = clock();
	ret = copy_nzcoo_and_nzdata_to_ans(h5dset->h5type->Rtype,
					   nzcoo_bufs, nzdata_buf,
					   ans);
	UNPROTECT(1);
	if (ret < 0)
		return R_NilValue;
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("copy_nzcoo_and_nzdata_to_ans(): %2.3f s\n", dt);
	return ans;
}

