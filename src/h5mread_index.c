/****************************************************************************
 *                 Workhorses behind h5mread methods 4, 5, 6                *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_index.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "h5mread_helpers.h"
#include "ChunkIterator.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */
//#include <time.h>


/****************************************************************************
 * copy_selected_chunk_data_to_Rarray()
 */

static void init_in_offset_and_out_offset(int ndim, SEXP index,
			const size_t *out_dim, const H5Viewport *out_vp,
			const H5Viewport *h5dset_vp,
			const hsize_t *h5chunkdim,
			size_t *in_offset, size_t *out_offset)
{
	size_t in_off = 0, out_off = 0;
	int along, h5along;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunkdim[h5along];
		out_off *= out_dim[along];
		size_t i = out_vp->off[along];
		SEXP start = GET_LIST_ELT(index, along);
		if (start != R_NilValue)
			in_off += get_trusted_elt(start, (R_xlen_t) i) - 1 -
				  h5dset_vp->h5off[h5along];
		out_off += i;
	}
	*in_offset = in_off;
	*out_offset = out_off;
	//printf("# in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static inline void update_in_offset_and_out_offset(int ndim,
		SEXP index,
		const hsize_t *h5chunkdim,
		const H5Viewport *out_vp,
		const size_t *inner_midx, int inner_moved_along,
		const size_t *out_dim,
		size_t *in_offset, size_t *out_offset)
{
	SEXP start = GET_LIST_ELT(index, inner_moved_along);
	long long int in_off_inc;
	if (start != R_NilValue) {
		R_xlen_t i1 = out_vp->off[inner_moved_along] +
			      inner_midx[inner_moved_along];
		R_xlen_t i0 = i1 - 1;
		in_off_inc = get_trusted_elt(start, i1) -
			     get_trusted_elt(start, i0);
	} else {
		in_off_inc = 1;
	}
	long long int out_off_inc = 1;
	if (inner_moved_along >= 1) {
		int along = inner_moved_along - 1;
		int h5along = ndim - inner_moved_along;
		do {
			in_off_inc  *= h5chunkdim[h5along];
			out_off_inc *= out_dim[along];
			long long int di = 1 - out_vp->dim[along];
			SEXP start = GET_LIST_ELT(index, along);
			if (start != R_NilValue) {
				R_xlen_t i1 = out_vp->off[along];
				R_xlen_t i0 = i1 - di;
				in_off_inc += get_trusted_elt(start, i1) -
					      get_trusted_elt(start, i0);
			} else {
				in_off_inc += di;
			}
			out_off_inc += di;
			along--;
			h5along++;
		} while (along >= 0);
	}
	*in_offset += in_off_inc;
	*out_offset += out_off_inc;
	//printf("## in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static inline void copy_vlen_string_to_character_Rarray(
		const H5DSetDescriptor *h5dset,
		const char * const *in, size_t in_offset,
		SEXP Rarray, size_t Rarray_offset)
{
	/* Variable length strings are always null-terminated. */
	const char *s = in[in_offset];
	int is_na = h5dset->as_na_attr &&
		    s[0] == 'N' && s[1] == 'A' && s[2] == 0;
	if (is_na) {
		SET_STRING_ELT(Rarray, Rarray_offset, NA_STRING);
	} else {
		SEXP Rarray_elt = PROTECT(mkChar(s));
		SET_STRING_ELT(Rarray, Rarray_offset, Rarray_elt);
		UNPROTECT(1);
	}
	return;
}

static inline void copy_string_to_character_Rarray(
		const H5DSetDescriptor *h5dset,
		const char *in, size_t in_offset,
		SEXP Rarray, size_t Rarray_offset)
{
	size_t h5type_size = h5dset->h5type->h5type_size;
	const char *s = in + in_offset * h5type_size;
	size_t s_len;
	for (s_len = 0; s_len < h5type_size; s_len++)
		if (s[s_len] == 0)
			break;
	int is_na = h5dset->as_na_attr && s_len == 2 &&
		    s[0] == 'N' && s[1] == 'A';
	if (is_na) {
		SET_STRING_ELT(Rarray, Rarray_offset, NA_STRING);
	} else {
		SEXP Rarray_elt = PROTECT(mkCharLen(s, (int) s_len));
		SET_STRING_ELT(Rarray, Rarray_offset, Rarray_elt);
		UNPROTECT(1);
	}
	return;
}

static long long int copy_selected_string_chunk_data_to_character_Rarray(
		const ChunkIterator *chunk_iter, size_t *inner_midx_buf,
		const void *in, size_t in_offset,
		const size_t *Rarray_dim, SEXP Rarray, size_t Rarray_offset)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	long long int nvals = 0;
	while (1) {
		if (h5dset->h5type->is_variable_str)
			copy_vlen_string_to_character_Rarray(h5dset,
						in, in_offset,
						Rarray, Rarray_offset);
		else
			copy_string_to_character_Rarray(h5dset,
						in, in_offset,
						Rarray, Rarray_offset);
		nvals++;
		int inner_moved_along = next_midx(ndim,
						chunk_iter->mem_vp.dim,
						inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset_and_out_offset(ndim,
				chunk_iter->index,
				h5dset->h5chunkdim,
				&chunk_iter->mem_vp,
				inner_midx_buf,
				inner_moved_along,
				Rarray_dim,
				&in_offset, &Rarray_offset);
	};
	return nvals;
}

#define	ARGS_AND_BODY_OF_COPY_FUNCTION(in_type, out_type)(		  \
		const ChunkIterator *chunk_iter, size_t *inner_midx_buf,  \
		const in_type *in, size_t in_offset,			  \
		const size_t *out_dim, out_type *out, size_t out_offset)  \
{									  \
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;		  \
	int ndim = h5dset->ndim;					  \
	long long int nvals = 0;					  \
	while (1) {							  \
		out[out_offset] = in[in_offset];			  \
		nvals++;						  \
		int inner_moved_along = next_midx(ndim,			  \
						  chunk_iter->mem_vp.dim, \
						  inner_midx_buf);	  \
		if (inner_moved_along == ndim)				  \
			break;						  \
		update_in_offset_and_out_offset(ndim,			  \
				chunk_iter->index,			  \
				h5dset->h5chunkdim,			  \
				&chunk_iter->mem_vp,			  \
				inner_midx_buf,				  \
				inner_moved_along,			  \
				out_dim,				  \
				&in_offset, &out_offset);		  \
	};								  \
	return nvals;							  \
}

/* copy_selected_XXX_chunk_data_to_int_array() functions: copy ints and
   any smaller standard native type to an array of ints. */
static long long int copy_selected_int_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(int, int)
static long long int copy_selected_char_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(char, int)
static long long int copy_selected_schar_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(signed char, int)
static long long int copy_selected_uchar_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, int)
static long long int copy_selected_short_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(short, int)
static long long int copy_selected_ushort_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned short, int)

/* copy_selected_XXX_chunk_data_to_double_array() functions: copy doubles
   and any smaller standard native type to an array of doubles. */
static long long int copy_selected_double_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(double, double)
static long long int copy_selected_char_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(char, double)
static long long int copy_selected_schar_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(signed char, double)
static long long int copy_selected_uchar_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, double)
static long long int copy_selected_short_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(short, double)
static long long int copy_selected_ushort_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned short, double)
static long long int copy_selected_int_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(int, double)
static long long int copy_selected_uint_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned int, double)
static long long int copy_selected_long_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(long, double)
static long long int copy_selected_ulong_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned long, double)
static long long int copy_selected_llong_chunk_data_to_double_array // be safe
	ARGS_AND_BODY_OF_COPY_FUNCTION(long long, double)
static long long int copy_selected_ullong_chunk_data_to_double_array // be safe
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned long long, double)
static long long int copy_selected_float_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(float, double)

/* Copy unsigned chars (a.k.a Rbytes) to an array of unsigned chars. */
static long long int copy_selected_uchar_chunk_data_to_uchar_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, unsigned char)

static long long int copy_selected_chunk_data_to_Rarray(
		const ChunkIterator *chunk_iter,
		ChunkDataBuffer *chunk_data_buf,
		size_t *inner_midx_buf,
		const size_t *Rarray_dim, SEXP Rarray)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	size_t in_offset, out_offset;
	init_in_offset_and_out_offset(h5dset->ndim, chunk_iter->index,
			Rarray_dim, &chunk_iter->mem_vp,
			&chunk_iter->h5dset_vp, h5dset->h5chunkdim,
			&in_offset, &out_offset);
	SEXPTYPE Rtype = h5dset->h5type->Rtype;
	if (Rtype == STRSXP)
		return copy_selected_string_chunk_data_to_character_Rarray(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, Rarray, out_offset);
	int copy_without_type_casting =
				chunk_data_buf->data_type_id ==
				h5dset->h5type->native_type_id_for_Rtype;
	switch (Rtype) {
	    case INTSXP: case LGLSXP: {
		int *out = Rtype == INTSXP ? INTEGER(Rarray) : LOGICAL(Rarray);
		if (copy_without_type_casting)
			return copy_selected_int_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_CHAR)
			return copy_selected_char_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SCHAR)
			return copy_selected_schar_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UCHAR)
			return copy_selected_uchar_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SHORT)
			return copy_selected_short_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_USHORT)
			return copy_selected_ushort_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		break;
	    }
	    case REALSXP: {
		double *out = REAL(Rarray);
		if (copy_without_type_casting)
			return copy_selected_double_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_CHAR)
			return copy_selected_char_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SCHAR)
			return copy_selected_schar_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UCHAR)
			return copy_selected_uchar_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SHORT)
			return copy_selected_short_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_USHORT)
			return copy_selected_ushort_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_INT)
			return copy_selected_int_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UINT)
			return copy_selected_uint_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_LONG)
			return copy_selected_long_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_ULONG)
			return copy_selected_ulong_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_LLONG)
			return copy_selected_llong_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_ULLONG)
			return copy_selected_ullong_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		if (chunk_data_buf->data_type_id == H5T_NATIVE_FLOAT)
			return copy_selected_float_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		break;
	    }
	    case RAWSXP: {
		Rbyte *out = RAW(Rarray);
		if (copy_without_type_casting)
			return copy_selected_uchar_chunk_data_to_uchar_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
		break;
	    }
	}
	PRINT_TO_ERRMSG_BUF("unsupported dataset type");
	return -1;
}


/****************************************************************************
 * read_data_4_5()
 *
 * method 4: One call to _read_H5Viewport() or _read_h5chunk() (wrappers for
 * H5Dread() or H5Dread_chunk(), respectively) per chunk touched by the
 * user-supplied array selection.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Make one call to _read_H5Viewport() or _read_h5chunk() to load the
 *     **entire** chunk data to an intermediate buffer.
 *   - Copy the user-selected data from the intermediate buffer to 'Rarray'.
 *
 * method 5: Like method 4 but bypasses the intermediate buffer if a
 * chunk is fully selected.
 */

static int read_data_4_5(ChunkIterator *chunk_iter,
		const size_t *Rarray_dim, SEXP Rarray,
		int method, int use_H5Dread_chunk)
{
	if (use_H5Dread_chunk)
		warning("using 'use.H5Dread_chunk=TRUE' is still "
			"experimental, use at your own risk");

	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	size_t *inner_midx_buf = R_alloc0_size_t_array(ndim);

	void *out = DATAPTR(Rarray);
	if (out == NULL)
		return -1;

	hid_t out_space_id = _create_mem_space(ndim, Rarray_dim);
	if (out_space_id < 0)
		return -1;

	ChunkDataBuffer chunk_data_buf;
	int ret = _init_ChunkDataBuffer(&chunk_data_buf, h5dset, 0);
	if (ret < 0) {
		H5Sclose(out_space_id);
		return ret;
	}
	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		//_print_tchunk_info(chunk_iter);
		int direct_load = method == 5 && _tchunk_is_fully_selected(ndim,
							&chunk_iter->h5dset_vp,
							&chunk_iter->mem_vp);
		if (direct_load) {
			/* Load the chunk **directly** into 'Rarray' (no
			   intermediate buffer). */
			ret = _read_H5Viewport(h5dset,
				&chunk_iter->h5dset_vp,
				h5dset->h5type->native_type_id_for_Rtype,
				out_space_id, out,
				&chunk_iter->mem_vp);
		} else {
			/* Load the **entire** chunk to an intermediate
			   buffer then copy the user-selected chunk data
			   from the intermediate buffer to 'Rarray'. */
			//clock_t t0 = clock();
			ret = _load_chunk(chunk_iter,
					&chunk_data_buf,
					use_H5Dread_chunk);
			if (ret < 0)
				break;
			//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
			//printf("- load chunk: %3.3f ms\n", dt);

			long long int nvals =
				copy_selected_chunk_data_to_Rarray(
						chunk_iter,
						&chunk_data_buf,
						inner_midx_buf,
						Rarray_dim, Rarray);
			if (nvals < 0)
				break;
			if (h5dset->h5type->is_variable_str)
				_reclaim_vlen_bufs(&chunk_data_buf);
		}
		if (ret < 0)
			break;
	}
	_destroy_ChunkDataBuffer(&chunk_data_buf);
	H5Sclose(out_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to _read_h5selection() (wrapper for H5Dread()) per chunk touched
 * by the user-supplied array selection. No intermediate buffer.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Select the hyperslabs obtained by intersecting the user-supplied
 *     array selection with the current chunk.
 *   - Call _read_h5selection(). This loads the selected data **directly**
 *     to the final R array.
 */

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP index,
		const H5Viewport *out_vp,
		LLongAEAE *inner_breakpoint_bufs, size_t *inner_nchip_buf)
{
	for (int along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		LLongAE *inner_breakpoint_buf =
					inner_breakpoint_bufs->elts[along];
		LLongAE_set_nelt(inner_breakpoint_buf, 0);
		long long int d = out_vp->dim[along];
		SEXP start = GET_LIST_ELT(index, along);
		if (start == R_NilValue) {
			LLongAE_insert_at(inner_breakpoint_buf, 0, d);
			inner_nchip_buf[along] = 1;
			continue;
		}
		R_xlen_t off = out_vp->off[along];
		long long int s1 = get_trusted_elt(start, off);
		size_t nchip = 0;
		for (long long int i = 1; i < d; i++) {
			long long int s0 = s1;
			s1 = get_trusted_elt(start, off + (R_xlen_t) i);
			if (s1 != s0 + 1)
				LLongAE_insert_at(inner_breakpoint_buf,
						  nchip++, i);
		}
		LLongAE_insert_at(inner_breakpoint_buf, nchip++, d);
		inner_nchip_buf[along] = nchip;
	}
	return;
}

static void init_inner_vp(int ndim, SEXP index,
		const H5Viewport *h5dset_vp,
		H5Viewport *inner_vp)
{
	int along, h5along;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		SEXP start = GET_LIST_ELT(index, along);
		hsize_t d;
		if (start == R_NilValue) {
			inner_vp->h5off[h5along] = h5dset_vp->h5off[h5along];
			d = h5dset_vp->h5dim[h5along];
		} else {
			d = 1;
		}
		inner_vp->h5dim[h5along] = d;
	}
	return;
}

static void update_inner_vp(int ndim,
		SEXP index, const H5Viewport *out_vp,
		const size_t *inner_midx, int inner_moved_along,
		const LLongAEAE *inner_breakpoint_bufs,
		H5Viewport *inner_vp)
{
	for (int along = 0; along < ndim; along++) {
		if (along > inner_moved_along)
			break;
		if (index == R_NilValue)
			continue;
		SEXP start = VECTOR_ELT(index, along);
		if (start == R_NilValue)
			continue;
		const long long int *inner_breakpoint =
				inner_breakpoint_bufs->elts[along]->elts;
		size_t idx = inner_midx[along];
		long long int off = idx == 0 ? 0 : inner_breakpoint[idx - 1];
		long long int d = inner_breakpoint[idx] - off;
		R_xlen_t i = out_vp->off[along] + off;
		int h5along = ndim - 1 - along;
		inner_vp->h5off[h5along] = (hsize_t)
						get_trusted_elt(start, i) - 1;
		inner_vp->h5dim[h5along] = (hsize_t) d;
	}
	return;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_intersection_of_chips_with_chunk(
		const H5DSetDescriptor *h5dset, SEXP index,
		const H5Viewport *out_vp, const H5Viewport *h5dset_vp,
		const LLongAEAE *inner_breakpoint_bufs,
		const size_t *inner_nchip,
		size_t *inner_midx_buf,
		H5Viewport *inner_vp)
{
	int ret = H5Sselect_none(h5dset->h5space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	int ndim = h5dset->ndim;

	init_inner_vp(ndim, index, h5dset_vp, inner_vp);

	/* Walk on the "inner chips" i.e. on the intersections between
	   the "chips" in the user-supplied array selection and the currrent
	   chunk. */
	long long int num_hyperslabs = 0;
	int inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_inner_vp(ndim, index, out_vp,
				inner_midx_buf, inner_moved_along,
				inner_breakpoint_bufs,
				inner_vp);
		ret = _add_H5Viewport_to_h5selection(h5dset->h5space_id,
						     inner_vp);
		if (ret < 0)
			return -1;
		inner_moved_along = next_midx(ndim, inner_nchip,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);
	return num_hyperslabs;
}

static int direct_load_selected_chunk_data(
		const ChunkIterator *chunk_iter,
		size_t *inner_midx_buf,
		H5Viewport *inner_vp,
		LLongAEAE *inner_breakpoint_bufs,
		size_t *inner_nchip_buf,
		hid_t out_space_id, void *out)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	update_inner_breakpoints(h5dset->ndim, chunk_iter->moved_along,
			chunk_iter->index, &chunk_iter->mem_vp,
			inner_breakpoint_bufs, inner_nchip_buf);
	int ret = select_intersection_of_chips_with_chunk(
			h5dset, chunk_iter->index,
			&chunk_iter->mem_vp, &chunk_iter->h5dset_vp,
			inner_breakpoint_bufs, inner_nchip_buf,
			inner_midx_buf, inner_vp);
	if (ret < 0)
		return ret;
	return _read_h5selection(h5dset,
				h5dset->h5type->native_type_id_for_Rtype,
				out_space_id, out,
				&chunk_iter->mem_vp);
}

static int read_data_6(ChunkIterator *chunk_iter,
		const size_t *Rarray_dim, SEXP Rarray)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	size_t *inner_midx_buf  = R_alloc0_size_t_array(ndim);
	size_t *inner_nchip_buf = R_alloc0_size_t_array(ndim);
	LLongAEAE *inner_breakpoint_bufs = new_LLongAEAE(ndim, ndim);

	void *out = DATAPTR(Rarray);
	if (out == NULL)
		return -1;

	hid_t out_space_id = _create_mem_space(ndim, Rarray_dim);
	if (out_space_id < 0)
		return -1;

	H5Viewport inner_vp;
	int ret = _alloc_H5Viewport(&inner_vp, ndim, ALLOC_H5OFF_AND_H5DIM);
	if (ret < 0) {
		H5Sclose(out_space_id);
		return ret;
	}

	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		ret = direct_load_selected_chunk_data(
			chunk_iter,
			inner_midx_buf,
			&inner_vp,
			inner_breakpoint_bufs, inner_nchip_buf,
			out_space_id, out);
		if (ret < 0)
			break;
	}
	_free_H5Viewport(&inner_vp);
	H5Sclose(out_space_id);
	return ret;
}


/****************************************************************************
 * _h5mread_index()
 *
 * Implements methods 4 to 6.
 * Return an ordinary array or R_NilValue if an error occured.
 */

SEXP _h5mread_index(ChunkIterator *chunk_iter, int method,
		int use_H5Dread_chunk, const size_t *ans_dim)
{
	const H5DSetDescriptor *h5dset = chunk_iter->h5dset;
	int ndim = h5dset->ndim;
	R_xlen_t ans_len = 1;
	for (int along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	SEXP ans = PROTECT(allocVector(h5dset->h5type->Rtype, ans_len));

	int ret;
	if (method <= 5) {
		/* methods 4 and 5 */
		ret = read_data_4_5(chunk_iter, ans_dim, ans,
				    method, use_H5Dread_chunk);
	} else {
		/* method 6 */
		ret = read_data_6(chunk_iter, ans_dim, ans);
	}

	UNPROTECT(1);
	return ret < 0 ? R_NilValue : ans;
}

