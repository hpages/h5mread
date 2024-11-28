/****************************************************************************
 *                   h5mread global error message buffer                    *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "global_errmsg_buf.h"


char * _h5mread_global_errmsg_buf()
{
	static char buf[ERRMSG_BUF_LENGTH];

	return buf;
}

