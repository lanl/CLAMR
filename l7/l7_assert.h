#ifndef L7_ASSERT_H
#define L7_ASSERT_H

/*
 * remove typesafe linkage if compiling under c++
 */
#ifdef __cpluscplus
extern "C"
{
#endif

	/*
	 * define L7_ASSERT (mainly used for checking return status)
	 */
	
#define L7_PRINT(                                                    \
	EXPRESSION,                                                      \
	ERROR_MESSAGE,                                                   \
	RETURN_VALUE                                                     \
	)                                                                \
	fflush(l7.assert_out_file ? l7.assert_out_file:stderr);          \
	fprintf(l7.assert_out_file ? l7.assert_out_file:stderr,          \
			" [pe:%d] **L7 Error** %s (%s:%d[!(%s)]) (%s) (%d)\n",   \
			l7.penum,                                                \
			(L7_LOCATION),                                           \
			__FILE__,                                                \
			__LINE__,                                                \
            #EXPRESSION,                                             \
			(ERROR_MESSAGE),                                         \
			(RETURN_VALUE)                                           \
			);                                                       \
	fflush(l7.assert_out_file ? l7.assert_out_file:stderr);
	
#if defined(L7_NO_ASSERT_PRINT) /* L7_NO_ASSERT_PRINT */
#define L7_ASSERT(                                                   \
		EXPRESSION,                                                  \
	    ERROR_MESSAGE,                                               \
	    RETURN_VALUE                                                 \
	    )                                                            \
	    if ( ! (EXPRESSION) ){                                       \
	      return (RETURN_VALUE);                                     \
        }

#define L7_ASSERTN(                                                  \
		EXPRESSION,                                                  \
	    ERROR_MESSAGE,                                               \
	    RETURN_VALUE                                                 \
	    )                                                            \
	    if ( ! (EXPRESSION) ){                                       \
	      return;                                                    \
        }
#else /* ! L7_NO_ASSERT_PRINT */
#define L7_ASSERT(                                                   \
		EXPRESSION,                                                  \
	    ERROR_MESSAGE,                                               \
	    RETURN_VALUE                                                 \
	    )                                                            \
	    if ( ! (EXPRESSION) ){                                       \
	      L7_PRINT( EXPRESSION, ERROR_MESSAGE, RETURN_VALUE)         \
	      return (RETURN_VALUE);                                     \
        }

#define L7_ASSERTN(                                                  \
		EXPRESSION,                                                  \
	    ERROR_MESSAGE,                                               \
	    RETURN_VALUE                                                 \
	    )                                                            \
	    if ( ! (EXPRESSION) ){                                       \
	      L7_PRINT( EXPRESSION, ERROR_MESSAGE, RETURN_VALUE)         \
	      return;                                                    \
        }
#endif /* ! L7_NO_ASSERT_PRINT */
	
	/*
	 * remove typesafe linkage if compiling under c++
	 */
	
#ifdef __cpluscplus
}
#endif

#endif /* L7_ASSERT_H */
