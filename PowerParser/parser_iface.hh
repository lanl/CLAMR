/* Copyright 2015.  Los Alamos National Security, LLC. This material was produced
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS
 * ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 * to produce derivative works, such modified software should be clearly marked,
 * so as not to confuse it with the version available from LANL.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 *
 * Under this license, it is required to include a reference to this work. We
 * request that each derivative work contain a reference to LANL Copyright
 * Disclosure C15076/LA-CC-15-054 so that this work's impact can be roughly
 * measured.
 *
 * This is LANL Copyright Disclosure C15076/LA-CC-15-054
 */

/*
 *  PowerParser is a general purpose input file parser for software applications.
 *
 *  Authors: Chuck Wingate   XCP-2   caw@lanl.gov
 *           Robert Robey    XCP-2   brobey@lanl.gov
 */

#ifndef PARSERIFACEHH
#define PARSERIFACEHH

#if !defined FC_FUNC
  #if defined aix || defined hpux || defined IBM
    #define    FC_FUNC(a,A)    a
  #elif defined _SQ && defined GNU
    #define    FC_FUNC(a,A)    a
  #else
    #define    FC_FUNC(a,A)    a ## _
  #endif
#endif

/* Function prototypes. */
void FC_FUNC(parser_init_comm,PARSER_INIT_COMM)();
void FC_FUNC(parser_comm_info,PARSER_COMM_INFO)(int *mype,
                                                int *npes, int *iope);

void FC_FUNC(parser_create,PARSER_CREATE)(char *oargs, int *len,
                                          char *fname, int *fname_len);
void FC_FUNC(parser_compile_buffer,PARSER_COMPILE_BUFFER)(int &return_value);
void FC_FUNC(parser_recreate,PARSER_RECREATE)();
void FC_FUNC(parser_destroy,PARSER_DESTROY)();


void FC_FUNC(get_logical0,GET_LOGICAL0)(char *cmdname, int *cmdvalue,
                                        int *len, int *iskip);
void FC_FUNC(get_logical1,GET_LOGICAL1)(char *cmdname, int *array_values,
                                        int *array_size, int *len,
                                        int *iskip);
void FC_FUNC(get_logical2,GET_LOGICAL2)(char *cmdname, int *array_values,
                                        int *size1, int *size2,
                                        int *len, int *iskip);
void FC_FUNC(get_logical3,GET_LOGICAL3)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *len, int *iskip);
void FC_FUNC(get_logical4,GET_LOGICAL4)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *size4, int *len, int *iskip);


void FC_FUNC(get_integer0,GET_INTEGER0)(char *cmdname, int *cmdvalue,
                                        int *len, int *iskip);
void FC_FUNC(get_integer1,GET_INTEGER1)(char *cmdname, int *array_values,
                                        int *array_size, int *len,
                                        int *iskip);
void FC_FUNC(get_integer2,GET_INTEGER2)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *len,
                                        int *iskip);
void FC_FUNC(get_integer3,GET_INTEGER3)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *len, int *iskip);
void FC_FUNC(get_integer4,GET_INTEGER4)(char *cmdname, int *array_values,
                                        int *size1, int *size2, int *size3,
                                        int *size4, int *len, int *iskip);

void FC_FUNC(get_int8_0,GET_INT8_0)(char *cmdname, int64_t *cmdvalue,
                                    int *len, int *iskip);

void FC_FUNC(get_real0,GET_REAL0)(char *cmdname, double *cmdvalue,
                                  int *len, int *iskip);
void FC_FUNC(get_real1,GET_REAL1)(char *cmdname, double *array_values,
                                  int *array_size, int *len, int *iskip);
void FC_FUNC(get_real2,GET_REAL2)(char *cmdname, double *array_values,
                                  int *size1, int *size2, int *len,
                                  int *iskip);
void FC_FUNC(get_real3,GET_REAL3)(char *cmdname, double *array_values,
                                  int *size1, int *size2, int *size3,
                                  int *len, int *iskip);
void FC_FUNC(get_real4,GET_REAL4)(char *cmdname, double *array_values,
                                  int *size1, int *size2, int *size3,
                                  int *size4, int *len, int *iskip);


void FC_FUNC(get_char00,GET_CHAR00)(char *cmdname, char *cmdvalue,
                                    int *len, int *iskip);
void FC_FUNC(get_char0,GET_CHAR0)(char *cmdname, char *cmdvalue,
                                  int *len, int *len_value, int *iskip);
void FC_FUNC(get_char1,GET_CHAR1)(char *cmdname, char *array_values,
                                  int *nchar, int *size1,
                                  int *len_name, int *iskip);
void FC_FUNC(get_char2,GET_CHAR2)(char *cmdname, char *array_values,
                                  int *nchar, int *size1, int *size2,
                                  int *len_name, int *iskip);
void FC_FUNC(get_char3,GET_CHAR3)(char *cmdname, char *array_values,
                                  int *nchar, int *size1, int *size2,
                                  int *size3, int *len_name, int *iskip);
void FC_FUNC(get_char4,GET_CHAR4)(char *cmdname, char *array_values,
                                  int *nchar, int *size1, int *size2,
                                  int *size3, int *size4, int *len_name,
                                  int *iskip);


void FC_FUNC(get_size1,GET_SIZE1)(char *cmdname, int *array_size, int *len);
void FC_FUNC(get_size2,GET_SIZE2)(char *cmdname, int *size1,
                                  int *size2, int *len);
void FC_FUNC(get_sizeb2,GET_SIZEB2)(char *cmdname, int *size1,
                                    int *size2, int *len);
void FC_FUNC(get_size3,GET_SIZE3)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *len);
void FC_FUNC(get_size4,GET_SIZE4)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *size4,
                                  int *len);
void FC_FUNC(get_size5,GET_SIZE5)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *size4,
                                  int *size5, int *len);
void FC_FUNC(get_size6,GET_SIZE6)(char *cmdname, int *size1,
                                  int *size2, int *size3, int *size4,
                                  int *size5, int *size6, int *len);


void FC_FUNC(process_error_final,PROCESS_ERROR_FINAL)(int &return_value);
void FC_FUNC(check_processed,CHECK_PROCESSED)(int *good_int);
void FC_FUNC(cmd_in_input,CMD_IN_INPUT)(char *cmdname, int *cmdlen,
                                        int *in_input_i, int *in_whenthen_i);
void FC_FUNC(cmd_set_processed,CMD_SET_PROCESSED)(char *cmdname, int *cmdlen,
                                                  int *bval_i);


void FC_FUNC(list_functions_start,LIST_FUNCTIONS_START)();
void FC_FUNC(list_variables_start,LIST_VARIABLES_START)();
void FC_FUNC(echo_ui_start,ECHO_UI_START)();
void FC_FUNC(echo_final_buffer,ECHO_FINAL_BUFFER)();
void FC_FUNC(echo_wt_final_buffer,ECHO_WT_FINAL_BUFFER)();
void FC_FUNC(get_output_line,GET_OUTPUT_LINE)(char *fline, int *nchar,
                                              int *bint);
void FC_FUNC(parser_chk_scalar_dup,PARSER_CHK_SCALAR_DUP)();

void FC_FUNC(parser_get_rbnum,PARSER_GET_RBNUM)(int *rbnum);
void FC_FUNC(parser_set_rbnum,PARSER_SET_RBNUM)(int *rbnum);
void FC_FUNC(parser_rb_check,PARSER_RB_CHECK)(char *varname_array,
    char *val_array, int *vv_active_a, int *nchar, int *nvv, int *rbci,
    int *rb_ntriggered, int *rb_triggered_indices);
void FC_FUNC(parser_set_rb_names,PARSER_SET_RB_NAMES)(char *rb_names_1d,
                                                      int *rb_num, int *ncstr_rb);
void FC_FUNC(parser_get_rb_names,PARSER_GET_RB_NAMES)(char *rb_names_1d,
                                                      int *rb_num, int *ncstr_rb);
void FC_FUNC(parser_set_rb_aflags,PARSER_SET_RB_AFLAGS)(int *rb_aflags,
                                                        int *rb_num);
void FC_FUNC(parser_get_rb_aflags,PARSER_GET_RB_AFLAGS)(int *rb_aflags);
void FC_FUNC(parser_get_rb_satsize,PARSER_GET_RB_SATSIZE)(int *rb_satsize);
void FC_FUNC(parser_set_rb_satsize,PARSER_SET_RB_SATSIZE)(int *rb_satsize);
void FC_FUNC(parser_get_rb_satprb,PARSER_GET_RB_SATPRB)(int *rb_satprb);
void FC_FUNC(parser_set_rb_satprb,PARSER_SET_RB_SATPRB)(int *rb_satprb,
                                                        int *rb_num);
void FC_FUNC(parser_get_rb_sat,PARSER_GET_RB_SAT)(int *rb_sat);
void FC_FUNC(parser_set_rb_sat,PARSER_SET_RB_SAT)(int *rb_sat, int *rb_satsize);
void FC_FUNC(echo_rb_info,ECHO_RB_INFO)();
void FC_FUNC(echo_rb1_info,ECHO_RB1_INFO)(int *rb);
void FC_FUNC(parser_get_rb_num_varnames,PARSER_GET_RB_NUM_VARNAMES)(int *vnum);
void FC_FUNC(parser_get_rb_varnames,PARSER_GET_RB_VARNAMES)(char *rb_varnames_1d,
    int *rb_vnum, int *ncstr_rb);


void FC_FUNC(parser_get_wtnum,PARSER_GET_WTNUM)(int *wtnum);
void FC_FUNC(parser_wt_check2,PARSER_WT_CHECK2)(int *wti, char *varname_array,
    char *val_array, int *vv_active_a, int *nchar, int *nvv, int *wtci);
void FC_FUNC(parser_wt_set_cmdsfp,PARSER_WT_SET_CMDSFP)(int *wti);
void FC_FUNC(parser_wt_reset,PARSER_WT_RESET)();

void FC_FUNC(parser_wt_casize,PARSER_WT_CASIZE)(int *wti, int *wt_casize);
void FC_FUNC(parser_wt_ca,PARSER_WT_CA)(int *wti, char *wt_ca, int *wt_casize);

void FC_FUNC(parser_wt_satsize,PARSER_WT_SATSIZE)(int *wti, int *wt_satsize);
void FC_FUNC(parser_wt_getsat,PARSER_WT_GETSAT)(int *wti, int *wt_sat, int *wt_satsize);
void FC_FUNC(parser_wt_setsat,PARSER_WT_SETSAT)(int *wti, int *wt_sat, int *wt_satsize);
void FC_FUNC(parser_wt_getprocessed,PARSER_WT_GETPROCESSED)(int *wti, int *wtp);
void FC_FUNC(parser_wt_setprocessed,PARSER_WT_SETPROCESSED)(int *wti, int *wtp);
void FC_FUNC(parser_wt_getseq,PARSER_WT_GETSEQ)(int *wti, int *wtseq);
void FC_FUNC(parser_wt_setseq,PARSER_WT_SETSEQ)(int *wti, int *wtseq);

void FC_FUNC(parser_dictionary_add,PARSER_DICTIONARY_ADD)(char *name, double value, bool pred, char *vdesc);
void FC_FUNC(parser_dictionary_env_add,PARSER_DICTIONARY_ENV_ADD)(char *name, bool pred, char *vdesc);

#endif /* PARSERIFACEHH */
