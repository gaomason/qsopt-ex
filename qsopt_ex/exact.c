/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2005 Daniel Espinoza.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * */
/* ========================================================================= */
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "exact.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "logging-private.h"

#include "util.h"
#include "eg_timer.h"
#include "eg_exutil.h"
#include "except.h"
#include "timing_log.h"

#include "basis_mpq.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "fct_mpq.h"
#include "lpdata_mpq.h"
#include "simplex_mpq.h"
#include "factor_mpq.h"
#include "factor_dbl.h"
#include "factor_mpf.h"
#include "basis_dbl.h"
#include "lpdata_dbl.h"
#include "fct_dbl.h"
#include "simplex_dbl.h"
#include "simplex_mpf.h"

/* ========================================================================= */
int QSexact_print_sol (mpq_QSdata * p,
											 EGioFile_t * out_f)
{
	int rval = 0,
	  status;
	const int ncols = mpq_QSget_colcount (p);
	const int nrows = mpq_QSget_rowcount (p);
	mpq_t *x = mpq_EGlpNumAllocArray (ncols);
	mpq_t *rc = mpq_EGlpNumAllocArray (ncols);
	mpq_t *slack = mpq_EGlpNumAllocArray (nrows);
	mpq_t *pi = mpq_EGlpNumAllocArray (nrows);
	mpq_t value;
	register int i;
	char *str1 = 0;
	mpq_init (value);
	EGcallD(mpq_QSget_status (p, &status));
	if(mpq_QSget_x_array (p, x)) mpq_EGlpNumFreeArray (x);
	if(mpq_QSget_slack_array (p, slack)) mpq_EGlpNumFreeArray (slack);
	if(mpq_QSget_pi_array (p, pi)) mpq_EGlpNumFreeArray (pi);
	if(mpq_QSget_rc_array (p, rc)) mpq_EGlpNumFreeArray (rc);

	switch (status)
	{
	case QS_LP_OPTIMAL:
		EGcallD(mpq_QSget_objval (p, &value));
		str1 = mpq_EGlpNumGetStr (value);
		EGioPrintf (out_f, "status OPTIMAL\n\tValue = %s\n", str1);
		free (str1);
		str1 = 0;
		break;
	case QS_LP_INFEASIBLE:
		EGioPrintf (out_f, "status INFEASIBLE\n");
		break;
	case QS_LP_UNBOUNDED:
		EGioPrintf (out_f, "status UNBOUNDED\n");
		break;
	case QS_LP_ITER_LIMIT:
	case QS_LP_TIME_LIMIT:
	case QS_LP_UNSOLVED:
	case QS_LP_ABORTED:
	case QS_LP_MODIFIED:
		EGioPrintf (out_f, "status NOT_SOLVED\n");
		break;
	}
	if (x)
	{
		EGioPrintf (out_f, "VARS:\n");
		for (i = 0; i < ncols; i++)
			if (!mpq_equal (x[i], __zeroLpNum_mpq__))
			{
				str1 = mpq_EGlpNumGetStr (x[i]);
				EGioPrintf (out_f, "%s = %s\n", p->qslp->colnames[i], str1);
				free (str1);
			}
	}
	if (rc)
	{
		EGioPrintf (out_f, "REDUCED COST:\n");
		for (i = 0; i < ncols; i++)
			if (!mpq_equal (rc[i], __zeroLpNum_mpq__))
			{
				str1 = mpq_EGlpNumGetStr (rc[i]);
				EGioPrintf (out_f, "%s = %s\n", p->qslp->colnames[i], str1);
				free (str1);
			}
	}
	if (pi)
	{
		EGioPrintf (out_f, "PI:\n");
		for (i = 0; i < nrows; i++)
			if (!mpq_equal (pi[i], __zeroLpNum_mpq__))
			{
				str1 = mpq_EGlpNumGetStr (pi[i]);
				EGioPrintf (out_f, "%s = %s\n", p->qslp->rownames[i], str1);
				free (str1);
			}
	}
	if (slack)
	{
		EGioPrintf (out_f, "SLACK:\n");
		for (i = 0; i < nrows; i++)
			if (!mpq_equal (slack[i], __zeroLpNum_mpq__))
			{
				str1 = mpq_EGlpNumGetStr (slack[i]);
				EGioPrintf (out_f, "%s = %s\n", p->qslp->rownames[i], str1);
				free (str1);
			}
	}

	/* ending */
CLEANUP:
	if (x)
		mpq_EGlpNumFreeArray (x);
	if (pi)
		mpq_EGlpNumFreeArray (pi);
	if (rc)
		mpq_EGlpNumFreeArray (rc);
	if (slack)
		mpq_EGlpNumFreeArray (slack);
	mpq_clear (value);
	return rval;
}

/* ========================================================================= */
dbl_QSdata *QScopy_prob_mpq_dbl (mpq_QSdata * p, const char *newname)
{
	// clock start for timing purposes
        clock_t start = clock();
	
	const int ncol = mpq_QSget_colcount(p);
	const int nrow = mpq_QSget_rowcount(p);
	char*sense=0;
	int*rowcnt=0;
	int*rowbeg=0;
	int*rowind=0;
	int objsense;
	mpq_t*mpq_lb=0;
	mpq_t*mpq_ub=0;
	mpq_t*mpq_obj=0;
	mpq_t*mpq_range=0;
	mpq_t*mpq_rhs=0;
	mpq_t*mpq_rowval=0;
	double*dbl_lb=0;
	double*dbl_ub=0;
	double*dbl_obj=0;
	double*dbl_range=0;
	double*dbl_rhs=0;
	double*dbl_rowval=0;
	dbl_QSdata *p2 = 0;
	int rval = 0;
	register int i;
	mpq_t mpq_val;
	double dbl_val;
	mpq_init(mpq_val); // AP: equivalent to mpq_EG1pNumInitVar(mpq_val)
	/* get all information */
	EGcallD(mpq_QSget_objsense(p,&objsense));
	mpq_lb = mpq_EGlpNumAllocArray(ncol);
	mpq_ub = mpq_EGlpNumAllocArray(ncol);
	EGcallD(mpq_QSget_bounds(p,mpq_lb,mpq_ub));
	dbl_lb = QScopy_array_mpq_dbl(mpq_lb);
	dbl_ub = QScopy_array_mpq_dbl(mpq_ub);
	mpq_EGlpNumFreeArray(mpq_ub);
	mpq_obj = mpq_lb;
	mpq_lb = 0;
	EGcallD(mpq_QSget_obj(p, mpq_obj));
	dbl_obj = QScopy_array_mpq_dbl(mpq_obj);
	mpq_EGlpNumFreeArray(mpq_obj);
	EGcallD(mpq_QSget_ranged_rows(p, &rowcnt, &rowbeg, &rowind, &mpq_rowval,
																&mpq_rhs, &sense, &mpq_range, 0));
	dbl_rowval = QScopy_array_mpq_dbl(mpq_rowval);
	mpq_EGlpNumFreeArray(mpq_rowval);
	dbl_range = QScopy_array_mpq_dbl(mpq_range);
	mpq_EGlpNumFreeArray(mpq_range);
	dbl_rhs = QScopy_array_mpq_dbl(mpq_rhs);
	mpq_EGlpNumFreeArray(mpq_rhs);
	/* create copy */
	p2 = dbl_QScreate_prob (newname, objsense);
	if (!p2) goto CLEANUP;
	for( i = 0 ; i < ncol; i++)
	{
		EGcallD(dbl_QSnew_col(p2, dbl_obj[i], dbl_lb[i], dbl_ub[i], 0));
	}
	dbl_EGlpNumFreeArray(dbl_lb);
	dbl_EGlpNumFreeArray(dbl_ub);
	dbl_EGlpNumFreeArray(dbl_obj);
	EGcallD(dbl_QSadd_ranged_rows(p2, nrow, rowcnt, rowbeg, rowind, dbl_rowval, dbl_rhs, sense, dbl_range, 0));
	/* set parameters */
	EGcallD(mpq_QSget_param(p, QS_PARAM_PRIMAL_PRICING, &objsense));
	EGcallD(dbl_QSset_param(p2, QS_PARAM_PRIMAL_PRICING, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_DUAL_PRICING, &objsense));
	EGcallD(dbl_QSset_param(p2, QS_PARAM_DUAL_PRICING, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_SIMPLEX_DISPLAY, &objsense));
	EGcallD(dbl_QSset_param(p2, QS_PARAM_SIMPLEX_DISPLAY, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_SIMPLEX_MAX_ITERATIONS, &objsense));
	EGcallD(dbl_QSset_param(p2, QS_PARAM_SIMPLEX_MAX_ITERATIONS, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_SIMPLEX_SCALING, &objsense));
	EGcallD(dbl_QSset_param(p2, QS_PARAM_SIMPLEX_SCALING, objsense));
	EGcallD(mpq_QSget_param_EGlpNum(p, QS_PARAM_SIMPLEX_MAX_TIME, &mpq_val));
	dbl_val = mpq_get_d(mpq_val);
	EGcallD(dbl_QSset_param_EGlpNum(p2, QS_PARAM_SIMPLEX_MAX_TIME, dbl_val));
	EGcallD(mpq_QSget_param_EGlpNum(p, QS_PARAM_OBJULIM, &mpq_val));
	dbl_val = mpq_get_d(mpq_val);
	EGcallD(dbl_QSset_param_EGlpNum(p2, QS_PARAM_OBJULIM, dbl_val));
	EGcallD(mpq_QSget_param_EGlpNum(p, QS_PARAM_OBJLLIM, &mpq_val));
	dbl_val = mpq_get_d(mpq_val);
	EGcallD(dbl_QSset_param_EGlpNum(p2, QS_PARAM_OBJLLIM, dbl_val));
	/* ending */
	CLEANUP:
	mpq_clear(mpq_val);
	dbl_EGlpNumFreeArray(dbl_rowval);
	dbl_EGlpNumFreeArray(dbl_range);
	dbl_EGlpNumFreeArray(dbl_rhs);
	dbl_EGlpNumFreeArray(dbl_lb);
	dbl_EGlpNumFreeArray(dbl_ub);
	dbl_EGlpNumFreeArray(dbl_obj);
	mpq_EGlpNumFreeArray(mpq_rowval);
	mpq_EGlpNumFreeArray(mpq_range);
	mpq_EGlpNumFreeArray(mpq_rhs);
	mpq_EGlpNumFreeArray(mpq_lb);
	mpq_EGlpNumFreeArray(mpq_ub);
	mpq_EGlpNumFreeArray(mpq_obj);
	EGfree(rowcnt);
	EGfree(rowbeg);
	EGfree(rowind);
	EGfree(sense);
	if (rval)
	{
		dbl_QSfree_prob (p2);
		p2 = 0;
	}
	// testing for log file
        clock_t end = clock();
        double duration = (double)(end - start) / CLOCKS_PER_SEC;
        log_timing("QScopy_prob_mpq_dbl took ", duration);

#if QSEXACT_SAVE_INT
	else
	{
		dbl_QSwrite_prob (p2, "prob.dbl.lp", "LP");
	}
#endif
	return p2;
}

/* ========================================================================= */
mpf_QSdata *QScopy_prob_mpq_mpf (mpq_QSdata * p, const char *newname)
{
	const int ncol = mpq_QSget_colcount(p);
	const int nrow = mpq_QSget_rowcount(p);
	char*sense=0;
	int*rowcnt=0;
	int*rowbeg=0;
	int*rowind=0;
	int objsense;
	mpq_t*mpq_lb=0;
	mpq_t*mpq_ub=0;
	mpq_t*mpq_obj=0;
	mpq_t*mpq_range=0;
	mpq_t*mpq_rhs=0;
	mpq_t*mpq_rowval=0;
	mpf_t*mpf_lb=0;
	mpf_t*mpf_ub=0;
	mpf_t*mpf_obj=0;
	mpf_t*mpf_range=0;
	mpf_t*mpf_rhs=0;
	mpf_t*mpf_rowval=0;
	mpf_QSdata *p2 = 0;
	int rval = 0;
	mpq_t mpq_val;
	mpf_t mpf_val;
	register int i;
	mpq_init(mpq_val);
	mpf_init(mpf_val);
	/* get all information */
	EGcallD(mpq_QSget_objsense(p,&objsense));
	mpq_lb = mpq_EGlpNumAllocArray(ncol);
	mpq_ub = mpq_EGlpNumAllocArray(ncol);
	EGcallD(mpq_QSget_bounds(p,mpq_lb,mpq_ub));
	mpf_lb = QScopy_array_mpq_mpf(mpq_lb);
	mpf_ub = QScopy_array_mpq_mpf(mpq_ub);
	mpq_EGlpNumFreeArray(mpq_ub);
	mpq_obj = mpq_lb;
	mpq_lb = 0;
	EGcallD(mpq_QSget_obj(p, mpq_obj));
	mpf_obj = QScopy_array_mpq_mpf(mpq_obj);
	mpq_EGlpNumFreeArray(mpq_obj);
	EGcallD(mpq_QSget_ranged_rows(p, &rowcnt, &rowbeg, &rowind, &mpq_rowval, &mpq_rhs, &sense, &mpq_range, 0));
	mpf_rowval = QScopy_array_mpq_mpf(mpq_rowval);
	mpq_EGlpNumFreeArray(mpq_rowval);
	mpf_range = QScopy_array_mpq_mpf(mpq_range);
	mpq_EGlpNumFreeArray(mpq_range);
	mpf_rhs = QScopy_array_mpq_mpf(mpq_rhs);
	mpq_EGlpNumFreeArray(mpq_rhs);
	/* create copy */
	p2 = mpf_QScreate_prob (newname, objsense);
	if (!p2) goto CLEANUP;
	for( i = 0 ; i < ncol; i++)
	{
		EGcallD(mpf_QSnew_col(p2, mpf_obj[i], mpf_lb[i], mpf_ub[i], 0));
	}
	mpf_EGlpNumFreeArray(mpf_lb);
	mpf_EGlpNumFreeArray(mpf_ub);
	mpf_EGlpNumFreeArray(mpf_obj);
	EGcallD(mpf_QSadd_ranged_rows(p2, nrow, rowcnt, rowbeg, rowind, (const mpf_t*)mpf_rowval,(const mpf_t*) mpf_rhs, sense,(const mpf_t*) mpf_range, 0));
	/* set parameters */
	EGcallD(mpq_QSget_param(p, QS_PARAM_PRIMAL_PRICING, &objsense));
	EGcallD(mpf_QSset_param(p2, QS_PARAM_PRIMAL_PRICING, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_DUAL_PRICING, &objsense));
	EGcallD(mpf_QSset_param(p2, QS_PARAM_DUAL_PRICING, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_SIMPLEX_DISPLAY, &objsense));
	EGcallD(mpf_QSset_param(p2, QS_PARAM_SIMPLEX_DISPLAY, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_SIMPLEX_MAX_ITERATIONS, &objsense));
	EGcallD(mpf_QSset_param(p2, QS_PARAM_SIMPLEX_MAX_ITERATIONS, objsense));
	EGcallD(mpq_QSget_param(p, QS_PARAM_SIMPLEX_SCALING, &objsense));
	EGcallD(mpf_QSset_param(p2, QS_PARAM_SIMPLEX_SCALING, objsense));
	EGcallD(mpq_QSget_param_EGlpNum(p, QS_PARAM_SIMPLEX_MAX_TIME, &mpq_val));
	mpf_set_q(mpf_val,mpq_val);
	EGcallD(mpf_QSset_param_EGlpNum(p2, QS_PARAM_SIMPLEX_MAX_TIME, mpf_val));
	EGcallD(mpq_QSget_param_EGlpNum(p, QS_PARAM_OBJULIM, &mpq_val));
	mpf_set_q(mpf_val,mpq_val);
	EGcallD(mpf_QSset_param_EGlpNum(p2, QS_PARAM_OBJULIM, mpf_val));
	EGcallD(mpq_QSget_param_EGlpNum(p, QS_PARAM_OBJLLIM, &mpq_val));
	mpf_set_q(mpf_val,mpq_val);
	EGcallD(mpf_QSset_param_EGlpNum(p2, QS_PARAM_OBJLLIM, mpf_val));
	/* ending */
	CLEANUP:
	mpq_clear(mpq_val);
	mpf_clear(mpf_val);
	mpf_EGlpNumFreeArray(mpf_rowval);
	mpf_EGlpNumFreeArray(mpf_range);
	mpf_EGlpNumFreeArray(mpf_rhs);
	mpf_EGlpNumFreeArray(mpf_lb);
	mpf_EGlpNumFreeArray(mpf_ub);
	mpf_EGlpNumFreeArray(mpf_obj);
	mpq_EGlpNumFreeArray(mpq_rowval);
	mpq_EGlpNumFreeArray(mpq_range);
	mpq_EGlpNumFreeArray(mpq_rhs);
	mpq_EGlpNumFreeArray(mpq_lb);
	mpq_EGlpNumFreeArray(mpq_ub);
	mpq_EGlpNumFreeArray(mpq_obj);
	EGfree(rowcnt);
	EGfree(rowbeg);
	EGfree(rowind);
	EGfree(sense);
	if (rval)
	{
		mpf_QSfree_prob (p2);
		p2 = 0;
	}
#if QSEXACT_SAVE_INT
	else
	{
		mpf_QSwrite_prob (p2, "prob.mpf.lp", "LP");
	}
#endif
	return p2;
}

#if QSEXACT_SAVE_OPTIMAL
/* ========================================================================= */
/** @brief used to enumerate the generated optimal tests */
static int QSEXACT_SAVE_OPTIMAL_IND = 0;
#endif

/* ========================================================================= */
int QSexact_optimal_test (mpq_QSdata * p,
													mpq_t * p_sol,
													mpq_t * d_sol,
													QSbasis * basis)
{
	// clock start for timing purposes
        clock_t start = clock();

	/* local variables */
	register int i,
	  j;
	mpq_ILLlpdata *qslp = p->lp->O;
	int *iarr1 = 0,
	 *rowmap = qslp->rowmap,
	 *structmap = qslp->structmap,
	  col;
	mpq_t *arr1 = 0,
	 *arr2 = 0,
	 *arr3 = 0,
	 *arr4 = 0,
	 *rhs_copy = 0;
	mpq_t *dz = 0;
	int objsense = (qslp->objsense == QS_MIN) ? 1 : -1;
	int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0 : 100000 * (1 - p->simplex_display);
	int rval = 1;									/* store whether or not the solution is optimal, we start 
																 * assuming it is. */
	mpq_t num1,
	  num2,
	  num3,
	  p_obj,
	  d_obj;
	mpq_init (num1);
	mpq_init (num2);
	mpq_init (num3);
	mpq_init (p_obj);
	mpq_init (d_obj);
	mpq_set_ui (p_obj, 0UL, 1UL);
	mpq_set_ui (d_obj, 0UL, 1UL);

	/* now check if the given basis is the optimal basis */
	arr3 = qslp->lower;
	arr4 = qslp->upper;
	if (mpq_QSload_basis (p, basis))
	{
		rval = 0;
		MESSAGE (msg_lvl, "QSload_basis failed");
		goto CLEANUP;
	}
	for (i = basis->nstruct; i--;)
	{
		/* check that the upper and lower bound define a non-empty space */
		if (mpq_cmp (arr3[structmap[i]], arr4[structmap[i]]) > 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "variable %s has empty feasible range [%lg,%lg]",
								 qslp->colnames[i], mpq_EGlpNumToLf(arr3[structmap[i]]), 
								 mpq_EGlpNumToLf(arr4[structmap[i]]));
			}
			goto CLEANUP;
		}
		/* set the variable to its apropiate values, depending its status */
		switch (basis->cstat[i])
		{
		case QS_COL_BSTAT_FREE:
		case QS_COL_BSTAT_BASIC:
			if (mpq_cmp (p_sol[i], arr4[structmap[i]]) > 0)
				mpq_set (p_sol[i], arr4[structmap[i]]);
			else if (mpq_cmp (p_sol[i], arr3[structmap[i]]) < 0)
				mpq_set (p_sol[i], arr3[structmap[i]]);
			break;
		case QS_COL_BSTAT_UPPER:
			mpq_set (p_sol[i], arr4[structmap[i]]);
			break;
		case QS_COL_BSTAT_LOWER:
			mpq_set (p_sol[i], arr3[structmap[i]]);
			break;
		default:
			rval = 0;
			MESSAGE (msg_lvl, "Unknown Variable basic status %d, for variable "
							 "(%s,%d)", basis->cstat[i], qslp->colnames[i], i);
			goto CLEANUP;
			break;
		}
	}
	for (i = basis->nrows; i--;)
	{
		/* check that the upper and lower bound define a non-empty space */
		if (mpq_cmp (arr3[rowmap[i]], arr4[rowmap[i]]) > 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "constraint %s logical has empty feasible range "
								 "[%lg,%lg]", qslp->rownames[i], 
								 mpq_EGlpNumToLf(arr3[rowmap[i]]), 
								 mpq_EGlpNumToLf(arr4[rowmap[i]]));
			}
			goto CLEANUP;
		}
		/* set the variable to its apropiate values, depending its status */
		switch (basis->rstat[i])
		{
		case QS_ROW_BSTAT_BASIC:
			if (mpq_cmp (p_sol[i + basis->nstruct], arr4[rowmap[i]]) > 0)
				mpq_set (p_sol[i + basis->nstruct], arr4[rowmap[i]]);
			else if (mpq_cmp (p_sol[i + basis->nstruct], arr3[rowmap[i]]) < 0)
				mpq_set (p_sol[i + basis->nstruct], arr3[rowmap[i]]);
			break;
		case QS_ROW_BSTAT_UPPER:
			mpq_set (p_sol[i + basis->nstruct], arr4[rowmap[i]]);
			break;
		case QS_ROW_BSTAT_LOWER:
			mpq_set (p_sol[i + basis->nstruct], arr3[rowmap[i]]);
			break;
		default:
			rval = 0;
			MESSAGE (msg_lvl, "Unknown Variable basic status %d, for constraint "
							 "(%s,%d)", basis->cstat[i], qslp->rownames[i], i);
			goto CLEANUP;
			break;
		}
	}

	/* compute the actual RHS */
	rhs_copy = mpq_EGlpNumAllocArray (qslp->nrows);
	for (i = qslp->nstruct; i--;)
	{
		if (!mpq_equal (p_sol[i], mpq_zeroLpNum))
		{
			arr1 = qslp->A.matval + qslp->A.matbeg[structmap[i]];
			iarr1 = qslp->A.matind + qslp->A.matbeg[structmap[i]];
			for (j = qslp->A.matcnt[structmap[i]]; j--;)
			{
				mpq_mul (num1, arr1[j], p_sol[i]);
				mpq_add (rhs_copy[iarr1[j]], rhs_copy[iarr1[j]], num1);
			}
		}
	}

	/* now check if both rhs and copy_rhs are equal */
	arr4 = qslp->upper;
	arr1 = qslp->rhs;
	arr2 = qslp->lower;
	for (i = qslp->nrows; i--;)
	{
		mpq_mul (num1, arr1[i], d_sol[i]);
		mpq_add (d_obj, d_obj, num1);
		mpq_sub (num2, arr1[i], rhs_copy[i]);
		EXIT (qslp->A.matcnt[rowmap[i]] != 1, "Imposible!");
		if (basis->rstat[i] == QS_ROW_BSTAT_BASIC)
			mpq_div (p_sol[qslp->nstruct + i], num2,
							 qslp->A.matval[qslp->A.matbeg[rowmap[i]]]);
		else
		{
			mpq_mul (num1, p_sol[qslp->nstruct + i],
							 qslp->A.matval[qslp->A.matbeg[rowmap[i]]]);
			if (!mpq_equal (num1, num2))
			{
				rval = 0;
				if(!msg_lvl)
				{
					MESSAGE(0, "solution is infeasible for constraint %s, violation"
								 " %lg", qslp->rownames[i],
								 mpq_get_d (num1) - mpq_get_d (num2));
				}
				goto CLEANUP;
			}
		}
		mpq_set (num2, p_sol[qslp->nstruct + i]);
		/* now we check the bounds on the logical variables */
		if (mpq_cmp (num2, arr2[rowmap[i]]) < 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "constraint %s artificial (%lg) bellow lower"
								 " bound (%lg), actual LHS (%lg), actual RHS (%lg)",
								 qslp->rownames[i], mpq_get_d (num2), 
								 mpq_get_d (arr2[rowmap[i]]), mpq_get_d (rhs_copy[i]), 
								 mpq_get_d (arr1[i]));
			}
			goto CLEANUP;
		}
		else if (mpq_cmp (num2, arr4[rowmap[i]]) > 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "constraint %s artificial (%lg) bellow lower bound"
								 " (%lg)", qslp->rownames[i], mpq_get_d (num2),
								 mpq_get_d (arr4[rowmap[i]]));
			}
			goto CLEANUP;
		}
	}

	/* compute the upper and lower bound dual variables, note that dl is the dual
	 * of the lower bounds, and du the dual of the upper bound, dl >= 0 and du <=
	 * 0 and A^t y + Idl + Idu = c, and the dual objective value is 
	 * max y*b + l*dl + u*du, we colapse both vector dl and du into dz, note that
	 * if we are maximizing, then dl <= 0 and du >=0 */
	dz = mpq_EGlpNumAllocArray (qslp->ncols);
	arr2 = qslp->obj;
	arr3 = qslp->lower;
	arr4 = qslp->upper;
	for (i = qslp->nstruct; i--;)
	{
		col = structmap[i];
		mpq_mul (num1, arr2[col], p_sol[i]);
		mpq_add (p_obj, p_obj, num1);
		arr1 = qslp->A.matval + qslp->A.matbeg[col];
		iarr1 = qslp->A.matind + qslp->A.matbeg[col];
		mpq_set (num1, arr2[col]);
		for (j = qslp->A.matcnt[col]; j--;)
		{
			mpq_mul (num2, arr1[j], d_sol[iarr1[j]]);
			mpq_sub (num1, num1, num2);
		}
		mpq_set (dz[col], num1);
		/* objective update */
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0)
		{
			mpq_mul (num3, dz[col], arr3[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		else
		{
			mpq_mul (num3, dz[col], arr4[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		/* now we check that only when the logical is tight then the dual
		 * variable may be non-zero, also check for primal feasibility with respect
		 * to lower/upper bounds. */
		mpq_set_ui (num2, 0UL, 1UL);
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0)
		{
			mpq_sub (num1, p_sol[i], arr3[col]);
			mpq_mul (num2, num1, dz[col]);
		}
		if (mpq_cmp_ui (num2, 0UL, 1UL) != 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "lower bound (%s,%d) slack (%lg) and dual variable (%lg)"
								 " don't satisfy complementary slacknes %s",
								 qslp->colnames[i], i, mpq_get_d(num1), mpq_get_d(dz[col]),
								 "(real)");
			}
			goto CLEANUP;
		}
		mpq_set_ui (num2, 0UL, 1UL);
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) < 0)
		{
			mpq_sub (num1, p_sol[i], arr4[col]);
			mpq_mul (num2, num1, dz[col]);
		}
		if (mpq_cmp_ui (num2, 0UL, 1UL) != 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "upper bound (%lg) variable (%lg) and dual variable"
									" (%lg) don't satisfy complementary slacknes for variable"
									" (%s,%d) %s", mpq_get_d(arr4[col]), mpq_get_d(p_sol[i]),
									mpq_get_d(dz[col]), qslp->colnames[i], i, "(real)");
			}
			goto CLEANUP;
		}
	}
	/* complenetary slackness checked, now update the same for the logical
	 * variables */
	for (i = qslp->nrows; i--;)
	{
		col = rowmap[i];
		mpq_mul (num1, arr2[col], p_sol[i + qslp->nstruct]);
		WARNING (mpq_cmp (arr2[col], mpq_zeroLpNum), "logical variable %s with "
						 "non-zero objective function %lf", qslp->rownames[i],
						 mpq_get_d (arr2[col]));
		mpq_add (p_obj, p_obj, num1);
		arr1 = qslp->A.matval + qslp->A.matbeg[col];
		iarr1 = qslp->A.matind + qslp->A.matbeg[col];
		mpq_set (num1, arr2[col]);
		for (j = qslp->A.matcnt[col]; j--;)
		{
			mpq_mul (num2, arr1[j], d_sol[iarr1[j]]);
			mpq_sub (num1, num1, num2);
		}
		mpq_set (dz[col], num1);
		/* objective update */
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0)
		{
			mpq_mul (num3, dz[col], arr3[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		else
		{
			mpq_mul (num3, dz[col], arr4[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		/* now we check that only when the primal variable is tight then the dual
		 * variable may be non-zero, also check for primal feasibility with respect
		 * to lower/upper bounds. */
		mpq_set_ui (num2, 0UL, 1UL);
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0)
		{
			mpq_sub (num1, p_sol[i + qslp->nstruct], arr3[col]);
			mpq_mul (num2, num1, dz[col]);
		}
		if (mpq_cmp_ui (num2, 0UL, 1UL) != 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "lower bound (%s,%d) slack (%lg) and dual variable (%lg)"
								 " don't satisfy complementary slacknes %s", 
								 qslp->colnames[col], i, mpq_get_d(num1), mpq_get_d(dz[col]), 
								 "(real)");
			}
			goto CLEANUP;
		}
		mpq_set_ui (num2, 0UL, 1UL);
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) < 0)
		{
			mpq_sub (num1, p_sol[i + qslp->nstruct], arr4[col]);
			mpq_mul (num2, num1, dz[col]);
		}
		if (mpq_cmp_ui (num2, 0UL, 1UL) != 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "upper bound (%lg) variable (%lg) and dual variable"
								" (%lg) don't satisfy complementary slacknes for variable "
								"(%s,%d) %s", mpq_get_d(arr4[col]), 
								mpq_get_d(p_sol[i+qslp->nstruct]), mpq_get_d(dz[col]), qslp->colnames[col], i,
								"(real)");
			}
			goto CLEANUP;
		}
	}

	/* now check the objective values */
	if (mpq_cmp (p_obj, d_obj) != 0)
	{
		rval = 0;
		if(!msg_lvl)
		{
			MESSAGE(0, "primal and dual objective value differ %lg %lg",
							 mpq_get_d(p_obj), mpq_get_d(d_obj));
		}
		goto CLEANUP;
	}
	/* now we report optimality */
	if(!msg_lvl)
	{
		MESSAGE(0, "Problem solved to optimality, LP value %lg", mpq_get_d(p_obj));
	}
	/* now we load into cache the solution */
	if (!p->cache)
	{
		p->cache = EGsMalloc (mpq_ILLlp_cache, 1);
		mpq_EGlpNumInitVar (p->cache->val);
		mpq_ILLlp_cache_init (p->cache);
	}
	if (qslp->nrows != p->cache->nrows || qslp->nstruct != p->cache->nstruct)
	{
		mpq_ILLlp_cache_free (p->cache);
		EGcallD(mpq_ILLlp_cache_alloc (p->cache, qslp->nstruct, qslp->nrows));
	}
	p->cache->status = QS_LP_OPTIMAL;
	p->qstatus = QS_LP_OPTIMAL;
	p->lp->basisstat.optimal = 1;
	mpq_set (p->cache->val, p_obj);
	for (i = qslp->nstruct; i--;)
	{
		mpq_set (p->cache->x[i], p_sol[i]);
		mpq_set (p->cache->rc[i], dz[structmap[i]]);
	}
	for (i = qslp->nrows; i--;)
	{
		mpq_set (p->cache->slack[i], p_sol[i + qslp->nstruct]);
		mpq_set (p->cache->pi[i], d_sol[i]);
	}

	/* save the problem and solution if enablred */
#if QSEXACT_SAVE_OPTIMAL
	{
		char stmp[1024];
		EGioFile_t *out_f = 0;
		snprintf (stmp, 1023, "%s-opt%03d.lp", p->name ? p->name : "UNNAMED",
							QSEXACT_SAVE_OPTIMAL_IND);
		if (mpq_QSwrite_prob (p, stmp, "LP"))
		{
			rval = 0;
			MESSAGE (0, "Couldn't write output problem %s", stmp);
			goto CLEANUP;
		}
		snprintf (stmp, 1023, "%s-opt%03d.sol.gz", p->name ? p->name : "UNNAMED",
							QSEXACT_SAVE_OPTIMAL_IND);
		if (!(out_f = EGioOpen (stmp, "w+")))
		{
			rval = 0;
			MESSAGE (0, "Couldn't open solution file %s", stmp);
			goto CLEANUP;
		}
		if (QSexact_print_sol (p, out_f))
		{
			rval = 0;
			MESSAGE (0, "Couldn't write output solution %s", stmp);
			goto CLEANUP;
		}
		EGioClose (out_f);
		QSEXACT_SAVE_OPTIMAL_IND++;
	}
#endif
	rval = 1;

	/* ending */
CLEANUP:
	// testing for log file
        clock_t end = clock();
        double duration = (double)(end - start) / CLOCKS_PER_SEC;
        log_timing("QSexact_optimal_test took ", duration);


	mpq_EGlpNumFreeArray (dz);
	mpq_EGlpNumFreeArray (rhs_copy);
	mpq_clear (num1);
	mpq_clear (num2);
	mpq_clear (num3);
	mpq_clear (p_obj);
	mpq_clear (d_obj);
	return rval;
}

/* ========================================================================= */
int QSexact_infeasible_test (mpq_QSdata * p, mpq_t * d_sol)
{
	// clock start for timing purposes
        clock_t start = clock();

	/* local variables */
	register int i,
	  j;
	int *iarr1;
	mpq_ILLlpdata *qslp = p->lp->O;
	mpq_t *arr1,
	 *arr2,
	 *arr3,
	 *arr4;
	mpq_t *dl = 0,
	 *du = 0;
	int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0 : 100000 * (1 - p->simplex_display);
	int rval = 1;									/* store whether or not the solution is optimal, we start 
																 * assuming it is. */
	mpq_t num1,
	  num2,
	  num3,
	  d_obj;
	mpq_init (num1);
	mpq_init (num2);
	mpq_init (num3);
	mpq_init (d_obj);
	mpq_set_ui (d_obj, 0UL, 1UL);

	/* compute the dual objective value */
	arr2 = qslp->rhs;
	for (i = qslp->nrows; i--;)
	{
		mpq_mul (num1, arr2[i], d_sol[i]);
		mpq_add (d_obj, d_obj, num1);
	}

	/* compute the upper and lower bound dual variables, note that dl is the dual
	 * of the lower bounds, and du the dual of the upper bound, dl <= 0 and du >=
	 * 0 and A^t y + Idl + Idu = c, and the dual objective value is 
	 * max y*b + l*dl + u*du */
	du = mpq_EGlpNumAllocArray (qslp->ncols);
	dl = mpq_EGlpNumAllocArray (qslp->ncols);
	arr3 = qslp->lower;
	arr4 = qslp->upper;
	for (i = qslp->ncols; i--;)
	{
		arr1 = qslp->A.matval + qslp->A.matbeg[i];
		iarr1 = qslp->A.matind + qslp->A.matbeg[i];
		mpq_set_ui (num1, 0UL, 1UL);
		mpq_set_ui (du[i], 0UL, 1UL);
		mpq_set_ui (dl[i], 0UL, 1UL);
		for (j = qslp->A.matcnt[i]; j--;)
		{
			mpq_mul (num2, arr1[j], d_sol[iarr1[j]]);
			mpq_sub (num1, num1, num2);
		}
		if (mpq_cmp_ui (num1, 0UL, 1UL) < 0)
			mpq_set (du[i], num1);
		else
			mpq_set (dl[i], num1);
		if (mpq_equal (arr4[i], mpq_ILL_MAXDOUBLE) &&
				mpq_cmp_ui (du[i], 0UL, 1UL) != 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "upper bound of variable is INFTY, and it's dual is "
								 "non-zero %lg", mpq_get_d(du[i]));
			}
			goto CLEANUP;
		}
		if (mpq_equal (arr3[i], mpq_ILL_MINDOUBLE) &&
				mpq_cmp_ui (dl[i], 0UL, 1UL) != 0)
		{
			rval = 0;
			if(!msg_lvl)
			{
				MESSAGE(0, "lower bound of variable is -INFTY, and it's dual is "
								 "non-zero %lg", mpq_get_d(dl[i]));
			}
			goto CLEANUP;
		}
		mpq_mul (num3, dl[i], arr3[i]);
		mpq_add (d_obj, d_obj, num3);
		mpq_mul (num3, du[i], arr4[i]);
		mpq_add (d_obj, d_obj, num3);
	}
	/* now check the objective values */
	if (mpq_cmp_ui (d_obj, 0UL, 1UL) <= 0)
	{
		rval = 0;
		if(!msg_lvl)
		{
			MESSAGE(0, "dual ray is feasible, but objective is non "
							"positive %lg", mpq_get_d(d_obj));
		}
		goto CLEANUP;
	}
	p->qstatus = QS_LP_INFEASIBLE;

	/* ending */
CLEANUP:
	// testing for log file
        clock_t end = clock();
        double duration = (double)(end - start) / CLOCKS_PER_SEC;
        log_timing("QSexact_infeasible_test took ", duration);

	mpq_EGlpNumFreeArray (dl);
	mpq_EGlpNumFreeArray (du);
	mpq_clear (num1);
	mpq_clear (num2);
	mpq_clear (num3);
	mpq_clear (d_obj);
	return rval;
}

/* ========================================================================= */
/** @brief Used as separator while printing output to the screen (controled by
 * enabling simplex_display in the mpq_QSdata */
/* ========================================================================= */
static const char __sp[81] =
	"================================================================================";

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully prove infeasibility, and save (if y is non
 * NULL ) the dual ray solution provided in y_mpq.
 * @param p_mpq the problem data.
 * @param y where to store the optimal dual solution (if not null).
 * @param y_mpq  the optimal dual solution.
 * */
/* ========================================================================= */
static void infeasible_output (mpq_QSdata * p_mpq,
															 mpq_t * const y,
															 mpq_t * y_mpq)
{
	if (p_mpq->simplex_display)
	{
		QSlog("Problem Is Infeasible");
	}
	if (y)
	{
		unsigned sz = __EGlpNumArraySize (y_mpq);
		while (sz--)
			mpq_set (y[sz], y_mpq[sz]);
	}
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully solved the problem at optimality, and save (if x and y are non
 * NULL respectivelly) the optimal primal/dual solution provided in x_mpq and
 * y_mpq. 
 * @param p_mpq the problem data.
 * @param x where to store the optimal primal solution (if not null).
 * @param y where to store the optimal dual solution (if not null).
 * @param x_mpq  the optimal primal solution.
 * @param y_mpq  the optimal dual solution.
 * */
/* ========================================================================= */
static void optimal_output (mpq_QSdata * p_mpq,
														mpq_t * const x,
														mpq_t * const y,
														mpq_t * x_mpq,
														mpq_t * y_mpq)
{
	if (p_mpq->simplex_display)
	{
		QSlog("Problem Solved Exactly");
	}
	if (y)
	{
		unsigned sz = __EGlpNumArraySize (y_mpq);
		while (sz--)
			mpq_set (y[sz], y_mpq[sz]);
	}
	if (x)
	{
		unsigned sz = __EGlpNumArraySize (x_mpq);
		while (sz--)
			mpq_set (x[sz], x_mpq[sz]);
	}
}

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
static int QSexact_basis_status (mpq_QSdata * p_mpq,
																 int *status,
																 QSbasis * const basis,
																 const int msg_lvl,
																 int *const simplexalgo)
{
	// clock start for timing purposes
        clock_t start = clock();

	int rval = 0,
	singular;
	mpq_feas_info fi;
	EGtimer_t local_timer;
	mpq_EGlpNumInitVar (fi.totinfeas);
	EGtimerReset (&local_timer);
	EGtimerStart (&local_timer);
	// load and check the basis
	EGcallD(mpq_QSload_basis (p_mpq, basis));
	if (p_mpq->cache) 
	{
		mpq_ILLlp_cache_free (p_mpq->cache);
		mpq_clear (p_mpq->cache->val);
		ILL_IFFREE(p_mpq->cache);
	}
	p_mpq->qstatus = QS_LP_MODIFIED;
	if(p_mpq->qslp->sinfo) 
	{
		mpq_ILLlp_sinfo_free(p_mpq->qslp->sinfo);
		ILL_IFFREE(p_mpq->qslp->sinfo);
	}
	if(p_mpq->qslp->rA)
	{
		mpq_ILLlp_rows_clear (p_mpq->qslp->rA);
		ILL_IFFREE(p_mpq->qslp->rA);
	}
	// rebuild internal lp 
	mpq_free_internal_lpinfo (p_mpq->lp);
	mpq_init_internal_lpinfo (p_mpq->lp);
	EGcallD(mpq_build_internal_lpinfo (p_mpq->lp));
	mpq_ILLfct_set_variable_type (p_mpq->lp);
	// factoring of basis
	EGcallD(mpq_ILLbasis_load (p_mpq->lp, p_mpq->basis, p_mpq->cached_baz));
	if (p_mpq->cached_lu == 0) 
	{
		EGcallD(mpq_ILLbasis_factor (p_mpq->lp, &singular));
		ILL_SAFE_MALLOC (p_mpq->cached_lu, 1, mpq_factor_work);
		mpq_EGlpNumInitVar (p_mpq->cached_lu->fzero_tol);
		mpq_EGlpNumInitVar (p_mpq->cached_lu->szero_tol); 
		mpq_EGlpNumInitVar (p_mpq->cached_lu->partial_tol);
		mpq_EGlpNumInitVar (p_mpq->cached_lu->maxelem_orig);
		mpq_EGlpNumInitVar (p_mpq->cached_lu->maxelem_factor);
		mpq_EGlpNumInitVar (p_mpq->cached_lu->maxelem_cur);
		mpq_EGlpNumInitVar (p_mpq->cached_lu->partial_cur);
		rval = mpq_ILLfactor_deep_copy(p_mpq->cached_lu, p_mpq->lp->f);
		if (rval) {
			QSlog("Failed to deep copy factor work");
			goto CLEANUP;
		}
		//cache the basis deep copy
		ILL_SAFE_MALLOC (p_mpq->cached_baz, p_mpq->lp->O->nrows, int);
		for (int i = 0; i < p_mpq->lp->O->nrows; ++i) {
			p_mpq->cached_baz[i] = p_mpq->lp->baz[i];
		}
	}
	else {
		int refactor = 0;
		// Collect all mismatch indices up front
		int* mismatch_indices = NULL;
		ILL_SAFE_MALLOC(mismatch_indices, p_mpq->lp->O->nrows, int);
		int mismatch_count = 0;
		for (int i = 0; i < p_mpq->lp->O->nrows; ++i) {
			if (p_mpq->cached_baz[i] != p_mpq->lp->baz[i]) {
				mismatch_indices[mismatch_count++] = i;
			}
		}
		if ((double)mismatch_count / p_mpq->lp->O->nrows > 0.05 ) {
			QSlog("Using refactorization");
			refactor = 1;
			free(mismatch_indices);
		}
		log_message("Mismatch's: %d/%d", mismatch_count, p_mpq->lp->O->nrows);
		while (!refactor) {
			int changed = 0;
			int update_pos = -1;
			int leaving_col = -1;

			// Find the next column that has changed
			for (int i = 0; i < mismatch_count; ++i) {
				if (p_mpq->cached_baz[mismatch_indices[i]] != p_mpq->lp->baz[mismatch_indices[i]]) {
					update_pos = mismatch_indices[i];
					leaving_col = p_mpq->cached_baz[update_pos];
					changed = 1;
					break;
				}
			}

			if (!changed) {
				break; // Bases are now in sync
			}
			else {
				int entering_col = p_mpq->lp->baz[update_pos];

				// Create a_s for the entering column
				mpq_svector a_s;
				a_s.nzcnt = p_mpq->lp->matcnt[entering_col];
				a_s.indx = &(p_mpq->lp->matind[p_mpq->lp->matbeg[entering_col]]);
				a_s.coef = &(p_mpq->lp->matval[p_mpq->lp->matbeg[entering_col]]);
		
				// Allocate and compute the spike vector
				mpq_svector spike, direction;
				mpq_ILLsvector_alloc(&spike, p_mpq->lp->nrows);
				mpq_ILLsvector_alloc(&direction, p_mpq->lp->nrows);
				
				// Use cached MPF 128-bit LU for faster and more accurate direction computation
				unsigned original_precision = EGLPNUM_PRECISION;
				QSexact_set_precision(128);  
				
				mpf_factor_work *mpf_cached_lu = NULL;
				mpf_svector mpf_a_s, mpf_spike, mpf_direction;
				
				ILL_SAFE_MALLOC(mpf_cached_lu, 1, mpf_factor_work);
				rval = mpq_factor_work_to_mpf_factor_work(mpf_cached_lu, p_mpq->cached_lu);
				if (rval) {
					QSlog("Failed to convert mpq_factor_work to mpf_factor_work");
					QSexact_set_precision(original_precision);
					refactor = 1;
					break;
				}	
				
				// Convert a_s to mpf precision
				mpf_ILLsvector_alloc(&mpf_a_s, p_mpq->lp->nrows);
				mpf_ILLsvector_alloc(&mpf_spike, p_mpq->lp->nrows);
				mpf_ILLsvector_alloc(&mpf_direction, p_mpq->lp->nrows);
				
				mpf_a_s.nzcnt = a_s.nzcnt;
				memcpy(mpf_a_s.indx, a_s.indx, a_s.nzcnt * sizeof(int));
				for (int j = 0; j < a_s.nzcnt; j++) {
					mpf_set_q(mpf_a_s.coef[j], a_s.coef[j]);
				}
				
				mpf_ILLfactor_ftran_update(mpf_cached_lu, &mpf_a_s, &mpf_spike, &mpf_direction);
				direction.nzcnt = mpf_direction.nzcnt;
				memcpy(direction.indx, mpf_direction.indx, mpf_direction.nzcnt * sizeof(int));
				for (int j = 0; j < mpf_direction.nzcnt; j++) {
					mpq_set_f(direction.coef[j], mpf_direction.coef[j]);
				}
				mpq_compute_spike(p_mpq->cached_lu, &a_s, &spike);
				
				// Clean up mpf precision structures
				mpf_ILLsvector_free(&mpf_a_s);
				mpf_ILLsvector_free(&mpf_spike);
				mpf_ILLsvector_free(&mpf_direction);
				mpf_ILLfactor_free_factor_work(mpf_cached_lu);
				ILL_IFFREE(mpf_cached_lu);
				
				// Restore original precision
				QSexact_set_precision(original_precision);

				// Look through mismatches for direction vector entry != 0
				int swap_pos = -1;
				char* is_mismatch = calloc(p_mpq->lp->O->nrows, sizeof(char));
				for (int j = 0; j < mismatch_count; ++j) {
					int pos = mismatch_indices[j];
					if (p_mpq->cached_baz[pos] != p_mpq->lp->baz[pos]) {
						is_mismatch[pos] = 1;
					}
				}
				
				// Find the mismatch position with maximum absolute value in direction vector
				double max_abs_val = 0.0;
				for (int k = 0; k < direction.nzcnt; ++k) {
					int pos = direction.indx[k];
					if (is_mismatch[pos]) {
						double abs_val = fabs(mpq_get_d(direction.coef[k]));
						if (abs_val > max_abs_val) {
							max_abs_val = abs_val;
								swap_pos = pos;
						}
					}
				}
				// If we found a position with direction > 0 within cached basis, we need to update that position in the current basis
				if (swap_pos != -1 && swap_pos != update_pos) {
				
					int temp_col = p_mpq->lp->baz[update_pos];
					p_mpq->lp->baz[update_pos] = p_mpq->lp->baz[swap_pos];
					p_mpq->lp->baz[swap_pos] = temp_col;
					update_pos = swap_pos;

				}
				if (swap_pos == -1) {
					QSlog("No swap found, increase copy precision");
					refactor = 1;
					break;
				}
				rval = mpq_ILLfactor_update(p_mpq->cached_lu, &spike, update_pos, &refactor);

				if (refactor || rval) {
					QSlog("LU update at position %d triggered refactorization (refactor=%d, rval=%d)\n", update_pos, refactor, rval);
				refactor = 1; 
				break; 
				}
				
				// Update was successful, update the cached basis for this position
				p_mpq->cached_baz[update_pos] = entering_col;

				mpq_ILLsvector_free(&spike);
				mpq_ILLsvector_free(&direction);
				free(is_mismatch);
				
				// Remove the mismatch from the list
				for (int idx = 0; idx < mismatch_count; ++idx) {
					if (mismatch_indices[idx] == update_pos) {
						/* swap current entry with the last one and shrink */
						mismatch_indices[idx] = mismatch_indices[mismatch_count - 1];
						mismatch_count--;
						break;
					}
				}

				// If no mismatches are left, we can stop looping
				if (mismatch_count == 0) {
					break;
				}
				
			}		
		}
		if (!refactor) {
			mpq_factor_work *temp_lu;
			ILL_SAFE_MALLOC (temp_lu, 1, mpq_factor_work);	
			mpq_EGlpNumInitVar (temp_lu->fzero_tol);
			mpq_EGlpNumInitVar (temp_lu->szero_tol); 
			mpq_EGlpNumInitVar (temp_lu->partial_tol);
			mpq_EGlpNumInitVar (temp_lu->maxelem_orig);
			mpq_EGlpNumInitVar (temp_lu->maxelem_factor);
			mpq_EGlpNumInitVar (temp_lu->maxelem_cur);
			mpq_EGlpNumInitVar (temp_lu->partial_cur);
			rval = mpq_ILLfactor_deep_copy(temp_lu, p_mpq->cached_lu);
			if (rval) {
				QSlog("Failed to deep copy factor work after refactorization");
				goto CLEANUP;
			}
			if (p_mpq->lp->f) {
				mpq_ILLfactor_free_factor_work(p_mpq->lp->f);
				ILL_IFFREE(p_mpq->lp->f);
			}
			p_mpq->lp->f = temp_lu;
			QSlog("Updated cached lu");
			free(mismatch_indices);	
		}
		if (refactor) {
		    int singular;
   			// Perform full refactorization
    		EGcallD(mpq_ILLbasis_factor(p_mpq->lp, &singular));
    		// Deep copy the new LU factorization into the cache
    		rval = mpq_ILLfactor_deep_copy(p_mpq->cached_lu, p_mpq->lp->f);
    		// Update cached_baz to match current baz
    		for (int i = 0; i < p_mpq->lp->O->nrows; ++i) {
        		p_mpq->cached_baz[i] = p_mpq->lp->baz[i];
    		}

		}
	}
	memset (&(p_mpq->lp->basisstat), 0, sizeof (mpq_lp_status_info));
	// feasibility check
	mpq_ILLfct_compute_piz (p_mpq->lp);
	mpq_ILLfct_compute_dz (p_mpq->lp);
	mpq_ILLfct_compute_xbz (p_mpq->lp);
	// primal and dual solution check
	mpq_ILLfct_check_pfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
	mpq_ILLfct_check_dfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
	mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus, PHASEII,
			PHASEII);
	if (p_mpq->lp->basisstat.optimal)
	{
		*status = QS_LP_OPTIMAL;
		EGcallD(mpq_QSgrab_cache (p_mpq, QS_LP_OPTIMAL));
	}
	else if (p_mpq->lp->basisstat.primal_infeasible
			|| p_mpq->lp->basisstat.dual_unbounded)
	{
		if (*status == QS_LP_INFEASIBLE)
			*simplexalgo = PRIMAL_SIMPLEX;
		*status = QS_LP_INFEASIBLE;
		p_mpq->lp->final_phase = PRIMAL_PHASEI;
		p_mpq->lp->pIpiz = mpq_EGlpNumAllocArray (p_mpq->lp->nrows);
		mpq_ILLfct_compute_phaseI_piz (p_mpq->lp);
	}
	else if (p_mpq->lp->basisstat.primal_unbounded)
		*status = QS_LP_UNBOUNDED;
	else
		*status = QS_LP_UNSOLVED;
	EGtimerStop (&local_timer);
	if(!msg_lvl)
	{
		MESSAGE(0, "Performing Rational Basic Solve on %s, %s, check"
				" done in %lg seconds, PS %s %lg, DS %s %lg", p_mpq->name, 
				(*status == QS_LP_OPTIMAL) ? "RAT_optimal" : 
				((*status == QS_LP_INFEASIBLE) ?  "RAT_infeasible" : 
				 ((*status == QS_LP_UNBOUNDED) ?  "RAT_unbounded" : "RAT_unsolved")),
				local_timer.time, p_mpq->lp->basisstat.primal_feasible ? 
				"F":(p_mpq->lp->basisstat.primal_infeasible ? "I" : "U"), 
				p_mpq->lp->basisstat.primal_feasible ?
				mpq_get_d(p_mpq->lp->objval) : 
				(p_mpq->lp->basisstat.primal_infeasible ?
				 mpq_get_d(p_mpq->lp->pinfeas) : mpq_get_d(p_mpq->lp->objbound)), 
				p_mpq->lp->basisstat.dual_feasible ? 
				"F":(p_mpq->lp->basisstat.dual_infeasible ? "I" : "U"), 
				p_mpq->lp->basisstat.dual_feasible ? mpq_get_d(p_mpq->lp->dobjval) 
				:(p_mpq->lp->basisstat.dual_infeasible ? 
					mpq_get_d(p_mpq->lp->dinfeas) : mpq_get_d(p_mpq->lp->objbound)) );
	}
CLEANUP:
	// testing for log file
        clock_t end = clock();
        double duration = (double)(end - start) / CLOCKS_PER_SEC;
        log_timing("QSexact_basis_status took ", duration);

	mpq_EGlpNumClearVar (fi.totinfeas);
	return rval;
}

/* ========================================================================= */
/** @brief test whether given basis is primal and dual feasible in rational arithmetic. 
 * @param p_mpq   the problem data.
 * @param basis   basis to be tested.
 * @param result  where to store whether given basis is primal and dual feasible.
 * @param msg_lvl message level.
 */
int QSexact_basis_optimalstatus(
   mpq_QSdata * p_mpq,
   QSbasis* basis,
   char* result,
   const int msg_lvl
   )
{
   int rval = 0,
      singular;
   mpq_feas_info fi;
   EGtimer_t local_timer;

   /* test primal and dual feasibility of basic solution */
   mpq_EGlpNumInitVar (fi.totinfeas);

   EGtimerReset (&local_timer);
   EGtimerStart (&local_timer);

   EGcallD(mpq_QSload_basis (p_mpq, basis));

   if (p_mpq->cache) 
   {
      mpq_ILLlp_cache_free (p_mpq->cache);
      mpq_clear (p_mpq->cache->val);
      ILL_IFFREE(p_mpq->cache);
   }
   p_mpq->qstatus = QS_LP_MODIFIED;

   if(p_mpq->qslp->sinfo) 
   {
      mpq_ILLlp_sinfo_free(p_mpq->qslp->sinfo);
      ILL_IFFREE(p_mpq->qslp->sinfo);
   }
   if(p_mpq->qslp->rA)
   {
      mpq_ILLlp_rows_clear (p_mpq->qslp->rA);
      ILL_IFFREE(p_mpq->qslp->rA);
   }

   mpq_free_internal_lpinfo (p_mpq->lp);
   mpq_init_internal_lpinfo (p_mpq->lp);
   EGcallD(mpq_build_internal_lpinfo (p_mpq->lp));

   mpq_ILLfct_set_variable_type (p_mpq->lp);

   EGcallD(mpq_ILLbasis_load (p_mpq->lp, p_mpq->basis, p_mpq->cached_baz));
   EGcallD(mpq_ILLbasis_factor (p_mpq->lp, &singular));

   memset (&(p_mpq->lp->basisstat), 0, sizeof (mpq_lp_status_info));
   mpq_ILLfct_compute_piz (p_mpq->lp); 
   mpq_ILLfct_compute_dz (p_mpq->lp);
   mpq_ILLfct_compute_xbz (p_mpq->lp);
   mpq_ILLfct_check_pfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
   mpq_ILLfct_check_dfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
   mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);
   
   if( p_mpq->lp->basisstat.optimal )
   {
      *result = 1;
   }
   else
   {
      *result = 0;
   }

   EGtimerStop (&local_timer);

   if( !msg_lvl )
   {
      MESSAGE(0, "Performing rational solution check for accuratelp on %s, sucess=%s", 
         p_mpq->name, 
         *result ? "YES" : "NO");
   }
   
 CLEANUP:
   mpq_EGlpNumClearVar (fi.totinfeas);
   return rval;
}

/* ========================================================================= */
/** @brief test whether given basis is dual feasible in rational arithmetic. 
 * @param p_mpq   the problem data.
 * @param basis   basis to be tested.
 * @param result  where to store whether given basis is dual feasible.
 * @param dobjval where to store dual solution value in case of dual feasibility (if not NULL).
 * @param msg_lvl message level.
 */
int QSexact_basis_dualstatus(
   mpq_QSdata * p_mpq,
   QSbasis* basis,
   char* result,
   mpq_t* dobjval,
   const int msg_lvl
   )
{
	int rval = 0,
	singular;
	mpq_feas_info fi;
	EGtimer_t local_timer;

	mpq_EGlpNumInitVar (fi.totinfeas);
	EGtimerReset (&local_timer);
	EGtimerStart (&local_timer);
	EGcallD(mpq_QSload_basis (p_mpq, basis));
	if (p_mpq->cache) 
	{
		mpq_ILLlp_cache_free (p_mpq->cache);
		mpq_clear (p_mpq->cache->val);
		ILL_IFFREE(p_mpq->cache);
	}
	p_mpq->qstatus = QS_LP_MODIFIED;

	if(p_mpq->qslp->sinfo) 
	{
		mpq_ILLlp_sinfo_free(p_mpq->qslp->sinfo);
		ILL_IFFREE(p_mpq->qslp->sinfo);
	}

	if(p_mpq->qslp->rA)
	{
		mpq_ILLlp_rows_clear (p_mpq->qslp->rA);
		ILL_IFFREE(p_mpq->qslp->rA);
	}

	mpq_free_internal_lpinfo (p_mpq->lp);
	mpq_init_internal_lpinfo (p_mpq->lp);
	EGcallD(mpq_build_internal_lpinfo (p_mpq->lp));
	mpq_ILLfct_set_variable_type (p_mpq->lp);
	EGcallD(mpq_ILLbasis_load (p_mpq->lp, p_mpq->basis, p_mpq->cached_baz));
	EGcallD(mpq_ILLbasis_factor (p_mpq->lp, &singular));

	memset (&(p_mpq->lp->basisstat), 0, sizeof (mpq_lp_status_info));
	mpq_ILLfct_compute_piz (p_mpq->lp); 
	mpq_ILLfct_compute_dz (p_mpq->lp);
	mpq_ILLfct_compute_dobj(p_mpq->lp); 
	mpq_ILLfct_check_dfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
	mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

	if( p_mpq->lp->basisstat.dual_feasible )
	{
		*result = 1;
		if( dobjval )
		{
			mpq_EGlpNumCopy(*dobjval, p_mpq->lp->dobjval);
		}
	}
	else if( p_mpq->lp->basisstat.dual_infeasible )
	{
		*result = 0;
	}
	else 
	{
		TESTG((rval=!(p_mpq->lp->basisstat.dual_unbounded)), CLEANUP, "Internal BUG, problem should be dual unbounded but is not");
		*result = 1;
		if( dobjval )
		{
			mpq_EGlpNumCopy(*dobjval, p_mpq->lp->objbound);
		}
	}

	EGtimerStop (&local_timer);

	if(!msg_lvl)
	{
		MESSAGE(0, "Performing Rational Basic Test on %s, check done in %lg seconds, DS %s %lg", 
				p_mpq->name, local_timer.time, 
				p_mpq->lp->basisstat.dual_feasible ? "F": (p_mpq->lp->basisstat.dual_infeasible ? "I" : "U"), 
				p_mpq->lp->basisstat.dual_feasible ? mpq_get_d(p_mpq->lp->dobjval) : (p_mpq->lp->basisstat.dual_infeasible ? mpq_get_d(p_mpq->lp->dinfeas) : mpq_get_d(p_mpq->lp->objbound)) );
	}

CLEANUP:
	mpq_EGlpNumClearVar (fi.totinfeas);
	return rval;
}

/* ========================================================================= */
/** @brief test whether given basis is dual feasible in rational arithmetic. 
 * if wanted it will first directly test the corresponding approximate dual and primal solution 
 * (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 * before performing the dual feasibility test on the more expensive exact basic solution. 
 * @param p_mpq   the problem data.
 * @param basis   basis to be tested.
 * @param useprestep whether to directly test approximate primal and dual solution first.
 * @param dbl_p_sol  approximate primal solution to use in prestep 
 *                   (NULL in order to compute it by dual simplex in double precision with given starting basis).
 * @param dbl_d_sol  approximate dual solution to use in prestep 
 *                   (NULL in order to compute it by dual simplex in double precision with given starting basis).
 * @param result  where to store whether given basis is dual feasible.
 * @param dobjval where to store dual solution value in case of dual feasibility (if not NULL).
 * @param msg_lvl message level.
 */
int QSexact_verify (
   mpq_QSdata * p_mpq,
   QSbasis* basis,
   int useprestep,
   double* dbl_p_sol,
   double* dbl_d_sol,
   char* result,
   mpq_t* dobjval,
   const int msg_lvl
)
{
   int rval = 0;

   //assert(basis);
   //assert(basis->nstruct);

   *result = 0;
            
   if( useprestep )
   {
      mpq_t *x_mpq = 0;
      mpq_t *y_mpq = 0;
      int status = 0;

      if( dbl_p_sol == NULL || dbl_d_sol == NULL  )
      {
         dbl_QSdata *p_dbl = 0;
         double *x_dbl = 0;
         double *y_dbl = 0;

         /* create double problem, warmstart with given basis and solve it using double precision 
          * this is only done to get approximate primal and dual solution corresponding to the given basis 
          */
         p_dbl = QScopy_prob_mpq_dbl(p_mpq, "dbl_problem");
   
         dbl_QSload_basis(p_dbl, basis);
		rval = dbl_ILLeditor_solve(p_dbl, DUAL_SIMPLEX);
         CHECKRVALG(rval, CLEANUP);
      
         rval = dbl_QSget_status(p_dbl, &status);
         CHECKRVALG(rval, CLEANUP);
      
         if( status == QS_LP_OPTIMAL )
         {
            /* get continued fraction approximation of approximate solution */
            x_dbl = dbl_EGlpNumAllocArray(p_dbl->qslp->ncols);
            y_dbl = dbl_EGlpNumAllocArray(p_dbl->qslp->nrows);

            rval = dbl_QSget_x_array(p_dbl, x_dbl);
            CHECKRVALG(rval, CLEANUP);
            rval = dbl_QSget_pi_array(p_dbl, y_dbl);
            CHECKRVALG(rval, CLEANUP);
            x_mpq = QScopy_array_dbl_mpq(x_dbl);
            y_mpq = QScopy_array_dbl_mpq(y_dbl);
            
            /* test optimality of constructed solution */
            basis = dbl_QSget_basis(p_dbl);
            rval = QSexact_optimal_test(p_mpq, x_mpq, y_mpq, basis);
            if( rval )
            {
               *result = 1;
               if( dobjval )
               {
                  rval = mpq_QSget_objval(p_mpq, dobjval);
                  if( rval )
                     *result = 0;
               }         
            }
            if( !msg_lvl )
            {
               MESSAGE(0, "Performing approximated solution check on %s, sucess=%s dobjval=%lg", 
                  p_mpq->name, 
                  *result ? "YES" : "NO",
                  *result ? mpq_get_d(*dobjval) : mpq_get_d(*dobjval));
            }
         }
      CLEANUP:
         dbl_EGlpNumFreeArray(x_dbl);
         dbl_EGlpNumFreeArray(y_dbl);
         mpq_EGlpNumFreeArray(x_mpq);
         mpq_EGlpNumFreeArray(y_mpq);
         dbl_QSfree_prob(p_dbl);
         rval = 0;
      }
      else
      {
         dbl_QSdata *p_dbl = 0;
         int i;

         /* for some reason, this help to avoid fails in QSexact_basis_dualstatus() after 
          * the test here fails, i.e., if we would not perform the test here, than QSexact_basis_dualstatus() would normally not fail
          * something happens with the basis... if we do not set up the dbl-prob (?) ????????????????????????
          */

         // AP: also confused by this
         p_dbl = QScopy_prob_mpq_dbl(p_mpq, "dbl_problem");
         dbl_QSload_basis(p_dbl, basis);

         x_mpq = mpq_EGlpNumAllocArray(p_mpq->qslp->ncols);
         y_mpq = mpq_EGlpNumAllocArray(p_mpq->qslp->nrows);

         /* get continued fraction approximation of approximate solution */
         for( i = 0; i < p_mpq->qslp->ncols; ++i )
            mpq_EGlpNumSet(x_mpq[i], dbl_p_sol[i]);

         for( i = 0; i < p_mpq->qslp->nrows; ++i )
            mpq_EGlpNumSet(y_mpq[i], dbl_d_sol[i]);
            
         /* test optimality of constructed solution */
         basis = dbl_QSget_basis(p_dbl);
         rval = QSexact_optimal_test(p_mpq, x_mpq, y_mpq, basis);
         if( rval )
         {
            *result = 1;
            if( dobjval )
            {
               rval = mpq_QSget_objval(p_mpq, dobjval);
               if( rval )
                  *result = 0;
            }         
         }
         if( !msg_lvl )
         {
            MESSAGE(0, "Performing approximated solution check on %s, sucess=%s dobjval=%lg", 
               p_mpq->name, 
               *result ? "YES" : "NO",
               *result ? mpq_get_d(*dobjval) : mpq_get_d(*dobjval));
         }
         mpq_EGlpNumFreeArray(x_mpq);
         mpq_EGlpNumFreeArray(y_mpq);
         dbl_QSfree_prob(p_dbl);
         rval = 0;
      }
   }

   if( !(*result) )
   {
      rval = QSexact_basis_dualstatus(p_mpq, basis, result, dobjval, msg_lvl);
      if( !msg_lvl )
      {
         MESSAGE(0, "Performing rational solution check on %s, sucess=%s dobjval=%lg", 
            p_mpq->name, 
            *result ? "YES" : "NO",
            *result ? mpq_get_d(*dobjval) : mpq_get_d(*dobjval));
      }
   }

   return rval;
}

/* ========================================================================= */
int QSexact_solver (mpq_QSdata * p_mpq, mpq_t * const x, mpq_t * const y, QSbasis * const ebasis, int simplexalgo, int *status)
{ 
	// clock start for timing purposes
        clock_t start = clock();     // full timer
	clock_t dbl_start = clock(); // double timer

	/* local variables */
	int last_status = 0, last_iter = 0;
	QSbasis *basis = 0;
	unsigned precision = EGLPNUM_PRECISION;
	int rval = 0,
	  it = QS_EXACT_MAX_ITER;
	dbl_QSdata *p_dbl = 0;
	mpf_QSdata *p_mpf = 0;
	double *x_dbl = 0,
	 *y_dbl = 0;
	mpq_t *x_mpq = 0,
	 *y_mpq = 0;
	mpf_t *x_mpf = 0,
	 *y_mpf = 0;
	int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0: (1 - p_mpq->simplex_display) * 10000;
	*status = 0;
	/* save the problem if we are really debugging */
	if(DEBUG >= __QS_SB_VERB)
	{
		EGcallD(mpq_QSwrite_prob(p_mpq, "qsxprob.lp","LP"));
	}
	/* try first with doubles */
	if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
	{
		QSlog("Trying double precision");

		// AP: save double precision to file
		EGioFile_t *out = 0;
		out = EGioOpen ("time_precision_data", "a");
		EGioPrintf (out, "64 ");
		EGioClose (out);
	}
	p_dbl = QScopy_prob_mpq_dbl (p_mpq, "dbl_problem");
	if(__QS_SB_VERB <= DEBUG) p_dbl->simplex_display = 1;
	if (ebasis && ebasis->nstruct)
		dbl_QSload_basis (p_dbl, ebasis); // AP: EGLPNUM_TYPENAME_QSload_basis in qsopt.c
	if (dbl_ILLeditor_solve (p_dbl, simplexalgo))
	{
		MESSAGE(p_mpq->simplex_display ? 0: __QS_SB_VERB, 
						"double approximation failed, code %d, "
						"continuing in extended precision", rval);
		goto MPF_PRECISION;
	}

	// testing for log file
	clock_t dbl_end = clock();
        double duration = (double)(dbl_end - dbl_start) / CLOCKS_PER_SEC;
        log_timing("DBL solve took ", duration);
	log_message("------------------------------------------------------------");

	EGcallD(dbl_QSget_status (p_dbl, status));
	if ((*status == QS_LP_INFEASIBLE) && (p_dbl->lp->final_phase != PRIMAL_PHASEI) && (p_dbl->lp->final_phase != DUAL_PHASEII))
		dbl_QSopt_primal (p_dbl, status);
	EGcallD(dbl_QSget_status (p_dbl, status));
	last_status = *status;
	EGcallD(dbl_QSget_itcnt(p_dbl, 0, 0, 0, 0, &last_iter));
	/* deal with the problem depending on what status we got from our optimizer */
	switch (*status)
	{
	case QS_LP_OPTIMAL:
		x_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->ncols);
		y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
		EGcallD(dbl_QSget_x_array (p_dbl, x_dbl));
		EGcallD(dbl_QSget_pi_array (p_dbl, y_dbl));
		x_mpq = QScopy_array_dbl_mpq (x_dbl);
		y_mpq = QScopy_array_dbl_mpq (y_dbl);
		dbl_EGlpNumFreeArray (x_dbl);
		dbl_EGlpNumFreeArray (y_dbl);
		basis = dbl_QSget_basis (p_dbl);
		if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
		{
			optimal_output (p_mpq, x, y, x_mpq, y_mpq);
			goto CLEANUP;
		}
		else
		{
			EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
			if (*status == QS_LP_OPTIMAL)
			{
				if(!msg_lvl)
				{
					MESSAGE(0,"Retesting solution");
				}
				EGcallD(mpq_QSget_x_array (p_mpq, x_mpq));
				EGcallD(mpq_QSget_pi_array (p_mpq, y_mpq));
				if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
				{
					optimal_output (p_mpq, x, y, x_mpq, y_mpq);
					goto CLEANUP;
				}
				else
				{
					last_status = *status = QS_LP_UNSOLVED;
				}
			}
			else
			{
				if(!msg_lvl)
				{
					MESSAGE(0,"Status is not optimal, but %d", *status);
				}
			}
		}
		mpq_EGlpNumFreeArray (x_mpq);
		mpq_EGlpNumFreeArray (y_mpq);
		break;
	case QS_LP_INFEASIBLE:
		y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
		if (dbl_QSget_infeas_array (p_dbl, y_dbl))
		{
			MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "double approximation"
							" failed, code %d, continuing in extended precision\n", rval);
			goto MPF_PRECISION;
		}
		y_mpq = QScopy_array_dbl_mpq (y_dbl);
		dbl_EGlpNumFreeArray (y_dbl);
		if (QSexact_infeasible_test (p_mpq, y_mpq))
		{
			infeasible_output (p_mpq, y, y_mpq);
			goto CLEANUP;
		}
		else
		{
			MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
			basis = dbl_QSget_basis (p_dbl);
			EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
			#if 0
			mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_MAX_ITERATIONS, 1);
			mpq_QSload_basis (p_mpq, basis);
			mpq_QSfree_basis (basis);
			EGcallD(mpq_ILLeditor_solve (p_mpq, simplexalgo));
			EGcallD(mpq_QSget_status (p_mpq, status));
			#endif
			if (*status == QS_LP_INFEASIBLE)
			{
				mpq_EGlpNumFreeArray (y_mpq);
				y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
				EGcallD(mpq_QSget_infeas_array (p_mpq, y_mpq));
				if (QSexact_infeasible_test (p_mpq, y_mpq))
				{
					infeasible_output (p_mpq, y, y_mpq);
					goto CLEANUP;
				}
				else
				{
					last_status = *status = QS_LP_UNSOLVED;
				}
			}
		}
		mpq_EGlpNumFreeArray (y_mpq);
		break;
	case QS_LP_UNBOUNDED:
		MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "%s\n\tUnbounded "
						"Problem found, not implemented to deal with this\n%s\n",__sp,__sp);
		break;
	case QS_LP_OBJ_LIMIT:
		rval=1;
		IFMESSAGE(p_mpq->simplex_display,"Objective limit reached (in floating point) ending now");
		goto CLEANUP;
		break;
	default:
		IFMESSAGE(p_mpq->simplex_display,"Re-trying inextended precision");
		break;
	}
	/* if we reach this point, then we have to keep going, we use the previous
	 * basis ONLY if the previous precision thinks that it has the optimal
	 * solution, otherwise we start from scratch. */
	precision = 128;
	MPF_PRECISION:
		dbl_QSfree_prob (p_dbl);
	p_dbl = 0;
	/* try with multiple precision floating points */
	for (; it--; precision = (unsigned) (precision * 1.5))
	{
		clock_t mpf_start = clock();

		// AP: save precision to file
		EGioFile_t *out = 0;
		out = EGioOpen ("time_precision_data", "a");
		EGioPrintf (out, "%d ", precision);
		EGioClose (out);

		QSexact_set_precision (precision);
		if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
		{
			QSlog("Trying mpf with %u bits", precision);
		}
		p_mpf = QScopy_prob_mpq_mpf (p_mpq, "mpf_problem");
		if(DEBUG >= __QS_SB_VERB)
		{
			EGcallD(mpf_QSwrite_prob(p_mpf, "qsxprob.mpf.lp","LP"));
		}
		if(__QS_SB_VERB <= DEBUG) p_mpf->simplex_display = 1;
		simplexalgo = PRIMAL_SIMPLEX;
		if(!last_iter) last_status = QS_LP_UNSOLVED;

		//AP : tester code
		QSlog("HUZAHHHHHHHHHH %d %d", last_iter, last_status);
		if(last_status == QS_LP_OPTIMAL || last_status == QS_LP_INFEASIBLE)
		{
			if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
			{
				QSlog("Re-using previous basis");

			}
			if (basis)
			{
				EGcallD(mpf_QSload_basis (p_mpf, basis));
				mpf_QSfree_basis (basis);
				simplexalgo = DUAL_SIMPLEX;
				basis = 0;
			}
			else if (ebasis && ebasis->nstruct)
			{
				mpf_QSload_basis (p_mpf, ebasis);
				simplexalgo = DUAL_SIMPLEX;
			}
		}
		else
		{
			if(p_mpf->basis)
			{
				mpf_ILLlp_basis_free(p_mpf->basis);
				p_mpf->lp->basisid = -1;
				p_mpf->factorok = 0;
			}
			if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
			{
				QSlog("Not-using previous basis");
			}
		}

		if (mpf_ILLeditor_solve (p_mpf, simplexalgo))
		{
			if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
			{
				QSlog("mpf_%u precision falied, error code %d, continuing with "
										"next precision", precision, rval);
			 }
			// stop timing before exiting iteration
            		clock_t mpf_end = clock();
            		double elapsed_mpf = (double)(mpf_end - mpf_start) / CLOCKS_PER_SEC;
            		char label[128];
            		snprintf(label, sizeof(label), "MPF solve at %u bits took ", precision);
            		log_timing(label, elapsed_mpf);
			goto NEXT_PRECISION;
		}
		EGcallD(mpf_QSget_status (p_mpf, status));
		if ((*status == QS_LP_INFEASIBLE) &&
				(p_mpf->lp->final_phase != PRIMAL_PHASEI) &&
				(p_mpf->lp->final_phase != DUAL_PHASEII))
			mpf_QSopt_primal (p_mpf, status);
		EGcallD(mpf_QSget_status (p_mpf, status));
		last_status = *status;
		EGcallD(mpf_QSget_itcnt(p_mpf, 0, 0, 0, 0, &last_iter));
		/* deal with the problem depending on status we got from our optimizer */
		switch (*status)
		{
		case QS_LP_OPTIMAL:
			basis = mpf_QSget_basis (p_mpf);
			x_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->ncols);
			y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
			EGcallD(mpf_QSget_x_array (p_mpf, x_mpf));
			EGcallD(mpf_QSget_pi_array (p_mpf, y_mpf));
			x_mpq = QScopy_array_mpf_mpq (x_mpf);
			y_mpq = QScopy_array_mpf_mpq (y_mpf);
			mpf_EGlpNumFreeArray (x_mpf);
			mpf_EGlpNumFreeArray (y_mpf);
			if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
			{
				optimal_output (p_mpq, x, y, x_mpq, y_mpq);
            			clock_t mpf_end = clock();
            			double elapsed_mpf = (double)(mpf_end - mpf_start) / CLOCKS_PER_SEC;
            			char label[128];
            			snprintf(label, sizeof(label), "MPF solve at %u bits took ", precision);
            			log_timing(label, elapsed_mpf);
				goto CLEANUP;
			}
			else
			{
				EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
				if (*status == QS_LP_OPTIMAL)
				{
					MESSAGE (msg_lvl, "Retesting solution");
					EGcallD(mpq_QSget_x_array (p_mpq, x_mpq));
					EGcallD(mpq_QSget_pi_array (p_mpq, y_mpq));
					if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
					{
						optimal_output (p_mpq, x, y, x_mpq, y_mpq);
            					clock_t mpf_end = clock();
            					double elapsed_mpf = (double)(mpf_end - mpf_start) / CLOCKS_PER_SEC;
            					char label[128];
            					snprintf(label, sizeof(label), "MPF solve at %u bits took ", precision);
            					log_timing(label, elapsed_mpf);
						goto CLEANUP;
					}
					else
					{
						last_status = *status = QS_LP_UNSOLVED;
					}
				}
				else
					MESSAGE (msg_lvl, "Status is not optimal, but %d", *status);
			}
			mpq_EGlpNumFreeArray (x_mpq);
			mpq_EGlpNumFreeArray (y_mpq);
			break;
		case QS_LP_INFEASIBLE:
			y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
			EGcallD(mpf_QSget_infeas_array (p_mpf, y_mpf));
			y_mpq = QScopy_array_mpf_mpq (y_mpf);
			mpf_EGlpNumFreeArray (y_mpf);
			if (QSexact_infeasible_test (p_mpq, y_mpq))
			{
				infeasible_output (p_mpq, y, y_mpq);
            			clock_t mpf_end = clock();
            			double elapsed_mpf = (double)(mpf_end - mpf_start) / CLOCKS_PER_SEC;
            			char label[128];
            			snprintf(label, sizeof(label), "MPF solve at %u bits took ", precision);
            			log_timing(label, elapsed_mpf);
				goto CLEANUP;
			}
			else
			{
				MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
				basis = mpf_QSget_basis (p_mpf);
				EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
#if 0
				mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_MAX_ITERATIONS, 1);
				mpq_QSload_basis (p_mpq, basis);
				mpq_QSfree_basis (basis);
				EGcallD(mpq_ILLeditor_solve (p_mpq, simplexalgo));
				EGcallD(mpq_QSget_status (p_mpq, status));
#endif
				if (*status == QS_LP_INFEASIBLE)
				{
					mpq_EGlpNumFreeArray (y_mpq);
					y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
					EGcallD(mpq_QSget_infeas_array (p_mpq, y_mpq));
					if (QSexact_infeasible_test (p_mpq, y_mpq))
					{
						infeasible_output (p_mpq, y, y_mpq);
						goto CLEANUP;
					}
					else
					{
						last_status = *status = QS_LP_UNSOLVED;
					}
				}
			}
			mpq_EGlpNumFreeArray (y_mpq);
			break;
			break;
		case QS_LP_OBJ_LIMIT:
			rval=1;
			IFMESSAGE(p_mpq->simplex_display,"Objective limit reached (in floating point) ending now");
			goto CLEANUP;
			break;
		case QS_LP_UNBOUNDED:
		default:
			MESSAGE(__QS_SB_VERB,"Re-trying inextended precision");
			break;
		}
	NEXT_PRECISION:
		clock_t mpf_end = clock();
        	double elapsed_mpf = (double)(mpf_end - mpf_start) / CLOCKS_PER_SEC;
        	char label[128];
        	snprintf(label, sizeof(label), "MPF solve at %u bits took ", precision);
        	log_timing(label, elapsed_mpf);	
	
		mpf_QSfree_prob (p_mpf);
		p_mpf = 0;
	}
	/* ending */
CLEANUP:
	dbl_EGlpNumFreeArray (x_dbl);
	dbl_EGlpNumFreeArray (y_dbl);
	mpq_EGlpNumFreeArray (x_mpq);
	mpq_EGlpNumFreeArray (y_mpq);
	mpf_EGlpNumFreeArray (x_mpf);
	mpf_EGlpNumFreeArray (y_mpf);
	if (ebasis && basis)
	{
		ILL_IFFREE(ebasis->cstat);
		ILL_IFFREE(ebasis->rstat);
		ebasis->nstruct = basis->nstruct;
		ebasis->nrows = basis->nrows;
		ebasis->cstat = basis->cstat;
		ebasis->rstat = basis->rstat;
		basis->cstat = basis->rstat = 0;
	}
	mpq_QSfree_basis (basis);
	dbl_QSfree_prob (p_dbl);
	mpf_QSfree_prob (p_mpf);

	// testing for log file
        clock_t end = clock();
        double total_duration = (double)(end - start) / CLOCKS_PER_SEC;
        log_timing("QSexact Solver took ", total_duration);

	return rval;
}

/* ========================================================================= */
int __QSexact_setup = 0;
/* ========================================================================= */
void QSexactStart(void)
{
	/* if we have been initialized before, do nothing */
	if(__QSexact_setup) return;
	/* we should call EGlpNumStart() */
	EGlpNumStart();
	
	/* now we call all setups */
	EXutilDoInit();
	dbl_ILLstart();
	mpf_ILLstart();
	mpq_ILLstart();
	/* ending */
	__QSexact_setup = 1;
}
/* ========================================================================= */
void QSexactClear(void)
{
	if(!__QSexact_setup) return;
	/* now we call all ends */
	dbl_ILLend();
	mpf_ILLend();
	mpq_ILLend();
	EXutilDoClear();
	/* ending */
	EGlpNumClear();
	__QSexact_setup = 0;
}
/* ========================================================================= */
/** @} */
/* end of exact.c */

/* ========================================================================= */
int mpq_factor_work_to_dbl_factor_work(dbl_factor_work *dest, const mpq_factor_work *src)
{
	int i;
	int rval = 0;
	int dsize;

	// Initialize the destination structure with defaults
	dbl_ILLfactor_init_factor_work(dest);

	// Copy scalar fields from src to dest, converting mpq_t to double where necessary
	dest->max_k = src->max_k;
	dest->fzero_tol = mpq_get_d(src->fzero_tol);
	dest->szero_tol = mpq_get_d(src->szero_tol);
	dest->partial_tol = mpq_get_d(src->partial_tol);
	dest->ur_space_mul = src->ur_space_mul;
	dest->uc_space_mul = src->uc_space_mul;
	dest->lc_space_mul = src->lc_space_mul;
	dest->lr_space_mul = src->lr_space_mul;
	dest->er_space_mul = src->er_space_mul;
	dest->grow_mul = src->grow_mul;
	dest->p = src->p;
	dest->etamax = src->etamax;
	dest->minmult = src->minmult;
	dest->maxmult = src->maxmult;
	dest->updmaxmult = src->updmaxmult;
	dest->dense_fract = src->dense_fract;
	dest->dense_min = src->dense_min;
	dest->maxelem_orig = mpq_get_d(src->maxelem_orig);
	dest->nzcnt_orig = src->nzcnt_orig;
	dest->maxelem_factor = mpq_get_d(src->maxelem_factor);
	dest->nzcnt_factor = src->nzcnt_factor;
	dest->maxelem_cur = mpq_get_d(src->maxelem_cur);
	dest->nzcnt_cur = src->nzcnt_cur;
	dest->partial_cur = mpq_get_d(src->partial_cur);
	dest->dim = src->dim;
	dest->stage = src->stage;
	dest->nstages = src->nstages;
	dest->etacnt = src->etacnt;
	dest->ur_space = src->ur_space;
	dest->uc_space = src->uc_space;
	dest->lc_space = src->lc_space;
	dest->lr_space = src->lr_space;
	dest->er_space = src->er_space;
	dest->ur_freebeg = src->ur_freebeg;
	dest->uc_freebeg = src->uc_freebeg;
	dest->lc_freebeg = src->lc_freebeg;
	dest->lr_freebeg = src->lr_freebeg;
	dest->er_freebeg = src->er_freebeg;
	dest->drows = src->drows;
	dest->dcols = src->dcols;
	dest->dense_base = src->dense_base;

	// Shallow copy for these pointers (they point to shared data)
	dest->p_nsing = src->p_nsing;
	dest->p_singr = src->p_singr;
	dest->p_singc = src->p_singc;

	// Convert work arrays
	if (src->work_coef) {
		dest->work_coef = dbl_EGlpNumAllocArray(src->dim);
		if (!dest->work_coef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->dim; i++) {
			dest->work_coef[i] = mpq_get_d(src->work_coef[i]);
		}
	}

	if (src->work_indx) {
		ILL_SAFE_MALLOC(dest->work_indx, src->dim, int);
		memcpy(dest->work_indx, src->work_indx, src->dim * sizeof(int));
	}

	// Convert info structs
	if (src->uc_inf) {
		ILL_SAFE_MALLOC(dest->uc_inf, src->dim + src->max_k + 1, dbl_uc_info);
		for (i = 0; i < src->dim + src->max_k + 1; i++) {
			dest->uc_inf[i].cbeg = src->uc_inf[i].cbeg;
			dest->uc_inf[i].nzcnt = src->uc_inf[i].nzcnt;
			dest->uc_inf[i].next = src->uc_inf[i].next;
			dest->uc_inf[i].prev = src->uc_inf[i].prev;
			dest->uc_inf[i].delay = src->uc_inf[i].delay;
		}
	}

	if (src->ur_inf) {
		ILL_SAFE_MALLOC(dest->ur_inf, src->dim + src->max_k + 1, dbl_ur_info);
		for (i = 0; i < src->dim + src->max_k + 1; i++) {
			dest->ur_inf[i].max = mpq_get_d(src->ur_inf[i].max);
			dest->ur_inf[i].rbeg = src->ur_inf[i].rbeg;
			dest->ur_inf[i].nzcnt = src->ur_inf[i].nzcnt;
			dest->ur_inf[i].pivcnt = src->ur_inf[i].pivcnt;
			dest->ur_inf[i].next = src->ur_inf[i].next;
			dest->ur_inf[i].prev = src->ur_inf[i].prev;
			dest->ur_inf[i].delay = src->ur_inf[i].delay;
		}
	}

	if (src->lc_inf) {
		ILL_SAFE_MALLOC(dest->lc_inf, src->dim, dbl_lc_info);
		for (i = 0; i < src->dim; i++) {
			dest->lc_inf[i].cbeg = src->lc_inf[i].cbeg;
			dest->lc_inf[i].nzcnt = src->lc_inf[i].nzcnt;
			dest->lc_inf[i].c = src->lc_inf[i].c;
			dest->lc_inf[i].crank = src->lc_inf[i].crank;
			dest->lc_inf[i].delay = src->lc_inf[i].delay;
		}
	}

	if (src->lr_inf) {
		ILL_SAFE_MALLOC(dest->lr_inf, src->dim, dbl_lr_info);
		for (i = 0; i < src->dim; i++) {
			dest->lr_inf[i].rbeg = src->lr_inf[i].rbeg;
			dest->lr_inf[i].nzcnt = src->lr_inf[i].nzcnt;
			dest->lr_inf[i].r = src->lr_inf[i].r;
			dest->lr_inf[i].rrank = src->lr_inf[i].rrank;
			dest->lr_inf[i].delay = src->lr_inf[i].delay;
		}
	}

	if (src->er_inf) {
		ILL_SAFE_MALLOC(dest->er_inf, src->etamax, dbl_er_info);
		for (i = 0; i < src->etamax; i++) {
			dest->er_inf[i].rbeg = src->er_inf[i].rbeg;
			dest->er_inf[i].nzcnt = src->er_inf[i].nzcnt;
			dest->er_inf[i].r = src->er_inf[i].r;
		}
	}

	// Convert U matrix data
	if (src->ucindx) {
		ILL_SAFE_MALLOC(dest->ucindx, src->uc_space + 1, int);
		memcpy(dest->ucindx, src->ucindx, (src->uc_space + 1) * sizeof(int));
	}

	if (src->ucrind) {
		ILL_SAFE_MALLOC(dest->ucrind, src->uc_space, int);
		memcpy(dest->ucrind, src->ucrind, src->uc_space * sizeof(int));
	}

	if (src->uccoef) {
		dest->uccoef = dbl_EGlpNumAllocArray(src->uc_space);
		if (!dest->uccoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->uc_space; i++) {
			dest->uccoef[i] = mpq_get_d(src->uccoef[i]);
		}
	}

	if (src->urindx) {
		ILL_SAFE_MALLOC(dest->urindx, src->ur_space + 1, int);
		memcpy(dest->urindx, src->urindx, (src->ur_space + 1) * sizeof(int));
	}

	if (src->urcind) {
		ILL_SAFE_MALLOC(dest->urcind, src->ur_space, int);
		memcpy(dest->urcind, src->urcind, src->ur_space * sizeof(int));
	}

	if (src->urcoef) {
		dest->urcoef = dbl_EGlpNumAllocArray(src->ur_space);
		if (!dest->urcoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->ur_space; i++) {
			dest->urcoef[i] = mpq_get_d(src->urcoef[i]);
		}
	}

	// Convert L matrix data
	if (src->lcindx) {
		ILL_SAFE_MALLOC(dest->lcindx, src->lc_space, int);
		memcpy(dest->lcindx, src->lcindx, src->lc_space * sizeof(int));
	}

	if (src->lccoef) {
		dest->lccoef = dbl_EGlpNumAllocArray(src->lc_space);
		if (!dest->lccoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->lc_space; i++) {
			dest->lccoef[i] = mpq_get_d(src->lccoef[i]);
		}
	}

	if (src->lrindx) {
		int lr_nzcnt = 0;
		for (i = 0; i < src->dim; i++) {
			lr_nzcnt += src->lr_inf[i].nzcnt;
		}
		ILL_SAFE_MALLOC(dest->lrindx, lr_nzcnt + 1, int);
		memcpy(dest->lrindx, src->lrindx, (lr_nzcnt + 1) * sizeof(int));
	}

	if (src->lrcoef) {
		int lr_nzcnt = 0;
		for (i = 0; i < src->dim; i++) {
			lr_nzcnt += src->lr_inf[i].nzcnt;
		}
		dest->lrcoef = dbl_EGlpNumAllocArray(lr_nzcnt);
		if (!dest->lrcoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < lr_nzcnt; i++) {
			dest->lrcoef[i] = mpq_get_d(src->lrcoef[i]);
		}
	}

	// Convert Eta data
	if (src->erindx) {
		ILL_SAFE_MALLOC(dest->erindx, src->er_space, int);
		memcpy(dest->erindx, src->erindx, src->er_space * sizeof(int));
	}

	if (src->ercoef) {
		dest->ercoef = dbl_EGlpNumAllocArray(src->er_space);
		if (!dest->ercoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->er_space; i++) {
			dest->ercoef[i] = mpq_get_d(src->ercoef[i]);
		}
	}

	// Convert permutations
	if (src->rperm) {
		ILL_SAFE_MALLOC(dest->rperm, src->dim, int);
		memcpy(dest->rperm, src->rperm, src->dim * sizeof(int));
	}

	if (src->rrank) {
		ILL_SAFE_MALLOC(dest->rrank, src->dim, int);
		memcpy(dest->rrank, src->rrank, src->dim * sizeof(int));
	}

	if (src->cperm) {
		ILL_SAFE_MALLOC(dest->cperm, src->dim, int);
		memcpy(dest->cperm, src->cperm, src->dim * sizeof(int));
	}

	if (src->crank) {
		ILL_SAFE_MALLOC(dest->crank, src->dim, int);
		memcpy(dest->crank, src->crank, src->dim * sizeof(int));
	}

	// Convert dense matrix
	if (src->dmat) {
		dsize = src->drows * src->dcols;
		dest->dmat = dbl_EGlpNumAllocArray(dsize);
		if (!dest->dmat) { rval = 1; goto CLEANUP; }
		for (i = 0; i < dsize; i++) {
			dest->dmat[i] = mpq_get_d(src->dmat[i]);
		}
	}

	// Convert svector xtmp
	rval = dbl_ILLsvector_alloc(&dest->xtmp, src->dim);
	if (rval) goto CLEANUP;
	dest->xtmp.nzcnt = src->xtmp.nzcnt;
	if (src->xtmp.nzcnt > 0) {
		memcpy(dest->xtmp.indx, src->xtmp.indx, src->xtmp.nzcnt * sizeof(int));
		for (i = 0; i < src->xtmp.nzcnt; i++) {
			dest->xtmp.coef[i] = mpq_get_d(src->xtmp.coef[i]);
		}
	}

	goto FINAL;

CLEANUP:
	// If any allocation fails, clean up
	printf("Error allocating memory for mpq_factor_work to dbl_factor_work conversion.\n");
	dbl_ILLfactor_free_factor_work(dest);

FINAL:
	return rval;
}

/* ========================================================================= */
int mpq_factor_work_to_mpf_factor_work(mpf_factor_work *dest, const mpq_factor_work *src)
{
	int i;
	int rval = 0;
	int dsize;
	// Initialize all mpf_t variables in the structure before calling init_factor_work
	mpf_EGlpNumInitVar(dest->fzero_tol);
	mpf_EGlpNumInitVar(dest->szero_tol);
	mpf_EGlpNumInitVar(dest->partial_tol);
	mpf_EGlpNumInitVar(dest->maxelem_orig);
	mpf_EGlpNumInitVar(dest->maxelem_factor);
	mpf_EGlpNumInitVar(dest->maxelem_cur);
	mpf_EGlpNumInitVar(dest->partial_cur);
	
	// Initialize the destination structure with defaults
	mpf_ILLfactor_init_factor_work(dest);
	// Copy scalar fields from src to dest, converting mpq_t to mpf_t where necessary
	dest->max_k = src->max_k;
	mpf_set_q(dest->fzero_tol, src->fzero_tol);
	mpf_set_q(dest->szero_tol, src->szero_tol);
	mpf_set_q(dest->partial_tol, src->partial_tol);
	dest->ur_space_mul = src->ur_space_mul;
	dest->uc_space_mul = src->uc_space_mul;
	dest->lc_space_mul = src->lc_space_mul;
	dest->lr_space_mul = src->lr_space_mul;
	dest->er_space_mul = src->er_space_mul;
	dest->grow_mul = src->grow_mul;
	dest->p = src->p;
	dest->etamax = src->etamax;
	dest->minmult = src->minmult;
	dest->maxmult = src->maxmult;
	dest->updmaxmult = src->updmaxmult;
	dest->dense_fract = src->dense_fract;
	dest->dense_min = src->dense_min;
	mpf_set_q(dest->maxelem_orig, src->maxelem_orig);
	dest->nzcnt_orig = src->nzcnt_orig;
	mpf_set_q(dest->maxelem_factor, src->maxelem_factor);
	dest->nzcnt_factor = src->nzcnt_factor;
	mpf_set_q(dest->maxelem_cur, src->maxelem_cur);
	dest->nzcnt_cur = src->nzcnt_cur;
	mpf_set_q(dest->partial_cur, src->partial_cur);
	dest->dim = src->dim;
	dest->stage = src->stage;
	dest->nstages = src->nstages;
	dest->etacnt = src->etacnt;
	dest->ur_space = src->ur_space;
	dest->uc_space = src->uc_space;
	dest->lc_space = src->lc_space;
	dest->lr_space = src->lr_space;
	dest->er_space = src->er_space;
	dest->ur_freebeg = src->ur_freebeg;
	dest->uc_freebeg = src->uc_freebeg;
	dest->lc_freebeg = src->lc_freebeg;
	dest->lr_freebeg = src->lr_freebeg;
	dest->er_freebeg = src->er_freebeg;
	dest->drows = src->drows;
	dest->dcols = src->dcols;
	dest->dense_base = src->dense_base;

	// Shallow copy for these pointers (they point to shared data)
	dest->p_nsing = src->p_nsing;
	dest->p_singr = src->p_singr;
	dest->p_singc = src->p_singc;
	// Convert work arrays
	if (src->work_coef) {
		dest->work_coef = mpf_EGlpNumAllocArray(src->dim);
		if (!dest->work_coef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->dim; i++) {
			mpf_set_q(dest->work_coef[i], src->work_coef[i]);
		}
	}
	if (src->work_indx) {
		ILL_SAFE_MALLOC(dest->work_indx, src->dim, int);
		memcpy(dest->work_indx, src->work_indx, src->dim * sizeof(int));
	}
	// Convert info structs
	if (src->uc_inf) {
		ILL_SAFE_MALLOC(dest->uc_inf, src->dim + src->max_k + 1, mpf_uc_info);
		for (i = 0; i < src->dim + src->max_k + 1; i++) {
			dest->uc_inf[i].cbeg = src->uc_inf[i].cbeg;
			dest->uc_inf[i].nzcnt = src->uc_inf[i].nzcnt;
			dest->uc_inf[i].next = src->uc_inf[i].next;
			dest->uc_inf[i].prev = src->uc_inf[i].prev;
			dest->uc_inf[i].delay = src->uc_inf[i].delay;
		}
	}
	if (src->ur_inf) {
		ILL_SAFE_MALLOC(dest->ur_inf, src->dim + src->max_k + 1, mpf_ur_info);
		for (i = 0; i < src->dim + src->max_k + 1; i++) {
			mpf_EGlpNumInitVar(dest->ur_inf[i].max);  // Initialize before setting!
			mpf_set_q(dest->ur_inf[i].max, src->ur_inf[i].max);
			dest->ur_inf[i].rbeg = src->ur_inf[i].rbeg;
			dest->ur_inf[i].nzcnt = src->ur_inf[i].nzcnt;
			dest->ur_inf[i].pivcnt = src->ur_inf[i].pivcnt;
			dest->ur_inf[i].next = src->ur_inf[i].next;
			dest->ur_inf[i].prev = src->ur_inf[i].prev;
			dest->ur_inf[i].delay = src->ur_inf[i].delay;
		}
	}
	if (src->lc_inf) {
		ILL_SAFE_MALLOC(dest->lc_inf, src->dim, mpf_lc_info);
		for (i = 0; i < src->dim; i++) {
			dest->lc_inf[i].cbeg = src->lc_inf[i].cbeg;
			dest->lc_inf[i].nzcnt = src->lc_inf[i].nzcnt;
			dest->lc_inf[i].c = src->lc_inf[i].c;
			dest->lc_inf[i].crank = src->lc_inf[i].crank;
			dest->lc_inf[i].delay = src->lc_inf[i].delay;
		}
	}	
	if (src->lr_inf) {
		ILL_SAFE_MALLOC(dest->lr_inf, src->dim, mpf_lr_info);
		for (i = 0; i < src->dim; i++) {
			dest->lr_inf[i].rbeg = src->lr_inf[i].rbeg;
			dest->lr_inf[i].nzcnt = src->lr_inf[i].nzcnt;
			dest->lr_inf[i].r = src->lr_inf[i].r;
			dest->lr_inf[i].rrank = src->lr_inf[i].rrank;
			dest->lr_inf[i].delay = src->lr_inf[i].delay;
		}
	}
	if (src->er_inf) {
		ILL_SAFE_MALLOC(dest->er_inf, src->etamax, mpf_er_info);
		for (i = 0; i < src->etamax; i++) {
			dest->er_inf[i].rbeg = src->er_inf[i].rbeg;
			dest->er_inf[i].nzcnt = src->er_inf[i].nzcnt;
			dest->er_inf[i].r = src->er_inf[i].r;
		}
	}
	// Convert U matrix data
	if (src->ucindx) {
		ILL_SAFE_MALLOC(dest->ucindx, src->uc_space + 1, int);
		memcpy(dest->ucindx, src->ucindx, (src->uc_space + 1) * sizeof(int));
	}
	if (src->ucrind) {
		ILL_SAFE_MALLOC(dest->ucrind, src->uc_space, int);
		memcpy(dest->ucrind, src->ucrind, src->uc_space * sizeof(int));
	}
	if (src->uccoef) {
		dest->uccoef = mpf_EGlpNumAllocArray(src->uc_space);
		if (!dest->uccoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->uc_space; i++) {
			mpf_set_q(dest->uccoef[i], src->uccoef[i]);
		}
	}
	if (src->urindx) {
		ILL_SAFE_MALLOC(dest->urindx, src->ur_space + 1, int);
		memcpy(dest->urindx, src->urindx, (src->ur_space + 1) * sizeof(int));
	}
	if (src->urcind) {
		ILL_SAFE_MALLOC(dest->urcind, src->ur_space, int);
		memcpy(dest->urcind, src->urcind, src->ur_space * sizeof(int));
	}
	if (src->urcoef) {
		dest->urcoef = mpf_EGlpNumAllocArray(src->ur_space);
		if (!dest->urcoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->ur_space; i++) {
			mpf_set_q(dest->urcoef[i], src->urcoef[i]);
		}
	}
	// Convert L matrix data
	if (src->lcindx) {
		ILL_SAFE_MALLOC(dest->lcindx, src->lc_space, int);
		memcpy(dest->lcindx, src->lcindx, src->lc_space * sizeof(int));
	}
	if (src->lccoef) {
		dest->lccoef = mpf_EGlpNumAllocArray(src->lc_space);
		if (!dest->lccoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->lc_space; i++) {
			mpf_set_q(dest->lccoef[i], src->lccoef[i]);
		}
	}
	if (src->lrindx) {
		int lr_nzcnt = 0;
		for (i = 0; i < src->dim; i++) {
			lr_nzcnt += src->lr_inf[i].nzcnt;
		}
		ILL_SAFE_MALLOC(dest->lrindx, lr_nzcnt + 1, int);
		memcpy(dest->lrindx, src->lrindx, (lr_nzcnt + 1) * sizeof(int));
	}
	if (src->lrcoef) {
		int lr_nzcnt = 0;
		for (i = 0; i < src->dim; i++) {
			lr_nzcnt += src->lr_inf[i].nzcnt;
		}
		dest->lrcoef = mpf_EGlpNumAllocArray(lr_nzcnt);
		if (!dest->lrcoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < lr_nzcnt; i++) {
			mpf_set_q(dest->lrcoef[i], src->lrcoef[i]);
		}
	}	
	// Convert Eta data
	if (src->erindx) {
		ILL_SAFE_MALLOC(dest->erindx, src->er_space, int);
		memcpy(dest->erindx, src->erindx, src->er_space * sizeof(int));
	}
	if (src->ercoef) {
		dest->ercoef = mpf_EGlpNumAllocArray(src->er_space);
		if (!dest->ercoef) { rval = 1; goto CLEANUP; }
		for (i = 0; i < src->er_space; i++) {
			mpf_set_q(dest->ercoef[i], src->ercoef[i]);
		}
	}
	// Convert permutations
	if (src->rperm) {
		ILL_SAFE_MALLOC(dest->rperm, src->dim, int);
		memcpy(dest->rperm, src->rperm, src->dim * sizeof(int));
	}
	if (src->rrank) {
		ILL_SAFE_MALLOC(dest->rrank, src->dim, int);
		memcpy(dest->rrank, src->rrank, src->dim * sizeof(int));
	}
	if (src->cperm) {
		ILL_SAFE_MALLOC(dest->cperm, src->dim, int);
		memcpy(dest->cperm, src->cperm, src->dim * sizeof(int));
	}
	if (src->crank) {
		ILL_SAFE_MALLOC(dest->crank, src->dim, int);
		memcpy(dest->crank, src->crank, src->dim * sizeof(int));
	}

	// Convert dense matrix
	if (src->dmat) {
		dsize = src->drows * src->dcols;
		dest->dmat = mpf_EGlpNumAllocArray(dsize);
		if (!dest->dmat) { rval = 1; goto CLEANUP; }
		for (i = 0; i < dsize; i++) {
			mpf_set_q(dest->dmat[i], src->dmat[i]);
		}
	}

	// Convert svector xtmp
	rval = mpf_ILLsvector_alloc(&dest->xtmp, src->dim);
	if (rval) goto CLEANUP;
	dest->xtmp.nzcnt = src->xtmp.nzcnt;
	if (src->xtmp.nzcnt > 0) {
		memcpy(dest->xtmp.indx, src->xtmp.indx, src->xtmp.nzcnt * sizeof(int));
		for (i = 0; i < src->xtmp.nzcnt; i++) {
			mpf_set_q(dest->xtmp.coef[i], src->xtmp.coef[i]);
		}
	}

	goto FINAL;

CLEANUP:
	// If any allocation fails, clean up
	printf("Error allocating memory for mpq_factor_work to mpf_factor_work conversion.\n");
	mpf_ILLfactor_free_factor_work(dest);

FINAL:
	return rval;
}

/* ========================================================================= */

