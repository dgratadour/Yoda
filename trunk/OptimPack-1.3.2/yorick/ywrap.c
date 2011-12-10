/* codger-generated yorick package wrapper file */
#include "play.h"
#include "ydata.h"

/*----------------begin OptimPack1.i */
extern BuiltIn Y___op_csrch;

extern int op_csrch(double , double , double *, double , double , 
  double , double , double , int *, char *, long *, double *);
void
Y___op_csrch(int n)
{
  if (n!=12) YError("__op_csrch takes exactly 12 arguments");
  PushIntValue(op_csrch(yarg_sd(11), yarg_sd(10), yarg_d(9,0), 
    yarg_sd(8), yarg_sd(7), yarg_sd(6), yarg_sd(5), yarg_sd(4), 
    yarg_i(3,0), yarg_c(2,0), yarg_l(1,0), yarg_d(0,0)));
}

extern BuiltIn Y___op_vmlmb_next;

extern int op_vmlmb_next(double *, double *, double *, void *, 
  void *, char *, long *, double *);
void
Y___op_vmlmb_next(int n)
{
  if (n!=8) YError("__op_vmlmb_next takes exactly 8 arguments");
  PushIntValue(op_vmlmb_next(yarg_d(7,0), yarg_d(6,0), yarg_d(5,0), 
    yarg_sp(4), yarg_sp(3), yarg_c(2,0), yarg_l(1,0), yarg_d(0,0)));
}

extern BuiltIn Y___op_vmlmb_first;

extern int op_vmlmb_first(long , long , double , double , double , 
  double , double , double , double , char *, long *, double *);
void
Y___op_vmlmb_first(int n)
{
  if (n!=12) YError("__op_vmlmb_first takes exactly 12 arguments");
  PushIntValue(op_vmlmb_first(yarg_sl(11), yarg_sl(10), yarg_sd(9), 
    yarg_sd(8), yarg_sd(7), yarg_sd(6), yarg_sd(5), yarg_sd(4), 
    yarg_sd(3), yarg_c(2,0), yarg_l(1,0), yarg_d(0,0)));
}

extern BuiltIn Y___op_vmlmb_set_fmin;

extern int op_vmlmb_set_fmin(char *, long *, double *, double , 
  double *);
void
Y___op_vmlmb_set_fmin(int n)
{
  if (n!=5) YError("__op_vmlmb_set_fmin takes exactly 5 arguments");
  PushIntValue(op_vmlmb_set_fmin(yarg_c(4,0), yarg_l(3,0), yarg_d(2,0), 
    yarg_sd(1), yarg_d(0,0)));
}

extern BuiltIn Y___op_vmlmb_get_fmin;

extern int op_vmlmb_get_fmin(char *, long *, double *, double *);
void
Y___op_vmlmb_get_fmin(int n)
{
  if (n!=4) YError("__op_vmlmb_get_fmin takes exactly 4 arguments");
  PushIntValue(op_vmlmb_get_fmin(yarg_c(3,0), yarg_l(2,0), yarg_d(1,0), 
    yarg_d(0,0)));
}


/*----------------list include files */

static char *y0_includes[] = {
  "OptimPack1.i",
  0,
  0
};

/*----------------collect pointers and names */

static BuiltIn *y0_routines[] = {
  &Y___op_csrch,
  &Y___op_vmlmb_next,
  &Y___op_vmlmb_first,
  &Y___op_vmlmb_set_fmin,
  &Y___op_vmlmb_get_fmin,
  0
};

static void *y0_values[] = {
  0
};

static char *y0_names[] = {
  "__op_csrch",
  "__op_vmlmb_next",
  "__op_vmlmb_first",
  "__op_vmlmb_set_fmin",
  "__op_vmlmb_get_fmin",
  0
};

/*----------------define package initialization function */

PLUG_EXPORT char *yk_OptimPack1(char ***,
                         BuiltIn ***, void ***, char ***);
static char *y0_pkgname = "OptimPack1";

char *
yk_OptimPack1(char ***ifiles,
       BuiltIn ***code, void ***data, char ***varname)
{
  *ifiles = y0_includes;
  *code = y0_routines;
  *data = y0_values;
  *varname = y0_names;
  return y0_pkgname;
}
