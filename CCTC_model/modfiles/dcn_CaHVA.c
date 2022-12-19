/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__dcnCaHVA
#define _nrn_initial _nrn_initial__dcnCaHVA
#define nrn_cur _nrn_cur__dcnCaHVA
#define _nrn_current _nrn_current__dcnCaHVA
#define nrn_jacob _nrn_jacob__dcnCaHVA
#define nrn_state _nrn_state__dcnCaHVA
#define _net_receive _net_receive__dcnCaHVA 
#define _f_rate _f_rate__dcnCaHVA 
#define rate rate__dcnCaHVA 
#define states states__dcnCaHVA 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define perm _p[0]
#define ica _p[1]
#define m _p[2]
#define cai _p[3]
#define cao _p[4]
#define minf _p[5]
#define taum _p[6]
#define T _p[7]
#define A _p[8]
#define Dm _p[9]
#define v _p[10]
#define _g _p[11]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_getGHKexp(void);
 static void _hoc_rate(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_dcnCaHVA", _hoc_setdata,
 "getGHKexp_dcnCaHVA", _hoc_getGHKexp,
 "rate_dcnCaHVA", _hoc_rate,
 0, 0
};
#define _f_getGHKexp _f_getGHKexp_dcnCaHVA
#define getGHKexp getGHKexp_dcnCaHVA
 extern double _f_getGHKexp( _threadargsprotocomma_ double );
 extern double getGHKexp( _threadargsprotocomma_ double );
 
static void _check_rate(double*, Datum*, Datum*, _NrnThread*); 
static void _check_getGHKexp(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_rate(_p, _ppvar, _thread, _nt);
   _check_getGHKexp(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define qdeltat qdeltat_dcnCaHVA
 double qdeltat = 1;
#define usetable usetable_dcnCaHVA
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_dcnCaHVA", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "perm_dcnCaHVA", "cm/seconds",
 "ica_dcnCaHVA", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "qdeltat_dcnCaHVA", &qdeltat_dcnCaHVA,
 "usetable_dcnCaHVA", &usetable_dcnCaHVA,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"dcnCaHVA",
 "perm_dcnCaHVA",
 0,
 "ica_dcnCaHVA",
 0,
 "m_dcnCaHVA",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	perm = 7.5e-006;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _dcn_CaHVA_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 dcnCaHVA D:/Projects/SchreglmannEtAl2020/CCTC_model/modfiles/dcn_CaHVA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_minf;
 static double *_t_taum;
 static double *_t_getGHKexp;
static int _reset;
static char *modelname = "High voltage activated calcium current (CaHVA) of deep cerebellar nucleus (DCN) neuron";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rate(_threadargsprotocomma_ double);
static int rate(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static double _n_getGHKexp(_threadargsprotocomma_ double _lv);
 static void _n_rate(_threadargsprotocomma_ double _lv);
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rate ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / taum ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rate ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taum )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rate ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taum)))*(- ( ( ( minf ) ) / taum ) / ( ( ( ( - 1.0 ) ) ) / taum ) - m) ;
   }
  return 0;
}
 static double _mfac_rate, _tmin_rate;
  static void _check_rate(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rate =  - 150.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rate)/300.; _mfac_rate = 1./_dx;
   for (_i=0, _x=_tmin_rate; _i < 301; _x += _dx, _i++) {
    _f_rate(_p, _ppvar, _thread, _nt, _x);
    _t_minf[_i] = minf;
    _t_taum[_i] = taum;
   }
  }
 }

 static int rate(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_rate(_p, _ppvar, _thread, _nt);
#endif
 _n_rate(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_rate(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rate(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_rate * (_lv - _tmin_rate);
 if (isnan(_xi)) {
  minf = _xi;
  taum = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 taum = _t_taum[0];
 return; }
 if (_xi >= 300.) {
 minf = _t_minf[300];
 taum = _t_taum[300];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 taum = _t_taum[_i] + _theta*(_t_taum[_i+1] - _t_taum[_i]);
 }

 
static int  _f_rate ( _threadargsprotocomma_ double _lv ) {
   if ( fabs ( _lv + 34.5 ) < 1e-6 ) {
     minf = 1.0 / ( 1.0 + exp ( ( _lv + 34.50001 ) / - 9.0 ) ) ;
     }
   else {
     minf = 1.0 / ( 1.0 + exp ( ( _lv + 34.5 ) / - 9.0 ) ) ;
     }
   taum = 1.0 / ( ( 31.746 * ( pow( ( exp ( ( _lv - 5.0 ) / - 13.89 ) + 1.0 ) , - 1.0 ) ) ) + ( 3.97e-4 * ( _lv + 8.9 ) ) * ( pow( ( exp ( ( _lv + 8.9 ) / 5.0 ) - 1.0 ) , - 1.0 ) ) ) ;
   taum = taum / qdeltat ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_rate(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 rate ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_getGHKexp, _tmin_getGHKexp;
  static void _check_getGHKexp(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_T;
  if (!usetable) {return;}
  if (_sav_T != T) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_getGHKexp =  - 150.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_getGHKexp)/300.; _mfac_getGHKexp = 1./_dx;
   for (_i=0, _x=_tmin_getGHKexp; _i < 301; _x += _dx, _i++) {
    _t_getGHKexp[_i] = _f_getGHKexp(_p, _ppvar, _thread, _nt, _x);
   }
   _sav_T = T;
  }
 }

 double getGHKexp(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_getGHKexp(_p, _ppvar, _thread, _nt);
#endif
 return _n_getGHKexp(_p, _ppvar, _thread, _nt, _lv);
 }

 static double _n_getGHKexp(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_getGHKexp(_p, _ppvar, _thread, _nt, _lv); 
}
 _xi = _mfac_getGHKexp * (_lv - _tmin_getGHKexp);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_getGHKexp[0];
 }
 if (_xi >= 300.) {
 return _t_getGHKexp[300];
 }
 _i = (int) _xi;
 return _t_getGHKexp[_i] + (_xi - (double)_i)*(_t_getGHKexp[_i+1] - _t_getGHKexp[_i]);
 }

 
double _f_getGHKexp ( _threadargsprotocomma_ double _lv ) {
   double _lgetGHKexp;
 _lgetGHKexp = exp ( - 23.20764929 * _lv / T ) ;
   
return _lgetGHKexp;
 }
 
static void _hoc_getGHKexp(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_getGHKexp(_p, _ppvar, _thread, _nt);
#endif
 _r =  getGHKexp ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  m = m0;
 {
   T = 273.15 + celsius ;
   rate ( _threadargscomma_ v ) ;
   m = minf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_rate(_p, _ppvar, _thread, _nt);
 _check_getGHKexp(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   A = getGHKexp ( _threadargscomma_ v ) ;
   ica = perm * m * m * m * ( 4.47814e6 * v / T ) * ( ( cai / 1000.0 ) - ( cao / 1000.0 ) * A ) / ( 1.0 - A ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
  cao = _ion_cao;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
   _t_minf = makevector(301*sizeof(double));
   _t_taum = makevector(301*sizeof(double));
   _t_getGHKexp = makevector(301*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "dcn_CaHVA.mod";
static const char* nmodl_file_text = 
  "TITLE High voltage activated calcium current (CaHVA) of deep cerebellar nucleus (DCN) neuron\n"
  "COMMENT\n"
  "    Translated from GENESIS by Johannes Luthman and Volker Steuber. \n"
  "    This mechanism and the other calcium channel (CaLVA.mod) are the only channel\n"
  "    mechanisms of the DCN model that use the GHK mechanism to calculate reversal\n"
  "    potential. Thus, extracellular Ca concentration is of importance and shall be \n"
  "    set from hoc.\n"
  "\n"
  "    Calcium entering the neuron through this channel is kept track of by CaConc.mod\n"
  "    and affects the conductance of the SK.mod channel mechanism.\n"
  "ENDCOMMENT \n"
  "\n"
  "NEURON { \n"
  "	SUFFIX dcnCaHVA \n"
  "	USEION ca READ cai, cao WRITE ica\n"
  "	RANGE perm, ica, m, cai\n"
  "	GLOBAL qdeltat\n"
  "} \n"
  " \n"
  "UNITS { \n"
  "	(mA) = (milliamp) \n"
  "	(mV) = (millivolt) \n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "} \n"
  " \n"
  "PARAMETER { \n"
  "    qdeltat = 1\n"
  "    perm = 7.5e-6 (cm/seconds)\n"
  "} \n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "    cai (mM)\n"
  "    cao (mM)     \n"
  "	ica (mA/cm2) \n"
  "	minf (1)\n"
  "	taum (ms) \n"
  "	celsius (degC)\n"
  "	T (kelvin)\n"
  "    A (1)\n"
  "} \n"
  "\n"
  "STATE {\n"
  "	m\n"
  "} \n"
  "\n"
  "INITIAL { \n"
  "    T = 273.15 + celsius\n"
  "    rate(v)\n"
  "    m = minf \n"
  "} \n"
  " \n"
  "BREAKPOINT { \n"
  "    SOLVE states METHOD cnexp \n"
  "    A = getGHKexp(v)\n"
  "    : \"4.47814e6 * v / T\" in the following is the simplification of the GHK\n"
  "    : current equation's (z^2 * F^2 * (0.001) * v) / (R * T). [*(0.001) is to get \n"
  "    : volt from NEURON's mV.] Together with the simplification in getGHKexp() \n"
  "    : (below), this speeds up the whole DCN simulation (without synapses) by 8%.\n"
  "    : The division of the calcium concentrations (mM) by 1000 gives molar as \n"
  "    : required by the GHK current equation.\n"
  "    ica = perm * m*m*m * (4.47814e6 * v / T) * ((cai/1000) - (cao/1000) * A) / (1 - A)\n"
  "} \n"
  " \n"
  "DERIVATIVE states { \n"
  "	rate(v) \n"
  "	m' = (minf - m)/taum \n"
  "} \n"
  "\n"
  "PROCEDURE rate(v(mV)) {\n"
  "	TABLE minf, taum FROM -150 TO 100 WITH 300\n"
  "	if(fabs(v + 34.5) < 1e-6) {\n"
  "		minf = 1 / (1 + exp((v + 34.50001) / -9))\n"
  "	} else {\n"
  "	minf = 1 / (1 + exp((v + 34.5) / -9))\n"
  "	}\n"
  "    taum = 1 / ((31.746 * ((exp((v - 5) / -13.89) + 1) ^ -1)) + (3.97e-4 * (v + 8.9)) * ((exp((v + 8.9) / 5) - 1) ^ -1))\n"
  "    taum = taum / qdeltat\n"
  "} \n"
  "\n"
  "FUNCTION getGHKexp(v(mV)) {\n"
  "    TABLE DEPEND T FROM -150 TO 100 WITH 300 \n"
  "    getGHKexp = exp(-23.20764929 * v / T): =the calculated values of\n"
  "            : getGHKexp = exp((-z * F * (0.001) * v) / (R * T)).\n"
  "}\n"
  ;
#endif
