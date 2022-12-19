/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
 
#define nrn_init _nrn_init__pcCaBK
#define _nrn_initial _nrn_initial__pcCaBK
#define nrn_cur _nrn_cur__pcCaBK
#define _nrn_current _nrn_current__pcCaBK
#define nrn_jacob _nrn_jacob__pcCaBK
#define nrn_state _nrn_state__pcCaBK
#define _net_receive _net_receive__pcCaBK 
#define rates rates__pcCaBK 
#define states states__pcCaBK 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkbar _p[0]
#define ik _p[1]
#define gk _p[2]
#define m _p[3]
#define z _p[4]
#define h _p[5]
#define ek _p[6]
#define cai _p[7]
#define qt _p[8]
#define Dm _p[9]
#define Dz _p[10]
#define Dh _p[11]
#define _g _p[12]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
#define _ion_cai	*_ppvar[3]._pval
 
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
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
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
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_pcCaBK", _hoc_setdata,
 "rates_pcCaBK", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define hinf hinf_pcCaBK
 double hinf = 0;
#define minf minf_pcCaBK
 double minf = 0;
#define tauz tauz_pcCaBK
 double tauz = 0;
#define tauh tauh_pcCaBK
 double tauh = 0;
#define taum taum_pcCaBK
 double taum = 0;
#define zhalf zhalf_pcCaBK
 double zhalf = 0.001;
#define zinf zinf_pcCaBK
 double zinf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "zhalf_pcCaBK", "mM",
 "taum_pcCaBK", "ms",
 "tauh_pcCaBK", "ms",
 "tauz_pcCaBK", "ms",
 "gkbar_pcCaBK", "S/cm2",
 "ik_pcCaBK", "mA/cm2",
 "gk_pcCaBK", "S/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 static double z0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "zhalf_pcCaBK", &zhalf_pcCaBK,
 "minf_pcCaBK", &minf_pcCaBK,
 "taum_pcCaBK", &taum_pcCaBK,
 "hinf_pcCaBK", &hinf_pcCaBK,
 "tauh_pcCaBK", &tauh_pcCaBK,
 "zinf_pcCaBK", &zinf_pcCaBK,
 "tauz_pcCaBK", &tauz_pcCaBK,
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
"pcCaBK",
 "gkbar_pcCaBK",
 0,
 "ik_pcCaBK",
 "gk_pcCaBK",
 0,
 "m_pcCaBK",
 "z_pcCaBK",
 "h_pcCaBK",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.014;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* cai */
 
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

 void _pc_CaBK_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 pcCaBK D:/Projects/SchreglmannEtAl2020/CCTC_model/modfiles/pc_CaBK.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double q10 = 2.2;
 static double cvm = 28.9;
 static double ckm = 6.2;
 static double ctm = 0.000505;
 static double cvtm1 = 86.4;
 static double cktm1 = -10.1;
 static double cvtm2 = -33.3;
 static double cktm2 = 10;
 static double ctauz = 1;
 static double ch = 0.085;
 static double cvh = 32;
 static double ckh = -5.8;
 static double cth = 0.0019;
 static double cvth1 = 48.5;
 static double ckth1 = -5.2;
 static double cvth2 = -54.2;
 static double ckth2 = 12.9;
static int _reset;
static char *modelname = "BK-type Purkinje calcium-activated potassium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / taum ;
   Dz = ( zinf - z ) / tauz ;
   Dh = ( hinf - h ) / tauh ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taum )) ;
 Dz = Dz  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauz )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauh )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taum)))*(- ( ( ( minf ) ) / taum ) / ( ( ( ( - 1.0 ) ) ) / taum ) - m) ;
    z = z + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauz)))*(- ( ( ( zinf ) ) / tauz ) / ( ( ( ( - 1.0 ) ) ) / tauz ) - z) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauh)))*(- ( ( ( hinf ) ) / tauh ) / ( ( ( ( - 1.0 ) ) ) / tauh ) - h) ;
   }
  return 0;
}
 
static int  rates (  double _lv ) {
   _lv = _lv + 5.0 ;
   minf = 1.0 / ( 1.0 + exp ( - ( _lv + cvm ) / ckm ) ) ;
   taum = ( 1e3 ) * ( ctm + 1.0 / ( exp ( - ( _lv + cvtm1 ) / cktm1 ) + exp ( - ( _lv + cvtm2 ) / cktm2 ) ) ) / qt ;
   zinf = 1.0 / ( 1.0 + zhalf / cai ) ;
   tauz = ctauz / qt ;
   hinf = ch + ( 1.0 - ch ) / ( 1.0 + exp ( - ( _lv + cvh ) / ckh ) ) ;
   tauh = ( 1e3 ) * ( cth + 1.0 / ( exp ( - ( _lv + cvth1 ) / ckth1 ) + exp ( - ( _lv + cvth2 ) / ckth2 ) ) ) / qt ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  z = z0;
 {
   qt = pow( q10 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   z = zinf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 v = _v;
  ek = _ion_ek;
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gk = gkbar * pow( m , 3.0 ) * pow( z , 2.0 ) * h ;
   ik = gk * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  ek = _ion_ek;
  cai = _ion_cai;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  ek = _ion_ek;
  cai = _ion_cai;
 { error =  states();
 if(error){fprintf(stderr,"at line 105 in file pc_CaBK.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(z) - _p;  _dlist1[1] = &(Dz) - _p;
 _slist1[2] = &(h) - _p;  _dlist1[2] = &(Dh) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "pc_CaBK.mod";
static const char* nmodl_file_text = 
  "TITLE BK-type Purkinje calcium-activated potassium current\n"
  "\n"
  "COMMENT\n"
  "\n"
  "	NEURON implementation of a BK-channel in Purkinje cells\n"
  "	Kinetical Scheme: Hodgkin-Huxley (m^3*z^2*h)\n"
  "\n"
  "	Modified from Khaliq et al., J.Neurosci. 23(2003)4899\n"
  "	 \n"
  "	Laboratory for Neuronal Circuit Dynamics\n"
  "	RIKEN Brain Science Institute, Wako City, Japan\n"
  "	http://www.neurodynamics.brain.riken.jp\n"
  "\n"
  "	Reference: Akemann and Knoepfel, J.Neurosci. 26 (2006) 4602\n"
  "	Date of Implementation: May 2005\n"
  "	Contact: akemann@brain.riken.jp\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "       SUFFIX pcCaBK\n"
  "       USEION k READ ek WRITE ik\n"
  "       USEION ca READ cai\n"
  "       RANGE gkbar, gk,  ik\n"
  "       GLOBAL minf, taum, hinf, tauh, zinf, tauz\n"
  "       GLOBAL zhalf\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "	(nA) = (nanoamp)\n"
  "	(pA) = (picoamp)\n"
  "	(S)  = (siemens)\n"
  "	(nS) = (nanosiemens)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)		\n"
  "}\n"
  "\n"
  "CONSTANT {\n"
  "	q10 = 2.2\n"
  "	\n"
  "	cvm = 28.9 (mV)\n"
  "	ckm = 6.2 (mV)\n"
  "\n"
  "	ctm = 0.000505 (s)\n"
  "	cvtm1 = 86.4 (mV)\n"
  "	cktm1 = -10.1 (mV)\n"
  "	cvtm2 = -33.3 (mV)\n"
  "	cktm2 = 10 (mV)\n"
  "\n"
  "	ctauz = 1 (ms)\n"
  "\n"
  "	ch = 0.085\n"
  "	cvh = 32 (mV)\n"
  "	ckh = -5.8 (mV)\n"
  "	cth = 0.0019 (s)\n"
  "	cvth1 = 48.5 (mV)\n"
  "	ckth1 = -5.2 (mV)\n"
  "	cvth2 = -54.2 (mV)\n"
  "	ckth2 = 12.9 (mV)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	v (mV)\n"
  "	celsius (degC)\n"
  "\n"
  "	gkbar = 0.014 (S/cm2)\n"
  "\n"
  "	ek (mV)\n"
  "	cai (mM)\n"
  "\n"
  "	zhalf = 0.001 (mM)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ik (mA/cm2)\n"
  "	qt\n"
  "      gk (S/cm2)   \n"
  "	minf\n"
  "	taum (ms)\n"
  "	hinf\n"
  "	tauh (ms)\n"
  "	zinf\n"
  "      tauz (ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m   FROM 0 TO 1\n"
  "	z   FROM 0 TO 1\n"
  "	h   FROM 0 TO 1\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	qt = q10^((celsius-22 (degC))/10 (degC))\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	z = zinf\n"
  "	h = hinf\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "      gk = gkbar * m^3 * z^2 * h      \n"
  "	ik = gk * (v - ek)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v)\n"
  "	m' = (minf-m)/taum\n"
  "	z' = (zinf-z)/tauz\n"
  "	h' = (hinf-h)/tauh\n"
  "}\n"
  "\n"
  "PROCEDURE rates( v (mV) ) {\n"
  "	v = v + 5 (mV)\n"
  "	minf = 1 / ( 1+exp(-(v+cvm)/ckm) )\n"
  "	taum = (1e3) * ( ctm + 1 (s) / ( exp(-(v+cvtm1)/cktm1) + exp(-(v+cvtm2)/cktm2) ) ) / qt\n"
  "	\n"
  "	zinf = 1 /(1 + zhalf/cai)\n"
  "      tauz = ctauz/qt\n"
  "\n"
  "	hinf = ch + (1-ch) / ( 1+exp(-(v+cvh)/ckh) )\n"
  "	tauh = (1e3) * ( cth + 1 (s) / ( exp(-(v+cvth1)/ckth1) + exp(-(v+cvth2)/ckth2) ) ) / qt\n"
  "}\n"
  ;
#endif
