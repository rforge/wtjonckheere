#include "wtJonckheere.hpp"

// delivers a jonckheere terpstra statistic

SEXP rcpp_wt_jonckheere(SEXP values, SEXP groups, SEXP weights)
{
    using namespace Rcpp ;
    
    NumericVector cWerte(values), cGewichte(weights);
    IntegerVector cGruppe(groups);
  
    double summe = 0;
  
    NumericVector::iterator iWert = cWerte.begin(), jWert,  iGewichte = cGewichte.begin(), jGewichte;
    IntegerVector::iterator iGruppe = cGruppe.begin(), jGruppe;
    
    for (; iWert != cWerte.end(); iWert++, iGruppe++, iGewichte++)
	    for (jWert = cWerte.begin(), jGruppe = cGruppe.begin(), jGewichte = cGewichte.begin(); jWert != cWerte.end(); jWert++, jGruppe++, jGewichte++)
    {
	    if (*iGruppe < *jGruppe) {
			double refGewicht = (*iGewichte) * (*jGewichte);
			if (*iWert < * jWert) summe += refGewicht; else if (*iWert == *jWert) summe += .5*refGewicht;
	    }
    }
    return wrap(summe);     
}
