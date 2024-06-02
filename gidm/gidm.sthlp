{smcl}
{* *! version 1.0.0  18jan2019}{...}
{viewerjumpto "Title" "gidm##title"}{...}
{viewerjumpto "Syntax" "gidm##syntax"}{...}
{viewerjumpto "Description" "gidm##description"}{...}
{viewerjumpto "Options" "gidm##options"}{...}
{viewerjumpto "Examples" "gidm##examples"}{...}
{viewerjumpto "Stored results" "gidm##results"}{...}
{viewerjumpto "Authors" "gidm##authors"}{...}
{viewerjumpto "Also see" "gidm##alsosee"}{...}
{cmd:help gidm}{right: ({browse "https://doi.org/10.1177/1536867X19874246":SJ19-3: st0574})}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 13 15 2}{...}
{p2col :{cmd:gidm} {hline 2}}Generalized inflated discrete model


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:gidm}
{cmd:(}{depvar} {indepvars}{cmd:)} 
{cmd:(}{varlist:1}{cmd:)}
[...
{cmd:(}{varlist:N}{cmd:)}]
{ifin}
[{it:{help gidm##weight:weight}}]{cmd:,} 
{opt inf:lation}{cmd:(}{it:numlist}{cmd:)} 
{opt lin:k}{cmd:(}{it:string}{cmd:)} 
[{cmdab:noinit:ial}
{opt vce(vcetype)} 
{opt l:evel(#)} {it:display_options} {it:maximize_options}]

{p 4 6 2}
{it:indepvars} and {it:varlist} may contain factor variables; see 
{help fvvarlist}.{p_end}
{p 4 6 2}
{cmd:bayes}, {cmd:bootstrap}, {cmd:by}, {cmd:fp}, {cmd:jackknife},
{cmd:rolling}, {cmd:statsby}, and {cmd:svy} are allowed; see {help prefix}.{p_end}
{p 4 6 2}
Weights are not allowed with the {helpb bootstrap} prefix.{p_end}
{p 4 6 2}
{opt vce()} and weights are not allowed with the {helpb svy}
prefix.{p_end}
{marker weight}{...}
{p 4 6 2}
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are allowed; see
{help weight}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
The {cmd:gidm} command fits a generalized inflated discrete model of
{it:depvar} on several sets of {it:indepvars} and {it:varlistN}.  The
{it:depvar} is a nonnegative integer of the response variable.  {it:indepvars}
is a set of explanatory variables for {it:depvar}, whereas {it:varlist1} to
{it:varlistN} are sets of explanatory variables for modeling the probabilities
of inflation at each of the points corresponding to the values specified in
the {it:numlist} in the option {opt inflation(numlist)}.  Specifically,
{cmd:(_con)} corresponds to an intercept-only model.


{marker options}{...}
{title:Options}

{phang}
{opt inflation(numlist)} specifies the list of values at which the inflations
are assumed.  The number of elements in {it:numlist} must be the same as the
number of equations specified by {indepvars} and {varlist:1} ... {varlist:N}.
{cmd:inflation()} is required.

{phang}
{opt link(string)} defines the distribution for both of the noninflated and
the inflated parts.  We use a four-letter combination to represent each model.
The first two letters, for example, {cmd:lg} for logit and {cmd:pb} for
probit, indicate the functional form used for the inflated part, and the last
two letters refer to the distribution of outcome.  The supported distributions
for outcomes are Poisson ({cmd:po}), negative binomial ({cmd:nb}), multinomial
({cmd:ml}), cumulative logit ({cmd:cl}), and cumulative probit ({cmd:cp}).
For instance, the keyword {cmd:lgpo} refers to a logit-inflated Poisson, and
{cmd:pbcp} is a probit-inflated cumulative probit.  {cmd:link()} is required.
A summary of the keywords of models supported by the {cmd:gidm} command is
given below.
		
{p 8 8 2}
{c TLC}{hline 11}{c TT}{hline 20}{c TT}{hline 21}{c -}{hline 21}{c TRC}{p_end}
{p 8 8 2}
{c |}Outcome{space 3} {c |}Model{space 15}{c |}{space 15} Option {opt link:}{cmd:(}{it:string}{cmd:)} {space 7}{c |}{p_end}
{p 8 8 2}
{c |}{space 11}{c |}{space 20}{c LT}{hline 21}{c TT}{hline 21}{c RT}{p_end}
{p 8 8 2}
{c |}{space 10} {c |}{space 20}{c |}Logit inflation{space 6}{c |}Probit inflation{space 5}{c |}{p_end}
{p 8 8 2}
{c LT}{hline 11}{c +}{hline 20}{c +}{hline 21}{c +}{hline 21}{c RT}{p_end}
{p 8 8 2}
{c |}Count{space 5} {c |}Poisson{space 13}{c |}{cmd:lgpo}{space 16} {c |}{cmd:pbpo}{space 17}{c |}{p_end}
{p 8 8 2}
{c |}Count{space 5} {c |}Negative binomial{space 3}{c |}{cmd:lgnb}{space 16}
{c |}{cmd:pbnb}{space 17}{c |}{p_end}
{p 8 8 2}
{c |}Category{space 2} {c |}Multinomial{space 9}{c |}{cmd:lgml}{space 16} {c |}{cmd:pbml}{space 17}{c |}{p_end}
{p 8 8 2}
{c |}Category{space 2} {c |}Ordered logit{space 7}{c |}{cmd:lgcl}{space 16} {c |}{cmd:pbcl}{space 17}{c |}{p_end}
{p 8 8 2}
{c |}Category{space 2} {c |}Ordered probit{space 6}{c |}{cmd:lgcp}{space 16} {c |}{cmd:pbcp}{space 17}{c |}{p_end}
{p 8 8 2}
{c BLC}{hline 11}{c BT}{hline 20}{c BT}{hline 21}{c BT}{hline 21}{c BRC}{p_end}

{phang}
{opt noinitial} suppresses the default initial values that are from results of
the separately fit model parts.  For example, with {cmd:link(lgpo)}, the
default initial values are obtained from a separately fit Poisson model for
the main part and logistic regressions for the inflated parts.

INCLUDE help vce_asymptall

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

INCLUDE help displayopts_list

{phang}
{it:maximize_options}: {opt dif:ficult},
{opth tech:nique(maximize##algorithm_spec:algorithm_spec)}, {opt iter:ate(#)},
[{cmd:{ul:no}}]{cmd:{ul:lo}}{cmd:g}, {opt tr:ace}, {opt grad:ient},
{opt showstep}, {opt hess:ian}, {opt showtol:erance}, {opt tol:erance(#)},
{opt ltol:erance(#)}, {opt nrtol:erance(#)}, {opt nonrtol:erance}, and
{opt from(init_specs)}; see {manhelp maximize R}.  These options are seldom
used.

{pmore}
Setting the optimization type to {cmd:technique(bhhh)} resets the default
{it:vcetype} to {cmd:vce(opg)}.


{marker examples}{...}
{title:Examples}

{pstd}
1) Inflated models for count {it:depvar}{p_end}

{pstd}
Setup{p_end}
{phang2}{cmd:. webuse fish}{p_end}

{pstd}
Zero-inflated Poisson regression with logit link for the inflation{p_end}
{phang2}{cmd:. gidm (count child camper) (persons), inflation(0) link(lgpo)}{p_end}

{pstd}
Zero-inflated Poisson regression with probit link for the inflation{p_end}
{phang2}{cmd:. gidm (count child camper) (persons), inflation(0) link(pbpo)}{p_end}

{pstd}
Zero-inflated negative binomial regression with logit link for the
inflation{p_end}
{phang2}{cmd:. gidm (count child camper) (persons), inflation(0) link(lgnb)}{p_end}

{pstd}
Zero-inflated negative binomial regression with probit link for the
inflation{p_end}
{phang2}{cmd:. gidm (count child camper) (persons), inflation(0) link(pbnb)}{p_end}

{pstd}
Generalized inflated Poisson regression with logit link for the inflated
points{p_end}
{phang2}{cmd:. gidm (count child camper) (persons) (persons), inflation(0 1) link(lgpo)}{p_end}

{pstd}
Generalized inflated Poisson regression with probit link for the inflated
points{p_end}
{phang2}{cmd:. gidm (count child camper) (persons) (persons), inflation(0 1) link(pbpo)}{p_end}

{pstd}
Generalized inflated negative binomial regression with logit link for the
inflated points{p_end}
{phang2}{cmd:. gidm (count child camper) (persons) (persons), inflation(0 1) link(lgnb)}{p_end}

{pstd}
Generalized inflated negative binomial regression with probit link for the
inflated points{p_end}
{phang2}{cmd:. gidm (count child camper) (persons) (persons), inflation(0 1) link(pbnb)}{p_end}

{pstd}
2) Inflated models for ordered {it:depvar}{p_end}

{pstd}
Setup{p_end}
{phang2}{cmd:. use "https://stats.idre.ucla.edu/stat/data/ologit.dta", clear}{p_end}

{pstd}
Zero-inflated ordered logistic model with logit link for the inflation{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared), inflation(0) link(lgcl)}{p_end}

{pstd}
Zero-inflated ordered logistic model with probit link for the inflation{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared), inflation(0) link(pbcl)}{p_end}

{pstd}
Zero-inflated ordered probit model with logit link for the inflation{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared), inflation(0) link(lgcp)}{p_end}

{pstd}
Zero-inflated ordered probit model with probit link for the inflation{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared), inflation(0) link(pbcp)}{p_end}

{pstd}
Generalized inflated ordered logistic model with logit link for the inflated
points{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared) (pared), inflation(0 1) link(lgcl)}{p_end}

{pstd}
Generalized inflated ordered logistic model with probit link for the inflated
points{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared) (pared), inflation(0 1) link(pbcl)}{p_end}

{pstd}
Generalized inflated ordered probit model with logit link for the inflated
points{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared) (pared), inflation(0 1) link(lgcp)}{p_end}

{pstd}
Generalized inflated ordered logistic model with probit link for the inflated
points{p_end}
{phang2}{cmd:. gidm (apply public gpa) (pared) (pared), inflation(0 1) link(pbcp)}{p_end}

{pstd} 
3) Inflated models for categorical {it:depvar}{p_end}

{pstd}
Setup{p_end}
{phang2}{cmd:. use "https://stats.idre.ucla.edu/stat/data/hsbdemo", clear}{p_end}

{pstd}
Inflated multinomial model with logit link for the inflation{p_end}
{phang2}{cmd:. gidm (prog  read) (math), inflation(3) link(lgml)}{p_end}

{pstd}
Inflated multinomial model with logit link for the inflation{p_end}
{phang2}{cmd:. gidm (prog  read) (math), inflation(3) link(pbml)}{p_end}

{pstd}
Generalized inflated multinomial model with logit link for the inflations{p_end}
{phang2}{cmd:. gidm (prog  read) (math) (math), inflation(2 3) link(lgml)}{p_end}

{pstd}
Generalized inflated multinomial model with probit link for the
inflations{p_end}
{phang2}{cmd:. gidm (prog  read) (math) (math), inflation(2 3) link(pbml)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:gidm} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(p)}}significance of model test{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:gidm}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared test{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(which)}}{cmd:max} or {cmd:min}; whether optimizer is to perform maximization or minimization{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker authors}{...}
{title:Authors}

{pstd}
Yiwei Xia{break}
Southwestern University of Finance and Economics{break}
Chengdu, China {break}

{pstd}
Yisu Zhou{break}
Faculty of Education{break}
University of Macau{break}
Macau, China

{pstd}
Tianji Cai{break}
Department of Sociology{break}
University of Macau{break}
Macau, China{break}
tjcai@um.edu.mo


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 19, number 3: {browse "https://doi.org/10.1177/1536867X19874246":st0574}{p_end}
