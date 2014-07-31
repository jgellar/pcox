Notes 2014-07-31
========================================================

## IDEAS/SUGGESTIONS

#### non-invertible penalties into `coxme`
Two possible approaches that might work with `coxme` _as it is_ come to mind: 

1. replace the zero eigenvalues of the penalties with very small non-zeroes $\rightarrow$ that instantly gives you positive definiteness, at the cost of changing the pen. slightly. E.g. for a first-differences penalty, this adds a **small** ridge penalty. This is also the idea behind `mgcv`'s `cs` and `ts` bases (c.f. Marra/Wood "Practical variable selection for GAMs") and the `pss`-constructor defined in `refund()` for `ff` -terms. Those should work with `coxme`. If you make the replacement eigenvalues just barely large enough to ensure numerical full rank (i.e., very very large variance versus infinite variance on the unpenalized function space), the change shouldn't really affect the fits for non-pathological cases.
* do the classic mixed model decomposition as in Section 6.6.1 of Wood (2006) and split off the unpenalized parts. 
    * messier implementation because you now have a single smooth term that has both     parametric/fixed effects and penalized/random effects, and for visualization/testing    etc you have to put them back together.
    * for tensor product splines this involves some additional design choices: do you do this split for each marginal basis (--> `t2()` in `mgcv`, c.f. Wood/Scheipl/Faraway(2013)) or do you do this after merging the marginal bases & penalties (c.f. Section 6 in Currie/Durban/Eilers (2006))

--> Good suggestions. I still would like to push Therneau to modify coxme to allow non-invertable penalties - he said he was willing to do this.
     
#### use `tv(...)` to indicate time-varying terms in the `pcox`-formula 
* no name clash/confusion with `survival::tt()`
* more intuitive naming

--> I agree


#### don't repeat all my pffr()-mistakes
some things I would do differently if I were to start over:
* coding style: consistent naming schemes, no use of periods as separators in names, more functional encapsulation especially for the pre- and post-processing, ...
* write more & better comments: explain how & why, not what. Right now you have way too much code that is not commented or documented at all, IMO. 
* stronger modification of the return object instead of simply adding more stuff to
the object returned by `mgcv`: easier to write methods, more memory efficient (`pffr` objects are often huge because of lots of duplication)
* debatable: doing most pre-processing etc. yourself and calling the low-level fitting functions (`mgcv:::gam.fit`, in yor case: `survival:::coxpenal.fit`) directly instead of calling `gam` or `coxph`. Advantages: not limited to functionality envisioned by package creator as long as you can define suitable design matrices etc., more control, less workarounds/jumping through hoops. Disadvantages: more code to maintain, partially re-inventing the wheel, calling a non-exported function (not sure if that is even still admissible on CRAN these days)  

--> For point #3 above regarding the return object, I was hoping to discuss this with you at some point precisely because of your experience with pffr. This is related to how to best design methods such as coef, plot, etc.

--> As for using `survival:::coxpenal.fit`, I don't think I want to go down that road unless necessary. Seems like more work overall. `coxph` has some nice functionality built into it that I would like to maintain, see e.g. strata() and cluster() terms.

## ISSUES
* `lf.vd.cox.R` defines the same function twice -- I guess the first one should be 
`lf.vd()` and the 2nd one `lf.vd.cox()`?
* `s.cox()` is defined twice as well: in `pterm.R` and `s.cox.R`
* there isn't a single example of using `pcox()` in the whole package 
 (or in Testing, for that matter)...
* deal with R CMD check errors/warnings
    * make folder structure of pcox more standard R-package like (move `Document` to `inst/doc` and `Testing` to `inst/tests` (for proper tests) or `inst/misc` for "sandbox"-scripts, e.g.)
    * don't use "<<-" !!
    * NAMESPACE imports and other Rd issues fixed (see pull request #??)
* turn `Document/pcox.tex` into a proper vignette (Rnw file) so examples etc can be included easily?

--> Actually, my intention is to get rid of `lf.vd.cox`, `s.cox`, etc. and just use `pterm` for generating all of the penalized terms in the format that `coxph` recognizes. These were in there because they were my first attempts - they can probably be deleted now.

--> As far as examples/tests for `pcox`, I haven't done any of that yet, and have mostly been working directly with `coxph` to fit all the models. This is mostly because `pcox` is far from finished now. I really haven't done anything in there related to the tt() terms.

--> Agreed on the package structure, things like the `Testing` folder weren't meant to be proper tests in the final version, they were more just me messing around and trying to figure out how to fit different models.

--> I only included `<<-` as a temporary workaround, so I am able to extract the smooth that comes from smoothCon(), and use it to recreate the functional estiamtes. In the final version, this will all be contained within `pcox`. I was thinking that `pcox` would return the list of all smooth terms as the element `smooth`, just like mgcv.

#### `simSurvTVC.R`
* `predict(fit2)` gives a 17254-vector for a 500 row dataset?!? more generally, all fits with `tt` terms return weird predictions in looong vectors -- is that a (`coxph`) bug or am I misunderstanding something?

--> `predict.coxph` doesn't work with models that contain tt() terms. We will eventually have to write a `predict.pcox` method, which will be able to do this. If you step through the code for coxph() with tt() terms, you'll see that including tt() terms causes coxph to duplicate each row, for all event times for which it is in the "risk set" of the Cox model. It looks like this results in 17254 rows in your example.

#### time varying terms via `pcox()`:
1. `pcox(Surv(time,event) ~ X1 + tt(X2), data=data1)` gives a time constant term for
`X2` unless you also specify `tt=function(x,t,...) x*s.cox(t)`, but then why have
the pcox wrapper at all -- that's the same syntax as for coxph?
(also see 2nd idea above)
* `lf.vd()` & `af.vd()` don't exist (yet?) in refund or refundDevel, AFAICS.
* special term type names don't correspond to those in `pcox.tex`
* for the formula interface to work like the spec in `pcox.tex`, I think we'll
need two steps following l. 87: one to get all the time-varying terms and to check
what kind they are (linear/nonlinear scalar, linear/nonlinear concurrent, linear
/nonlinear functional (vd yes/no), linear/nonlinear historical functional 
(vd yes/no)), and then a second one to gather all of them and the time-constant
special terms for further processing (`where.ttlf`, `where.lf` etc.). ATM, 
`tt(x)` and a `tt(s(x))` terms would be processed the same way before being sent 
on to the fitter, so that's probably not going to work. What we would want for the
time-varying terms is to be replaced by `tt`-terms that `coxph` understands with the appropriate `tt`-function calling `s.cox` & friends ... right? Alternatively, you
could allow special terms `tts(x)` (for `tt(s(x))`)  `ttlf.vd(x)` etc. in the
formula. I think the first approach as given in `pcox.tex` is cleaner.

--> I have not updated `pcox` to allow for time-varying terms yet. You can see that I have figured out how to fit these models in coxph(), but I just haven't written the code to parse the formula and correctly create the appropriate term to pass to coxph.

--> I agree with your long point above, I was thinking along these lines, but haven't worked out the details yet. I also prefer `tv(s(x))` to `tts(x)`. One thing I was thinking would be nice is if there was a way to add another dimension to an existing smooth.spec object. e.g., `tv(s(x))` would know to add an interaction with time to the s(x) object, just like `tv(af(x))` would do to the existing bivariate smooth. We could have an option as to whether this is done with a tensor product approach or with a bi/tri-variate smoother such as TPRS. The only way I have come up with to implement this is to build in logic that does something separately for each of the "internal" smooths. If you can think of a way to re-use the same code to do this regardless of what the internal smooth.spec is, that would be fantastic.

--> You're right that `lf.vd` and `af.vd` are not in refund or refundDevel. I've written `lf.vd` for use in `fgam`, I just haven't gotten around to pushing it to refundDevel yet.

  
## QUESTIONS 

* l.14 in `s.cox.R`: `lambda <- ifelse(theta<=0, 0, theta/(1-theta))`: why is `theta` guaranteed to be in $[-\infty, 1)$ -- what kind of weird transform of `lambda` is this?
* do we really want all the simulation code in the package itself? maybe put the `genBeta` stuff into the `demo` folder or similar?

--> I blatently ripped off the code for fitting penalized Cox models from Therneau's functions for ridge regression and frailty models. See his 1998 technical report titled "Penalized Cox models and frailty" for details (not to be confused with his JCGS paper on "Penalized Survival models and frailty"). His code uses a "theta" smoothing parameter that is scaled from 0 to 1, and then converts it to a "lambda" that goes from 0 to infinity. Since all his optimization routines use the "theta", I found it easier to just re-use his code here.

--> The genBeta stuff was added because I am running some simulations on the historical Cox model as we speak, and it was easiest/cleanest to insert the code in there, and just install pcox on our computing cluster here. They won't end up in the final version.
